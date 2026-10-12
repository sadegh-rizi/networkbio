"""Stage 04: preregistered CORNETO runs and stage-04 readouts.

This script performs the real-data stage only after constructing the TF and
metabolite tables from stages 01-03.  It does not add biological controls by
hand and it does not cap path length before CORNETO preprocessing.
"""

from __future__ import annotations

import hashlib
import inspect
import json
import re
import subprocess
import sys
from collections import deque
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from activities import (  # noqa: E402
    leave_one_line_out_ulm,
    load_regulon,
    moderated_t_statistics,
    select_stable_top_tfs,
)
from inference import (  # noqa: E402
    CANONICAL_CHOLESTEROL_GENES,
    CHOLESTEROL_HMDB,
    EXCLUDED_INPUTS,
    LAMBDA_GRID,
    COSMOS_TIME_LIMIT,
    OMNIPATH_TIME_LIMIT,
    SOLVER,
    SOLVER_RELATIVE_GAP,
    SOLVER_SEED,
    SOLVER_THREADS,
    choose_lambda,
    is_excluded_input,
    make_corneto_data,
    lambda_grid_table,
    lambda_label,
    read_pkn,
    sample_edges,
    solve_carnival,
)
from network_plot import draw_signed_network, network_legend_handles  # noqa: E402
from plot_config import apply_plot_style, save_figure  # noqa: E402
from preprocess import stage_provenance, write_stage_provenance  # noqa: E402

PLAN = "doc/decisions/2026-10-10-stage04-inference-plan.md"
ANALYSIS = "results/ionescu_corneto/04_inference"
ACTIVITIES = "results/ionescu_corneto/02_activities"
PREPROCESSED = "results/ionescu_corneto/01_preprocessing"
PKN = "results/ionescu_corneto/03_pkn"


def smoke_corneto() -> None:
    """Print and check the pinned API before any real MILP is started."""
    from corneto.backend._base import ProblemDef
    from corneto.methods import CarnivalFlow
    from corneto.methods import sampler

    print("CarnivalFlow signature:", inspect.signature(CarnivalFlow))
    print("ProblemDef.solve signature:", inspect.signature(ProblemDef.solve))
    print("HiGHS option mapping:", {
        "time_limit": "max_seconds -> time_limit",
        "relative_gap": "mip_rel_gap",
        "threads": "threads",
        "seed": "random_seed",
    })
    print("sample_alternative_solutions:", hasattr(sampler, "sample_alternative_solutions"))
    if not hasattr(sampler, "sample_alternative_solutions"):
        raise RuntimeError("Pinned CORNETO lacks sample_alternative_solutions")
    import highspy

    options = highspy.Highs()
    for name in ("time_limit", "mip_rel_gap", "threads", "random_seed"):
        status, _ = options.getOptionValue(name)
        if status.name != "kOk":
            raise RuntimeError(f"Pinned HiGHS lacks option {name}")


def _groups(columns: list[str]) -> pd.Series:
    groups = pd.Series(
        ["Ctrl" if column.startswith("Ctrl_") else "PMS" for column in columns],
        index=columns,
    )
    if groups.value_counts().to_dict() != {"Ctrl": 3, "PMS": 4}:
        raise ValueError(f"Expected 3 Ctrl and 4 PMS line columns, got {groups.value_counts().to_dict()}")
    return groups


def _pkn_paths(root: Path) -> dict[str, Path]:
    return {
        "omnipath": root / PKN / "omnipath/primary/pkn.tsv",
        "cosmos_tf": root / PKN / "cosmos/primary/pkn.tsv",
        "cosmos_multi": root / PKN / "cosmos/primary/pkn.tsv",
        "cosmos_repair": root / PKN / "cosmos/variants/canonical_repair/pkn.tsv",
        "omnipath_lines": root / PKN / "omnipath/primary/pkn.tsv",
    }


def _node_degrees(nodes_path: Path) -> dict[str, int]:
    nodes = pd.read_csv(nodes_path, sep="\t", dtype={"node": str})
    required = {"node", "in_degree", "out_degree"}
    if not required.issubset(nodes.columns):
        raise ValueError(f"Missing degree columns in {nodes_path}")
    return {
        str(row.node): int(row.in_degree) + int(row.out_degree)
        for row in nodes.itertuples(index=False)
    }


def _choose_metabolite_node(rows: pd.DataFrame, degrees: dict[str, int], preferred: str) -> tuple[str, str]:
    rows = rows[rows["mapped_to_pkn"].astype(bool)].copy()
    if rows.empty:
        raise ValueError("Measured metabolite has no COSMOS node mapping")
    rows = rows.drop_duplicates("pkn_node")
    preferred_rows = rows[rows["compartment"].eq(preferred)]
    if len(preferred_rows) == 1:
        return str(preferred_rows.iloc[0]["pkn_node"]), f"preferred_{preferred}"
    if len(preferred_rows) > 1:
        preferred_ranked = preferred_rows.assign(
            degree=preferred_rows["pkn_node"].map(degrees).fillna(-1)
        )
        maximum = preferred_ranked["degree"].max()
        candidates = sorted(preferred_ranked.loc[preferred_ranked["degree"].eq(maximum), "pkn_node"])
        if len(candidates) != 1:
            raise ValueError(f"Degree tie requires an explicit metabolite mapping decision: {candidates}")
        return str(candidates[0]), f"preferred_{preferred}_highest_degree"
    ranked = rows.assign(degree=rows["pkn_node"].map(degrees).fillna(-1))
    maximum = ranked["degree"].max()
    candidates = sorted(ranked.loc[ranked["degree"].eq(maximum), "pkn_node"])
    if len(candidates) != 1:
        raise ValueError(f"Degree tie requires an explicit metabolite mapping decision: {candidates}")
    return str(candidates[0]), "highest_degree"


def _metabolite_table(root: Path, *, sensitivity: bool = False) -> pd.DataFrame:
    """Compute primary or sensitivity metabolite contrasts and node choices."""
    mapping = pd.read_csv(root / PKN / "mapping/metabolite_mapping.tsv", sep="\t", dtype=str)
    mapping["mapped_to_pkn"] = mapping["mapped_to_pkn"].eq("True")
    degrees = _node_degrees(root / PKN / "cosmos/primary/nodes.tsv")
    rows: list[dict[str, Any]] = []

    for study, preferred, classes, values_name in (
        (
            "ST003331", "c",
            {"m0_unconfounded", "isotopologue_sum_partial", "m0_confounded"} if sensitivity
            else {"m0_unconfounded"},
            "line_log2.tsv",
        ),
        (
            "ST003332", "e",
            {"m0_unconfounded", "isotopologue_sum_partial", "m0_confounded"} if sensitivity
            else {"m0_unconfounded"},
            "line_log2_vs_blank.tsv",
        ),
    ):
        directory = root / PREPROCESSED / "metabolomics" / study
        features = pd.read_csv(directory / "features.tsv", sep="\t", dtype=str).set_index("feature_id")
        values = pd.read_csv(directory / values_name, sep="\t", index_col=0)
        features = features.loc[features.index.intersection(values.index)]
        features = features[features["labelling_class"].isin(classes)]
        features = features.loc[values.loc[features.index].notna().all(axis=1)]
        if features.empty:
            raise ValueError(f"No approved metabolite features remain for {study}")
        selected_nodes = {}
        for feature_id, feature in features.iterrows():
            candidates = mapping[(mapping["study"].eq(study)) & mapping["feature_id"].eq(feature_id)]
            if not candidates["mapped_to_pkn"].any():
                continue
            node, choice = _choose_metabolite_node(candidates, degrees, preferred)
            selected_nodes[feature_id] = (node, choice)
        if len(set(node for node, _ in selected_nodes.values())) != len(selected_nodes):
            raise ValueError(f"Multiple measured features map to the same COSMOS node in {study}")
        line_groups = _groups(values.columns.tolist())
        # Empirical-Bayes prior from every complete feature of the study, then
        # subset to the measured features (audit 2026-10-10, finding 2).
        complete_values = values.dropna()
        contrasts, prior = moderated_t_statistics(complete_values, line_groups)
        for feature_id, (node, choice) in selected_nodes.items():
            rows.append({
                "study": study,
                "feature_id": feature_id,
                "feature_name": features.loc[feature_id, "metabolite_name"],
                "labelling_class": features.loc[feature_id, "labelling_class"],
                "pkn_node": node,
                "mapping_choice": choice,
                "value": float(contrasts.loc[feature_id, "moderated_t"]),
                "df_g": prior["df_g"],
                "d0": prior["d0"],
                "s0_squared": prior["s0_squared"],
                "n_features_prior": len(complete_values),
            })

    directory = root / PREPROCESSED / "metabolomics/ST003328"
    features = pd.read_csv(directory / "features.tsv", sep="\t", dtype=str).set_index("feature_id")
    values = pd.read_csv(directory / "line_log2.tsv", sep="\t", index_col=0)
    cholesterol = features[features["metabolite_name"].str.lower().eq("cholesterol")]
    if len(cholesterol) != 1:
        raise ValueError("Expected exactly one free-cholesterol feature in ST003328")
    feature_id = str(cholesterol.index[0])
    untreated = [column for column in values.columns if column.endswith("_untreated")]
    line_groups = _groups(untreated)
    # Prior from all complete untreated lipid species, not from one feature.
    complete_lipids = values.loc[:, untreated].dropna()
    contrast, prior = moderated_t_statistics(complete_lipids, line_groups)
    if feature_id not in contrast.index:
        raise ValueError("Free cholesterol has missing untreated line values")
    candidates = mapping[(mapping["study"].eq("ST003328")) & mapping["feature_id"].eq(feature_id)]
    node, choice = _choose_metabolite_node(candidates, degrees, "c")
    rows.append({
        "study": "ST003328",
        "feature_id": feature_id,
        "feature_name": "cholesterol",
        "labelling_class": "not_labelled",
        "pkn_node": node,
        "mapping_choice": choice + "; untreated_only",
        "value": float(contrast.loc[feature_id, "moderated_t"]),
        "df_g": prior["df_g"],
        "d0": prior["d0"],
        "s0_squared": prior["s0_squared"],
        "n_features_prior": len(complete_lipids),
    })
    return pd.DataFrame(rows)


def _tf_inputs(root: Path, pkns: dict[str, list[tuple[str, int, str]]]) -> tuple[dict[str, dict[str, dict[str, float]]], pd.DataFrame, pd.DataFrame, Path]:
    expression = pd.read_csv(root / PREPROCESSED / "rnaseq/log2cpm.tsv", sep="\t", index_col=0)
    samples = pd.read_csv(root / PREPROCESSED / "rnaseq/samples.tsv", sep="\t").set_index("column")
    groups = samples["group"].reindex(expression.columns)
    regulon, regulon_path, resource_metadata = load_regulon(root / "data/resources/regulons", "collectri")
    primary = pd.read_csv(root / ACTIVITIES / "collectri/tf_activity_contrast.tsv", sep="\t", index_col=0)
    per_line = pd.read_csv(root / ACTIVITIES / "collectri/tf_activity_per_line.tsv", sep="\t", index_col=0)
    loo = leave_one_line_out_ulm(expression, groups, regulon, tmin=5)

    selections = []
    measured_rows = []
    grouped: dict[str, dict[str, dict[str, float]]] = {}
    for pkn_name in pkns:
        selected, selection = select_stable_top_tfs(
            primary["score"], loo, {node for edge in pkns[pkn_name] for node in (edge[0], edge[2])}, k=50,
        )
        selection.insert(0, "pkn", pkn_name)
        selections.append(selection)
        measured_rows.extend({"pkn": pkn_name, "TF": tf, "score": float(primary.loc[tf, "score"])} for tf in selected)
        contrast = {"PMS_vs_Ctrl": {tf: float(primary.loc[tf, "score"]) for tf in selected}}
        grouped[pkn_name] = contrast
        if pkn_name == "omnipath_lines":
            grouped[pkn_name] = {
                str(line): {tf: float(per_line.loc[tf, line]) for tf in selected}
                for line in per_line.columns
            }
    measured = pd.DataFrame(measured_rows)
    measured["resource"] = "collectri"
    measured["tmin"] = 5
    measured["statistic"] = "constant-prior moderated t; ULM score"
    selection_table = pd.concat(selections, ignore_index=True)
    return grouped, measured, selection_table, regulon_path


def _metabolite_measurements(table: pd.DataFrame) -> dict[str, float]:
    return {str(row.pkn_node): float(row.value) for row in table.itertuples(index=False)}


def _write_lambda_result(directory: Path, result: dict[str, Any]) -> None:
    """Write edges, node states (from CORNETO vertex_value) and fit for one lambda."""
    directory.mkdir(parents=True, exist_ok=True)
    result["selected_edges"].to_csv(directory / "edges.tsv", sep="\t", index=False)
    result["nodes"].to_csv(directory / "nodes.tsv", sep="\t", index=False)
    result["fit"].to_csv(directory / "fit.tsv", sep="\t", index=False)


def _draw_result(out: Path, run_name: str, edges: list[tuple[str, int, str]], result: dict[str, Any], measurements: dict[str, dict[str, float]]) -> None:
    selected = result["selected_edges"]
    if selected.empty:
        return
    selected_pkn = list(dict.fromkeys((row.source, int(row.pkn_sign), row.target) for row in selected.itertuples(index=False)))
    fig, ax = plt.subplots(figsize=(max(8, len(selected_pkn) * 0.4), 7))
    draw_signed_network(
        ax,
        selected_pkn,
        selected[["sample", "source", "target", "edge_value"]],
        inputs=result["inputs"],
        measured={node: int(np.sign(value)) for values in measurements.values() for node, value in values.items()},
        title=f"{run_name}: selected CORNETO subnetwork; lambda={result['lambda_reg']:g}",
    )
    ax.legend(handles=network_legend_handles(), loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=7)
    save_figure(fig, out / "figures" / f"{run_name}_network")
    plt.close(fig)


def _parse_backfill_log(path: Path) -> dict[tuple[str, float], dict[str, Any]]:
    """Read per-lambda solver lines ("<run>: lambda=<x>" then "status=... error=...")."""
    records: dict[tuple[str, float], dict[str, Any]] = {}
    current: tuple[str, float] | None = None
    header = re.compile(r"^(\w+): lambda=([0-9.eE+-]+)$")
    status = re.compile(r"^\s*status=(\S+) accepted=(\S+) gap=(\S+) seconds=(\S+) error=(\S+)$")
    for line in Path(path).read_text(errors="replace").splitlines():
        if match := header.match(line.strip()):
            current = (match.group(1), float(match.group(2)))
        elif current and (match := status.match(line)):
            records[current] = {
                "status": match.group(1),
                "accepted": match.group(2) == "True",
                "mip_gap": float(match.group(3)),
                "solve_seconds": float(match.group(4)),
                "weighted_fitting_error": float(match.group(5)),
            }
            current = None
    return records


def _check_resume_inputs(out: Path, measured_tfs: pd.DataFrame, metabolites: pd.DataFrame) -> None:
    """Refuse to resume if the measurements differ from the interrupted run's."""
    for name, new, keys in (
        ("measured_tfs.tsv", measured_tfs, ["pkn", "TF", "score"]),
        ("measured_metabolites.tsv", metabolites, ["pkn_node", "value"]),
    ):
        path = out / "inputs" / name
        if not path.exists():
            raise ValueError(f"Cannot resume: {path} is missing")
        old = pd.read_csv(path, sep="\t")[keys].sort_values(keys[:-1]).reset_index(drop=True)
        cur = new[keys].sort_values(keys[:-1]).reset_index(drop=True)
        same = old[keys[:-1]].astype(str).equals(cur[keys[:-1]].astype(str)) and np.allclose(
            old[keys[-1]].to_numpy(float), cur[keys[-1]].to_numpy(float), rtol=0, atol=1e-12)
        if not same:
            raise ValueError(f"Cannot resume: {name} differs from the interrupted run; archive the folder instead")


def _task_signature(task: dict[str, Any]) -> str:
    """Stable hash of everything that defines one (run, lambda) problem."""
    payload = json.dumps(task, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _solve_worker(task_path: Path) -> int:
    """Solve one (run, lambda) problem in this process and write its tables.

    Called as a separate Python process for every problem, so all solver
    memory is returned to the operating system when it exits (two stage-04
    runs on 2026-10-11 were killed when problems were solved in one process).
    """
    task = json.loads(task_path.read_text())
    edges = read_pkn(ROOT / task["pkn_path"])
    result = solve_carnival(
        edges,
        task["measurements"],
        float(task["lambda_reg"]),
        time_limit=int(task["time_limit"]),
        solver=SOLVER,
        solver_relative_gap=SOLVER_RELATIVE_GAP,
        threads=SOLVER_THREADS,
        seed=SOLVER_SEED,
    )
    directory = task_path.parent
    _write_lambda_result(directory, result)
    summary = {
        key: result[key] for key in (
            "lambda_reg", "status", "accepted", "mip_gap", "objectives", "weighted_fitting_error",
            "inputs", "missing_measurements", "n_missing_measurements",
            "preprocessed_node_count", "preprocessed_edge_count", "solve_seconds",
        )
    }
    summary["accepted"] = bool(summary["accepted"])
    summary["task_signature"] = task["signature"]
    (directory / "result.json").write_text(json.dumps(summary, default=float))
    print(f"  status={result['status']} accepted={result['accepted']} gap={result['mip_gap']:.4g} "
          f"seconds={result['solve_seconds']:.0f} error={result['weighted_fitting_error']:.3f}", flush=True)
    return 0


def _read_table(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t", dtype={"sample": str, "source": str, "target": str, "node": str})
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return pd.DataFrame(columns=columns or [])


def _load_lambda_result(directory: Path, lambda_reg: float) -> dict[str, Any]:
    """Rebuild the in-memory result dict from a worker's written tables."""
    summary = json.loads((directory / "result.json").read_text())
    edges = _read_table(directory / "edges.tsv", ["sample", "source", "target", "edge_value", "pkn_sign"])
    nodes = _read_table(directory / "nodes.tsv", ["sample", "node", "state", "is_measured", "observed_value", "is_input", "reachability"])
    fit = _read_table(directory / "fit.tsv") if summary["accepted"] else pd.DataFrame()
    states: dict[str, dict[str, int]] = {}
    for row in nodes.itertuples(index=False):
        if row.reachability == "reachable" or int(row.state) != 0:
            states.setdefault(str(row.sample), {})[str(row.node)] = int(row.state)
    summary.update({"lambda_reg": float(lambda_reg), "selected_edges": edges, "nodes": nodes, "fit": fit, "states": states})
    return summary


def _failed_result(lambda_reg: float, status: str) -> dict[str, Any]:
    return {
        "lambda_reg": float(lambda_reg), "status": status, "accepted": False, "mip_gap": float("nan"),
        "objectives": {}, "weighted_fitting_error": float("nan"), "inputs": [], "missing_measurements": {},
        "n_missing_measurements": 0, "preprocessed_node_count": 0, "preprocessed_edge_count": 0,
        "solve_seconds": float("nan"),
        "selected_edges": pd.DataFrame(columns=["sample", "source", "target", "edge_value", "pkn_sign"]),
        "nodes": pd.DataFrame(columns=["sample", "node", "state", "is_measured", "observed_value", "is_input", "reachability"]),
        "fit": pd.DataFrame(), "states": {},
    }


def _run_one(
    out: Path,
    run_name: str,
    pkn_path: Path,
    edges: list[tuple[str, int, str]],
    measurements: dict[str, dict[str, float]],
    time_limit: int,
    lambdas: tuple[float, ...] = LAMBDA_GRID,
    *,
    resume: bool = False,
    backfill: dict[tuple[str, float], dict[str, Any]] | None = None,
) -> tuple[list[dict[str, Any]], float | None]:
    """Solve every lambda of one run, each in a fresh Python process."""
    output = out / run_name
    results = []
    backfill_inputs: dict[str, tuple[list[str], dict[str, list[str]]]] = {}
    for lambda_reg in lambdas:
        directory = output / f"lambda_{lambda_label(lambda_reg)}"
        task = {
            "run": run_name, "pkn_path": _portable_path(pkn_path),
            "pkn_sha256": _file_sha256(pkn_path),
            "measurements": measurements, "lambda_reg": float(lambda_reg), "time_limit": int(time_limit),
            "solver": [SOLVER, SOLVER_RELATIVE_GAP, SOLVER_THREADS, SOLVER_SEED],
        }
        task["signature"] = _task_signature(task)
        done = directory / "result.json"
        tables_exist = all((directory / name).exists() for name in ("edges.tsv", "nodes.tsv", "fit.tsv"))
        logged = (backfill or {}).get((run_name, float(lambda_reg)))
        if resume and not done.exists() and tables_exist and logged is not None:
            # One-time backfill for results written before result.json existed
            # (interrupted run of 2026-10-11): solver fields come from the run log.
            if run_name not in backfill_inputs:
                _, used, missing = make_corneto_data(edges, measurements)
                backfill_inputs[run_name] = (used, missing)
            used, missing = backfill_inputs[run_name]
            done.write_text(json.dumps({
                "lambda_reg": float(lambda_reg), **logged, "objectives": {},
                # full-precision error from fit.tsv (the log prints 3 decimals)
                "weighted_fitting_error": float(pd.read_csv(directory / "fit.tsv", sep="\t")["weighted_fitting_error"].sum()),
                "inputs": used, "missing_measurements": missing,
                "n_missing_measurements": sum(len(nodes) for nodes in missing.values()),
                "preprocessed_node_count": float("nan"), "preprocessed_edge_count": float("nan"),
                "backfilled_from_log": True, "task_signature": task["signature"],
            }))
            print(f"{run_name}: lambda={lambda_reg:g} (backfilled from log)", flush=True)
            results.append(_load_lambda_result(directory, lambda_reg))
            continue
        if resume and done.exists() and json.loads(done.read_text()).get("task_signature") == task["signature"]:
            print(f"{run_name}: lambda={lambda_reg:g} (resumed from {done.relative_to(out)})", flush=True)
            results.append(_load_lambda_result(directory, lambda_reg))
            continue
        directory.mkdir(parents=True, exist_ok=True)
        done.unlink(missing_ok=True)
        task_path = directory / "task.json"
        task_path.write_text(json.dumps(task))
        print(f"{run_name}: lambda={lambda_reg:g}", flush=True)
        completed = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", str(task_path)], check=False)
        if completed.returncode != 0 or not done.exists():
            status = f"worker_failed_exit_{completed.returncode}"
            print(f"  {status}", flush=True)
            results.append(_failed_result(lambda_reg, status))
            continue
        results.append(_load_lambda_result(directory, lambda_reg))
    # The lambda rule needs the lambda = 0 baseline; a pilot without it selects nothing.
    chosen = choose_lambda(results)
    if chosen is not None:
        chosen_result = next(result for result in results if result["lambda_reg"] == chosen)
        _draw_result(out, run_name, edges, chosen_result, measurements)
    return results, chosen


def _portable_path(path: Path) -> str:
    """Repository-relative path (identical under Windows and WSL mounts), else absolute."""
    try:
        return Path(path).resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return Path(path).resolve().as_posix()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _overlap_rows(selected: dict[str, dict[str, Any]]) -> pd.DataFrame:
    rows = []
    for left, right in (("omnipath", "cosmos_tf"), ("omnipath", "cosmos_multi"), ("cosmos_tf", "cosmos_multi")):
        if left not in selected or right not in selected:
            continue
        left_edges = selected[left]["selected_edges"]
        right_edges = selected[right]["selected_edges"]
        left_nodes = {node for column in ("source", "target") for node in left_edges[column] if not str(node).startswith("Metab__")}
        right_nodes = {node for column in ("source", "target") for node in right_edges[column] if not str(node).startswith("Metab__")}
        left_pairs = {(row.source, row.target) for row in left_edges.itertuples(index=False)}
        right_pairs = {(row.source, row.target) for row in right_edges.itertuples(index=False)}
        rows.extend([
            {"comparison": f"{left}_vs_{right}", "entity": "protein_node", "jaccard": len(left_nodes & right_nodes) / len(left_nodes | right_nodes) if left_nodes | right_nodes else np.nan},
            {"comparison": f"{left}_vs_{right}", "entity": "edge", "jaccard": len(left_pairs & right_pairs) / len(left_pairs | right_pairs) if left_pairs | right_pairs else np.nan},
        ])
    return pd.DataFrame(rows)


def _reachable(edges: pd.DataFrame, sources: set[str], target: str) -> bool:
    adjacency: dict[str, list[str]] = {}
    for row in edges.itertuples(index=False):
        adjacency.setdefault(str(row.source), []).append(str(row.target))
    queue = deque(sources)
    seen = set(sources)
    while queue:
        node = queue.popleft()
        if node == target:
            return True
        for child in adjacency.get(node, []):
            if child not in seen:
                seen.add(child)
                queue.append(child)
    return False


def _reaction_gene_map(nodes_path: Path) -> dict[str, set[str]]:
    """Map COSMOS reaction nodes to their recognised HGNC symbols (stage-03 node table)."""
    nodes = pd.read_csv(nodes_path, sep="\t", dtype=str)
    reactions = nodes[nodes["type"].eq("enzyme_or_reaction")]
    return {
        str(row.node): {symbol for symbol in str(row.approved_symbol).split(";") if symbol and symbol != "nan"}
        for row in reactions.itertuples(index=False)
    }


def _aim5_rows(selected: dict[str, dict[str, Any]], reaction_genes: dict[str, dict[str, set[str]]]) -> pd.DataFrame:
    """Aim-5 readout per COSMOS run; active means CORNETO state != 0."""
    rows = []
    for run_name in ("cosmos_multi", "cosmos_repair"):
        result = selected.get(run_name)
        if result is None:
            continue
        edges = result["selected_edges"]
        nodes = result["nodes"]
        active = nodes[nodes["state"].ne(0)]
        active_nodes = set(active["node"].astype(str))
        genes = reaction_genes.get(run_name, {})
        canonical_active = sorted(
            node for node in active_nodes if genes.get(node, set()) & CANONICAL_CHOLESTEROL_GENES
        )
        cholesterol = sorted(node for node in active_nodes if CHOLESTEROL_HMDB in node)
        hmgcr_state = int(active.loc[active["node"].eq("HMGCR"), "state"].iloc[0]) if "HMGCR" in active_nodes else 0
        for target in ("E2F1", "EGR1", "SP1", "JUN", "JUNB"):
            rows.append({
                "run": run_name,
                "hmgcr_state": hmgcr_state,
                "n_canonical_reaction_nodes_active": len(canonical_active),
                "canonical_reaction_nodes_active": ",".join(canonical_active),
                "cholesterol_nodes_active": ",".join(
                    f"{node}:{int(active.loc[active['node'].eq(node), 'state'].iloc[0])}" for node in cholesterol
                ),
                "cholesterol_to_target": target,
                "target_state": int(active.loc[active["node"].eq(target), "state"].iloc[0]) if target in active_nodes else 0,
                "selected_route_exists": _reachable(edges, set(cholesterol), target) if cholesterol else False,
            })
    return pd.DataFrame(rows)


RUN_NAMES = ("omnipath", "cosmos_tf", "cosmos_multi", "cosmos_repair", "omnipath_lines")


def _parse_args(argv: list[str] | None) -> Any:
    """Optional run/lambda subsets for a timing pilot; registered values only."""
    import argparse

    parser = argparse.ArgumentParser(description="Stage 04 CORNETO inference")
    parser.add_argument("--runs", nargs="+", choices=RUN_NAMES, default=list(RUN_NAMES))
    parser.add_argument("--lambdas", nargs="+", type=float, default=list(LAMBDA_GRID))
    parser.add_argument("--resume", action="store_true",
                        help="reuse finished (run, lambda) results in an incomplete output folder")
    parser.add_argument("--backfill-log", type=Path, default=None,
                        help="with --resume: stage-04 log of an interrupted run whose finished lambdas lack result.json")
    args = parser.parse_args(argv)
    unregistered = [value for value in args.lambdas if value not in LAMBDA_GRID]
    if unregistered:
        parser.error(f"lambda values must come from the registered grid {LAMBDA_GRID}: {unregistered}")
    args.lambdas = tuple(sorted(set(args.lambdas)))
    args.runs = [name for name in RUN_NAMES if name in set(args.runs)]
    return args


def _inputs_rows(run_name: str, edges: list[tuple[str, int, str]], used_inputs: list[str], measured: set[str]) -> list[dict[str, Any]]:
    """Describe every PKN source node and whether it was passed as an input."""
    nodes = {node for edge in edges for node in (edge[0], edge[2])}
    targets = {edge[2] for edge in edges}
    used = set(used_inputs)
    return [
        {
            "run": run_name,
            "node": node,
            "excluded": is_excluded_input(node),
            "measured": node in measured,
            "used_as_input": node in used,
        }
        for node in sorted(nodes - targets)
    ]


def main(root: Path = ROOT, argv: list[str] | None = None) -> int:
    if argv and argv[0] == "--worker":
        return _solve_worker(Path(argv[1]))
    args = _parse_args(argv)
    pilot = args.runs != list(RUN_NAMES) or args.lambdas != tuple(LAMBDA_GRID)
    smoke_corneto()
    apply_plot_style()
    pkn_paths = _pkn_paths(root)
    statuses: list[dict[str, Any]] = []
    runs = list(args.runs)
    # The repaired COSMOS PKN comes from the Aim-5 route plan; until that stage
    # has been run, cosmos_repair is recorded as not run (audit finding 1).
    if "cosmos_repair" in runs and not pkn_paths["cosmos_repair"].exists():
        runs.remove("cosmos_repair")
        statuses.append({"run": "cosmos_repair", "status": "not_run_missing_pkn",
                         "selected_lambda": np.nan, "n_missing_measurements": np.nan})
        print("cosmos_repair: repaired PKN missing; run recorded as not_run_missing_pkn")
    missing = [pkn_paths[name] for name in runs if not pkn_paths[name].exists()]
    if missing:
        raise FileNotFoundError("Required stage-03 PKN is missing: " + ", ".join(map(str, missing)))
    if not runs:
        raise ValueError("No runnable stage-04 run selected")
    pkns = {name: read_pkn(pkn_paths[name]) for name in runs}
    tf_measurements, measured_tfs, selection, regulon_path = _tf_inputs(root, pkns)
    metabolites = _metabolite_table(root)
    metabolites_sensitivity = _metabolite_table(root, sensitivity=True)

    out = root / (ANALYSIS + "_pilot" if pilot else ANALYSIS)
    if out.exists() and any(out.iterdir()) and not (out / "04_inference.provenance.json").exists() and not args.resume:
        raise ValueError(f"Incomplete stage directory {out}; archive it or rerun with --resume")
    if args.backfill_log is not None and not args.resume:
        raise ValueError("--backfill-log requires --resume")
    backfill = _parse_backfill_log(args.backfill_log) if args.backfill_log is not None else {}
    if args.resume and (out / "inputs").exists():
        _check_resume_inputs(out, measured_tfs, metabolites)
    (out / "inputs").mkdir(parents=True, exist_ok=True)
    measured_tfs.to_csv(out / "inputs/measured_tfs.tsv", sep="\t", index=False)
    selection.to_csv(out / "inputs/tf_selection_loo.tsv", sep="\t", index=False)
    metabolites.to_csv(out / "inputs/measured_metabolites.tsv", sep="\t", index=False)
    metabolites_sensitivity.to_csv(out / "inputs/measured_metabolites_sensitivity.tsv", sep="\t", index=False)

    metabolite_values = _metabolite_measurements(metabolites)
    run_specs: dict[str, tuple[list[tuple[str, int, str]], dict[str, dict[str, float]], int]] = {}
    for run_name in runs:
        measurements = tf_measurements[run_name]
        if run_name in ("cosmos_multi", "cosmos_repair"):
            measurements = {"PMS_vs_Ctrl": {**measurements["PMS_vs_Ctrl"], **metabolite_values}}
        limit = COSMOS_TIME_LIMIT if run_name.startswith("cosmos") else OMNIPATH_TIME_LIMIT
        run_specs[run_name] = (pkns[run_name], measurements, limit)

    reaction_genes: dict[str, dict[str, set[str]]] = {}
    node_tables = {
        "cosmos_multi": root / PKN / "cosmos/primary/nodes.tsv",
        "cosmos_repair": root / PKN / "cosmos/variants/canonical_repair/nodes.tsv",
    }
    for run_name, nodes_path in node_tables.items():
        if run_name in runs:
            reaction_genes[run_name] = _reaction_gene_map(nodes_path)

    all_grid = []
    input_rows: list[dict[str, Any]] = []
    chosen: dict[str, dict[str, Any]] = {}
    for run_name, (edges, measurements, time_limit) in run_specs.items():
        results, selected_lambda = _run_one(
            out, run_name, pkn_paths[run_name], edges, measurements, time_limit, args.lambdas,
            resume=args.resume, backfill=backfill,
        )
        measured_nodes = {node for values in measurements.values() for node in values}
        used_inputs = next((result["inputs"] for result in results if result["inputs"]), None)
        if used_inputs is None:
            used_inputs = make_corneto_data(edges, measurements)[1]
        input_rows.extend(_inputs_rows(run_name, edges, used_inputs, measured_nodes))
        all_grid.append(lambda_grid_table(results).assign(run=run_name))
        if selected_lambda is None:
            statuses.append({"run": run_name, "status": "no_accepted_lambda", "selected_lambda": np.nan,
                             "n_missing_measurements": int(results[0]["n_missing_measurements"])})
            continue
        result = next(result for result in results if result["lambda_reg"] == selected_lambda)
        chosen[run_name] = result
        statuses.append({"run": run_name, "status": "accepted", "selected_lambda": selected_lambda,
                         "n_missing_measurements": result["n_missing_measurements"]})
    pd.DataFrame(input_rows).to_csv(out / "inputs/inputs.tsv", sep="\t", index=False)
    (out / "summary").mkdir(parents=True, exist_ok=True)
    pd.concat(all_grid, ignore_index=True).to_csv(out / "summary/lambda_grid.tsv", sep="\t", index=False)
    pd.DataFrame(statuses).to_csv(out / "summary/run_status.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"run": run_name, "selected_lambda": result["lambda_reg"], "weighted_fitting_error": result["weighted_fitting_error"]}
        for run_name, result in chosen.items()
    ]).to_csv(out / "summary/selected_lambda.tsv", sep="\t", index=False)
    _overlap_rows(chosen).to_csv(out / "summary/aim4_overlap.tsv", sep="\t", index=False)
    _aim5_rows(chosen, reaction_genes).to_csv(out / "summary/aim5_route.tsv", sep="\t", index=False)

    metabolomics = root / PREPROCESSED / "metabolomics"
    stage_inputs = [
        root / PLAN,
        root / PREPROCESSED / "rnaseq/log2cpm.tsv",
        root / PREPROCESSED / "rnaseq/samples.tsv",
        root / ACTIVITIES / "collectri/tf_activity_contrast.tsv",
        root / ACTIVITIES / "collectri/tf_activity_per_line.tsv",
        root / PKN / "mapping/metabolite_mapping.tsv",
        metabolomics / "ST003331/features.tsv",
        metabolomics / "ST003331/line_log2.tsv",
        metabolomics / "ST003332/features.tsv",
        metabolomics / "ST003332/line_log2_vs_blank.tsv",
        metabolomics / "ST003328/features.tsv",
        metabolomics / "ST003328/line_log2.tsv",
        root / PKN / "cosmos/primary/nodes.tsv",
        Path(regulon_path),
    ] + [pkn_paths[name] for name in runs] + [
        node_tables[name] for name in ("cosmos_repair",) if name in runs
    ]
    stage_inputs = list(dict.fromkeys(stage_inputs))
    resources = {
        "collectri": str(Path(regulon_path).relative_to(root)),
        "pkns": {name: str(pkn_paths[name].relative_to(root)) for name in runs},
    }
    parameters = {
        "runs": runs,
        "pilot": pilot,
        "lambda_grid": list(args.lambdas),
        "registered_lambda_grid": list(LAMBDA_GRID),
        "lambda_selection": "largest accepted lambda with weighted error <= 1.10 * lambda=0",
        "solver": SOLVER,
        "solver_relative_gap": SOLVER_RELATIVE_GAP,
        "solver_threads": SOLVER_THREADS,
        "solver_seed": SOLVER_SEED,
        "time_limits_seconds": {"omnipath": OMNIPATH_TIME_LIMIT, "cosmos": COSMOS_TIME_LIMIT},
        "acceptance_gap": 0.05,
        "excluded_inputs": sorted(EXCLUDED_INPUTS),
        "measured_tfs": "CollecTRI tmin=5, top 50 | primary moderated t, six-LOO sign stability",
        "metabolites": "ST003331/ST003332 m0_unconfounded plus untreated ST003328 free cholesterol; "
                       "moderated-t prior from all complete features of each study",
        "node_states": "CORNETO vertex_value",
        "canonical_cholesterol_genes": sorted(CANONICAL_CHOLESTEROL_GENES),
        "e1": "constant-prior moderated t; limma-trend deferred",
        "execution": "one Python process per (run, lambda); resume=" + str(args.resume),
        "backfilled_from_log": sorted(f"{run}:{lam:g}" for run, lam in backfill),
    }
    code = [
        root / "src/inference.py", root / "src/activities.py", root / "src/network_plot.py",
        root / "src/plot_config.py", root / "src/preprocess.py", root / "scripts/analysis/05_corneto_inference.py",
        root / PLAN,
    ]
    provenance = stage_provenance(root, "04_inference", stage_inputs, resources, parameters, code)
    provenance["accessions"] = {
        "GSE297192": "Park et al. 2025, Neuron; DOI 10.1016/j.neuron.2025.09.022",
        "ST003331": "Metabolomics Workbench; Ionescu et al. 2024",
        "ST003332": "Metabolomics Workbench; Ionescu et al. 2024",
        "ST003328": "Metabolomics Workbench; Ionescu et al. 2024",
    }
    write_stage_provenance(out, provenance)
    print(f"Stage 04 complete: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(argv=sys.argv[1:]))
