"""Stage 05: stability, analysis-choice and topology controls for stage 04."""

from __future__ import annotations

import itertools
import random
import sys
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from activities import (  # noqa: E402
    leave_one_line_out_ulm,
    load_regulon,
    moderated_t_statistics,
    run_ulm,
    select_stable_top_tfs,
)
from inference import (  # noqa: E402
    LAMBDA_GRID,
    OMNIPATH_TIME_LIMIT,
    COSMOS_TIME_LIMIT,
    SOLVER,
    SOLVER_RELATIVE_GAP,
    SOLVER_SEED,
    SOLVER_THREADS,
    lambda_label,
    read_pkn,
    sample_edges,
    solve_carnival,
)
from pkn import rewire_signed_graph  # noqa: E402
from plot_config import save_figure  # noqa: E402
from preprocess import stage_provenance, write_stage_provenance  # noqa: E402

ANALYSIS = "results/ionescu_corneto/05_robustness"
STAGE04 = "results/ionescu_corneto/04_inference"
ACTIVITIES = "results/ionescu_corneto/02_activities"
PREPROCESSED = "results/ionescu_corneto/01_preprocessing"
PKN = "results/ionescu_corneto/03_pkn"
PLAN = "doc/decisions/2026-10-10-stage04-inference-plan.md"
SEED_BASE = 20261010


def _edges(root: Path, variant: str = "primary") -> list[tuple[str, int, str]]:
    if variant == "primary":
        path = root / PKN / "omnipath/primary/pkn.tsv"
    else:
        path = root / PKN / f"omnipath/variants/{variant}/pkn.tsv"
    return read_pkn(path)


def _selected_lambda(root: Path, run: str) -> float:
    table = pd.read_csv(root / STAGE04 / "summary/selected_lambda.tsv", sep="\t")
    row = table[table["run"].eq(run)]
    if len(row) != 1:
        raise ValueError(f"Stage 04 has no selected lambda for {run}")
    return float(row.iloc[0]["selected_lambda"])


def _expression_inputs(root: Path, regulon_name: str = "collectri") -> tuple[pd.DataFrame, pd.Series, pd.DataFrame, pd.DataFrame]:
    expression = pd.read_csv(root / PREPROCESSED / "rnaseq/log2cpm.tsv", sep="\t", index_col=0)
    samples = pd.read_csv(root / PREPROCESSED / "rnaseq/samples.tsv", sep="\t").set_index("column")
    groups = samples["group"].reindex(expression.columns)
    regulon, _, _ = load_regulon(root / "data/resources/regulons", regulon_name)
    loo = leave_one_line_out_ulm(expression, groups, regulon, tmin=5)
    return expression, groups, regulon, loo


def _primary_score(root: Path, statistic: str, regulon_name: str) -> pd.Series:
    file_name = {
        "moderated_t": "tf_activity_contrast.tsv",
        "welch": "tf_activity_contrast_welch.tsv",
        "logFC": "tf_activity_contrast_logfc.tsv",
    }[statistic]
    table = pd.read_csv(root / ACTIVITIES / regulon_name / file_name, sep="\t", index_col=0)
    return table["score"]


def _measurements_for_selection(
    root: Path,
    edges: list[tuple[str, int, str]],
    *,
    statistic: str = "moderated_t",
    regulon_name: str = "collectri",
    k: int = 50,
) -> tuple[dict[str, dict[str, float]], list[str]]:
    expression, groups, regulon, loo = _expression_inputs(root, regulon_name)
    primary = _primary_score(root, statistic, regulon_name)
    selected, _ = select_stable_top_tfs(
        primary,
        loo,
        {node for source, _, target in edges for node in (source, target)},
        k=k,
    )
    return {"PMS_vs_Ctrl": {tf: float(primary[tf]) for tf in selected}}, selected


def _edge_set(edges: pd.DataFrame) -> set[tuple[str, str]]:
    return {(str(row.source), str(row.target)) for row in edges.itertuples(index=False)}


def _jaccard(left: set[Any], right: set[Any]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else np.nan


def _save_edges(directory: Path, edges: pd.DataFrame) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    edges.to_csv(directory / "edges.tsv", sep="\t", index=False)


def _node_type(node: str) -> str:
    if node.startswith("Metab__"):
        return "metabolite"
    if node.startswith("Gene"):
        return "reaction"
    return "protein"


def _rewire_type_preserving(
    edges: list[tuple[str, int, str]], *, seed: int, accepted_swaps: int | None = None
) -> list[tuple[str, int, str]]:
    """Swap targets only within source-type/target-type/sign classes."""
    current = list(edges)
    target_swaps = accepted_swaps if accepted_swaps is not None else 3 * len(current)
    rng = random.Random(seed)
    edge_keys = {(source, target) for source, _, target in current}
    groups: dict[tuple[str, str, int], list[int]] = {}
    for index, (source, sign, target) in enumerate(current):
        groups.setdefault((_node_type(source), _node_type(target), sign), []).append(index)
    accepted = 0
    attempts = 0
    max_attempts = max(1000, target_swaps * 200)
    while accepted < target_swaps and attempts < max_attempts:
        attempts += 1
        candidates = [indices for indices in groups.values() if len(indices) >= 2]
        if not candidates:
            break
        first, second = rng.sample(rng.choice(candidates), 2)
        source_a, sign_a, target_b = current[first]
        source_c, sign_c, target_d = current[second]
        proposed = [(source_a, sign_a, target_d), (source_c, sign_c, target_b)]
        proposed_keys = [(source, target) for source, _, target in proposed]
        old_keys = {(source_a, target_b), (source_c, target_d)}
        if (source_a == target_d or source_c == target_b or
                proposed_keys[0] == proposed_keys[1] or
                any(key in edge_keys and key not in old_keys for key in proposed_keys)):
            continue
        current[first], current[second] = proposed
        edge_keys.difference_update(old_keys)
        edge_keys.update(proposed_keys)
        accepted += 1
    if accepted < target_swaps:
        raise RuntimeError(f"Could not complete {target_swaps} type-preserving swaps after {attempts} attempts")
    return sorted(current)


def _sample_stage04_runs(root: Path, out: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sample alternatives by rebuilding each accepted stage-04 problem."""
    rows = []
    frequencies = []
    run_settings = {
        "omnipath": (_edges(root), OMNIPATH_TIME_LIMIT),
        "cosmos_tf": (read_pkn(root / PKN / "cosmos/primary/pkn.tsv"), COSMOS_TIME_LIMIT),
        "cosmos_multi": (read_pkn(root / PKN / "cosmos/primary/pkn.tsv"), COSMOS_TIME_LIMIT),
        "cosmos_repair": (read_pkn(root / PKN / "cosmos/variants/canonical_repair/pkn.tsv"), COSMOS_TIME_LIMIT),
        "omnipath_lines": (_edges(root), OMNIPATH_TIME_LIMIT),
    }
    measured = pd.read_csv(root / STAGE04 / "inputs/measured_tfs.tsv", sep="\t")
    for run, (edges, time_limit) in run_settings.items():
        selected_table = pd.read_csv(root / STAGE04 / "summary/selected_lambda.tsv", sep="\t")
        selected_row = selected_table[selected_table["run"].eq(run)]
        if len(selected_row) != 1:
            rows.append({"run": run, "status": "stage04_unsolved"})
            continue
        lambda_reg = float(selected_row.iloc[0]["selected_lambda"])
        tf = measured[measured["pkn"].eq(run)]
        measurement = {"PMS_vs_Ctrl": dict(zip(tf["TF"], tf["score"]))}
        if run == "omnipath_lines":
            per_line = pd.read_csv(root / ACTIVITIES / "collectri/tf_activity_per_line.tsv", sep="\t", index_col=0)
            measurement = {
                line: {tf_name: float(per_line.loc[tf_name, line]) for tf_name in tf["TF"]}
                for line in per_line.columns
            }
        result = solve_carnival(edges, measurement, lambda_reg, time_limit=time_limit)
        if not result["accepted"]:
            rows.append({"run": run, "status": "unsolved"})
            continue
        alternatives = sample_edges(result["problem"], result["method"], time_limit=time_limit)
        n_alternatives = int(alternatives["alternative"].nunique()) if not alternatives.empty else 0
        rows.append({"run": run, "status": "accepted", "n_alternatives": n_alternatives, "lambda_reg": lambda_reg})
        if not alternatives.empty:
            run_frequency = (
                alternatives.groupby(["source", "target"], as_index=False)
                .agg(n_alternatives=("alternative", "nunique"))
                .assign(run=run, n_total_alternatives=n_alternatives,
                        frequency=lambda x: x["n_alternatives"] / max(n_alternatives, 1))
            )
            (out / "sampling" / run).mkdir(parents=True, exist_ok=True)
            run_frequency.to_csv(out / "sampling" / run / "edge_frequency.tsv", sep="\t", index=False)
            frequencies.append(run_frequency)
    frequency_table = pd.concat(frequencies, ignore_index=True) if frequencies else pd.DataFrame()
    return pd.DataFrame(rows), frequency_table


def _loo_controls(
    root: Path,
    out: Path,
    lambda_reg: float,
    *,
    pkn_edges: list[tuple[str, int, str]] | None = None,
    directory_prefix: str = "",
    time_limit: int = OMNIPATH_TIME_LIMIT,
    extra_measurements: dict[str, float] | None = None,
) -> tuple[pd.DataFrame, list[pd.DataFrame]]:
    expression, groups, regulon, full_loo = _expression_inputs(root)
    edges = pkn_edges if pkn_edges is not None else _edges(root)
    pkn_nodes = {node for source, _, target in edges for node in (source, target)}
    rows = []
    edge_tables = []
    for left_out in expression.columns:
        primary = full_loo[f"without_{left_out}"]
        selected, selection = select_stable_top_tfs(primary, full_loo, pkn_nodes, k=50)
        measurements = {"without_" + left_out: {tf: float(primary[tf]) for tf in selected}}
        if extra_measurements:
            measurements["without_" + left_out].update(extra_measurements)
        result = solve_carnival(edges, measurements, lambda_reg, time_limit=time_limit)
        directory = out / "loo" / directory_prefix / left_out if directory_prefix else out / "loo" / left_out
        directory.mkdir(parents=True, exist_ok=True)
        if result["accepted"]:
            _save_edges(directory, result["selected_edges"])
            edge_tables.append(result["selected_edges"].assign(
                loo=left_out, pkn=directory_prefix or "omnipath"
            ))
        selection.to_csv(directory / "tf_selection.tsv", sep="\t", index=False)
        rows.append({
            "left_out": left_out,
            "status": result["status"],
            "accepted": result["accepted"],
            "n_selected_tfs": len(selected),
            "n_edges": len(result["selected_edges"]),
            "weighted_fitting_error": result["weighted_fitting_error"],
        })
    return pd.DataFrame(rows), edge_tables


def _choice_grid(root: Path, out: Path, primary_edges: pd.DataFrame, primary_measurements: dict[str, dict[str, float]], lambda_reg: float) -> pd.DataFrame:
    rows = []
    factors = [
        ("statistic", "welch"), ("statistic", "logFC"),
        ("regulon", "dorothea_ABC"), ("regulon", "collectri_nocomplex"),
        ("k", 25), ("k", 100),
        ("pkn", "nc_rule"), ("pkn", "signor_only"), ("pkn", "unpruned"),
    ]
    lambda_index = list(LAMBDA_GRID).index(lambda_reg)
    factors.extend(
        ("lambda", value)
        for value in LAMBDA_GRID[max(0, lambda_index - 1):min(len(LAMBDA_GRID), lambda_index + 2)]
    )
    for factor, level in factors:
        edges = _edges(root, level) if factor == "pkn" else _edges(root)
        kwargs = {
            "statistic": level if factor == "statistic" else "moderated_t",
            "regulon_name": level if factor == "regulon" else "collectri",
            "k": level if factor == "k" else 50,
        }
        measurements, selected = _measurements_for_selection(root, edges, **kwargs)
        solve_lambda = float(level) if factor == "lambda" else lambda_reg
        result = solve_carnival(edges, measurements, solve_lambda, time_limit=OMNIPATH_TIME_LIMIT)
        edge_set = _edge_set(result["selected_edges"])
        baseline_edges = _edge_set(primary_edges)
        rows.append({
            "factor": factor,
            "level": str(level),
            "status": result["status"],
            "accepted": result["accepted"],
            "n_selected_tfs": len(selected),
            "edge_jaccard": _jaccard(edge_set, baseline_edges),
            "node_jaccard": _jaccard(
                {node for edge in edge_set for node in edge},
                {node for edge in baseline_edges for node in edge},
            ),
            "weighted_fitting_error": result["weighted_fitting_error"],
        })
        _save_edges(out / "grid" / factor / str(level), result["selected_edges"])
    return pd.DataFrame(rows)


def _permutations(root: Path, out: Path, lambda_reg: float) -> pd.DataFrame:
    expression = pd.read_csv(root / PREPROCESSED / "rnaseq/log2cpm.tsv", sep="\t", index_col=0)
    edges = _edges(root)
    regulon, _, _ = load_regulon(root / "data/resources/regulons", "collectri")
    rows = []
    true_edges: set[tuple[str, str]] | None = None
    for assignment_index, ctrl_tuple in enumerate(itertools.combinations(expression.columns, 3)):
        ctrl = set(ctrl_tuple)
        groups = pd.Series(["Ctrl" if line in ctrl else "PMS" for line in expression.columns], index=expression.columns)
        moderated, _ = moderated_t_statistics(expression, groups)
        vector = pd.DataFrame([moderated["moderated_t"]], index=["PMS_vs_Ctrl"])
        scores, _, _, _ = run_ulm(vector, regulon, tmin=5)
        loo = leave_one_line_out_ulm(expression, groups, regulon, tmin=5)
        selected, _ = select_stable_top_tfs(
            scores.loc["PMS_vs_Ctrl"], loo, {node for source, _, target in edges for node in (source, target)}, k=50,
        )
        measurements = {"PMS_vs_Ctrl": {tf: float(scores.loc["PMS_vs_Ctrl", tf]) for tf in selected}}
        result = solve_carnival(edges, measurements, lambda_reg, time_limit=OMNIPATH_TIME_LIMIT)
        edge_set = _edge_set(result["selected_edges"])
        if ctrl == {"C1", "C2", "C3"}:
            true_edges = edge_set
        rows.append({
            "assignment": assignment_index,
            "ctrl_lines": ",".join(sorted(ctrl)),
            "pms_lines": ",".join(sorted(set(expression.columns) - ctrl)),
            "status": result["status"],
            "accepted": result["accepted"],
            "weighted_fitting_error": result["weighted_fitting_error"],
            "n_edges": len(edge_set),
            "n_selected_tfs": len(selected),
            "jaccard_with_true": np.nan if true_edges is None else _jaccard(edge_set, true_edges),
        })
        _save_edges(out / "permutation" / str(assignment_index), result["selected_edges"])
    table = pd.DataFrame(rows)
    if true_edges is not None:
        table.loc[table["ctrl_lines"].eq("C1,C2,C3"), "jaccard_with_true"] = 1.0
    table["rank_by_error"] = table["weighted_fitting_error"].rank(method="min")
    return table


def _rewired(
    root: Path,
    out: Path,
    measurements: dict[str, dict[str, float]],
    lambda_reg: float,
    cosmos_measurements: dict[str, dict[str, float]] | None = None,
    cosmos_lambda: float | None = None,
) -> pd.DataFrame:
    edges = _edges(root)
    real_edge_set = {(source, target) for source, _, target in edges}
    rows = []
    for replicate in range(10):
        rewired = rewire_signed_graph(edges, seed=SEED_BASE + replicate)
        result = solve_carnival(rewired, measurements, lambda_reg, time_limit=OMNIPATH_TIME_LIMIT)
        selected = _edge_set(result["selected_edges"])
        rows.append({
            "pkn": "omnipath",
            "replicate": replicate,
            "status": result["status"],
            "accepted": result["accepted"],
            "weighted_fitting_error": result["weighted_fitting_error"],
            "n_edges": len(selected),
            "jaccard_with_real": _jaccard(selected, real_edge_set),
        })
        _save_edges(out / "rewired" / "omnipath" / str(replicate), result["selected_edges"])
    if cosmos_measurements is not None and cosmos_lambda is not None:
        cosmos_edges = read_pkn(root / PKN / "cosmos/primary/pkn.tsv")
        for replicate in range(5):
            rewired = _rewire_type_preserving(cosmos_edges, seed=SEED_BASE + replicate)
            result = solve_carnival(
                rewired, cosmos_measurements, cosmos_lambda, time_limit=COSMOS_TIME_LIMIT
            )
            selected = _edge_set(result["selected_edges"])
            rows.append({
                "pkn": "cosmos",
                "replicate": replicate,
                "status": result["status"],
                "accepted": result["accepted"],
                "weighted_fitting_error": result["weighted_fitting_error"],
                "n_edges": len(selected),
                "jaccard_with_real": np.nan,
            })
            _save_edges(out / "rewired" / "cosmos" / str(replicate), result["selected_edges"])
    return pd.DataFrame(rows)


def _hub_bias(root: Path, frequencies: pd.DataFrame, out: Path, rewired_nodes: set[str] | None = None) -> pd.DataFrame:
    frequencies = frequencies[frequencies["run"].eq("omnipath")] if not frequencies.empty else frequencies
    if frequencies.empty and not rewired_nodes:
        return pd.DataFrame(columns=["network", "odds_ratio_per_degree_doubling", "pvalue"])
    edges = _edges(root)
    degree: dict[str, int] = {}
    for source, _, target in edges:
        degree[source] = degree.get(source, 0) + 1
        degree[target] = degree.get(target, 0) + 1
    selected_nodes = set(frequencies.loc[frequencies["frequency"].gt(0), "source"]) | set(
        frequencies.loc[frequencies["frequency"].gt(0), "target"]
    ) if not frequencies.empty else set()
    nodes = sorted(degree)
    model_data = pd.DataFrame({
        "selected": [int(node in selected_nodes) for node in nodes],
        "log2_degree": [np.log2(degree[node]) for node in nodes],
    })
    rows = []
    for network, network_nodes in (("omnipath", selected_nodes), ("rewired_omnipath", rewired_nodes or set())):
        if not network_nodes:
            continue
        model_data["selected"] = [int(node in network_nodes) for node in nodes]
        fit = sm.Logit(model_data["selected"], sm.add_constant(model_data["log2_degree"])).fit(disp=False)
        rows.append({
            "network": network,
            "odds_ratio_per_degree_doubling": float(np.exp(fit.params["log2_degree"])),
            "pvalue": float(fit.pvalues["log2_degree"]),
            "n_nodes": len(nodes),
        })
    table = pd.DataFrame(rows)
    table.to_csv(out / "summary/hub_bias.tsv", sep="\t", index=False)
    return table


def main(root: Path = ROOT) -> int:
    stage04 = root / STAGE04 / "04_inference.provenance.json"
    if not stage04.exists():
        raise FileNotFoundError("Run stage 04 before stage 05")
    out = root / ANALYSIS
    (out / "summary").mkdir(parents=True, exist_ok=True)
    status, frequencies = _sample_stage04_runs(root, out)
    status.to_csv(out / "summary/sampling_status.tsv", sep="\t", index=False)
    if not frequencies.empty:
        frequencies.to_csv(out / "summary/edge_frequencies.tsv", sep="\t", index=False)
    lambda_reg = _selected_lambda(root, "omnipath")
    loo_status, loo_edges = _loo_controls(root, out, lambda_reg)
    loo_status.to_csv(out / "summary/loo_status.tsv", sep="\t", index=False)

    primary_edges = pd.read_csv(root / STAGE04 / f"omnipath/lambda_{lambda_label(lambda_reg)}/edges.tsv", sep="\t")
    measured = pd.read_csv(root / STAGE04 / "inputs/measured_tfs.tsv", sep="\t")
    omni = measured[measured["pkn"].eq("omnipath")]
    measurements = {"PMS_vs_Ctrl": dict(zip(omni["TF"], omni["score"]))}
    grid = _choice_grid(root, out, primary_edges, measurements, lambda_reg)
    grid.to_csv(out / "summary/choice_grid.tsv", sep="\t", index=False)
    permutations = _permutations(root, out, lambda_reg)
    permutations.to_csv(out / "summary/permutation_rank.tsv", sep="\t", index=False)
    cosmos_measurements = None
    cosmos_lambda = None
    lambda_grid = pd.read_csv(root / STAGE04 / "summary/lambda_grid.tsv", sep="\t")
    cosmos_selected = lambda_grid[lambda_grid["run"].eq("cosmos_multi")]
    selected_status = pd.read_csv(root / STAGE04 / "summary/selected_lambda.tsv", sep="\t")
    cosmos_status = selected_status[selected_status["run"].eq("cosmos_multi")]
    cosmos_selected_lambda = float(cosmos_status.iloc[0]["selected_lambda"]) if len(cosmos_status) == 1 else None
    if cosmos_selected_lambda is not None:
        cosmos_selected = cosmos_selected[cosmos_selected["lambda_reg"].eq(cosmos_selected_lambda)]
    if cosmos_selected_lambda is not None and len(cosmos_selected) == 1 and bool(cosmos_selected.iloc[0]["accepted"]) and float(cosmos_selected.iloc[0]["solve_seconds"]) < 1800:
        cosmos_tf = measured[measured["pkn"].eq("cosmos_multi")]
        metabolites = pd.read_csv(root / STAGE04 / "inputs/measured_metabolites.tsv", sep="\t")
        cosmos_measurements = {
            "PMS_vs_Ctrl": {
                **dict(zip(cosmos_tf["TF"], cosmos_tf["score"])),
                **dict(zip(metabolites["pkn_node"], metabolites["value"])),
            }
        }
        cosmos_lambda = cosmos_selected_lambda
        cosmos_loo_status, cosmos_loo_edges = _loo_controls(
            root,
            out,
            cosmos_lambda,
            pkn_edges=read_pkn(root / PKN / "cosmos/primary/pkn.tsv"),
            directory_prefix="cosmos",
            time_limit=COSMOS_TIME_LIMIT,
            extra_measurements={
                str(node): float(value)
                for node, value in zip(metabolites["pkn_node"], metabolites["value"])
            },
        )
        loo_status = pd.concat([loo_status, cosmos_loo_status], ignore_index=True)
        loo_edges.extend(cosmos_loo_edges)
        loo_status.to_csv(out / "summary/loo_status.tsv", sep="\t", index=False)
    rewired = _rewired(root, out, measurements, lambda_reg, cosmos_measurements, cosmos_lambda)
    rewired.to_csv(out / "summary/rewired.tsv", sep="\t", index=False)
    jaccard_rows = []
    primary_edge_set = _edge_set(primary_edges)
    for loo_edges_table in loo_edges:
        left_out = str(loo_edges_table["loo"].iloc[0])
        jaccard_rows.append({
            "analysis": "leave_one_line_out",
            "level": left_out,
            "edge_jaccard": _jaccard(_edge_set(loo_edges_table), primary_edge_set),
            "node_jaccard": _jaccard(
                {node for edge in _edge_set(loo_edges_table) for node in edge},
                {node for edge in primary_edge_set for node in edge},
            ),
        })
    for row in grid.itertuples(index=False):
        jaccard_rows.append({
            "analysis": "choice_grid",
            "level": f"{row.factor}={row.level}",
            "edge_jaccard": row.edge_jaccard,
            "node_jaccard": row.node_jaccard,
        })
    pd.DataFrame(jaccard_rows).to_csv(out / "summary/jaccard.tsv", sep="\t", index=False)

    if not frequencies.empty or loo_edges:
        core_records = []
        for row in frequencies[frequencies["run"].eq("omnipath")].itertuples(index=False) if not frequencies.empty else []:
            core_records.append({
                "source": row.source,
                "target": row.target,
                "n_present": row.n_alternatives,
                "n_steps": row.n_total_alternatives,
            })
        for loo_table in loo_edges:
            if "pkn" in loo_table.columns and loo_table["pkn"].iloc[0] != "omnipath":
                continue
            for source, target in sorted(_edge_set(loo_table)):
                core_records.append({"source": source, "target": target, "n_present": 1, "n_steps": 1})
        if core_records:
            core = pd.DataFrame(core_records).groupby(["source", "target"], as_index=False).agg(
                n_present=("n_present", "sum"), n_steps=("n_steps", "sum")
            )
            core["frequency"] = core["n_present"] / core["n_steps"]
            core["core"] = core["frequency"].ge(0.8)
            core.to_csv(out / "summary/core_network.tsv", sep="\t", index=False)
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.hist(core["frequency"], bins=np.linspace(0, 1, 11), color="#0072B2", edgecolor="white")
            ax.set(xlabel="Combined edge frequency", ylabel="Edges", title="Stage-05 edge stability")
            save_figure(fig, out / "figures/edge_frequency")
            plt.close(fig)
    rewired_nodes = set()
    for edge_file in (out / "rewired/omnipath").glob("*/edges.tsv"):
        rewired_edges = pd.read_csv(edge_file, sep="\t")
        if not rewired_edges.empty:
            rewired_nodes.update(rewired_edges["source"])
            rewired_nodes.update(rewired_edges["target"])
    _hub_bias(root, frequencies, out, rewired_nodes)

    inputs = [
        stage04,
        root / STAGE04 / "inputs/measured_tfs.tsv",
        root / STAGE04 / "summary/selected_lambda.tsv",
        root / PREPROCESSED / "rnaseq/log2cpm.tsv",
        root / PREPROCESSED / "rnaseq/samples.tsv",
        root / PLAN,
    ]
    code = [
        root / "src/inference.py", root / "src/activities.py", root / "src/pkn.py",
        root / "scripts/analysis/06_corneto_robustness.py", root / PLAN,
    ]
    parameters = {
        "alternative_solutions": {"max_samples": 30, "percentage": 0.03, "scale": 0.03, "rel_opt_tol": 0.05},
        "seed_base": SEED_BASE,
        "permutation_assignments": 20,
        "omnipath_rewired_replicates": 10,
        "solver": SOLVER,
        "solver_threads": SOLVER_THREADS,
        "solver_seed": SOLVER_SEED,
    }
    provenance = stage_provenance(root, "05_robustness", inputs, {}, parameters, code)
    write_stage_provenance(out, provenance)
    print(f"Stage 05 complete: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
