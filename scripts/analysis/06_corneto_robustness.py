"""Stage 05: stability, analysis-choice and topology controls for stage 04.

Implements steps 16-21 of doc/decisions/2026-10-10-stage04-inference-plan.md
and the controls added in doc/decisions/2026-10-11-stage05-robustness-addendum.md
(tighter-gap re-solve, simvastatin positive control, receptor-only inputs).

Every CORNETO solve runs in its own Python process (src/corneto_worker.py),
because solving several problems in one process exhausted memory twice in
stage 04.  Finished solves are reused with ``--resume``.

The OmniPath controls (about 3 h on the course laptop) run by default.  COSMOS
controls take 45-95 min per solve and run only with ``--cosmos``:

    python scripts/analysis/06_corneto_robustness.py                 # OmniPath + SV + receptors
    python scripts/analysis/06_corneto_robustness.py --resume --cosmos
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from functools import lru_cache
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from activities import (  # noqa: E402
    leave_one_line_out_ulm,
    load_regulon,
    moderated_t_statistics,
    one_sample_moderated_t,
    run_ulm,
    select_stable_top_tfs,
)
from corneto_worker import make_task, solve_in_subprocess  # noqa: E402
from inference import (  # noqa: E402
    CANONICAL_CHOLESTEROL_GENES,
    COSMOS_TIME_LIMIT,
    LAMBDA_GRID,
    MEVALONATE_HMDB,
    OMNIPATH_TIME_LIMIT,
    SOLVER,
    SOLVER_RELATIVE_GAP,
    SOLVER_SEED,
    SOLVER_THREADS,
    choose_lambda,
    is_excluded_input,
    lambda_label,
    read_pkn,
)
from pkn import _fetch_bytes, rewire_signed_graph, rewire_type_preserving_signed_graph  # noqa: E402
from plot_config import save_figure  # noqa: E402
from preprocess import stage_provenance, write_stage_provenance  # noqa: E402

ANALYSIS = "results/ionescu_corneto/05_robustness"
STAGE04 = "results/ionescu_corneto/04_inference"
ACTIVITIES = "results/ionescu_corneto/02_activities"
PREPROCESSED = "results/ionescu_corneto/01_preprocessing"
PKN = "results/ionescu_corneto/03_pkn"
PLAN = "doc/decisions/2026-10-10-stage04-inference-plan.md"
ADDENDUM = "doc/decisions/2026-10-11-stage05-robustness-addendum.md"
SEED_BASE = 20261010

PARTS = ("sampling", "loo", "grid", "permutations", "rewired", "receptors", "sv_control", "tight_gap",
         "no_drug_inputs", "late_integration", "context_cosmos")
STAGE04_RUNS = ("omnipath", "cosmos_tf", "cosmos_multi", "cosmos_repair", "omnipath_lines")
COSMOS_RUNS = frozenset({"cosmos_tf", "cosmos_multi", "cosmos_repair"})
RUN_PKN = {
    "omnipath": "omnipath", "omnipath_lines": "omnipath",
    "cosmos_tf": "cosmos", "cosmos_multi": "cosmos", "cosmos_repair": "cosmos_repair",
}
N_SAMPLES = 30
DEFAULT_COSMOS_SAMPLES = 10
N_OMNIPATH_REWIRED = 10
N_COSMOS_REWIRED = 5
CORE_THRESHOLD = 0.8
TIGHT_GAP = 1e-4
TIGHT_LAMBDAS = (0.0, 0.01, 0.03, 0.1, 0.3)
SV_LAMBDA = 0.01
RECEPTOR_URL = "https://omnipathdb.org/intercell?categories=receptor&format=tsv"
RECEPTOR_CACHE = "data/resources/pkn/omnipath_intercell_receptor.tsv"
RECEPTOR_MIN_CONSENSUS = 3


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------

def _pkn_path(root: Path, name: str) -> Path:
    paths = {
        "omnipath": PKN + "/omnipath/primary/pkn.tsv",
        "nc_rule": PKN + "/omnipath/variants/nc_rule/pkn.tsv",
        "signor_only": PKN + "/omnipath/variants/signor_only/pkn.tsv",
        "unpruned": PKN + "/omnipath/variants/unpruned/pkn.tsv",
        "cosmos": PKN + "/cosmos/primary/pkn.tsv",
        "cosmos_repair": PKN + "/cosmos/variants/canonical_repair/pkn.tsv",
    }
    return root / paths[name]


def _time_limit(pkn_name: str) -> int:
    return COSMOS_TIME_LIMIT if pkn_name.startswith("cosmos") else OMNIPATH_TIME_LIMIT


def _nodes(edges: list[tuple[str, int, str]]) -> set[str]:
    return {node for source, _, target in edges for node in (source, target)}


def _stage04_selected(root: Path) -> dict[str, float]:
    table = pd.read_csv(root / STAGE04 / "summary/selected_lambda.tsv", sep="\t")
    return {str(row.run): float(row.selected_lambda) for row in table.itertuples(index=False)}


def _stage04_edges(root: Path, run: str, lambda_reg: float) -> pd.DataFrame:
    return pd.read_csv(root / STAGE04 / run / f"lambda_{lambda_label(lambda_reg)}/edges.tsv", sep="\t",
                       dtype={"source": str, "target": str})


def _stage04_measurements(root: Path, run: str) -> dict[str, dict[str, float]]:
    """Rebuild the measurements stage 04 used for ``run`` from its input tables."""
    measured = pd.read_csv(root / STAGE04 / "inputs/measured_tfs.tsv", sep="\t")
    tfs = measured[measured["pkn"].eq(run)]
    values = {str(tf): float(score) for tf, score in zip(tfs["TF"], tfs["score"])}
    if run == "omnipath_lines":
        per_line = pd.read_csv(root / ACTIVITIES / "collectri/tf_activity_per_line.tsv", sep="\t", index_col=0)
        return {str(line): {tf: float(per_line.loc[tf, line]) for tf in values} for line in per_line.columns}
    if run in ("cosmos_multi", "cosmos_repair"):
        metabolites = pd.read_csv(root / STAGE04 / "inputs/measured_metabolites.tsv", sep="\t")
        values.update({str(node): float(value) for node, value in zip(metabolites["pkn_node"], metabolites["value"])})
    return {"PMS_vs_Ctrl": values}


def _check_against_stage04_task(root: Path, run: str, lambda_reg: float, measurements: dict) -> None:
    """Where stage 04 left a task.json, the rebuilt measurements must match it."""
    path = root / STAGE04 / run / f"lambda_{lambda_label(lambda_reg)}/task.json"
    if not path.exists():
        return
    stored = json.loads(path.read_text())["measurements"]
    for sample, values in measurements.items():
        other = stored.get(sample, {})
        if set(other) != set(values) or any(abs(float(other[node]) - values[node]) > 1e-9 for node in values):
            raise ValueError(f"Rebuilt measurements for {run} differ from {path}")


@lru_cache(maxsize=None)
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
    return pd.read_csv(root / ACTIVITIES / regulon_name / file_name, sep="\t", index_col=0)["score"]


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _edge_set(edges: pd.DataFrame) -> set[tuple[str, str]]:
    if edges is None or edges.empty:
        return set()
    return {(str(row.source), str(row.target)) for row in edges.itertuples(index=False)}


def _node_set(edges: set[tuple[str, str]]) -> set[str]:
    return {node for edge in edges for node in edge}


def _jaccard(left: set[Any], right: set[Any]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else np.nan


def _solve(
    root: Path, out: Path, subdir: str, label: str, pkn_path: Path,
    measurements: dict[str, dict[str, float]], lambda_reg: float, time_limit: int, *,
    resume: bool, gap: float = SOLVER_RELATIVE_GAP, inputs: dict[str, float] | None = None,
    max_samples: int = 0,
) -> dict[str, Any]:
    task = make_task(
        label=label, pkn_path=pkn_path, measurements=measurements, lambda_reg=lambda_reg,
        time_limit=time_limit, solver_relative_gap=gap, inputs=inputs, max_samples=max_samples, root=root,
    )
    return solve_in_subprocess(task, out / subdir, resume=resume, root=root)


def _fit_columns(result: dict[str, Any]) -> dict[str, Any]:
    """Fit summary shared by every control table (NaN when unsolved)."""
    fit = result.get("fit")
    has_fit = isinstance(fit, pd.DataFrame) and not fit.empty
    total = (lambda column: float(fit[column].sum()) if has_fit and column in fit else np.nan)
    edges = _edge_set(result["selected_edges"])
    return {
        "status": result["status"],
        "accepted": bool(result["accepted"]),
        "mip_gap": result["mip_gap"],
        "solve_seconds": result["solve_seconds"],
        "weighted_fitting_error": result["weighted_fitting_error"],
        "n_edges": len(edges) if result["accepted"] else np.nan,
        "n_nodes": len(_node_set(edges)) if result["accepted"] else np.nan,
        "n_measurements": total("n_measurements"),
        "n_reachable": total("n_reachable"),
        "n_measurements_correct_sign": total("n_measurements_correct_sign"),
        "n_inputs_offered": len(result.get("inputs", [])),
    }


def _comparison(result: dict[str, Any], reference: set[tuple[str, str]]) -> dict[str, float]:
    if not result["accepted"]:
        return {"edge_jaccard": np.nan, "node_jaccard": np.nan}
    edges = _edge_set(result["selected_edges"])
    return {"edge_jaccard": _jaccard(edges, reference), "node_jaccard": _jaccard(_node_set(edges), _node_set(reference))}


def _write(table: pd.DataFrame, out: Path, name: str) -> pd.DataFrame:
    (out / "summary").mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "summary" / name, sep="\t", index=False)
    return table


def _write_pkn(path: Path, edges: list[tuple[str, int, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(edges, columns=["source", "sign", "target"])
    path.write_bytes(frame.to_csv(sep="\t", index=False, lineterminator="\n").encode("utf-8"))
    return path


# --------------------------------------------------------------------------
# step 16: alternative solutions
# --------------------------------------------------------------------------

def _sampling(root: Path, out: Path, runs: list[str], samples: dict[str, int], resume: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Re-solve each selected stage-04 problem and sample alternatives in the same worker.

    The ensemble for edge frequencies is the re-solved optimum plus every
    alternative within the 5 % objective tolerance (addendum, item 1).
    """
    selected = _stage04_selected(root)
    rows, frequencies = [], []
    for run in runs:
        if run not in selected:
            rows.append({"run": run, "status": "stage04_no_selected_lambda"})
            continue
        lambda_reg = selected[run]
        if lambda_reg == 0.0:
            # Without the sparsity penalty every sign-consistent edge set is
            # optimal, so alternatives say nothing about stability (omnipath_lines:
            # no incumbent for any lambda > 0 within 1,800 s in stage 04).
            rows.append({"run": run, "lambda_reg": lambda_reg, "status": "skipped_selected_lambda_0"})
            continue
        measurements = _stage04_measurements(root, run)
        _check_against_stage04_task(root, run, lambda_reg, measurements)
        pkn_name = RUN_PKN[run]
        result = _solve(root, out, f"sampling/{run}", f"sampling {run}", _pkn_path(root, pkn_name),
                        measurements, lambda_reg, _time_limit(pkn_name), resume=resume, max_samples=samples[run])
        stage04 = _edge_set(_stage04_edges(root, run, lambda_reg))
        rows.append({
            "run": run, "lambda_reg": lambda_reg, "requested_alternatives": samples[run],
            **_fit_columns(result), "sampling_status": result.get("sampling_status"),
            "n_alternatives": result.get("n_alternatives", 0),
            "resolve_edge_jaccard_with_stage04": _comparison(result, stage04)["edge_jaccard"],
        })
        if not result["accepted"]:
            continue
        networks = [_edge_set(result["selected_edges"])]
        alternatives = result["alternatives"]
        if not alternatives.empty:
            networks += [
                {(str(s), str(t)) for s, t in group[["source", "target"]].itertuples(index=False, name=None)}
                for _, group in alternatives.groupby("alternative")
            ]
        counts: dict[tuple[str, str], int] = {}
        for network in networks:
            for edge in network:
                counts[edge] = counts.get(edge, 0) + 1
        frequencies.append(pd.DataFrame(
            [{"run": run, "source": s, "target": t, "n_present": n, "n_networks": len(networks),
              "frequency": n / len(networks), "in_stage04_network": (s, t) in stage04}
             for (s, t), n in sorted(counts.items())]
        ))
    frequency_table = pd.concat(frequencies, ignore_index=True) if frequencies else pd.DataFrame(
        columns=["run", "source", "target", "n_present", "n_networks", "frequency", "in_stage04_network"])
    return pd.DataFrame(rows), frequency_table


# --------------------------------------------------------------------------
# step 17: leave one line out
# --------------------------------------------------------------------------

def _loo(root: Path, out: Path, run: str, lambda_reg: float, resume: bool) -> tuple[pd.DataFrame, list[pd.DataFrame]]:
    """Redo TF selection and the solve without each RNA-seq line.

    With one line removed the LOO-of-LOO stability filter would need 1-vs-3
    contrasts, which the moderated t cannot estimate; the stage-04 six-LOO
    sign filter is therefore applied to the LOO contrast (addendum, item 2).
    Metabolite measurements (cosmos_multi) come from different cell lines and
    are kept unchanged.
    """
    expression, _, _, full_loo = _expression_inputs(root)
    pkn_name = RUN_PKN[run]
    pkn_path = _pkn_path(root, pkn_name)
    pkn_nodes = _nodes(read_pkn(pkn_path))
    extra = {}
    if run == "cosmos_multi":
        extra = {node: value for node, value in _stage04_measurements(root, run)["PMS_vs_Ctrl"].items()
                 if node.startswith("Metab__")}
    reference = _edge_set(_stage04_edges(root, run, lambda_reg))
    rows, tables = [], []
    for left_out in expression.columns:
        contrast = full_loo[f"without_{left_out}"]
        selected, selection = select_stable_top_tfs(contrast, full_loo, pkn_nodes, k=50)
        values = {tf: float(contrast[tf]) for tf in selected}
        values.update(extra)
        subdir = f"loo/{run}/{left_out}"
        result = _solve(root, out, subdir, f"loo {run} without {left_out}", pkn_path,
                        {f"without_{left_out}": values}, lambda_reg, _time_limit(pkn_name), resume=resume)
        selection.to_csv(out / subdir / "tf_selection.tsv", sep="\t", index=False)
        rows.append({"run": run, "left_out": left_out, "n_selected_tfs": len(selected),
                     **_fit_columns(result), **_comparison(result, reference)})
        if result["accepted"]:
            tables.append(result["selected_edges"].assign(loo=left_out, run=run))
    return pd.DataFrame(rows), tables


# --------------------------------------------------------------------------
# step 18: analysis-choice grid (OmniPath)
# --------------------------------------------------------------------------

def _choice_grid(root: Path, out: Path, lambda_reg: float, resume: bool) -> pd.DataFrame:
    """One factor at a time around the stage-04 OmniPath run.

    For the Welch and logFC levels the LOO sign filter still uses the
    moderated-t LOO contrasts (stage 02 stores no per-statistic LOO).
    """
    reference = _edge_set(_stage04_edges(root, "omnipath", lambda_reg))
    factors: list[tuple[str, Any]] = [
        ("statistic", "welch"), ("statistic", "logFC"),
        ("regulon", "dorothea_ABC"), ("regulon", "collectri_nocomplex"),
        ("k", 25), ("k", 100),
        ("pkn", "nc_rule"), ("pkn", "signor_only"), ("pkn", "unpruned"),
    ]
    index = list(LAMBDA_GRID).index(lambda_reg)
    factors += [("lambda", value) for value in LAMBDA_GRID[max(0, index - 1):index + 2]]
    rows = []
    for factor, level in factors:
        pkn_name = str(level) if factor == "pkn" else "omnipath"
        pkn_path = _pkn_path(root, pkn_name)
        statistic = str(level) if factor == "statistic" else "moderated_t"
        regulon_name = str(level) if factor == "regulon" else "collectri"
        k = int(level) if factor == "k" else 50
        _, _, _, loo = _expression_inputs(root, regulon_name)
        primary = _primary_score(root, statistic, regulon_name)
        selected, _ = select_stable_top_tfs(primary, loo, _nodes(read_pkn(pkn_path)), k=k)
        solve_lambda = float(level) if factor == "lambda" else lambda_reg
        result = _solve(root, out, f"grid/{factor}/{level}", f"grid {factor}={level}", pkn_path,
                        {"PMS_vs_Ctrl": {tf: float(primary[tf]) for tf in selected}},
                        solve_lambda, OMNIPATH_TIME_LIMIT, resume=resume)
        rows.append({
            "factor": factor, "level": str(level),
            "is_stage04_setting": factor == "lambda" and solve_lambda == lambda_reg,
            "n_selected_tfs": len(selected), **_fit_columns(result), **_comparison(result, reference),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# step 19: group-label permutations (OmniPath)
# --------------------------------------------------------------------------

def _permutations(root: Path, out: Path, lambda_reg: float, resume: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    """All 20 three-vs-three assignments of the six RNA-seq lines.

    Inputs are free-sign sources, so most assignments are fitted perfectly
    (error 0) and the fit error cannot rank them.  The table therefore also
    reports signal statistics computed before CORNETO (number of LOO-stable
    TFs, summed |score| of the selected TFs) and the rank of the true split
    for each statistic (addendum, item 3).
    """
    expression, _, regulon, _ = _expression_inputs(root)
    pkn_path = _pkn_path(root, "omnipath")
    pkn_nodes = _nodes(read_pkn(pkn_path))
    stage04_tfs = set(pd.read_csv(root / STAGE04 / "inputs/measured_tfs.tsv", sep="\t").query("pkn == 'omnipath'")["TF"])
    rows, networks = [], {}
    true_ctrl, mirror_ctrl = "C1,C2,C3", "C4,C5,C6"
    for index, ctrl_tuple in enumerate(itertools.combinations(expression.columns, 3)):
        ctrl = set(ctrl_tuple)
        groups = pd.Series(["Ctrl" if line in ctrl else "PMS" for line in expression.columns], index=expression.columns)
        moderated, _ = moderated_t_statistics(expression, groups)
        scores, _, _, _ = run_ulm(pd.DataFrame([moderated["moderated_t"]], index=["PMS_vs_Ctrl"]), regulon, tmin=5)
        loo = leave_one_line_out_ulm(expression, groups, regulon, tmin=5)
        contrast = scores["PMS_vs_Ctrl"]  # run_ulm returns TF x observation
        selected, _ = select_stable_top_tfs(contrast, loo, pkn_nodes, k=50)
        values = {tf: float(contrast[tf]) for tf in selected}
        result = _solve(root, out, f"permutation/{index}", f"permutation {index}", pkn_path,
                        {"PMS_vs_Ctrl": values}, lambda_reg, OMNIPATH_TIME_LIMIT, resume=resume)
        ctrl_label = ",".join(sorted(ctrl))
        networks[index] = _edge_set(result["selected_edges"]) if result["accepted"] else None
        rows.append({
            "assignment": index, "ctrl_lines": ctrl_label,
            "pms_lines": ",".join(sorted(set(expression.columns) - ctrl)),
            "role": "true" if ctrl_label == true_ctrl else "mirror_of_true" if ctrl_label == mirror_ctrl else "permuted",
            "n_stable_tfs": len(selected),
            "sum_abs_score_selected": float(sum(abs(value) for value in values.values())),
            "matches_stage04_selection": set(selected) == stage04_tfs if ctrl_label == true_ctrl else np.nan,
            **_fit_columns(result),
        })
    table = pd.DataFrame(rows)
    true_index = int(table.loc[table["role"].eq("true"), "assignment"].iloc[0])
    true_edges = networks[true_index]
    table["edge_jaccard_with_true"] = [
        _jaccard(networks[i], true_edges) if networks[i] is not None and true_edges is not None else np.nan
        for i in table["assignment"]
    ]
    summary_rows = []
    for statistic, ascending in (("n_stable_tfs", False), ("sum_abs_score_selected", False),
                                 ("weighted_fitting_error", True), ("n_edges", True)):
        ranks = table[statistic].rank(method="min", ascending=ascending)
        summary_rows.append({
            "statistic": statistic, "better_is": "lower" if ascending else "higher",
            "true_value": table.loc[table["role"].eq("true"), statistic].iloc[0],
            "true_rank_of_20": (lambda r: int(r) if pd.notna(r) else np.nan)(ranks[table["role"].eq("true")].iloc[0]),
            "n_ties_with_true": int((table[statistic] == table.loc[table["role"].eq("true"), statistic].iloc[0]).sum() - 1),
            "median_permuted": float(table.loc[table["role"].eq("permuted"), statistic].median()),
        })
    return table, pd.DataFrame(summary_rows)


# --------------------------------------------------------------------------
# step 20: rewired PKNs
# --------------------------------------------------------------------------

def _rewired(root: Path, out: Path, pkn_name: str, run: str, lambda_reg: float, replicates: int, resume: bool) -> pd.DataFrame:
    """Solve the stage-04 measurements on degree/sign-preserving rewired PKNs.

    Overlap is measured against the real stage-04 *selected* network (nodes
    and edges), not against the whole PKN (addendum, item 4).
    """
    edges = read_pkn(_pkn_path(root, pkn_name))
    measurements = _stage04_measurements(root, run)
    reference = _edge_set(_stage04_edges(root, run, lambda_reg))
    rows = []
    for replicate in range(replicates):
        seed = SEED_BASE + replicate
        rewired = (rewire_type_preserving_signed_graph if pkn_name.startswith("cosmos") else rewire_signed_graph)(edges, seed=seed)
        pkn_path = _write_pkn(out / "rewired" / pkn_name / str(replicate) / "pkn.tsv", rewired)
        result = _solve(root, out, f"rewired/{pkn_name}/{replicate}", f"rewired {pkn_name} {replicate}", pkn_path,
                        measurements, lambda_reg, _time_limit(pkn_name), resume=resume)
        rows.append({"pkn": pkn_name, "run": run, "replicate": replicate, "seed": seed,
                     **_fit_columns(result), **_comparison(result, reference)})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# step 21: hub bias
# --------------------------------------------------------------------------

def _hub_bias(root: Path, frequencies: pd.DataFrame, out: Path) -> pd.DataFrame:
    """Logistic regression of node selection on log2 PKN degree.

    Measured TFs are forced into every network as outputs and are removed
    from the regression (addendum, item 5).
    """
    edges = read_pkn(_pkn_path(root, "omnipath"))
    degree: dict[str, int] = {}
    for source, _, target in edges:
        degree[source] = degree.get(source, 0) + 1
        degree[target] = degree.get(target, 0) + 1
    measured = set(pd.read_csv(root / STAGE04 / "inputs/measured_tfs.tsv", sep="\t").query("pkn == 'omnipath'")["TF"])
    nodes = sorted(set(degree) - measured)
    networks = {}
    real = frequencies[frequencies["run"].eq("omnipath")] if not frequencies.empty else frequencies
    if not real.empty:
        networks["omnipath_ensemble"] = set(real["source"]) | set(real["target"])
    rewired_nodes: set[str] = set()
    for path in sorted((out / "rewired/omnipath").glob("*/edges.tsv")):
        table = pd.read_csv(path, sep="\t", dtype={"source": str, "target": str})
        rewired_nodes.update(table["source"]); rewired_nodes.update(table["target"])
    if rewired_nodes:
        networks["rewired_omnipath_union"] = rewired_nodes
    rows = []
    log_degree = sm.add_constant(pd.Series([np.log2(degree[node]) for node in nodes], name="log2_degree"))
    for network, selected in networks.items():
        outcome = pd.Series([int(node in selected) for node in nodes], name="selected")
        try:
            fit = sm.Logit(outcome, log_degree).fit(disp=False)
            odds, pvalue, note = float(np.exp(fit.params["log2_degree"])), float(fit.pvalues["log2_degree"]), ""
        except Exception as error:  # e.g. perfect separation
            odds, pvalue, note = np.nan, np.nan, f"{type(error).__name__}: {error}"
        rows.append({"network": network, "n_nodes_tested": len(nodes), "n_selected": int(outcome.sum()),
                     "odds_ratio_per_degree_doubling": odds, "pvalue": pvalue, "note": note})
    return pd.DataFrame(rows, columns=["network", "n_nodes_tested", "n_selected",
                                       "odds_ratio_per_degree_doubling", "pvalue", "note"])


def _core_network(frequencies: pd.DataFrame, loo_tables: list[pd.DataFrame], out: Path) -> pd.DataFrame:
    """Edge frequency over the OmniPath ensemble plus accepted LOO networks.

    Every edge shares one denominator: networks in the ensemble plus accepted
    LOO networks.
    """
    networks: list[set[tuple[str, str]]] = []
    real = frequencies[frequencies["run"].eq("omnipath")] if not frequencies.empty else frequencies
    counts: dict[tuple[str, str], int] = {}
    n_ensemble = int(real["n_networks"].iloc[0]) if not real.empty else 0
    for row in real.itertuples(index=False):
        counts[(row.source, row.target)] = int(row.n_present)
    for table in loo_tables:
        if table["run"].iloc[0] == "omnipath":
            networks.append(_edge_set(table))
    for network in networks:
        for edge in network:
            counts[edge] = counts.get(edge, 0) + 1
    total = n_ensemble + len(networks)
    if total == 0:
        return pd.DataFrame(columns=["source", "target", "n_present", "n_networks", "frequency", "core"])
    core = pd.DataFrame([{"source": s, "target": t, "n_present": n} for (s, t), n in sorted(counts.items())])
    core["n_networks"] = total
    core["frequency"] = core["n_present"] / total
    core["core"] = core["frequency"].ge(CORE_THRESHOLD)
    fig, ax = plt.subplots(figsize=(4.5, 3))
    ax.hist(core["frequency"], bins=np.linspace(0, 1, 11), color="#0072B2", edgecolor="white")
    ax.axvline(CORE_THRESHOLD, color="0.3", linestyle="--", linewidth=1)
    ax.set(xlabel=f"Edge frequency ({n_ensemble} sampled + {len(networks)} LOO networks)", ylabel="Edges")
    save_figure(fig, out / "figures/edge_frequency")
    plt.close(fig)
    return core


# --------------------------------------------------------------------------
# addendum controls
# --------------------------------------------------------------------------

def _receptor_list(root: Path) -> tuple[set[str], dict[str, Any]]:
    """OmniPath intercell receptors (cached with a checksum sidecar).

    Only the OmniPath composite (consensus) annotation of transmembrane
    proteins with consensus score >= 3 is used.  The resource-specific rows
    also label intracellular kinases (MAPK1, MAPK14, GSK3B) as receptors.
    """
    path, metadata = _fetch_bytes(root / RECEPTOR_CACHE, RECEPTOR_URL)
    table = pd.read_csv(path, sep="\t", dtype=str)
    required = {"category", "genesymbol", "source", "entity_type", "consensus_score", "plasma_membrane_transmembrane"}
    if not required.issubset(table.columns):
        raise ValueError(f"OmniPath intercell table lacks {sorted(required - set(table.columns))}: {path}")
    table = table[
        table["category"].eq("receptor") & table["source"].eq("composite") & table["entity_type"].eq("protein")
        & table["plasma_membrane_transmembrane"].eq("True")
        & pd.to_numeric(table["consensus_score"], errors="coerce").ge(RECEPTOR_MIN_CONSENSUS)
    ]
    receptors = {symbol for symbol in table["genesymbol"].dropna() if symbol}
    if not receptors:
        raise ValueError("No receptors left after filtering the OmniPath intercell table")
    return receptors, metadata


def _receptor_inputs(root: Path, out: Path, lambda_reg: float, resume: bool) -> pd.DataFrame:
    """OmniPath run with inputs restricted to expressed receptors in the PKN."""
    receptors, metadata = _receptor_list(root)
    pkn_path = _pkn_path(root, "omnipath")
    edges = read_pkn(pkn_path)
    measurements = _stage04_measurements(root, "omnipath")
    measured = set(measurements["PMS_vs_Ctrl"])
    inputs = {node: 0.0 for node in sorted(_nodes(edges) & receptors)
              if node not in measured and not is_excluded_input(node)}
    sources = {source for source, _, _ in edges} - {target for _, _, target in edges}
    (out / "receptors").mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"node": sorted(inputs), "is_pkn_source": [node in sources for node in sorted(inputs)]}).to_csv(
        out / "receptors/inputs.tsv", sep="\t", index=False)
    result = _solve(root, out, "receptors/omnipath", "receptor-only inputs", pkn_path, measurements,
                    lambda_reg, OMNIPATH_TIME_LIMIT, resume=resume, inputs=inputs)
    reference = _edge_set(_stage04_edges(root, "omnipath", lambda_reg))
    return pd.DataFrame([{
        "run": "omnipath", "lambda_reg": lambda_reg, "n_receptors_in_resource": len(receptors),
        "n_receptor_inputs": len(inputs), "n_receptor_inputs_that_are_pkn_sources": len(set(inputs) & sources),
        "receptor_resource_sha256": metadata.get("sha256"), **_fit_columns(result), **_comparison(result, reference),
    }])


def _sv_measurement(root: Path) -> tuple[str, pd.DataFrame]:
    """Free cholesterol, simvastatin vs untreated, paired within each line."""
    directory = root / PREPROCESSED / "metabolomics/ST003328"
    features = pd.read_csv(directory / "features.tsv", sep="\t", dtype=str).set_index("feature_id")
    values = pd.read_csv(directory / "line_log2.tsv", sep="\t", index_col=0)
    lines = sorted({column.rsplit("_", 1)[0] for column in values.columns})
    if any(f"{line}_SV" not in values or f"{line}_untreated" not in values for line in lines):
        raise ValueError("Every ST003328 line needs an SV and an untreated column")
    differences = pd.DataFrame({line: values[f"{line}_SV"] - values[f"{line}_untreated"] for line in lines})
    contrast, prior = one_sample_moderated_t(differences)
    cholesterol = features.index[features["metabolite_name"].str.lower().eq("cholesterol")]
    if len(cholesterol) != 1 or cholesterol[0] not in contrast.index:
        raise ValueError("Expected one complete free-cholesterol feature in ST003328")
    metabolites = pd.read_csv(root / STAGE04 / "inputs/measured_metabolites.tsv", sep="\t")
    node = str(metabolites.loc[metabolites["study"].eq("ST003328"), "pkn_node"].iloc[0])
    row = contrast.loc[[cholesterol[0]]].assign(feature_id=cholesterol[0], pkn_node=node, n_lines=len(lines), **prior)
    return node, row.reset_index(drop=True)


def _reaction_genes(nodes_path: Path) -> dict[str, set[str]]:
    nodes = pd.read_csv(nodes_path, sep="\t", dtype=str)
    reactions = nodes[nodes["type"].eq("enzyme_or_reaction")]
    return {str(row.node): {s for s in str(row.approved_symbol).split(";") if s and s != "nan"}
            for row in reactions.itertuples(index=False)}


def _metabolic_pkn(root: Path, out: Path, pkn_name: str) -> Path:
    """Metabolism-only, cofactor-free signed PKN (the Aim-5 route graph).

    Keeps metabolite-reaction edges and each enzyme's edges to its own
    reactions (src/metabolic_route.metabolism_graph) and drops the registered
    non-protected cofactor nodes of doc/decisions/2026-10-10-aim5-metabolic-route-plan.md.
    """
    from metabolic_route import metabolism_graph, read_pkn as read_pkn_table

    table = read_pkn_table(_pkn_path(root, pkn_name))
    cofactors = pd.read_csv(root / PKN / "aim5_route/cofactor_nodes.tsv", sep="\t")
    drop = set(cofactors.loc[cofactors["mapped"].astype(bool) & ~cofactors["protected"].astype(bool), "pkn_node"])
    keep = {(u, v) for u, v in metabolism_graph(table).edges() if u not in drop and v not in drop}
    edges = [(str(s), int(sign), str(t)) for s, sign, t in table[["source", "sign", "target"]].itertuples(index=False, name=None)
             if (str(s), str(t)) in keep]
    return _write_pkn(out / "sv_control" / f"{pkn_name}_metabolic" / "pkn.tsv", edges)


def _sv_control(root: Path, out: Path, resume: bool, full: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Positive control: known input HMGCR = -1, measured SV cholesterol drop.

    Asks whether CORNETO connects a known perturbation to its known metabolic
    consequence, and through which reactions, on the full and the
    metabolism-only (cofactor-free) PKN, each primary and canonically
    repaired (addendum, item B).
    """
    node, contrast = _sv_measurement(root)
    rows = []
    for pkn_name in ("cosmos", "cosmos_repair"):
        if not _pkn_path(root, pkn_name).exists():
            rows.append({"pkn": pkn_name, "status": "missing_pkn"})
            continue
        genes = _reaction_genes(_pkn_path(root, pkn_name).with_name("nodes.tsv"))
        # Full-PKN scopes found no incumbent in 600 s trials; they run only with --sv-full.
        for scope in (("full", "metabolic") if full else ("metabolic",)):
            pkn_path = _pkn_path(root, pkn_name) if scope == "full" else _metabolic_pkn(root, out, pkn_name)
            result = _solve(root, out, f"sv_control/{pkn_name}_{scope}", f"SV control {pkn_name} {scope}", pkn_path,
                            {"SV_vs_untreated": {node: float(contrast["moderated_t"].iloc[0])}},
                            SV_LAMBDA, COSMOS_TIME_LIMIT, resume=resume, inputs={"HMGCR": -1.0})
            active = result["nodes"][result["nodes"]["state"].ne(0)] if result["accepted"] else result["nodes"].iloc[0:0]
            active_nodes = set(active["node"].astype(str))
            canonical = sorted({g for n in active_nodes for g in genes.get(n, set()) & CANONICAL_CHOLESTEROL_GENES})
            rows.append({
                "pkn": pkn_name, "scope": scope, "lambda_reg": SV_LAMBDA, "input": "HMGCR=-1", "measured_node": node,
                "measured_t": float(contrast["moderated_t"].iloc[0]), **_fit_columns(result),
                "cholesterol_state": int(active.loc[active["node"].eq(node), "state"].iloc[0]) if node in active_nodes else 0,
                "mevalonate_active": any(MEVALONATE_HMDB in n for n in active_nodes),
                "n_reaction_nodes_active": sum(n in genes for n in active_nodes),
                "n_metabolite_nodes_active": sum(n.startswith("Metab__") for n in active_nodes),
                "canonical_genes_used": ",".join(canonical),
                "n_canonical_genes_used": len(canonical),
                "n_canonical_genes_total": len(CANONICAL_CHOLESTEROL_GENES),
            })
    return pd.DataFrame(rows), contrast


def _tight_gap(root: Path, out: Path, resume: bool) -> pd.DataFrame:
    """Re-solve cosmos_multi at mip_rel_gap 1e-4 to test the lambda choice."""
    measurements = _stage04_measurements(root, "cosmos_multi")
    grid = pd.read_csv(root / STAGE04 / "summary/lambda_grid.tsv", sep="\t").query("run == 'cosmos_multi'")
    selected = _stage04_selected(root).get("cosmos_multi", np.nan)
    rows, results = [], []
    for lambda_reg in TIGHT_LAMBDAS:
        result = _solve(root, out, f"tight_gap/cosmos_multi/lambda_{lambda_label(lambda_reg)}",
                        f"tight gap cosmos_multi lambda={lambda_reg:g}", _pkn_path(root, "cosmos"), measurements, lambda_reg,
                        COSMOS_TIME_LIMIT, resume=resume, gap=TIGHT_GAP)
        results.append(result)
        stage04 = grid[grid["lambda_reg"].eq(lambda_reg)]
        reference = _edge_set(_stage04_edges(root, "cosmos_multi", lambda_reg))
        rows.append({
            "lambda_reg": lambda_reg, "requested_gap": TIGHT_GAP, **_fit_columns(result),
            "stage04_error": float(stage04["weighted_fitting_error"].iloc[0]) if len(stage04) else np.nan,
            "stage04_gap": float(stage04["mip_gap"].iloc[0]) if len(stage04) else np.nan,
            **_comparison(result, reference),
        })
    table = pd.DataFrame(rows)
    # lambda = 1 is not re-solved: its stage-04 error (28.8) is far outside the
    # 10 % band of any lambda = 0 error in this run.
    chosen = choose_lambda(results)
    table["tight_selected_lambda"] = chosen if chosen is not None else np.nan
    table["stage04_selected_lambda"] = selected
    return table


def _compound_inputs(edges: list[tuple[str, int, str]], candidates: list[str]) -> dict[str, list[str]]:
    """Split metabolite inputs into drug/compound-layer nodes and the rest.

    A metabolite source whose out-edges all go to protein nodes takes part in
    no PKN reaction: it comes from the COSMOS compound-target layer (e.g.
    tamoxifen Metab__5376_c, cinnamaldehyde Metab__HMDB0003441_c).
    """
    targets: dict[str, set[str]] = {}
    for source, _, target in edges:
        targets.setdefault(source, set()).add(target)
    compound = [n for n in candidates if n.startswith("Metab__") and targets.get(n)
                and all(not t.startswith(("Metab__", "Gene")) for t in targets[n])]
    return {"compound_layer": compound, "all_metabolites": [n for n in candidates if n.startswith("Metab__")]}


def _no_drug_inputs(root: Path, out: Path, lambda_reg: float, resume: bool) -> pd.DataFrame:
    """cosmos_multi with metabolite free inputs removed (addendum, item D).

    In stage 04 one unmeasured input carries about half the fit: tamoxifen
    (Metab__5376_c) explains 39 of 79 reachable measurements in cosmos_multi,
    and cinnamaldehyde (Metab__HMDB0003441_c) 44 of 79 in cosmos_repair.  Two
    variants: drop compound-layer metabolite inputs (only protein targets), or
    drop every metabolite input.  The nodes stay in the PKN.
    """
    pkn_path = _pkn_path(root, "cosmos")
    edges = read_pkn(pkn_path)
    measurements = _stage04_measurements(root, "cosmos_multi")
    measured = set(measurements["PMS_vs_Ctrl"])
    sources = sorted({s for s, _, _ in edges} - {t for _, _, t in edges})
    default = [n for n in sources if n not in measured and not is_excluded_input(n)]
    reference = _edge_set(_stage04_edges(root, "cosmos_multi", lambda_reg))
    stage04_fit = pd.read_csv(root / STAGE04 / f"cosmos_multi/lambda_{lambda_label(lambda_reg)}/fit.tsv", sep="\t")
    rows = []
    for variant, removed in _compound_inputs(edges, default).items():
        inputs = {n: 0.0 for n in default if n not in set(removed)}
        result = _solve(root, out, f"no_drug_inputs/{variant}", f"cosmos_multi without {variant} inputs", pkn_path,
                        measurements, lambda_reg, COSMOS_TIME_LIMIT, resume=resume, inputs=inputs)
        rows.append({
            "run": "cosmos_multi", "variant": variant, "lambda_reg": lambda_reg, "n_default_inputs": len(default),
            "n_inputs_removed": len(removed),
            "stage04_error": float(stage04_fit["weighted_fitting_error"].sum()),
            "stage04_n_correct_sign": int(stage04_fit["n_measurements_correct_sign"].sum()),
            **_fit_columns(result), **_comparison(result, reference),
        })
    return pd.DataFrame(rows)

def _layer(node: str) -> str:
    return "metabolic" if node.startswith(("Metab__", "Gene")) else "protein"


def _network_summary(name: str, edges: pd.DataFrame, nodes: pd.DataFrame, tfs: set[str], metabolites: set[str]) -> dict[str, Any]:
    """Topology of one selected network; states are CORNETO vertex values (0 = inactive)."""
    pairs = {(str(s), str(t)) for s, t in zip(edges["source"], edges["target"])}
    graph = nx.DiGraph(list(pairs))
    reachable = {n: nx.descendants(graph, n) for n in metabolites if n in graph}
    linked = {(m, tf) for m, down in reachable.items() for tf in down & tfs}
    return {
        "network": name, "n_edges": len(pairs), "n_nodes": graph.number_of_nodes(),
        "n_weak_components": nx.number_weakly_connected_components(graph) if graph.number_of_nodes() else 0,
        "n_cross_layer_edges": sum(_layer(s) != _layer(t) for s, t in pairs),
        "n_metabolite_to_tf_pairs": len(linked),
        "n_measured_metabolites_upstream_of_a_tf": len({m for m, _ in linked}),
        "n_tfs_downstream_of_a_metabolite": len({tf for _, tf in linked}),
    }


def _late_integration(root: Path, out: Path, lambda_reg: float, resume: bool) -> pd.DataFrame:
    """Late versus early integration on COSMOS at the joint run's lambda.

    Late: a metabolite-only network (solved here) merged with the stage-04
    TF-only network (cosmos_tf) at the same lambda.  Early: the stage-04
    joint network (cosmos_multi).  Compares topology, cross-layer edges,
    metabolite-to-TF connections and node-state conflicts in the merge.
    """
    metabolites = pd.read_csv(root / STAGE04 / "inputs/measured_metabolites.tsv", sep="\t")
    values = {str(n): float(v) for n, v in zip(metabolites["pkn_node"], metabolites["value"])}
    result = _solve(root, out, "late_integration/metabolites_only", "metabolites only (COSMOS)", _pkn_path(root, "cosmos"),
                    {"PMS_vs_Ctrl": values}, lambda_reg, COSMOS_TIME_LIMIT, resume=resume)
    label = lambda_label(lambda_reg)
    tf_dir, joint_dir = root / STAGE04 / "cosmos_tf" / f"lambda_{label}", root / STAGE04 / "cosmos_multi" / f"lambda_{label}"
    read = lambda d, f: pd.read_csv(d / f, sep="\t", dtype={"source": str, "target": str, "node": str})
    tf_edges, tf_nodes = read(tf_dir, "edges.tsv"), read(tf_dir, "nodes.tsv")
    joint_edges, joint_nodes = read(joint_dir, "edges.tsv"), read(joint_dir, "nodes.tsv")
    met_edges, met_nodes = result["selected_edges"], result["nodes"]
    tfs = set(tf_nodes.loc[tf_nodes["is_measured"], "node"])
    mets = set(values)
    late_edges = pd.concat([tf_edges[["source", "target"]], met_edges[["source", "target"]]]).drop_duplicates()
    tf_state = dict(zip(tf_nodes["node"], tf_nodes["state"]))
    met_state = dict(zip(met_nodes["node"], met_nodes["state"])) if result["accepted"] else {}
    shared = {n for n in set(tf_state) & set(met_state) if tf_state[n] != 0 and met_state[n] != 0}
    conflicts = sorted(n for n in shared if np.sign(tf_state[n]) != np.sign(met_state[n]))
    rows = [
        _network_summary("tf_only", tf_edges, tf_nodes, tfs, mets),
        _network_summary("metabolites_only", met_edges, met_nodes, tfs, mets),
        _network_summary("late_union", late_edges, pd.DataFrame(), tfs, mets),
        _network_summary("early_joint", joint_edges, joint_nodes, tfs, mets),
    ]
    joint_set = _edge_set(joint_edges)
    late_set = {(str(s), str(t)) for s, t in zip(late_edges["source"], late_edges["target"])}
    single = _edge_set(tf_edges) | _edge_set(met_edges)
    for row in rows:
        row.update({"lambda_reg": lambda_reg, "metabolite_solve_status": result["status"],
                    "metabolite_solve_gap": result["mip_gap"], "metabolite_solve_error": result["weighted_fitting_error"]})
    rows[2].update({"n_shared_active_nodes": len(shared), "n_state_conflicts": len(conflicts),
                    "state_conflicts": ",".join(conflicts),
                    "edge_jaccard_with_joint": _jaccard(late_set, joint_set)})
    rows[3].update({"n_edges_in_neither_single_layer_network": len(joint_set - single),
                    "edge_jaccard_with_late": _jaccard(late_set, joint_set)})
    return pd.DataFrame(rows)


def _context_cosmos(root: Path, out: Path, resume: bool) -> pd.DataFrame:
    """Refit the stage-04 joint and stage-05 SV controls on the iMAT PKN."""
    pkn_path = root / "results/ionescu_corneto/06_validation/metabolic_context/context_cosmos/pkn.tsv"
    if not pkn_path.exists():
        return pd.DataFrame([{"run": "context_cosmos", "status": "missing_context_pkn"}])
    rows = []
    measurements = _stage04_measurements(root, "cosmos_multi")
    joint = _solve(root, out, "context_cosmos/joint", "context COSMOS joint", pkn_path, measurements,
                   0.1, COSMOS_TIME_LIMIT, resume=resume)
    rows.append({"run": "cosmos_multi", "scope": "context", "lambda_reg": 0.1,
                 **_fit_columns(joint)})
    node, contrast = _sv_measurement(root)
    # The pruned PKN may contain no directed path from the fixed input to the
    # measured node; CarnivalFlow then prunes the graph to nothing and the worker
    # fails (2026-10-15). Record that as a result instead of a crash.
    pkn = pd.read_csv(pkn_path, sep="\t", dtype={"source": str, "target": str})
    graph = nx.DiGraph(list(zip(pkn["source"], pkn["target"])))
    if not ("HMGCR" in graph and node in graph and nx.has_path(graph, "HMGCR", node)):
        rows.append({"run": "cosmos_multi", "scope": "context_sv", "lambda_reg": SV_LAMBDA,
                     "measured_node": node, "measured_t": float(contrast["moderated_t"].iloc[0]),
                     "status": "no_directed_path_from_HMGCR", "accepted": False})
        return pd.DataFrame(rows)
    sv = _solve(root, out, "context_cosmos/sv_control", "context COSMOS simvastatin", pkn_path,
                {"SV_vs_untreated": {node: float(contrast["moderated_t"].iloc[0])}}, SV_LAMBDA,
                COSMOS_TIME_LIMIT, resume=resume, inputs={"HMGCR": -1.0})
    active = sv["nodes"]
    active = active[active["state"].ne(0)] if sv["accepted"] else active.iloc[0:0]
    rows.append({"run": "cosmos_multi", "scope": "context_sv", "lambda_reg": SV_LAMBDA,
                 "measured_node": node, "measured_t": float(contrast["moderated_t"].iloc[0]),
                 "cholesterol_state": int(active.loc[active["node"].eq(node), "state"].iloc[0]) if node in set(active["node"]) else 0,
                 **_fit_columns(sv)})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stage 05 robustness controls")
    parser.add_argument("--parts", nargs="+", choices=PARTS, default=list(PARTS))
    parser.add_argument("--cosmos", action="store_true",
                        help="also run COSMOS controls, including context-specific refits")
    parser.add_argument("--cosmos-samples", type=int, default=DEFAULT_COSMOS_SAMPLES)
    parser.add_argument("--sv-full", action="store_true",
                        help="also run the SV control on the full COSMOS PKNs (up to 2 x 3,600 s)")
    parser.add_argument("--resume", action="store_true", help="reuse finished solves in the output folder")
    args = parser.parse_args(argv)
    args.parts = [part for part in PARTS if part in set(args.parts)]
    return args


def main(root: Path = ROOT, argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    if not (root / STAGE04 / "04_inference.provenance.json").exists():
        raise FileNotFoundError("Run stage 04 to completion before stage 05")
    out = root / ANALYSIS
    if out.exists() and any(out.iterdir()) and not args.resume:
        raise ValueError(f"{out} is not empty; archive it or rerun with --resume")
    out.mkdir(parents=True, exist_ok=True)
    selected = _stage04_selected(root)
    omnipath_lambda = selected["omnipath"]
    cosmos_lambda = selected.get("cosmos_multi")
    done: list[str] = []

    frequencies = pd.DataFrame(columns=["run", "source", "target", "n_present", "n_networks", "frequency"])
    loo_tables: list[pd.DataFrame] = []
    if "sampling" in args.parts:
        runs = ["omnipath", "omnipath_lines"] + (["cosmos_multi"] if args.cosmos else [])
        samples = {run: N_SAMPLES if run not in COSMOS_RUNS else args.cosmos_samples for run in runs}
        status, frequencies = _sampling(root, out, runs, samples, args.resume)
        _write(status, out, "sampling_status.tsv")
        _write(frequencies, out, "edge_frequencies.tsv")
        done.append("sampling")
    if "loo" in args.parts:
        tables = []
        for run in ["omnipath"] + (["cosmos_multi"] if args.cosmos and cosmos_lambda is not None else []):
            status, edges = _loo(root, out, run, selected[run], args.resume)
            tables.append(status)
            loo_tables += edges
        _write(pd.concat(tables, ignore_index=True), out, "loo_status.tsv")
        done.append("loo")
    if "grid" in args.parts:
        _write(_choice_grid(root, out, omnipath_lambda, args.resume), out, "choice_grid.tsv")
        done.append("grid")
    if "permutations" in args.parts:
        table, summary = _permutations(root, out, omnipath_lambda, args.resume)
        _write(table, out, "permutations.tsv")
        _write(summary, out, "permutation_rank.tsv")
        done.append("permutations")
    if "rewired" in args.parts:
        tables = [_rewired(root, out, "omnipath", "omnipath", omnipath_lambda, N_OMNIPATH_REWIRED, args.resume)]
        if args.cosmos and cosmos_lambda is not None:
            tables.append(_rewired(root, out, "cosmos", "cosmos_multi", cosmos_lambda, N_COSMOS_REWIRED, args.resume))
        _write(pd.concat(tables, ignore_index=True), out, "rewired.tsv")
        done.append("rewired")
    if "receptors" in args.parts:
        try:
            _write(_receptor_inputs(root, out, omnipath_lambda, args.resume), out, "receptor_inputs.tsv")
            done.append("receptors")
        except Exception as error:  # network or resource-format failure must not stop the other controls
            print(f"receptors: skipped ({type(error).__name__}: {error})", flush=True)
            _write(pd.DataFrame([{"run": "omnipath", "status": f"failed: {type(error).__name__}: {error}"}]),
                   out, "receptor_inputs.tsv")
    if "sv_control" in args.parts:
        table, contrast = _sv_control(root, out, args.resume, full=args.sv_full)
        _write(contrast, out, "sv_cholesterol_contrast.tsv")
        _write(table, out, "sv_control.tsv")
        done.append("sv_control")
    if "tight_gap" in args.parts and args.cosmos:
        _write(_tight_gap(root, out, args.resume), out, "tight_gap.tsv")
        done.append("tight_gap")
    if "late_integration" in args.parts and args.cosmos and cosmos_lambda is not None:
        _write(_late_integration(root, out, cosmos_lambda, args.resume), out, "late_integration.tsv")
        done.append("late_integration")
    if "no_drug_inputs" in args.parts and args.cosmos and cosmos_lambda is not None:
        _write(_no_drug_inputs(root, out, cosmos_lambda, args.resume), out, "no_drug_inputs.tsv")
        done.append("no_drug_inputs")
    if "context_cosmos" in args.parts and args.cosmos:
        _write(_context_cosmos(root, out, args.resume), out, "context_cosmos.tsv")
        done.append("context_cosmos")
    if "sampling" in args.parts:
        _write(_hub_bias(root, frequencies, out), out, "hub_bias.tsv")
        if "loo" in args.parts:
            _write(_core_network(frequencies, loo_tables, out), out, "core_network.tsv")

    inputs = [
        root / STAGE04 / "04_inference.provenance.json",
        root / STAGE04 / "inputs/measured_tfs.tsv",
        root / STAGE04 / "inputs/measured_metabolites.tsv",
        root / STAGE04 / "summary/selected_lambda.tsv",
        root / STAGE04 / "summary/lambda_grid.tsv",
        root / PREPROCESSED / "rnaseq/log2cpm.tsv",
        root / PREPROCESSED / "rnaseq/samples.tsv",
        root / PREPROCESSED / "metabolomics/ST003328/line_log2.tsv",
        root / ACTIVITIES / "collectri/tf_activity_per_line.tsv",
        root / PLAN, root / ADDENDUM,
        root / "doc/decisions/2026-10-13-imat-fba-plan.md",
    ] + [_pkn_path(root, name) for name in ("omnipath", "nc_rule", "signor_only", "unpruned", "cosmos", "cosmos_repair")]
    if "receptors" in done:
        inputs.append(root / RECEPTOR_CACHE)
    if "sv_control" in done:
        inputs.append(root / PKN / "aim5_route/cofactor_nodes.tsv")
    code = [root / "src" / name for name in ("inference.py", "corneto_worker.py", "activities.py", "pkn.py",
                                             "metabolic_route.py", "gem.py", "gem_worker.py")] + [
        root / "scripts/analysis/06_corneto_robustness.py", root / "scripts/analysis/08_metabolic_context.py"]
    parameters = {
        "parts_completed_this_invocation": done,
        "cosmos": args.cosmos,
        "alternatives": {"omnipath": N_SAMPLES, "omnipath_lines": N_SAMPLES, "cosmos_multi": args.cosmos_samples,
                         "percentage": 0.03, "scale": 0.03, "rel_opt_tol": 0.05},
        "rewired_replicates": {"omnipath": N_OMNIPATH_REWIRED, "cosmos": N_COSMOS_REWIRED},
        "core_threshold": CORE_THRESHOLD,
        "tight_gap": {"mip_rel_gap": TIGHT_GAP, "lambdas": list(TIGHT_LAMBDAS)},
        "sv_control": {"full_pkn": args.sv_full, "input": "HMGCR=-1", "lambda": SV_LAMBDA, "measurement": "paired moderated t, SV - untreated",
                       "pkns": "cosmos/cosmos_repair x full/metabolism-only cofactor-free"},
        "receptor_resource": RECEPTOR_URL,
        "receptor_filter": f"composite, transmembrane protein, consensus_score >= {RECEPTOR_MIN_CONSENSUS}",
        "seed_base": SEED_BASE, "solver": SOLVER, "solver_relative_gap": SOLVER_RELATIVE_GAP,
        "solver_threads": SOLVER_THREADS, "solver_seed": SOLVER_SEED,
        "execution": "one Python process per solve (src/corneto_worker.py); resume=" + str(args.resume),
    }
    existing = [path for path in inputs if path.exists()]
    provenance = stage_provenance(root, "05_robustness", existing, {}, parameters, code)
    write_stage_provenance(out, provenance)
    print(f"Stage 05 parts completed: {', '.join(done) or 'none'} -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(argv=sys.argv[1:]))
