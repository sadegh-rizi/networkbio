"""CORNETO/CarnivalFlow data construction, solving and readouts for stage 04.

The functions here keep the scientific choices visible: every source node is
an unknown-sign input unless it is explicitly excluded, measured values are
passed as signed values without scaling, and graph pruning is left to
CORNETO's own input-to-measurement reachability preprocessing.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import corneto as cn
import numpy as np
import pandas as pd
from corneto.io import load_graph_from_sif_tuples
from corneto.methods import CarnivalFlow
from corneto.methods.sampler import sample_alternative_solutions

Edge = tuple[str, int, str]
Measurements = Mapping[str, Mapping[str, float]]

LAMBDA_GRID = (0.0, 0.01, 0.03, 0.1, 0.3, 1.0)
SOLVER = "highs"
SOLVER_SEED = 20261010
SOLVER_THREADS = 4
SOLVER_RELATIVE_GAP = 0.01
ACCEPTANCE_GAP = 0.05
OMNIPATH_TIME_LIMIT = 1800
COSMOS_TIME_LIMIT = 3600
EXCLUDED_INPUTS = frozenset({
    "HMGCR", "SREBF1", "SREBF2", "SCAP", "INSIG1", "E2F1", "EGR1",
    "WT1", "SP1", "JUN", "JUNB",
})
MEVALONATE_HMDB = "HMDB0000227"
CHOLESTEROL_HMDB = "HMDB0000067"
# Pre-registered canonical cholesterol-biosynthesis genes
# (doc/decisions/2026-10-10-aim5-metabolic-route-plan.md, method step 1).
CANONICAL_CHOLESTEROL_GENES = frozenset({
    "HMGCR", "MVK", "PMVK", "MVD", "IDI1", "IDI2", "FDPS", "FDFT1", "SQLE",
    "LSS", "CYP51A1", "TM7SF2", "LBR", "MSMO1", "NSDHL", "HSD17B7", "EBP",
    "SC5D", "DHCR7", "DHCR24",
})


def read_pkn(path: Path) -> list[Edge]:
    """Read a signed SIF-style stage-03 table and reject duplicate edges."""
    table = pd.read_csv(path, sep="\t", dtype={"source": str, "target": str})
    if not {"source", "sign", "target"}.issubset(table.columns):
        raise ValueError(f"PKN lacks source/sign/target columns: {path}")
    if table[["source", "sign", "target"]].isna().any().any():
        raise ValueError(f"PKN contains missing edge fields: {path}")
    edges = [(str(source), int(sign), str(target)) for source, sign, target in
             table[["source", "sign", "target"]].itertuples(index=False, name=None)]
    if any(sign not in (-1, 1) for _, sign, _ in edges):
        raise ValueError("PKN signs must be -1 or +1")
    if len({(source, target) for source, _, target in edges}) != len(edges):
        raise ValueError(f"PKN has duplicate source-target pairs: {path}")
    return edges


def pkn_nodes(edges: Iterable[Edge]) -> set[str]:
    """Return the vertices of a signed directed edge list."""
    return {node for source, _, target in edges for node in (source, target)}


def source_nodes(edges: Iterable[Edge]) -> list[str]:
    """Return all zero-in-degree nodes before CORNETO preprocessing."""
    edges = list(edges)
    nodes = pkn_nodes(edges)
    targets = {target for _, _, target in edges}
    return sorted(nodes - targets)


def is_excluded_input(node: str, excluded: Iterable[str] = EXCLUDED_INPUTS) -> bool:
    """Apply the fixed protein and mevalonate/cholesterol input exclusions."""
    excluded = set(excluded)
    if node in excluded:
        return True
    return any(identifier in node for identifier in (MEVALONATE_HMDB, CHOLESTEROL_HMDB))


def make_corneto_data(
    edges: Iterable[Edge],
    measurements: Measurements,
    *,
    excluded_inputs: Iterable[str] = EXCLUDED_INPUTS,
    inputs: Mapping[str, float] | None = None,
) -> tuple[cn.Data, list[str], dict[str, list[str]]]:
    """Build CORNETO data and return the inputs used and missing measurements.

    A measured source node is kept as an output rather than assigned both
    roles, so it is not returned as an input.  The returned input list is the
    union of the per-sample input sets actually written to the CORNETO data.
    Missing measured nodes are reported per sample and never silently passed
    to CORNETO.

    ``inputs`` replaces the default input set (every non-excluded PKN source
    node with value 0, i.e. unknown sign).  It maps node -> value, where 0 is
    an unknown-sign input and +1/-1 a known perturbation.  Explicit inputs may
    be internal nodes and are not filtered by ``excluded_inputs``; nodes
    absent from the PKN are dropped (the caller compares with the returned
    input list).  Used by the stage-05 positive control and receptor-only
    sensitivity.
    """
    edges = list(edges)
    nodes = pkn_nodes(edges)
    if inputs is None:
        input_values = {
            node: 0.0 for node in source_nodes(edges) if not is_excluded_input(node, excluded_inputs)
        }
    else:
        input_values = {str(node): float(value) for node, value in inputs.items() if str(node) in nodes}
    data: dict[str, dict[str, dict[str, Any]]] = {}
    missing: dict[str, list[str]] = {}
    used_inputs: set[str] = set()
    for sample, values in measurements.items():
        sample_values = {str(node): float(value) for node, value in values.items()}
        if not all(np.isfinite(value) for value in sample_values.values()):
            raise ValueError(f"Measurements for {sample} contain non-finite values")
        absent = sorted(set(sample_values) - nodes)
        missing[str(sample)] = absent
        if absent:
            sample_values = {node: value for node, value in sample_values.items() if node in nodes}
        features: dict[str, dict[str, Any]] = {
            node: {"value": value, "role": "input", "mapping": "vertex"}
            for node, value in input_values.items() if node not in sample_values
        }
        used_inputs.update(features)
        features.update({
            node: {"value": value, "role": "output", "mapping": "vertex"}
            for node, value in sample_values.items()
        })
        data[str(sample)] = features
    if not data:
        raise ValueError("At least one CORNETO sample is required")
    return cn.Data.from_cdict(data), sorted(used_inputs), missing


def _objective_value(value: Any) -> float:
    array = np.asarray(value)
    if array.size != 1:
        raise ValueError(f"Expected scalar objective value, got shape {array.shape}")
    return float(array.item())


def _mip_gap(problem: Any) -> float:
    extra = getattr(getattr(problem, "solver_stats", None), "extra_stats", None)
    value = getattr(extra, "mip_gap", np.nan)
    if isinstance(extra, Mapping):
        value = extra.get("mip_gap", np.nan)
    try:
        value = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return value if np.isfinite(value) else float("nan")


def solution_accepted(status: str, mip_gap: float, acceptance_gap: float = ACCEPTANCE_GAP) -> bool:
    """Apply the preregistered optimality/time-limit acceptance rule."""
    return status == "optimal" or (
        status == "user_limit" and np.isfinite(mip_gap) and mip_gap <= acceptance_gap
    )


def _selected_edges(
    method: CarnivalFlow,
    problem: Any,
    measurements: Measurements,
    pkn_sign: Mapping[tuple[str, str], int],
) -> pd.DataFrame:
    """Extract selected non-boundary edges for every CORNETO sample."""
    values = getattr(problem.expr, "edge_value", None)
    if values is None or values.value is None:
        return pd.DataFrame(columns=["sample", "source", "target", "edge_value", "pkn_sign"])
    edge_values = np.asarray(values.value).reshape(method.processed_graph.num_edges, -1)
    samples = list(measurements)
    rows: list[dict[str, Any]] = []
    for edge_index, (sources, targets) in enumerate(method.processed_graph.E):
        if len(sources) == 0 or len(targets) == 0:
            continue
        source, target = str(next(iter(sources))), str(next(iter(targets)))
        if (source, target) not in pkn_sign:
            raise ValueError(f"CORNETO returned an edge outside the PKN: {source}->{target}")
        for sample_index, sample in enumerate(samples):
            value = float(edge_values[edge_index, sample_index])
            if np.isfinite(value) and abs(value) > 0.5:
                rows.append({
                    "sample": sample,
                    "source": source,
                    "target": target,
                    "edge_value": int(round(value)),
                    "pkn_sign": pkn_sign[(source, target)],
                })
    return pd.DataFrame(rows, columns=["sample", "source", "target", "edge_value", "pkn_sign"])


def vertex_states(method: CarnivalFlow, problem: Any, samples: Iterable[str]) -> dict[str, dict[str, int]]:
    """Read node states per sample from CORNETO's ``vertex_value`` variable.

    This is the variable the CarnivalFlow error term uses, so fit readouts
    computed from it match the optimised objective.  Vertices outside the
    processed (reachability-pruned) graph have no state and are absent.
    """
    samples = [str(sample) for sample in samples]
    variable = getattr(problem.expr, "vertex_value", None)
    if variable is None or variable.value is None:
        return {sample: {} for sample in samples}
    vertices = [str(vertex) for vertex in method.processed_graph.V]
    values = np.asarray(variable.value, dtype=float).reshape(len(vertices), -1)
    if values.shape[1] != len(samples):
        raise ValueError(f"vertex_value has {values.shape[1]} columns for {len(samples)} samples")
    return {
        sample: {vertex: int(np.rint(values[i, j])) for i, vertex in enumerate(vertices)}
        for j, sample in enumerate(samples)
    }


def _edge_state_conflicts(selected: pd.DataFrame, states: Mapping[str, int]) -> int:
    """Count selected edges whose signal disagrees with either endpoint state."""
    conflicts = 0
    for row in selected.itertuples(index=False):
        value = int(row.edge_value)
        source_expected = states.get(str(row.source), 0) * int(row.pkn_sign)
        target_expected = states.get(str(row.target), 0)
        if value != source_expected or value != target_expected:
            conflicts += 1
    return conflicts


def _fit_table(
    selected_edges: pd.DataFrame,
    states: Mapping[str, Mapping[str, int]],
    measurements: Measurements,
    missing: Mapping[str, list[str]],
    objectives: Mapping[str, float],
) -> pd.DataFrame:
    """Return per-sample fit and selection readouts from ``vertex_value``.

    Measured nodes absent from the PKN are counted as missing; measured PKN
    nodes that CORNETO pruned (no input reaches them) are counted as
    unreachable.  Neither enters CORNETO's error term, so the sign-agreement
    ratio is computed over reachable measured nodes only.
    """
    rows = []
    for sample, measured in measurements.items():
        sample = str(sample)
        sample_states = states.get(sample, {})
        sample_edges = selected_edges[selected_edges["sample"].eq(sample)]
        absent = set(missing.get(sample, []))
        reachable = [node for node in measured if node not in absent and node in sample_states]
        unreachable = [node for node in measured if node not in absent and node not in sample_states]
        matches = [sample_states[node] == int(np.sign(measured[node])) for node in reachable]
        sample_error = sum(
            value for name, value in objectives.items()
            if name.startswith(f"error_{sample}_")
        )
        rows.append({
            "sample": sample,
            "weighted_fitting_error": sample_error,
            "n_measurements": len(measured),
            "n_missing_from_pkn": len(absent),
            "n_unreachable": len(unreachable),
            "n_reachable": len(reachable),
            "n_measurements_correct_sign": int(sum(matches)),
            "sign_agreement_ratio": float(np.mean(matches)) if matches else np.nan,
            "n_selected_edges": len(sample_edges),
            "n_selected_nodes": sum(1 for state in sample_states.values() if state != 0),
            "n_edge_state_conflicts": _edge_state_conflicts(sample_edges, sample_states),
        })
    return pd.DataFrame(rows)


def node_table(
    states: Mapping[str, Mapping[str, int]],
    measurements: Measurements,
    missing: Mapping[str, list[str]],
    inputs: Iterable[str],
) -> pd.DataFrame:
    """Active nodes (state != 0) plus every measured node, per sample."""
    inputs = set(inputs)
    rows = []
    for sample, measured in measurements.items():
        sample = str(sample)
        sample_states = states.get(sample, {})
        absent = set(missing.get(sample, []))
        listed = {node for node, state in sample_states.items() if state != 0} | set(measured)
        for node in sorted(listed):
            if node in absent:
                reach = "missing_from_pkn"
            elif node in sample_states:
                reach = "reachable"
            else:
                reach = "unreachable"
            rows.append({
                "sample": sample,
                "node": node,
                "state": int(sample_states.get(node, 0)),
                "is_measured": node in measured,
                "observed_value": measured.get(node, np.nan),
                "is_input": node in inputs,
                "reachability": reach,
            })
    return pd.DataFrame(
        rows,
        columns=["sample", "node", "state", "is_measured", "observed_value", "is_input", "reachability"],
    )


def solve_carnival(
    edges: Iterable[Edge],
    measurements: Measurements,
    lambda_reg: float,
    *,
    time_limit: int,
    solver: str = SOLVER,
    solver_relative_gap: float = SOLVER_RELATIVE_GAP,
    acceptance_gap: float = ACCEPTANCE_GAP,
    threads: int = SOLVER_THREADS,
    seed: int = SOLVER_SEED,
    excluded_inputs: Iterable[str] = EXCLUDED_INPUTS,
    inputs: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    """Build and solve one CarnivalFlow problem with the preregistered rule.

    ``inputs`` optionally replaces the default input set; see
    :func:`make_corneto_data`.
    """
    edges = list(edges)
    data, inputs, missing = make_corneto_data(
        edges, measurements, excluded_inputs=excluded_inputs, inputs=inputs
    )
    method = CarnivalFlow(lambda_reg=lambda_reg)
    problem = method.build_from_data(load_graph_from_sif_tuples(edges), data)
    cvx_problem = problem.solve(
        solver=solver,
        max_seconds=time_limit,
        verbosity=0,
        mip_rel_gap=solver_relative_gap,
        threads=threads,
        random_seed=seed,
    )
    status = str(getattr(cvx_problem, "status", "unknown"))
    mip_gap = _mip_gap(cvx_problem)
    edge_value = getattr(problem.expr, "edge_value", None)
    accepted = (
        solution_accepted(status, mip_gap, acceptance_gap)
        and edge_value is not None
        and edge_value.value is not None
    )
    objectives = {
        str(objective.name): _objective_value(objective.value)
        for objective in problem.objectives
        if objective.value is not None
    }
    signs = {(source, target): sign for source, sign, target in edges}
    if accepted:
        selected = _selected_edges(method, problem, measurements, signs)
        states = vertex_states(method, problem, measurements)
        fit = _fit_table(selected, states, measurements, missing, objectives)
        nodes = node_table(states, measurements, missing, inputs)
    else:
        selected = pd.DataFrame(columns=["sample", "source", "target", "edge_value", "pkn_sign"])
        states = {str(sample): {} for sample in measurements}
        fit = pd.DataFrame()
        nodes = node_table(states, measurements, missing, inputs).iloc[0:0]
    return {
        "lambda_reg": float(lambda_reg),
        "status": status,
        "accepted": accepted,
        "mip_gap": mip_gap,
        "objectives": objectives,
        "weighted_fitting_error": sum(value for name, value in objectives.items() if name.startswith("error_")),
        "selected_edges": selected,
        "states": states,
        "nodes": nodes,
        "fit": fit,
        "inputs": inputs,
        "missing_measurements": missing,
        "n_missing_measurements": sum(len(nodes) for nodes in missing.values()),
        "preprocessed_node_count": int(method.processed_graph.num_vertices),
        "preprocessed_edge_count": int(method.processed_graph.num_edges),
        "method": method,
        "problem": problem,
        "cvx_problem": cvx_problem,
        "solve_seconds": float(getattr(getattr(cvx_problem, "solver_stats", None), "solve_time", np.nan)),
    }


LAMBDA_ABSOLUTE_TOLERANCE = 1e-9


def _within_tolerance(error: float, baseline: float, tolerance: float) -> bool:
    """Registered rule error <= tolerance * baseline, with a 1e-9 allowance for solver rounding.

    Without it, a baseline of exactly 0 rejects errors such as 1.4e-15 that are
    floating-point noise, not misfit (stage-04 run of 2026-10-11, OmniPath lambda = 1).
    """
    return error <= tolerance * baseline + LAMBDA_ABSOLUTE_TOLERANCE


def choose_lambda(results: Iterable[Mapping[str, Any]], tolerance: float = 1.10) -> float | None:
    """Choose the largest accepted lambda within 10% of the lambda-zero fit."""
    results = list(results)
    baseline = next((result for result in results if result["lambda_reg"] == 0.0 and result["accepted"]), None)
    if baseline is None:
        return None
    baseline_error = float(baseline["weighted_fitting_error"])
    accepted = [
        result for result in results
        if result["accepted"] and _within_tolerance(float(result["weighted_fitting_error"]), baseline_error, tolerance)
    ]
    if not accepted:
        return None
    return max(float(result["lambda_reg"]) for result in accepted)


def lambda_label(value: float) -> str:
    """Return a stable directory label for a registered lambda value."""
    return f"{value:g}"


def lambda_grid_table(results: Iterable[Mapping[str, Any]], tolerance: float = 1.10) -> pd.DataFrame:
    """Summarise all lambda points before the selected lambda is reported."""
    results = list(results)
    baseline = next((result for result in results if result["lambda_reg"] == 0.0 and result["accepted"]), None)
    baseline_error = float(baseline["weighted_fitting_error"]) if baseline is not None else np.nan
    rows = []
    for result in results:
        fit = result["fit"]
        rows.append({
            "lambda_reg": result["lambda_reg"],
            "status": result["status"],
            "accepted": result["accepted"],
            "mip_gap": result["mip_gap"],
            "solve_seconds": result["solve_seconds"],
            "weighted_fitting_error": result["weighted_fitting_error"],
            "n_missing_measurements": result["n_missing_measurements"],
            "within_error_tolerance": (
                bool(result["accepted"] and np.isfinite(baseline_error) and
                     _within_tolerance(float(result["weighted_fitting_error"]), baseline_error, tolerance))
            ),
            "n_selected_edges": int(result["selected_edges"].drop_duplicates(["source", "target"]).shape[0]),
            "n_selected_nodes": len({
                node for sample_states in result.get("states", {}).values()
                for node, state in sample_states.items() if state != 0
            }),
            "sign_agreement_ratio": float(fit["sign_agreement_ratio"].mean()) if not fit.empty else np.nan,
            "n_unreachable_measurements": int(fit["n_unreachable"].sum()) if not fit.empty else np.nan,
            "n_edge_state_conflicts": int(fit["n_edge_state_conflicts"].sum()) if not fit.empty else np.nan,
            "preprocessed_nodes": result["preprocessed_node_count"],
            "preprocessed_edges": result["preprocessed_edge_count"],
        })
    return pd.DataFrame(rows)


def sample_edges(
    problem: Any,
    method: CarnivalFlow,
    *,
    time_limit: int,
    solver: str = SOLVER,
    solver_relative_gap: float = SOLVER_RELATIVE_GAP,
    threads: int = SOLVER_THREADS,
    seed: int = SOLVER_SEED,
    max_samples: int = 30,
) -> pd.DataFrame:
    """Sample near-optimal edge solutions using the pinned CORNETO helper."""
    if "edge_value" not in problem.expr:
        raise KeyError("CORNETO problem has no edge_value variable")
    sampled = sample_alternative_solutions(
        problem,
        "edge_value",
        percentage=0.03,
        scale=0.03,
        rel_opt_tol=0.05,
        max_samples=max_samples,
        solver_kwargs={
            "solver": solver,
            "max_seconds": time_limit,
            "mip_rel_gap": solver_relative_gap,
            "threads": threads,
            "random_seed": seed,
        },
        rng=seed,
        collect_vars=["edge_value"],
        verbose=0,
    )
    values = np.asarray(sampled.get("edge_value", []), dtype=float)
    if values.size == 0:
        # No perturbed solution stayed within rel_opt_tol: report zero alternatives.
        return pd.DataFrame(columns=["alternative", "source", "target", "edge_value", "n_samples_selected"])
    values = values.reshape(values.shape[0], method.processed_graph.num_edges, -1)
    rows: list[dict[str, Any]] = []
    for sample_index in range(values.shape[0]):
        for edge_index, (sources, targets) in enumerate(method.processed_graph.E):
            if len(sources) == 0 or len(targets) == 0:
                continue
            # Multi-sample problems (omnipath_lines) have one column per
            # sample; an edge belongs to the alternative if any sample uses it.
            per_sample = values[sample_index, edge_index, :]
            used = np.abs(per_sample) > 0.5
            if not used.any():
                continue
            rows.append({
                "alternative": sample_index,
                "source": str(next(iter(sources))),
                "target": str(next(iter(targets))),
                "edge_value": int(round(float(per_sample[used][0]))),
                "n_samples_selected": int(used.sum()),
            })
    return pd.DataFrame(rows, columns=["alternative", "source", "target", "edge_value", "n_samples_selected"])
