"""Run one isolated Recon3D FBA or iMAT solve.

The full Recon3D model is loaded and solved in a fresh Python process.  This
matches the stage-04 worker pattern and prevents successive MILP solves from
retaining CVXPY/HiGHS memory in the parent process.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parent
ROOT = SRC.parent
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gem import add_cholesterol_demand, load_recon3d, reaction_table  # noqa: E402

SOLVER = "highs"
SOLVER_THREADS = 4
SOLVER_SEED = 20261010
SOLVER_RELATIVE_GAP = 0.01
ACCEPTANCE_GAP = 0.05


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def portable_path(path: Path, root: Path = ROOT) -> str:
    path = Path(path).resolve()
    try:
        return path.relative_to(Path(root).resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def task_signature(task: Mapping[str, Any]) -> str:
    payload = {key: value for key, value in task.items() if key != "signature"}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode("utf-8")).hexdigest()


def _json_bounds(bounds: Mapping[str, Mapping[str, tuple[float, float]]]) -> dict[str, dict[str, list[float]]]:
    return {
        str(condition): {str(reaction): [float(pair[0]), float(pair[1])] for reaction, pair in values.items()}
        for condition, values in bounds.items()
    }


def make_task(
    *,
    label: str,
    model_path: Path,
    hgnc_path: Path,
    method: str,
    gene_expression: Mapping[str, Mapping[str, float]] | None = None,
    objectives: Mapping[str, Mapping[str, float]] | None = None,
    reaction_bounds: Mapping[str, Mapping[str, tuple[float, float]]] | None = None,
    time_limit: int = 3600,
    solver_relative_gap: float = SOLVER_RELATIVE_GAP,
    eps: float = 1e-2,
    high_expression_threshold: float | None = None,
    low_expression_threshold: float | None = None,
    add_demand: bool = False,
    root: Path = ROOT,
) -> dict[str, Any]:
    if method not in {"fba", "imat"}:
        raise ValueError(f"Unknown GEM method: {method}")
    task = {
        "label": label,
        "model_path": portable_path(model_path, root),
        "model_sha256": file_sha256(model_path),
        "hgnc_path": portable_path(hgnc_path, root),
        "hgnc_sha256": file_sha256(hgnc_path),
        "method": method,
        "gene_expression": None if gene_expression is None else {
            str(condition): {str(gene): float(value) for gene, value in values.items()}
            for condition, values in gene_expression.items()
        },
        "objectives": None if objectives is None else {
            str(condition): {str(reaction): float(value) for reaction, value in values.items()}
            for condition, values in objectives.items()
        },
        "reaction_bounds": None if reaction_bounds is None else _json_bounds(reaction_bounds),
        "time_limit": int(time_limit),
        "solver_relative_gap": float(solver_relative_gap),
        "eps": float(eps),
        "high_expression_threshold": high_expression_threshold,
        "low_expression_threshold": low_expression_threshold,
        "add_cholesterol_demand": bool(add_demand),
        "solver": [SOLVER, SOLVER_THREADS, SOLVER_SEED],
    }
    task["signature"] = task_signature(task)
    return task


def _bounds(task: Mapping[str, Any]) -> dict[str, dict[str, tuple[float, float]]]:
    return {
        str(condition): {str(reaction): (float(pair[0]), float(pair[1])) for reaction, pair in values.items()}
        for condition, values in (task.get("reaction_bounds") or {}).items()
    }


def _empty_fluxes() -> pd.DataFrame:
    return pd.DataFrame(columns=["condition", "reaction_number", "reaction", "flux"])


def _solve_task(task: Mapping[str, Any], root: Path) -> tuple[dict[str, Any], pd.DataFrame]:
    model_path = root / task["model_path"]
    hgnc_path = root / task["hgnc_path"]
    if file_sha256(model_path) != task["model_sha256"]:
        raise ValueError(f"Recon3D changed since the task was written: {model_path}")
    if file_sha256(hgnc_path) != task["hgnc_sha256"]:
        raise ValueError(f"HGNC resource changed since the task was written: {hgnc_path}")
    model, _ = load_recon3d(model_path, hgnc_path)
    if task.get("add_cholesterol_demand", False):
        add_cholesterol_demand(model)
    conditions = set((task.get("objectives") or {}).keys()) | set((task.get("reaction_bounds") or {}).keys()) | set((task.get("gene_expression") or {}).keys())
    if not conditions:
        raise ValueError("A GEM task needs at least one named condition")
    ordered_conditions = list((task.get("objectives") or {}).keys()) or list(conditions)
    ordered_conditions += [condition for condition in conditions if condition not in ordered_conditions]
    bounds = _bounds(task)
    objectives = task.get("objectives") or {condition: {} for condition in ordered_conditions}
    for condition in ordered_conditions:
        bounds.setdefault(condition, {})
        objectives.setdefault(condition, {})

    if task["method"] == "imat":
        from corneto.methods import MultiSampleIMAT

        method = MultiSampleIMAT(
            lambda_reg=0.0,
            eps=float(task["eps"]),
            high_expression_threshold=task.get("high_expression_threshold"),
            low_expression_threshold=task.get("low_expression_threshold"),
            use_bigm_constraints=True,
            use_mean_for_missing_reactions=False,
        )
        problem = method.build_many(
            model,
            gene_expression=task["gene_expression"],
            reaction_bounds=bounds,
            objectives=objectives,
        )
    else:
        from corneto.methods import MultiSampleFBA

        method = MultiSampleFBA(lambda_reg=0.0)
        problem = method.build_many(model, objectives=objectives, reaction_bounds=bounds)

    started = time.perf_counter()
    solved = problem.solve(
        solver=SOLVER,
        max_seconds=int(task["time_limit"]),
        mip_rel_gap=float(task["solver_relative_gap"]),
        threads=SOLVER_THREADS,
        random_seed=SOLVER_SEED,
    )
    elapsed = time.perf_counter() - started
    stats = getattr(getattr(solved, "solver_stats", None), "extra_stats", None)
    gap = float(getattr(stats, "mip_gap", 0.0) or 0.0)
    status = str(getattr(solved, "status", "unknown"))
    accepted = status in {"optimal", "optimal_inaccurate", "user_limit"} and np.isfinite(gap) and gap <= ACCEPTANCE_GAP
    flow = problem.expr.flow.value
    fluxes = _empty_fluxes()
    if flow is not None:
        values = np.asarray(flow, dtype=float)
        if values.ndim == 1:
            values = values[:, None]
        rows = []
        table = reaction_table(model)
        for column, condition in enumerate(ordered_conditions):
            for reaction_number, reaction, value in zip(table["reaction_number"], table["reaction"], values[:, column]):
                rows.append({"condition": condition, "reaction_number": int(reaction_number), "reaction": str(reaction), "flux": float(value)})
        fluxes = pd.DataFrame(rows)
    summary = {
        "label": task["label"],
        "method": task["method"],
        "status": status,
        "accepted": bool(accepted),
        "mip_gap": gap,
        "solve_seconds": float(getattr(getattr(solved, "solver_stats", None), "solve_time", elapsed) or elapsed),
        "wall_seconds": float(elapsed),
        "objective_value": float(solved.value) if solved.value is not None else float("nan"),
        "n_conditions": len(ordered_conditions),
        "n_reactions": model.num_edges,
        "n_metabolites": model.num_vertices,
        "task_signature": task["signature"],
    }
    return summary, fluxes


def _write_condition_outputs(directory: Path, task: Mapping[str, Any], summary: Mapping[str, Any], fluxes: pd.DataFrame) -> None:
    for condition in sorted(set(fluxes["condition"]) if not fluxes.empty else set((task.get("objectives") or {}).keys())):
        child = directory / str(condition)
        child.mkdir(parents=True, exist_ok=True)
        child.joinpath("task.json").write_text(json.dumps(task, indent=2))
        child.joinpath("result.json").write_text(json.dumps({**summary, "condition": condition}, indent=2, default=float))
        fluxes.loc[fluxes["condition"].eq(condition)].to_csv(child / "fluxes.tsv", sep="\t", index=False)


def run_task(task_path: Path, root: Path = ROOT) -> int:
    task_path = Path(task_path)
    task = json.loads(task_path.read_text())
    summary, fluxes = _solve_task(task, Path(root))
    directory = task_path.parent
    directory.mkdir(parents=True, exist_ok=True)
    task_path.write_text(json.dumps(task, indent=2))
    fluxes.to_csv(directory / "fluxes.tsv", sep="\t", index=False)
    (directory / "result.json").write_text(json.dumps(summary, indent=2, default=float))
    _write_condition_outputs(directory, task, summary, fluxes)
    print(f"  status={summary['status']} accepted={summary['accepted']} gap={summary['mip_gap']:.4g} seconds={summary['wall_seconds']:.0f}", flush=True)
    return 0


def load_result(directory: Path) -> dict[str, Any]:
    directory = Path(directory)
    result = json.loads((directory / "result.json").read_text())
    result["fluxes"] = pd.read_csv(directory / "fluxes.tsv", sep="\t") if (directory / "fluxes.tsv").exists() else _empty_fluxes()
    result["directory"] = directory
    return result


def failed_result(task: Mapping[str, Any], status: str, directory: Path | None = None) -> dict[str, Any]:
    return {
        "label": task["label"], "method": task["method"], "status": status, "accepted": False,
        "mip_gap": float("nan"), "solve_seconds": float("nan"), "wall_seconds": float("nan"),
        "objective_value": float("nan"), "task_signature": task["signature"],
        "fluxes": _empty_fluxes(), "directory": directory,
    }


def solve_in_subprocess(task: Mapping[str, Any], directory: Path, *, resume: bool = False, root: Path = ROOT) -> dict[str, Any]:
    directory = Path(directory)
    done = directory / "result.json"
    if resume and done.exists() and json.loads(done.read_text()).get("task_signature") == task["signature"]:
        print(f"{task['label']}: resumed", flush=True)
        return load_result(directory)
    directory.mkdir(parents=True, exist_ok=True)
    task_path = directory / "task.json"
    task_path.write_text(json.dumps(task, indent=2))
    print(f"{task['label']}: method={task['method']}", flush=True)
    completed = subprocess.run([sys.executable, str(Path(__file__).resolve()), str(task_path), str(Path(root).resolve())], check=False)
    if completed.returncode != 0 or not done.exists():
        status = f"worker_failed_exit_{completed.returncode}"
        print(f"  {status}", flush=True)
        return failed_result(task, status, directory)
    return load_result(directory)


if __name__ == "__main__":
    raise SystemExit(run_task(Path(sys.argv[1]), Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT))
