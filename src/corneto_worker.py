"""Solve one CORNETO CarnivalFlow problem per Python process.

Stage 04 was killed twice on 2026-10-11 when several COSMOS problems were
solved in one process: HiGHS/cvxpy memory is not returned to the operating
system between solves.  This module gives stage 05 the same isolation as the
stage-04 worker: the parent writes a ``task.json``, a fresh interpreter runs
:func:`run_task`, writes the tables plus ``result.json`` and exits.

Alternative-solution sampling needs the live CORNETO problem object, so a
sampling task solves and samples inside the same worker process.

A finished task is reused with ``resume=True`` only when the signature stored
in its ``result.json`` matches the new task (PKN checksum, measurements,
inputs, lambda, limits and solver settings).
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

SRC = Path(__file__).resolve().parent
ROOT = SRC.parent
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from inference import (  # noqa: E402
    SOLVER,
    SOLVER_RELATIVE_GAP,
    SOLVER_SEED,
    SOLVER_THREADS,
    read_pkn,
    sample_edges,
    solve_carnival,
)

EDGE_COLUMNS = ["sample", "source", "target", "edge_value", "pkn_sign"]
NODE_COLUMNS = ["sample", "node", "state", "is_measured", "observed_value", "is_input", "reachability"]
ALTERNATIVE_COLUMNS = ["alternative", "source", "target", "edge_value", "n_samples_selected"]
SUMMARY_KEYS = (
    "lambda_reg", "status", "accepted", "mip_gap", "objectives", "weighted_fitting_error",
    "inputs", "missing_measurements", "n_missing_measurements",
    "preprocessed_node_count", "preprocessed_edge_count", "solve_seconds",
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def portable_path(path: Path, root: Path = ROOT) -> str:
    """Repository-relative POSIX path, so signatures match across WSL/Windows."""
    path = Path(path).resolve()
    try:
        return path.relative_to(Path(root).resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def task_signature(task: Mapping[str, Any]) -> str:
    payload = {key: value for key, value in task.items() if key != "signature"}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode("utf-8")).hexdigest()


def make_task(
    *,
    label: str,
    pkn_path: Path,
    measurements: Mapping[str, Mapping[str, float]],
    lambda_reg: float,
    time_limit: int,
    solver_relative_gap: float = SOLVER_RELATIVE_GAP,
    inputs: Mapping[str, float] | None = None,
    max_samples: int = 0,
    root: Path = ROOT,
) -> dict[str, Any]:
    """Describe one solve (and optional sampling) completely and reproducibly."""
    task = {
        "label": label,
        "pkn_path": portable_path(pkn_path, root),
        "pkn_sha256": file_sha256(pkn_path),
        "measurements": {
            str(sample): {str(node): float(value) for node, value in values.items()}
            for sample, values in measurements.items()
        },
        "inputs": None if inputs is None else {str(node): float(value) for node, value in inputs.items()},
        "lambda_reg": float(lambda_reg),
        "time_limit": int(time_limit),
        "solver_relative_gap": float(solver_relative_gap),
        "max_samples": int(max_samples),
        "solver": [SOLVER, SOLVER_THREADS, SOLVER_SEED],
    }
    task["signature"] = task_signature(task)
    return task


def run_task(task_path: Path, root: Path = ROOT) -> int:
    """Worker entry point: solve (and sample) one task, write its tables."""
    task_path = Path(task_path)
    task = json.loads(task_path.read_text())
    pkn_path = Path(root) / task["pkn_path"]
    if file_sha256(pkn_path) != task["pkn_sha256"]:
        raise ValueError(f"PKN changed since the task was written: {pkn_path}")
    edges = read_pkn(pkn_path)
    result = solve_carnival(
        edges,
        task["measurements"],
        float(task["lambda_reg"]),
        time_limit=int(task["time_limit"]),
        solver=SOLVER,
        solver_relative_gap=float(task["solver_relative_gap"]),
        threads=SOLVER_THREADS,
        seed=SOLVER_SEED,
        inputs=task["inputs"],
    )
    directory = task_path.parent
    # Tables are written before sampling: sampling perturbs the problem.
    result["selected_edges"].to_csv(directory / "edges.tsv", sep="\t", index=False)
    result["nodes"].to_csv(directory / "nodes.tsv", sep="\t", index=False)
    result["fit"].to_csv(directory / "fit.tsv", sep="\t", index=False)
    summary = {key: result[key] for key in SUMMARY_KEYS}
    summary["accepted"] = bool(summary["accepted"])
    summary["requested_inputs_absent_from_pkn"] = sorted(
        set(task["inputs"] or {}) - set(result["inputs"]) - {
            node for values in task["measurements"].values() for node in values
        }
    )
    summary["n_alternatives"] = 0
    summary["sampling_status"] = "not_requested"
    if int(task["max_samples"]) > 0:
        if not result["accepted"]:
            summary["sampling_status"] = "skipped_unaccepted_solution"
        else:
            try:
                alternatives = sample_edges(
                    result["problem"], result["method"],
                    time_limit=int(task["time_limit"]),
                    solver_relative_gap=float(task["solver_relative_gap"]),
                    max_samples=int(task["max_samples"]),
                )
                alternatives.to_csv(directory / "alternatives.tsv", sep="\t", index=False)
                summary["n_alternatives"] = int(alternatives["alternative"].nunique()) if not alternatives.empty else 0
                summary["sampling_status"] = "sampled"
            except Exception as error:  # recorded, not hidden: the solve itself stays usable
                summary["sampling_status"] = f"sampling_failed: {type(error).__name__}: {error}"
    summary["task_signature"] = task["signature"]
    (directory / "result.json").write_text(json.dumps(summary, default=float))
    print(f"  status={result['status']} accepted={result['accepted']} gap={result['mip_gap']:.4g} "
          f"seconds={result['solve_seconds']:.0f} error={result['weighted_fitting_error']:.3f} "
          f"alternatives={summary['n_alternatives']}", flush=True)
    return 0


def read_table(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t", dtype={"sample": str, "source": str, "target": str, "node": str})
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return pd.DataFrame(columns=columns or [])


def load_result(directory: Path) -> dict[str, Any]:
    """Rebuild a result dict from the tables a worker wrote."""
    directory = Path(directory)
    summary = json.loads((directory / "result.json").read_text())
    nodes = read_table(directory / "nodes.tsv", NODE_COLUMNS)
    states: dict[str, dict[str, int]] = {}
    for row in nodes.itertuples(index=False):
        if row.reachability == "reachable" or int(row.state) != 0:
            states.setdefault(str(row.sample), {})[str(row.node)] = int(row.state)
    summary.update({
        "selected_edges": read_table(directory / "edges.tsv", EDGE_COLUMNS),
        "nodes": nodes,
        "fit": read_table(directory / "fit.tsv") if summary["accepted"] else pd.DataFrame(),
        "states": states,
        "alternatives": read_table(directory / "alternatives.tsv", ALTERNATIVE_COLUMNS),
        "directory": directory,
    })
    return summary


def failed_result(task: Mapping[str, Any], status: str, directory: Path | None = None) -> dict[str, Any]:
    return {
        "lambda_reg": float(task["lambda_reg"]), "status": status, "accepted": False,
        "mip_gap": float("nan"), "objectives": {}, "weighted_fitting_error": float("nan"),
        "inputs": [], "missing_measurements": {}, "n_missing_measurements": 0,
        "preprocessed_node_count": 0, "preprocessed_edge_count": 0, "solve_seconds": float("nan"),
        "n_alternatives": 0, "sampling_status": "not_run",
        "selected_edges": pd.DataFrame(columns=EDGE_COLUMNS),
        "nodes": pd.DataFrame(columns=NODE_COLUMNS), "fit": pd.DataFrame(), "states": {},
        "alternatives": pd.DataFrame(columns=ALTERNATIVE_COLUMNS), "directory": directory,
    }


def solve_in_subprocess(
    task: Mapping[str, Any], directory: Path, *, resume: bool = False, root: Path = ROOT
) -> dict[str, Any]:
    """Run ``task`` in a fresh interpreter (or reuse a matching finished result).

    ``root`` must be the root passed to :func:`make_task`; relative PKN paths
    in the task are resolved against it.
    """
    directory = Path(directory)
    done = directory / "result.json"
    if resume and done.exists() and json.loads(done.read_text()).get("task_signature") == task["signature"]:
        print(f"{task['label']}: resumed", flush=True)
        return load_result(directory)
    directory.mkdir(parents=True, exist_ok=True)
    done.unlink(missing_ok=True)
    (directory / "alternatives.tsv").unlink(missing_ok=True)
    task_path = directory / "task.json"
    task_path.write_text(json.dumps(task))
    print(f"{task['label']}: lambda={task['lambda_reg']:g}", flush=True)
    completed = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), str(task_path), str(Path(root).resolve())], check=False
    )
    if completed.returncode != 0 or not done.exists():
        status = f"worker_failed_exit_{completed.returncode}"
        print(f"  {status}", flush=True)
        return failed_result(task, status, directory)
    return load_result(directory)


if __name__ == "__main__":
    raise SystemExit(run_task(Path(sys.argv[1]), Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT))
