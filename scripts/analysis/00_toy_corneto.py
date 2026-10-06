"""Toy CORNETO problem with a known answer (plan step 1).

Checks, before any real data:
1. The solver (HiGHS) runs and returns a solution.
2. CarnivalFlow recovers a planted signed network that contains a
   metabolite node (COSMOS-style name): a signalling + metabolism graph
   is just another signed, directed graph to CORNETO.
3. An input given with unknown sign (value 0) has its sign inferred.
4. The multi-sample problem selects the same edges for two identical
   samples.

Writes results/corneto_toy/04_inference/{toy_edges.tsv,toy_checks.json}
and exits non-zero if any check fails.

Run from the repository root (WSL):
    .venv/bin/python scripts/analysis/00_toy_corneto.py
"""

from __future__ import annotations

import json
import platform
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd

import corneto as cn
from corneto.io import load_graph_from_sif_tuples
from corneto.methods import CarnivalFlow

OUT_DIR = Path("results/corneto_toy/04_inference")
METAB = "Metab__HMDB0000067_c"  # cholesterol, named as in the COSMOS meta-PKN
SOLVER = "highs"
MAX_SECONDS = 60
LAMBDA = 0.1

# Signed prior-knowledge network: (source, sign, target)
PKN = [
    ("R1", 1, "K1"),
    ("R1", 1, "K2"),
    ("K1", 1, "TF1"),
    ("K2", -1, "TF2"),
    ("R1", 1, "ENZ"),
    ("ENZ", 1, METAB),
    (METAB, 1, "TF3"),
    # decoys a correct solution must not use
    ("R2", 1, "K3"),  # R2 is not an input
    ("K3", 1, "TF1"),
    ("K1", 1, "TF2"),  # wrong sign for TF2 = -1
    ("K2", 1, "K4"),
    ("K4", 1, "TF4"),  # TF4 is not measured
]
EXPECTED = {
    ("R1", "K1"),
    ("R1", "K2"),
    ("K1", "TF1"),
    ("K2", "TF2"),
    ("R1", "ENZ"),
    ("ENZ", METAB),
    (METAB, "TF3"),
}
MEASURED = {"TF1": 1, "TF2": -1, "TF3": 1, METAB: 1}


def make_sample(input_value: float) -> dict:
    """One sample: R1 is the input, the TFs and the metabolite are measured."""
    feats = {"R1": {"value": input_value, "role": "input", "mapping": "vertex"}}
    for node, val in MEASURED.items():
        feats[node] = {"value": val, "role": "output", "mapping": "vertex"}
    return feats


def solve(samples: dict) -> tuple[list[dict], dict]:
    """Solve one CarnivalFlow problem; return the selected edges and objective values."""
    G = load_graph_from_sif_tuples(PKN)
    data = cn.Data.from_cdict(samples)
    method = CarnivalFlow(lambda_reg=LAMBDA)
    P = method.build_from_data(G, data)
    P.solve(solver=SOLVER, max_seconds=MAX_SECONDS, verbosity=0)
    objectives = {o.name: float(np.sum(o.value)) for o in P.objectives}

    edge_vals = np.asarray(P.expr.edge_value.value).reshape(method.processed_graph.num_edges, -1)
    rows = []
    for i, (src, tgt) in enumerate(method.processed_graph.E):
        if len(src) == 0 or len(tgt) == 0:
            continue  # boundary (flow) edges added by CORNETO
        s, t = next(iter(src)), next(iter(tgt))
        for j, name in enumerate(samples):
            v = float(edge_vals[i, j])
            if abs(v) > 0.5:
                rows.append({"sample": name, "source": s, "target": t, "edge_value": int(round(v))})
    return rows, objectives


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    checks: dict[str, bool] = {}
    all_rows: list[dict] = []

    rows, obj_known = solve({"known_sign": make_sample(1)})
    got = {(r["source"], r["target"]) for r in rows}
    checks["recovers_planted_network"] = got == EXPECTED
    checks["uses_metabolite_node"] = {("ENZ", METAB), (METAB, "TF3")} <= got
    checks["inhibition_sign_correct"] = any(
        r["source"] == "K2" and r["target"] == "TF2" and r["edge_value"] == -1 for r in rows
    )
    all_rows += rows

    rows, obj_unknown = solve({"unknown_sign": make_sample(0)})
    got0 = {(r["source"], r["target"]) for r in rows}
    checks["unknown_input_sign_inferred"] = got0 == EXPECTED
    all_rows += rows

    rows, obj_multi = solve({"s1": make_sample(1), "s2": make_sample(1)})
    per = {s: {(r["source"], r["target"]) for r in rows if r["sample"] == s} for s in ("s1", "s2")}
    checks["multisample_identical_samples_agree"] = per["s1"] == per["s2"] == EXPECTED
    all_rows += rows

    pd.DataFrame(all_rows).to_csv(OUT_DIR / "toy_edges.tsv", sep="\t", index=False)
    report = {
        "checks": checks,
        "all_passed": all(checks.values()),
        "expected_edges": sorted(map(list, EXPECTED)),
        "objectives": {"known_sign": obj_known, "unknown_sign": obj_unknown, "multisample": obj_multi},
        "settings": {"solver": SOLVER, "max_seconds": MAX_SECONDS, "lambda_reg": LAMBDA},
        "versions": {p: version(p) for p in ("corneto", "highspy", "cvxpy", "numpy", "pandas")},
        "python": platform.python_version(),
        "platform": platform.platform(),
    }
    (OUT_DIR / "toy_checks.json").write_text(json.dumps(report, indent=2))
    for name, ok in checks.items():
        print(f"{'PASS' if ok else 'FAIL'}  {name}")
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
