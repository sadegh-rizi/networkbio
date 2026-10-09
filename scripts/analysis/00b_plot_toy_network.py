"""Figure and Cytoscape tables for the toy CORNETO run (no solver needed).

Reads the outputs of `00_toy_corneto.py`:
    results/corneto_toy/04_inference/toy_edges.tsv
    results/corneto_toy/04_inference/toy_checks.json
and takes the toy PKN, inputs and measurements from the constants in
`00_toy_corneto.py` itself (parsed with `ast`, not imported, so CORNETO is
not needed and the PKN has a single source of truth).

Writes:
    results/corneto_toy/04_inference/toy_network_edges.tsv   (Cytoscape edge table)
    results/corneto_toy/04_inference/toy_network_nodes.tsv   (Cytoscape node table)
    results/corneto_toy/04_inference/figures/toy_network.{pdf,_presentation.png,_manuscript.png}

Run from the repository root (WSL):
    .venv/bin/python scripts/analysis/00b_plot_toy_network.py
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from network_plot import (  # noqa: E402
    METABOLITE_PREFIX,
    draw_signed_network,
    edge_table,
    layered_positions,
    network_legend_handles,
    node_states,
    pkn_signs,
)
from plot_config import apply_plot_style, save_figure  # noqa: E402

TOY_SCRIPT = ROOT / "scripts/analysis/00_toy_corneto.py"
STAGE = ROOT / "results/corneto_toy/04_inference"
INPUT_NODE = "R1"  # the only input in make_sample() of 00_toy_corneto.py
PANELS = [
    ("known_sign", "A  One sample, input sign given (R1 = +1)"),
    ("unknown_sign", "B  One sample, input sign inferred (R1 = 0)"),
    ("multisample", "C  Two identical samples, solved jointly"),
]


def load_toy_constants(path: Path = TOY_SCRIPT) -> dict:
    """Return PKN, MEASURED, EXPECTED and METAB from the toy script's source."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    wanted = {"METAB", "PKN", "MEASURED", "EXPECTED"}
    found: dict = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            name = getattr(node.targets[0], "id", None)
            if name in wanted:
                expr = ast.Expression(node.value)
                found[name] = eval(compile(expr, str(path), "eval"), {"__builtins__": {}}, dict(found))
    missing = wanted - set(found)
    if missing:
        raise ValueError(f"{path.name}: constants not found: {sorted(missing)}")
    return found


def split_runs(edges: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """toy_edges.tsv stacks three runs; the multi-sample run has samples s1, s2."""
    runs = {
        "known_sign": edges[edges["sample"] == "known_sign"],
        "unknown_sign": edges[edges["sample"] == "unknown_sign"],
        "multisample": edges[edges["sample"].isin(["s1", "s2"])],
    }
    if sum(len(v) for v in runs.values()) != len(edges):
        raise ValueError(f"unexpected sample names in toy_edges.tsv: {sorted(set(edges['sample']))}")
    return runs


def main() -> int:
    const = load_toy_constants()
    pkn, measured, metab = const["PKN"], const["MEASURED"], const["METAB"]
    signs = pkn_signs(pkn)
    edges = pd.read_csv(STAGE / "toy_edges.tsv", sep="\t")
    checks = json.loads((STAGE / "toy_checks.json").read_text())
    runs = split_runs(edges)
    settings = checks["settings"]

    nodes = list(dict.fromkeys([n for s, _, t in pkn for n in (s, t)]))
    labels = {metab: "Chol"}
    n_decoys = len(signs) - len(const["EXPECTED"])
    pos = layered_positions(nodes, signs)

    # Cytoscape tables: one edge row per PKN edge and run; one node row per node and run.
    etabs, ntabs = [], []
    for run, df in runs.items():
        et = edge_table(pkn, df)
        et.insert(0, "run", run)
        etabs.append(et)
        for sample, sdf in df.groupby("sample", sort=False):
            st = node_states(sdf, signs)
            for n in nodes:
                ntabs.append(
                    {
                        "run": run,
                        "sample": sample,
                        "node": n,
                        "node_type": "metabolite" if n.startswith(METABOLITE_PREFIX) else "protein",
                        "role": "input" if n == INPUT_NODE else ("measured" if n in measured else "other"),
                        "measured_sign": measured.get(n, ""),
                        "inferred_state": st.get(n, 0),
                        "in_solution": n in st,
                    }
                )
    edge_out, node_out = pd.concat(etabs, ignore_index=True), pd.DataFrame(ntabs)
    edge_out.to_csv(STAGE / "toy_network_edges.tsv", sep="\t", index=False)
    node_out.to_csv(STAGE / "toy_network_nodes.tsv", sep="\t", index=False)

    # Fit check written into the figure: does every measured node get its observed sign?
    sol = node_out[node_out["role"] == "measured"]
    mismatches = int((sol["inferred_state"] != sol["measured_sign"].astype(int)).sum())

    apply_plot_style()
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.3))
    for ax, (run, title) in zip(axes, PANELS):
        draw_signed_network(
            ax, pkn, runs[run], inputs=[INPUT_NODE], measured=measured, labels=labels, pos=pos, title=title
        )
    fig.legend(handles=network_legend_handles(), loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.5, -0.02))
    fig.text(
        0.5,
        -0.10,
        f"Toy PKN: {len(nodes)} nodes, {len(signs)} signed edges ({n_decoys} decoys); "
        f"Chol = {metab} (cholesterol, COSMOS-style name). CORNETO "
        f"{checks['versions']['corneto']} CarnivalFlow, \u03bb = {settings['lambda_reg']}, solver "
        f"{settings['solver']} (limit {settings['max_seconds']} s); one optimal solution per panel, no sampling. "
        f"Edge width and n/N in C = fraction of samples selecting the edge. "
        f"Measured nodes with sign mismatches across panels: {mismatches}.",
        ha="center",
        va="top",
        fontsize=7,
        wrap=True,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.20, wspace=0.05)
    written = save_figure(fig, STAGE / "figures" / "toy_network")
    plt.close(fig)

    for p in [STAGE / "toy_network_edges.tsv", STAGE / "toy_network_nodes.tsv", *written]:
        print(p.relative_to(ROOT))
    print(f"measured-sign mismatches: {mismatches}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
