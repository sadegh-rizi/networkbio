"""Network properties of the inferred networks (no new solves).

Reads stage-03 PKNs, stage-04 networks and stage-05 rewired networks and
writes, to results/ionescu_corneto/07_downstream/network_properties/:

- degree_bias.tsv: does node selection depend on PKN degree?  Logistic
  regression of "selected" on log2 total degree over the unmeasured PKN nodes,
  per network (measured TFs/metabolites are in every network by design).
- input_degree.tsv: out-degree percentile of the starting nodes CORNETO used,
  relative to every candidate starting node of that run.
- centrality.tsv: degree and betweenness percentiles in the PKN for the
  recurring kinase core and the hypothesis nodes.  Betweenness is exact for
  OmniPath and sampled (k pivots, fixed seed) for COSMOS.
- metabolite_tf_distance.tsv: shortest directed distance from every COSMOS
  metabolite node to the hypothesis TFs; how many metabolites are as close as
  cholesterol.

    python scripts/analysis/07_network_properties.py [--out DIR] [--betweenness-k 1000]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

PKN = ROOT / "results/ionescu_corneto/03_pkn"
STAGE04 = ROOT / "results/ionescu_corneto/04_inference"
STAGE05 = ROOT / "results/ionescu_corneto/05_robustness"
OUT = ROOT / "results/ionescu_corneto/07_downstream/network_properties"
RUN_PKN = {"omnipath": "omnipath/primary", "cosmos_tf": "cosmos/primary", "cosmos_multi": "cosmos/primary",
           "cosmos_repair": "cosmos/variants/canonical_repair"}
LAMBDAS = ("0.01", "0.03", "0.1", "0.3", "1")
CORE = ["MAPK1", "MAPK14", "PRKACA", "FYN", "AKT1"]
HYPOTHESIS = ["HMGCR", "SREBF1", "SREBF2", "E2F1", "EGR1", "WT1", "SP1", "JUN", "JUNB"]
CHOLESTEROL = "Metab__HMDB0000067_c"
TARGET_TFS = ["E2F1", "EGR1", "SP1", "JUN", "JUNB"]
SEED = 20261010


def _pkn(name: str) -> nx.DiGraph:
    table = pd.read_csv(PKN / RUN_PKN.get(name, name) / "pkn.tsv", sep="\t", dtype={"source": str, "target": str})
    return nx.DiGraph(list(zip(table["source"], table["target"])))


def _network(directory: Path) -> tuple[set[str], pd.DataFrame]:
    edges = pd.read_csv(directory / "edges.tsv", sep="\t", dtype={"source": str, "target": str})
    nodes = pd.read_csv(directory / "nodes.tsv", sep="\t", dtype={"node": str})
    return set(edges["source"]) | set(edges["target"]), nodes


def _logit(degree: dict[str, int], selected: set[str], universe: list[str]) -> tuple[float, float, float]:
    y = np.array([int(n in selected) for n in universe])
    x = sm.add_constant(np.log2(np.array([degree[n] for n in universe], dtype=float)))
    try:
        fit = sm.Logit(y, x).fit(disp=False)
        return float(np.exp(fit.params[1])), float(np.exp(fit.conf_int()[1][0])), float(fit.pvalues[1])
    except Exception:  # perfect separation or too few selected nodes
        return np.nan, np.nan, np.nan


def degree_bias() -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(STAGE04 / "inputs/inputs.tsv", sep="\t")
    rows, input_rows = [], []
    networks = [(run, lam, STAGE04 / run / f"lambda_{lam}", run) for run in RUN_PKN for lam in LAMBDAS]
    networks += [("omnipath_rewired", str(rep), STAGE05 / "rewired/omnipath" / str(rep), None)
                 for rep in range(10) if (STAGE05 / "rewired/omnipath" / str(rep) / "edges.tsv").exists()]
    graphs: dict[str, nx.DiGraph] = {}
    for run, level, directory, pkn_name in networks:
        if pkn_name is None:  # rewired: its own PKN, same degrees as OmniPath by construction
            table = pd.read_csv(directory / "pkn.tsv", sep="\t", dtype={"source": str, "target": str})
            graph = nx.DiGraph(list(zip(table["source"], table["target"])))
        else:
            graph = graphs.setdefault(pkn_name, _pkn(pkn_name))
        selected, nodes = _network(directory)
        measured = set(nodes.loc[nodes["is_measured"], "node"])
        degree = dict(graph.degree())
        universe = sorted(n for n in graph if n not in measured)
        odds, low, pvalue = _logit(degree, selected - measured, universe)
        chosen = [n for n in nodes.loc[nodes["is_input"] & nodes["state"].ne(0), "node"] if n in graph]
        rows.append({"network": run, "level": level, "n_universe": len(universe),
                     "n_selected_unmeasured": len(selected - measured), "odds_ratio_per_degree_doubling": odds,
                     "or_ci_low": low, "pvalue": pvalue,
                     "median_degree_selected": float(np.median([degree[n] for n in selected - measured if n in degree] or [np.nan])),
                     "median_degree_all": float(np.median([degree[n] for n in universe]))})
        base_run = "omnipath" if pkn_name is None else run
        pool = candidates.loc[candidates["run"].eq(base_run) & candidates["used_as_input"], "node"]
        out_degree = pd.Series({n: graph.out_degree(n) for n in pool if n in graph})
        if chosen and not out_degree.empty:
            reach = set(nodes.loc[nodes["is_measured"] & nodes["reachability"].eq("reachable"), "node"])
            sub = nx.DiGraph(list(zip(*[pd.read_csv(directory / "edges.tsv", sep="\t", dtype=str)[c] for c in ("source", "target")])))
            explained = {n: len(nx.descendants(sub, n) & reach) if n in sub else 0 for n in chosen}
            hub = max(explained, key=explained.get)
            for n in chosen:
                input_rows.append({"network": run, "level": level, "input": n, "is_dominant": n == hub,
                                   "n_explained": explained[n], "out_degree": int(graph.out_degree(n)),
                                   "out_degree_percentile": float((out_degree < graph.out_degree(n)).mean() * 100),
                                   "n_candidates": len(out_degree)})
    return pd.DataFrame(rows), pd.DataFrame(input_rows)


def centrality(k: int) -> pd.DataFrame:
    rows = []
    for name, graph_name in (("omnipath", "omnipath"), ("cosmos", "cosmos_multi")):
        graph = _pkn(graph_name)
        exact = name == "omnipath"
        between = nx.betweenness_centrality(graph, k=None if exact else min(k, graph.number_of_nodes()), seed=SEED)
        degree = dict(graph.degree())
        b = pd.Series(between); d = pd.Series(degree)
        proteins = [n for n in graph if not n.startswith(("Metab__", "Gene"))]
        for node in CORE + HYPOTHESIS + ([CHOLESTEROL] if name == "cosmos" else []):
            if node not in graph:
                rows.append({"pkn": name, "node": node, "in_pkn": False})
                continue
            peers = [n for n in graph if n.startswith("Metab__")] if node.startswith("Metab__") else proteins
            rows.append({
                "pkn": name, "node": node, "in_pkn": True, "degree": degree[node],
                "degree_percentile": float((d[peers] < degree[node]).mean() * 100),
                "betweenness": between[node],
                "betweenness_percentile": float((b[peers] < between[node]).mean() * 100),
                "peer_group": "metabolites" if node.startswith("Metab__") else "proteins",
                "betweenness_method": "exact" if exact else f"sampled k={k}, seed {SEED}",
            })
    return pd.DataFrame(rows)


def metabolite_tf_distance() -> pd.DataFrame:
    graph = _pkn("cosmos_multi")
    reverse = graph.reverse(copy=False)
    distance = {tf: nx.single_source_shortest_path_length(reverse, tf) for tf in TARGET_TFS if tf in graph}
    metabolites = [n for n in graph if n.startswith("Metab__")]
    rows = []
    for tf, lengths in distance.items():
        values = pd.Series({m: lengths.get(m, np.inf) for m in metabolites})
        chol = values.get(CHOLESTEROL, np.inf)
        rows.append({"tf": tf, "cholesterol_distance": chol,
                     "n_metabolites": len(values),
                     "n_metabolites_reaching": int(np.isfinite(values).sum()),
                     "n_metabolites_as_close_or_closer": int((values <= chol).sum()),
                     "fraction_as_close_or_closer": float((values <= chol).mean()),
                     "median_distance_reaching": float(values[np.isfinite(values)].median())})
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--betweenness-k", type=int, default=1000)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    bias, inputs = degree_bias()
    bias.to_csv(args.out / "degree_bias.tsv", sep="\t", index=False)
    inputs.to_csv(args.out / "input_degree.tsv", sep="\t", index=False)
    centrality(args.betweenness_k).to_csv(args.out / "centrality.tsv", sep="\t", index=False)
    metabolite_tf_distance().to_csv(args.out / "metabolite_tf_distance.tsv", sep="\t", index=False)
    print(f"Network properties written to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
