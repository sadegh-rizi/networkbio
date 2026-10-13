"""Synthetic checks for the Stage-06 model import and context helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(ROOT / "src"))

from gem import (  # noqa: E402
    context_reactions,
    gene_scores,
    load_recon3d,
    medium_bounds,
    prune_cosmos_to_context,
)
from gem_worker import make_task, solve_in_subprocess  # noqa: E402


def _toy_graph():
    try:
        from corneto.io import graph_from_vertex_incidence
    except ImportError:
        from corneto.io._metabolism import graph_from_vertex_incidence
    from corneto.methods import MultiSampleFBA, MultiSampleIMAT

    stoich = np.array([
        [-1, -1, -1, 0],
        [0, 1, 0, -1],
    ], dtype=float)
    graph = graph_from_vertex_incidence(stoich, ["a", "b"], ["EX_a", "R_HIGH", "R_LOW", "BIOMASS"])
    for reaction, gene, lower in (
        ("EX_a", "", -10.0),
        ("R_HIGH", "HIGH", 0.0),
        ("R_LOW", "LOW", 0.0),
        ("BIOMASS", "", 0.0),
    ):
        graph.get_attr_edge(list(graph.get_attr_from_edges("id")).index(reaction)).update(
            {"default_lb": lower, "default_ub": 1000.0, "GPR": gene}
        )
    return graph, MultiSampleFBA, MultiSampleIMAT


def test_toy_fba_has_known_optimum():
    graph, fba, _ = _toy_graph()
    problem = fba().build_many(graph, objectives={"toy": {"BIOMASS": -1}}, reaction_bounds={"toy": {"EX_a": (-2, 1000)}})
    solved = problem.solve(solver="highs", max_seconds=30, mip_rel_gap=0.01, threads=1, random_seed=1)
    assert solved.status == "optimal"
    assert problem.expr.flow.value[-1] == pytest.approx(2.0)


def test_imat_selects_the_high_expression_branch():
    graph, _, imat = _toy_graph()
    problem = imat(eps=1e-2, high_expression_threshold=1, low_expression_threshold=-1).build_many(
        graph,
        gene_expression={"toy": {"HIGH": 1.0, "LOW": -1.0}},
        reaction_bounds={"toy": {"EX_a": (-2, 1000)}},
        objectives={"toy": {}},
    )
    solved = problem.solve(solver="highs", max_seconds=30, mip_rel_gap=0.01, threads=1, random_seed=1)
    assert solved.status == "optimal"
    flux = dict(zip(graph.get_attr_from_edges("id"), np.asarray(problem.expr.flow.value)))
    assert flux["R_HIGH"] >= 0.01
    assert abs(flux["R_LOW"]) <= 0.01


def test_gpr_is_rewritten_to_unique_symbols(tmp_path):
    model = {
        "id": "toy",
        "version": "1",
        "metabolites": [{"id": "a", "name": "a", "annotation": {}}, {"id": "b", "name": "b", "annotation": {}}],
        "genes": [
            {"id": "1_AT1", "name": "AAA", "annotation": {"ncbigene": ["1"]}},
            {"id": "2_AT1", "name": "BBB", "annotation": {"ncbigene": ["2"]}},
        ],
        "reactions": [{"id": "R", "name": "r", "metabolites": {"a": -1, "b": 1}, "lower_bound": 0, "upper_bound": 1000,
                        "gene_reaction_rule": "1_AT1 and 2_AT1", "subsystem": "toy", "annotation": {}}],
    }
    path = tmp_path / "model.json"
    path.write_text(__import__("json").dumps(model))
    hgnc = pd.DataFrame({"symbol": ["AAA", "BBB"], "status": ["Approved", "Approved"], "entrez_id": ["1", "2"]})
    graph, mapping = load_recon3d(path, hgnc)
    assert graph.get_attr_edge(0)["GPR"] == "AAA and BBB"
    assert set(mapping["mapped_gene"]) == {"AAA", "BBB"}


def test_gpr_and_or_scores_follow_corneto_convention():
    from corneto.methods.metabolism._utils import evaluate_gpr

    assert evaluate_gpr("AAA and BBB", {"AAA": 1, "BBB": -1}) == -1
    assert evaluate_gpr("AAA or BBB", {"AAA": -1, "BBB": 1}) == 1


def test_gene_scores_add_unexpressed_model_genes_as_low():
    expression = pd.DataFrame({"C1": [1.0, 2.0, 3.0], "C2": [1.0, 2.0, 3.0], "P1": [3.0, 2.0, 1.0]}, index=["A", "B", "C"])
    scores = gene_scores(expression, {"C1": "Ctrl", "C2": "Ctrl", "P1": "PMS"}, model_genes={"A", "D"})
    assert scores["Ctrl"]["D"] == -1
    assert scores["Ctrl"]["C"] == 1


def test_blocked_uptake_makes_toy_demand_infeasible():
    graph, fba, _ = _toy_graph()
    problem = fba().build_many(
        graph,
        objectives={"toy": {"BIOMASS": -1}},
        reaction_bounds={"toy": {"EX_a": (0, 0), "BIOMASS": (1, 1000)}},
    )
    solved = problem.solve(solver="highs", max_seconds=30, mip_rel_gap=0.01, threads=1, random_seed=1)
    assert solved.status == "infeasible"


def test_context_pruning_removes_only_inactive_reaction_nodes():
    pkn = pd.DataFrame([
        ("A", 1, "Gene1__AAA"), ("Gene1__AAA", 1, "B"),
        ("A", 1, "Gene2__BBB"), ("Gene2__BBB", 1, "B"),
        ("B", 1, "C"),
    ], columns=["source", "sign", "target"])
    pruned = prune_cosmos_to_context(pkn, {"Gene1__AAA"})
    assert set(zip(pruned["source"], pruned["target"])) == {("A", "Gene1__AAA"), ("Gene1__AAA", "B"), ("B", "C")}


def test_context_reactions_marks_direction_and_activity():
    fluxes = pd.DataFrame({"condition": ["Ctrl", "PMS"], "reaction_number": [1, 1], "reaction": ["R", "R"], "flux": [0.0, -0.2]})
    result = context_reactions(fluxes)
    assert result.loc[0, "direction"] == "inactive"
    assert result.loc[1, "direction"] == "reverse"


def test_gem_worker_resumes_a_completed_fba_task(tmp_path):
    model = {
        "id": "toy",
        "version": "1",
        "metabolites": [{"id": "a", "name": "a", "annotation": {}}, {"id": "b", "name": "b", "annotation": {}}],
        "genes": [],
        "reactions": [
            {"id": "EX_a", "name": "uptake", "metabolites": {"a": -1}, "lower_bound": -2, "upper_bound": 1000, "gene_reaction_rule": "", "subsystem": "exchange", "annotation": {}},
            {"id": "R", "name": "conversion", "metabolites": {"a": -1, "b": 1}, "lower_bound": 0, "upper_bound": 1000, "gene_reaction_rule": "", "subsystem": "toy", "annotation": {}},
            {"id": "BIOMASS", "name": "biomass", "metabolites": {"b": -1}, "lower_bound": 0, "upper_bound": 1000, "gene_reaction_rule": "", "subsystem": "toy", "annotation": {}},
        ],
    }
    model_path = tmp_path / "model.json"
    hgnc_path = tmp_path / "hgnc.tsv"
    model_path.write_text(__import__("json").dumps(model))
    hgnc_path.write_text("symbol\tstatus\tentrez_id\n")
    task = make_task(
        label="toy worker", model_path=model_path, hgnc_path=hgnc_path, method="fba",
        objectives={"toy": {"BIOMASS": -1}}, reaction_bounds={"toy": {"EX_a": (-2, 1000)}},
        time_limit=30, root=tmp_path,
    )
    result = solve_in_subprocess(task, tmp_path / "solve", root=tmp_path)
    assert result["accepted"]
    assert result["fluxes"].loc[result["fluxes"]["reaction"].eq("BIOMASS"), "flux"].iloc[0] == pytest.approx(2)
    resumed = solve_in_subprocess(task, tmp_path / "solve", root=tmp_path, resume=True)
    assert resumed["task_signature"] == task["signature"]
