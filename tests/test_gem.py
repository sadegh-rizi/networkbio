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
    directional_reaction_nodes,
    is_active,
    gene_scores,
    gpr_token,
    load_recon3d,
    cap_synthesis_direction,
    medium_bounds,
    mevalonate_synthesis_sign,
    minimum_extra_exchanges,
    prune_cosmos_to_context,
)
from gem_worker import SOLVER_THREADS, clear_model_cache, make_task, solve_in_process, solve_in_subprocess  # noqa: E402

# All in-process HiGHS solves in this file use the worker's thread count: HiGHS keeps
# one thread pool per process, and a different count makes later solves fail.


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
    solved = problem.solve(solver="highs", max_seconds=30, mip_rel_gap=0.01, threads=SOLVER_THREADS, random_seed=1)
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
    solved = problem.solve(solver="highs", max_seconds=30, mip_rel_gap=0.01, threads=SOLVER_THREADS, random_seed=1)
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
    solved = problem.solve(solver="highs", max_seconds=30, mip_rel_gap=0.01, threads=SOLVER_THREADS, random_seed=1)
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


def test_gpr_tokens_are_safe_for_the_corneto_parser():
    assert gpr_token("MT-ND3") == "MT_ND3"
    assert gpr_token("NT5C1B-RDH14") == "NT5C1B_RDH14"
    assert gpr_token("2713") == "ENTREZ_2713"
    assert gpr_token("HMGCR") == "HMGCR"


def test_hyphenated_and_unmapped_genes_are_scored_correctly(tmp_path):
    """MT-ND3 must not parse as subtraction; an unmapped Entrez ID must not score as its number."""
    from corneto.methods.metabolism._utils import evaluate_gpr

    model = {
        "id": "toy",
        "version": "1",
        "metabolites": [{"id": "a", "name": "a", "annotation": {}}, {"id": "b", "name": "b", "annotation": {}}],
        "genes": [
            {"id": "4537_AT1", "name": "MT-ND3", "annotation": {"ncbigene": ["4537"]}},
            {"id": "2713_AT1", "name": "GK3P", "annotation": {"ncbigene": ["2713"]}},
        ],
        "reactions": [
            {"id": "R1", "name": "r1", "metabolites": {"a": -1, "b": 1}, "lower_bound": 0, "upper_bound": 1000,
             "gene_reaction_rule": "4537_AT1", "subsystem": "toy", "annotation": {}},
            {"id": "R2", "name": "r2", "metabolites": {"a": -1, "b": 1}, "lower_bound": 0, "upper_bound": 1000,
             "gene_reaction_rule": "2713_AT1", "subsystem": "toy", "annotation": {}},
        ],
    }
    path = tmp_path / "model.json"
    path.write_text(__import__("json").dumps(model))
    hgnc = pd.DataFrame({"symbol": ["MT-ND3"], "status": ["Approved"], "entrez_id": ["4537"]})
    graph, mapping = load_recon3d(path, hgnc)
    assert graph.get_attr_edge(0)["GPR"] == "MT_ND3"
    assert graph.get_attr_edge(1)["GPR"] == "ENTREZ_2713"
    unscored = set(mapping.loc[mapping["mapping_status"].ne("unique_approved_symbol"), "gpr_token"])
    assert unscored == {"ENTREZ_2713"}

    expression = pd.DataFrame({"C1": [10.0, 1.0], "P1": [10.0, 1.0]}, index=["MT-ND3", "OTHER"])
    expression.index = expression.index.map(gpr_token)
    scores = gene_scores(expression, {"C1": "Ctrl", "P1": "PMS"}, quantiles=(0.25, 0.75),
                         model_genes=set(mapping["gpr_token"]), unscored_genes=unscored)
    assert evaluate_gpr(graph.get_attr_edge(0)["GPR"], scores["Ctrl"]) == 1
    assert evaluate_gpr(graph.get_attr_edge(1)["GPR"], scores["Ctrl"]) == 0


def test_defined_medium_closes_sink_and_demand_supply(tmp_path):
    """SK_/DM_ reactions with a negative lower bound are uptakes and must be closed."""
    model = {
        "id": "toy",
        "version": "1",
        "metabolites": [{"id": m, "name": m, "annotation": {}} for m in ("a_e", "b_c", "c_c", "d_c")],
        "genes": [],
        "reactions": [
            {"id": "EX_a_e", "name": "ex", "metabolites": {"a_e": -1}, "lower_bound": -1000, "upper_bound": 1000,
             "gene_reaction_rule": "", "subsystem": "Exchange", "annotation": {}},
            {"id": "SK_b_c", "name": "sink", "metabolites": {"b_c": -1}, "lower_bound": -1000, "upper_bound": 1000,
             "gene_reaction_rule": "", "subsystem": "Sink", "annotation": {}},
            {"id": "DM_c_c", "name": "demand", "metabolites": {"c_c": -1}, "lower_bound": -1000, "upper_bound": 1000,
             "gene_reaction_rule": "", "subsystem": "Demand", "annotation": {}},
            {"id": "DM_d_c", "name": "demand out", "metabolites": {"d_c": -1}, "lower_bound": 0, "upper_bound": 1000,
             "gene_reaction_rule": "", "subsystem": "Demand", "annotation": {}},
        ],
    }
    path = tmp_path / "model.json"
    path.write_text(__import__("json").dumps(model))
    graph, _ = load_recon3d(path, pd.DataFrame({"symbol": [], "status": [], "entrez_id": []}))
    bounds = medium_bounds(graph, "defined_chol_minus")
    assert bounds["EX_a_e"] == (0.0, 1000.0)
    assert bounds["SK_b_c"] == (0.0, 1000.0)
    assert bounds["DM_c_c"] == (0.0, 1000.0)
    assert "DM_d_c" not in bounds
    assert medium_bounds(graph, "open") == {}


def test_minimum_extra_exchanges_reopens_the_needed_uptake(tmp_path):
    """Biomass needs a; with EX_a closed the MILP must select EX_a, not report infeasible."""
    model = {
        "id": "toy",
        "version": "1",
        "metabolites": [{"id": "a_e", "name": "a", "annotation": {}}, {"id": "b_e", "name": "b", "annotation": {}}],
        "genes": [],
        "reactions": [
            {"id": "EX_a_e", "name": "ex a", "metabolites": {"a_e": -1}, "lower_bound": -1000, "upper_bound": 1000,
             "gene_reaction_rule": "", "subsystem": "Exchange", "annotation": {}},
            {"id": "EX_b_e", "name": "ex b", "metabolites": {"b_e": -1}, "lower_bound": -1000, "upper_bound": 1000,
             "gene_reaction_rule": "", "subsystem": "Exchange", "annotation": {}},
            {"id": "BIOMASS", "name": "biomass", "metabolites": {"a_e": -1}, "lower_bound": 0, "upper_bound": 1000,
             "gene_reaction_rule": "", "subsystem": "Biomass", "annotation": {}},
        ],
    }
    path = tmp_path / "model.json"
    path.write_text(__import__("json").dumps(model))
    graph, _ = load_recon3d(path, pd.DataFrame({"symbol": [], "status": [], "entrez_id": []}))
    bounds = {"EX_a_e": (0.0, 1000.0), "EX_b_e": (0.0, 1000.0)}
    result = minimum_extra_exchanges(graph, bounds, "BIOMASS", 5.0, time_limit=30)
    assert result["status"] == "optimal"
    assert result["selected_exchanges"] == ["EX_a_e"]


def _hmg_toy(tmp_path):
    model = {
        "id": "toy",
        "version": "1",
        "metabolites": [{"id": m, "name": m, "annotation": {}} for m in
                        ("hmgcoa_c", "mev__R_c", "hmgcoa_x", "mev__R_x", "phe__L_e", "ile__L_e", "leu__L_e")],
        "genes": [],
        "reactions": [
            # Recon3D orientation: cytosolic HMGCOAR is written mevalonate -> HMG-CoA.
            {"id": "HMGCOAR", "name": "c", "metabolites": {"mev__R_c": -1, "hmgcoa_c": 1}, "lower_bound": -1000,
             "upper_bound": 1000, "gene_reaction_rule": "", "subsystem": "Sterol", "annotation": {}},
            {"id": "HMGCOARx", "name": "x", "metabolites": {"hmgcoa_x": -1, "mev__R_x": 1}, "lower_bound": 0,
             "upper_bound": 1000, "gene_reaction_rule": "", "subsystem": "Sterol", "annotation": {}},
        ] + [{"id": f"EX_{m}", "name": m, "metabolites": {m: -1}, "lower_bound": -1000, "upper_bound": 1000,
              "gene_reaction_rule": "", "subsystem": "Exchange", "annotation": {}} for m in ("phe__L_e", "ile__L_e", "leu__L_e")],
    }
    path = tmp_path / "model.json"
    path.write_text(__import__("json").dumps(model))
    graph, _ = load_recon3d(path, pd.DataFrame({"symbol": [], "status": [], "entrez_id": []}))
    return graph


def test_simvastatin_cap_limits_the_mevalonate_forming_direction(tmp_path):
    graph = _hmg_toy(tmp_path)
    assert mevalonate_synthesis_sign(graph, "HMGCOAR") == -1
    assert mevalonate_synthesis_sign(graph, "HMGCOARx") == 1
    # Synthesis capacity 1000 in each direction of writing; cap at 10 %.
    assert cap_synthesis_direction((-1000.0, 500.0), -1, 1000.0, 0.1) == (-100.0, 500.0)
    assert cap_synthesis_direction((0.0, 1000.0), 1, 1000.0, 0.1) == (0.0, 100.0)


def test_defined_medium_opens_unmeasured_essential_amino_acids(tmp_path):
    graph = _hmg_toy(tmp_path)
    bounds = medium_bounds(graph, "defined_chol_minus")
    assert bounds["EX_phe__L_e"][0] < 0
    assert bounds["EX_ile__L_e"][0] < 0
    assert bounds["EX_leu__L_e"][0] == 0  # not detected in the (empty) blank set here


# --------------------------------------------------------------------------
# Addendum D9: activity threshold at eps, and in-process FBA
# --------------------------------------------------------------------------

def test_context_reactions_tolerate_solver_values_at_eps():
    """iMAT 'on' reactions sit at eps; the solver returns 0.01 or 0.0099999999999909."""
    flux = [0.01, 0.0099999999999909, -0.01, -0.0099999999999909, 0.005, 2.0, 0.0]
    fluxes = pd.DataFrame({"condition": "Ctrl", "reaction_number": range(1, 8),
                           "reaction": [f"R{i}" for i in range(1, 8)], "flux": flux})
    result = context_reactions(fluxes, 1e-2)
    assert result["active"].tolist() == [True, True, True, True, False, True, False]
    assert result["direction"].tolist() == ["forward", "forward", "reverse", "reverse", "inactive", "forward", "inactive"]
    # every active reaction has a direction, every inactive one has none
    assert (result["active"] == result["direction"].ne("inactive")).all()
    assert is_active(flux, 1e-2).tolist() == result["active"].tolist()


def test_directional_nodes_keep_reactions_running_at_eps():
    fluxes = pd.DataFrame({"condition": ["Ctrl", "Ctrl", "PMS"], "reaction_number": [1, 2, 2],
                           "reaction": ["R1", "R2", "R2"], "flux": [0.0099999999999909, 0.0, -0.01]})
    nodes = {"Gene1__AAA", "Gene1__AAA_reverse", "Gene2__BBB", "Gene2__BBB_reverse", "X"}
    assert directional_reaction_nodes(fluxes, nodes, 1e-2) == {"Gene1__AAA", "Gene2__BBB_reverse"}


def _worker_toy(tmp_path):
    model = {
        "id": "toy", "version": "1",
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
    return model_path, hgnc_path


def test_in_process_fba_matches_subprocess_and_reuses_the_cached_model(tmp_path):
    model_path, hgnc_path = _worker_toy(tmp_path)
    clear_model_cache()
    results = []
    for bound in (-2, -3, -2):  # same cached model, different bounds per task
        task = make_task(label=f"toy {bound}", model_path=model_path, hgnc_path=hgnc_path, method="fba",
                         objectives={"toy": {"BIOMASS": -1}}, reaction_bounds={"toy": {"EX_a": (bound, 1000)}},
                         time_limit=30, root=tmp_path)
        result = solve_in_process(task, tmp_path / f"in_process_{len(results)}", root=tmp_path)
        assert result["accepted"]
        assert result.get("process") == "in_process"
        results.append(result["fluxes"].loc[result["fluxes"]["reaction"].eq("BIOMASS"), "flux"].iloc[0])
    assert results == pytest.approx([2.0, 3.0, 2.0])
    task = make_task(label="toy", model_path=model_path, hgnc_path=hgnc_path, method="fba",
                     objectives={"toy": {"BIOMASS": -1}}, reaction_bounds={"toy": {"EX_a": (-2, 1000)}},
                     time_limit=30, root=tmp_path)
    in_process = solve_in_process(task, tmp_path / "a", root=tmp_path)
    subprocess_result = solve_in_subprocess(task, tmp_path / "b", root=tmp_path)
    assert in_process["fluxes"]["flux"].tolist() == pytest.approx(subprocess_result["fluxes"]["flux"].tolist())
    for name in ("task.json", "result.json", "fluxes.tsv", "toy/result.json"):
        assert (tmp_path / "a" / name).exists()
    resumed = solve_in_process(task, tmp_path / "a", root=tmp_path, resume=True)
    assert resumed["task_signature"] == task["signature"]


def test_in_process_solve_refuses_imat(tmp_path):
    model_path, hgnc_path = _worker_toy(tmp_path)
    task = make_task(label="toy imat", model_path=model_path, hgnc_path=hgnc_path, method="imat",
                     gene_expression={"toy": {}}, reaction_bounds={"toy": {}}, time_limit=30, root=tmp_path)
    with pytest.raises(ValueError):
        solve_in_process(task, tmp_path / "imat", root=tmp_path)
