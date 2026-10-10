"""Synthetic stage-03 tests; no network or project data are used."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pkn import (  # noqa: E402
    clean_cosmos_edges,
    filter_cosmos,
    filter_omnipath,
    load_hgnc_maps,
    map_metabolite_features,
    parse_cosmos_node,
    rewire_signed_graph,
    signed_shortest_paths,
)


def _maps() -> dict:
    return {
        "approved": {"A", "B", "C", "D"},
        "prev_symbol": {"OLD": {"A"}},
        "alias_symbol": {"AMB": {"A", "B"}, "ALIAS": {"C"}},
        "entrez": {"1": {"A"}, "2": {"B"}},
        "expressed": {"A", "B", "C", "D"},
    }


def _cosmos_maps() -> dict:
    maps = _maps()
    maps["approved"] = set(maps["approved"]) | {"HMGCR", "SREBF1", "SLC7A6"}
    maps["expressed"] = set(maps["expressed"]) | {"HMGCR"}
    return maps


def _omni_row(source, target, *, direction=True, stimulation=True, inhibition=False, sources="x"):
    return {
        "source_genesymbol": source,
        "target_genesymbol": target,
        "consensus_direction": direction,
        "consensus_stimulation": stimulation,
        "consensus_inhibition": inhibition,
        "curation_effort": 2,
        "sources": sources,
    }


def test_consensus_sign_rule_and_networkcommons_sensitivity():
    raw = pd.DataFrame(
        [
            _omni_row("A", "B", stimulation=True, inhibition=True),
            _omni_row("B", "C", direction=False),
            _omni_row("C", "D", stimulation=False, inhibition=True),
        ]
    )
    primary, _, _ = filter_omnipath(raw, expressed={"A", "B", "C", "D"}, hgnc_maps=_maps())
    nc_rule, _, _ = filter_omnipath(raw, expressed={"A", "B", "C", "D"}, hgnc_maps=_maps(), variant="nc_rule")
    assert set(map(tuple, primary.itertuples(index=False, name=None))) == {("C", -1, "D")}
    assert set(map(tuple, nc_rule.itertuples(index=False, name=None))) == {("A", 1, "B"), ("C", -1, "D")}


def test_complex_loops_and_duplicate_sign_conflicts_are_logged():
    raw = pd.DataFrame(
        [
            _omni_row("COMPLEX:AP1", "B"),
            _omni_row("A_A", "B"),
            _omni_row("A", "A"),
            _omni_row("B", "C"),
            _omni_row("B", "C"),
            _omni_row("C", "D", stimulation=False, inhibition=True),
            _omni_row("C", "D"),
        ]
    )
    edges, log, _ = filter_omnipath(raw, expressed={"A", "B", "C", "D"}, hgnc_maps=_maps())
    assert set(map(tuple, edges.itertuples(index=False, name=None))) == {("B", 1, "C")}
    p4 = log.set_index("step")
    assert p4.loc["P4_complexes", "n_complex_rows"] == 2
    assert p4.loc["P4_self_loops", "n_self_loop_rows"] == 1
    assert p4.loc["P4_duplicates", "conflicting_pairs"] == 1
    assert p4.loc["P4_duplicates", "identical_rows_merged"] == 1


def test_symbol_harmonisation_prefers_unique_previous_or_alias_only():
    assert load_hgnc_maps
    raw = pd.DataFrame([_omni_row("OLD", "ALIAS"), _omni_row("AMB", "C")])
    edges, _, _ = filter_omnipath(raw, expressed={"A", "B", "C", "D", "AMB"}, hgnc_maps=_maps(), variant="unpruned")
    assert ("A", 1, "C") in set(map(tuple, edges.itertuples(index=False, name=None)))
    assert ("AMB", 1, "C") in set(map(tuple, edges.itertuples(index=False, name=None)))


def test_expression_pruning_keeps_metabolites_but_removes_unexpressed_genes():
    raw = pd.DataFrame(
        [
            ("HMGCR", 1, "Gene1600__HMGCR"),
            ("Gene1600__HMGCR_reverse", 1, "Metab__HMDB0000227_c"),
            ("Metab__2676_c", 1, "Gene1601__EX_glc__D_e"),
            ("Gene1602__HMGCR_SREBF1", 1, "Metab__2676_c"),
            ("Metab__HMDB0000227_c", 1, "Gene1603__B"),
        ],
        columns=["source", "sign", "target"],
    )
    edges, _, nodes = filter_cosmos(raw, hgnc_maps=_cosmos_maps(), expressed={"HMGCR"})
    edge_set = set(map(tuple, edges.itertuples(index=False, name=None)))
    assert ("HMGCR", 1, "Gene1600__HMGCR") in edge_set
    assert ("Gene1600__HMGCR_reverse", 1, "Metab__HMDB0000227_c") in edge_set
    assert ("Metab__2676_c", 1, "Gene1601__EX_glc__D_e") in edge_set
    assert ("Gene1602__HMGCR_SREBF1", 1, "Metab__2676_c") in edge_set
    assert "Metab__2676_c" in set(nodes["node"])
    assert "Gene1603__B" not in set(nodes["node"])


def test_cosmos_parser_uses_suffix_genes_and_preserves_numeric_ids():
    enzyme = parse_cosmos_node("Gene1600__HMGCR_reverse", _cosmos_maps())
    multi_gene = parse_cosmos_node("Gene1602__HMGCR_SREBF1", _cosmos_maps())
    transport_gene = parse_cosmos_node("Gene10001__SLC7A6_TRANSPORTER1", _cosmos_maps())
    pseudo = parse_cosmos_node("Gene1601__EX_glc__D_e", _cosmos_maps())
    numeric = parse_cosmos_node("Metab__2676_c", _cosmos_maps())

    assert enzyme["approved_symbol"] == "HMGCR"
    assert enzyme["type"] == "enzyme_or_reaction"
    assert multi_gene["approved_symbol"] == "HMGCR;SREBF1"
    assert transport_gene["approved_symbol"] == "SLC7A6"
    assert transport_gene["type"] == "enzyme_or_reaction"
    assert pseudo["type"] == "transport_or_pseudo"
    assert numeric["node"] == "Metab__2676_c"
    assert numeric["type"] == "metabolite_numeric_id"


def test_cosmos_parser_rejects_incompatible_x_grammar():
    for node in ("XMetab__227___c____", "X3156", "XGene1600__3156"):
        with pytest.raises(ValueError, match="X-prefixed"):
            parse_cosmos_node(node, _cosmos_maps())


def test_cosmos_parser_accepts_plain_symbols_starting_with_x():
    maps = _cosmos_maps()
    maps["approved"] = set(maps["approved"]) | {"XBP1", "XIAP"}
    maps["expressed"] = set(maps["expressed"]) | {"XBP1"}
    xbp1 = parse_cosmos_node("XBP1", maps)
    xiap = parse_cosmos_node("XIAP", maps)
    assert (xbp1["type"], xbp1["node"], xbp1["expressed"]) == ("gene", "XBP1", True)
    assert (xiap["type"], xiap["expressed"]) == ("gene", False)


def test_cosmos_parser_keeps_model_metabolites_with_underscores():
    for node, compartment in (("Metab__gd1b2_hs_g", "g"), ("Metab__2hibup_S_r", "r"), ("Metab__pa_hs_e", "e")):
        parsed = parse_cosmos_node(node, _cosmos_maps())
        assert parsed["type"] == "metabolite_model_id"
        assert parsed["node"] == node
        assert parsed["compartment"] == compartment
        assert parsed["hmdb_id"] is None
    hmdb = parse_cosmos_node("Metab__HMDB0000067_c", _cosmos_maps())
    assert (hmdb["type"], hmdb["node"], hmdb["compartment"]) == ("metabolite", "Metab__HMDB0000067_c", "c")
    padded = parse_cosmos_node("Metab__HMDB10384 _c", _cosmos_maps())
    assert (padded["type"], padded["node"]) == ("metabolite", "Metab__HMDB0010384_c")


def test_cosmos_pruning_keeps_model_metabolites_and_logs_removed_types():
    raw = pd.DataFrame(
        [
            ("Gene1600__HMGCR", 1, "Metab__gd1b2_hs_g"),
            ("Metab__gd1b2_hs_g", 1, "HMGCR"),
            ("HMGCR", 1, "Gene1600__HMGCR"),
            ("A_B", 1, "HMGCR"),
            ("B", 1, "HMGCR"),
        ],
        columns=["source", "sign", "target"],
    )
    edges, log, nodes = filter_cosmos(raw, hgnc_maps=_cosmos_maps(), expressed={"HMGCR"})
    assert "Metab__gd1b2_hs_g" in set(nodes["node"])
    assert {"A_B", "B"}.isdisjoint(set(nodes["node"]))
    steps = log.set_index("step")
    assert steps.loc["P5_hgnc", "n_metabolite_model_id_nodes"] == 1
    assert steps.loc["P5_hgnc", "n_other_nodes"] == 1
    assert steps.loc["P5_hgnc", "n_gene_nodes"] == 2
    assert steps.loc["P6_expression", "n_removed_other_nodes"] == 1
    assert steps.loc["P6_expression", "n_removed_gene_nodes"] == 1
    assert steps.loc["P6_expression", "n_removed_metabolite_model_id_nodes"] == 0


def test_cosmos_duplicate_signs_average_to_zero_and_drop():
    raw = pd.DataFrame(
        [("A", 1, "B"), ("A", -1, "B"), ("B", 1, "C")],
        columns=["source", "sign", "target"],
    )
    cleaned = clean_cosmos_edges(raw)
    assert set(map(tuple, cleaned.itertuples(index=False, name=None))) == {("B", 1, "C")}


def test_metabolite_mapping_uses_all_recon_ids_and_carries_labelling_class():
    features = {
        "ST003331": pd.DataFrame(
            [
                {
                    "feature_id": "f1",
                    "metabolite_name": "test",
                    "refmet_name": "test",
                    "kegg_id": "C00001",
                    "labelling_class": "m0_unconfounded",
                }
            ]
        )
    }
    mapping = map_metabolite_features(
        features,
        pkn_nodes={"Metab__HMDB0000001_c", "Metab__HMDB0000001_m", "Metab__HMDB0000002_c", "Metab__HMDB0000002_m"},
        recon_by_kegg={"C00001": ["HMDB0000001", "HMDB0000002"]},
        workbench=pd.DataFrame(columns=["kegg_id", "hmdb_id", "source", "error"]),
    )
    assert set(mapping["hmdb_id"]) == {"HMDB0000001", "HMDB0000002"}
    assert set(mapping["compartment"]) == {"c", "m"}
    assert set(mapping["mapping_source"]) == {"recon3d"}
    assert set(mapping["labelling_class"]) == {"m0_unconfounded"}


def test_signed_bfs_reports_both_signs_and_respects_cap():
    edges = [("A", 1, "B"), ("B", -1, "C"), ("A", -1, "D"), ("D", -1, "C"), ("C", 1, "E")]
    result = signed_shortest_paths(edges, ["A"], ["C", "E"], max_length=2).set_index("target")
    assert result.loc["C", "shortest_length"] == 2
    assert result.loc["C", "shortest_signs"] == "+1;-1"
    assert result.loc["C", "n_shortest_positive"] == 1
    assert result.loc["C", "n_shortest_negative"] == 1
    assert pd.isna(result.loc["E", "shortest_length"])
    assert result.loc["E", "shortest_signs"] == ""


def test_rewiring_preserves_directed_degrees_sign_counts_and_seed():
    edges = [(str(i), 1 if i % 2 else -1, str((i + 1) % 6)) for i in range(6)]
    rewired_a = rewire_signed_graph(edges, seed=20261009, accepted_swaps=3)
    rewired_b = rewire_signed_graph(edges, seed=20261009, accepted_swaps=3)
    assert rewired_a == rewired_b
    assert len({(source, target) for source, _, target in rewired_a}) == len(rewired_a)
    assert all(source != target for source, _, target in rewired_a)
    assert sorted((source, target) for source, _, target in edges) != []
    for sign in (-1, 1):
        assert sum(edge[1] == sign for edge in edges) == sum(edge[1] == sign for edge in rewired_a)
    assert {source: sum(edge[0] == source for edge in edges) for source in {edge[0] for edge in edges}} == {
        source: sum(edge[0] == source for edge in rewired_a) for source in {edge[0] for edge in rewired_a}
    }
    assert {target: sum(edge[2] == target for edge in edges) for target in {edge[2] for edge in edges}} == {
        target: sum(edge[2] == target for edge in rewired_a) for target in {edge[2] for edge in rewired_a}
    }
