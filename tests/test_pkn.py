"""Synthetic stage-03 tests; no network or project data are used."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pkn import (  # noqa: E402
    clean_cosmos_edges,
    filter_cosmos,
    filter_omnipath,
    load_hgnc_maps,
    map_metabolite_features,
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
            ("X1", 1, "XMetab__67___c____"),
            ("XMetab__67___c____", 1, "X2"),
        ],
        columns=["source", "sign", "target"],
    )
    edges, _, nodes = filter_cosmos(raw, hgnc_maps=_maps(), expressed={"A"})
    assert set(map(tuple, edges.itertuples(index=False, name=None))) == {("A", 1, "Metab__HMDB0000067_c")}
    assert "Metab__HMDB0000067_c" in set(nodes["node"])
    assert "B" not in set(nodes["node"])


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
