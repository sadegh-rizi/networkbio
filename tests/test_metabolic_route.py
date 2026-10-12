"""Synthetic checks for the Aim-5 structural route audit."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from metabolic_route import (  # noqa: E402
    canonical_step_table,
    gap_audit,
    path_variant_tables,
    repaired_pkn,
    typed_metabolite_tf_null,
)


def _maps() -> dict:
    return {
        "entrez": {
            "3156": {"HMGCR"},
            "2222": {"FDFT1"},
        }
    }


def _model() -> dict:
    return {
        "metabolites": [
            {"id": name, "compartment": name.rsplit("_", 1)[-1], "annotation": {}}
            for name in ("hmgcoa_c", "mev__R_c", "frdp_c", "HC01118_c", "sql_c")
        ],
        "reactions": [
            {
                "id": "HMGCOAR",
                "name": "HMG-CoA reductase",
                "subsystem": "cholesterol",
                "gene_reaction_rule": "3156_AT1",
                "metabolites": {"hmgcoa_c": 1.0, "mev__R_c": -1.0},
            },
            {
                "id": "r0170",
                "name": "FPP to presqualene",
                "subsystem": "cholesterol",
                "gene_reaction_rule": "2222_AT1",
                "metabolites": {"frdp_c": -2.0, "HC01118_c": 1.0},
            },
            {
                "id": "r0575",
                "name": "Presqualene to squalene",
                "subsystem": "cholesterol",
                "gene_reaction_rule": "2222_AT1",
                "metabolites": {"HC01118_c": -1.0, "sql_c": 1.0},
            },
        ],
    }


def _edges() -> pd.DataFrame:
    return pd.DataFrame([
        ("HMGCR", 1, "Gene1__HMGCR"),
        ("Gene1__HMGCR", 1, "Metab__hmgcoa_c"),
        ("Metab__hmgcoa_c", 1, "Gene1__HMGCR_reverse"),
        ("Gene1__HMGCR_reverse", 1, "Metab__mev__R_c"),
        ("Metab__HC01118_c", 1, "Gene3__FDFT1"),
        ("Gene3__FDFT1", 1, "Metab__sql_c"),
    ], columns=["source", "sign", "target"])


def test_gap_audit_identifies_missing_fdft1_and_repair_adds_only_it():
    pkn_nodes = set(_edges()["source"]) | set(_edges()["target"]) | {"Metab__frdp_c"}
    steps = canonical_step_table(_model(), _maps(), pkn_nodes=pkn_nodes)
    gaps = gap_audit(steps, _edges(), _edges())
    missing = gaps[gaps["recon_id"].eq("r0170")]
    assert len(missing) == 1
    assert missing.iloc[0]["status"] == "absent_from_cosmos"

    repaired, _, _ = repaired_pkn(
        _edges(),
        pd.DataFrame({"node": sorted(set(_edges()["source"]) | set(_edges()["target"])), "type": "gene"}),
        steps,
        gaps,
    )
    added = set(map(tuple, repaired.itertuples(index=False, name=None))) - set(map(tuple, _edges().itertuples(index=False, name=None)))
    assert added == {
        ("Metab__frdp_c", 1, "Gene2__FDFT1"),
        ("Gene2__FDFT1", 1, "Metab__HC01118_c"),
    }


def test_unannotated_intermediate_is_named_from_cosmos_neighbours():
    # Recon3D has no HMDB annotation for 7-dehydrocholesterol; COSMOS names it
    # by HMDB. The name must come from the neighbours of the COSMOS reaction node.
    model = {
        "metabolites": [
            {"id": "7dhchsterol_r", "compartment": "r", "annotation": {}},
            {"id": "chsterol_r", "compartment": "r", "annotation": {"hmdb": ["HMDB0000067"]}},
            {"id": "nadph_r", "compartment": "r", "annotation": {}},
        ],
        "reactions": [{
            "id": "DHCR72r", "name": "7-DHC reductase", "subsystem": "cholesterol",
            "gene_reaction_rule": "1717_AT1",
            "metabolites": {"7dhchsterol_r": -1.0, "nadph_r": -1.0, "chsterol_r": 1.0},
        }],
    }
    edges = pd.DataFrame([
        ("DHCR7", 1, "Gene1__DHCR7"),
        ("Metab__HMDB0000032_r", 1, "Gene1__DHCR7"),
        ("Gene1__DHCR7", 1, "Metab__HMDB0000067_r"),
    ], columns=["source", "sign", "target"])
    nodes = set(edges["source"]) | set(edges["target"])
    maps = {"entrez": {"1717": {"DHCR7"}}}
    without = canonical_step_table(model, maps, pkn_nodes=nodes)
    assert pd.isna(without.loc[0, "substrate"])
    steps = canonical_step_table(model, maps, pkn_nodes=nodes, edges=edges)
    assert steps.loc[0, "substrate"] == "Metab__HMDB0000032_r"
    assert steps.loc[0, "substrate_source"] == "cosmos_neighbour_of_Gene1"
    assert steps.loc[0, "product_source"] == "recon3d_annotation"
    gaps = gap_audit(steps, edges, edges)
    assert gaps.loc[0, "status"] == "present_forward"


def test_same_side_metabolites_are_not_a_conversion_step():
    # HMR_1532-like reaction: GGPP + PPi <=> DMAPP + FPP has DMAPP and FPP on
    # the same side, so it is not a DMAPP -> FPP step.
    model = {
        "metabolites": [
            {"id": name, "compartment": "c", "annotation": {}}
            for name in ("dmpp_c", "frdp_c", "ppi_c", "ttc_ggdp_c")
        ],
        "reactions": [
            {"id": "HMGCOAR", "name": "", "subsystem": "", "gene_reaction_rule": "3156_AT1",
             "metabolites": {"hmgcoa_c": 1.0, "mev__R_c": -1.0}},
            {"id": "HMR_1532", "name": "", "subsystem": "", "gene_reaction_rule": "2224_AT1",
             "metabolites": {"dmpp_c": 1.0, "frdp_c": 1.0, "ppi_c": -1.0, "ttc_ggdp_c": -1.0}},
        ],
    }
    maps = {"entrez": {"3156": {"HMGCR"}, "2224": {"FDPS"}}}
    steps = canonical_step_table(model, maps, pkn_nodes=set())
    assert "HMR_1532" not in set(steps["recon_id"])


def test_top_degree_sensitivity_keeps_protected_cholesterol_target():
    steps = pd.DataFrame([{
        "substrate": "Metab__HMDB0000227_c", "product": "Metab__HMDB0000067_c",
        "previous_intermediate": "Mevalonate", "next_intermediate": "Cholesterol",
        "reaction_node": "Gene2__X", "reaction_node_candidates": "Gene2__X", "branch": "core",
    }])
    edges = [("HMGCR", 1, "Gene1__HMGCR"), ("Gene1__HMGCR", 1, "Metab__HMDB0000227_c"),
             ("Metab__HMDB0000227_c", 1, "Gene2__X"), ("Gene2__X", 1, "Metab__HMDB0000067_c")]
    # Make cholesterol the highest-degree metabolite.
    edges += [("Metab__HMDB0000067_c", 1, f"Gene{n}__Y") for n in range(10, 20)]
    edges += [(f"Gene{n}__Y", 1, f"Metab__m{n}_c") for n in range(10, 20)]
    route_edges = pd.DataFrame(edges, columns=["source", "sign", "target"])
    cofactor_columns = ["model_id", "bigg_id", "compartment", "pkn_node", "degree", "protected", "mapped"]
    paths, _ = path_variant_tables({"repaired": route_edges}, steps, pd.DataFrame(columns=cofactor_columns))
    top = paths[(paths["path_variant"] == "metabolism_no_cofactors_weighted_top1pct")
                & (paths["target"] == "Metab__HMDB0000067_c")]
    assert bool(top["path_found"].iloc[0])
    assert top["fraction_waypoints_in_order"].iloc[0] > 0


def test_path_variants_report_missing_and_repaired_routes():
    steps = pd.DataFrame([{
        "substrate": "HMGCR",
        "product": "Metab__HMDB0000067_c",
        "reaction_node": "Gene1__HMGCR",
        "reaction_node_candidates": "Gene1__HMGCR",
        "branch": "core",
    }])
    route_edges = pd.DataFrame([
        ("HMGCR", 1, "Gene1__HMGCR"),
        ("Gene1__HMGCR", 1, "Metab__HMDB0000067_c"),
    ], columns=["source", "sign", "target"])
    cofactor_columns = ["model_id", "bigg_id", "compartment", "pkn_node", "degree", "protected", "mapped"]
    paths, _ = path_variant_tables(
        {"primary": route_edges, "repaired": route_edges},
        steps,
        pd.DataFrame(columns=cofactor_columns),
    )
    assert not paths[(paths["variant"] == "repaired") & paths["path_found"]].empty


def test_type_preserving_null_handles_no_observed_positive_route():
    edges = pd.DataFrame([
        ("Metab__HMDB0000067_c", 1, "A"),
        ("B", 1, "E2F1"),
    ], columns=["source", "sign", "target"])
    result = typed_metabolite_tf_null(
        edges,
        variant="synthetic",
        cholesterol_nodes=["Metab__HMDB0000067_c"],
        target_tfs=["E2F1"],
        n_replicates=0,
    )
    assert pd.isna(result.loc[0, "observed_positive_shortest_length"])
