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
