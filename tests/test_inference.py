"""Synthetic CORNETO checks for stage-04 construction and selection rules."""

from __future__ import annotations

import sys

import pandas as pd

ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from inference import (  # noqa: E402
    _edge_state_conflicts,
    choose_lambda,
    make_corneto_data,
    solve_carnival,
    solution_accepted,
)


EDGES = [
    ("R1", 1, "K1"),
    ("R1", 1, "K2"),
    ("K1", 1, "TF1"),
    ("K2", -1, "TF2"),
    ("R1", 1, "ENZ"),
    ("ENZ", 1, "Metab__HMDB0000067_c"),
    ("Metab__HMDB0000067_c", 1, "TF3"),
    ("R2", 1, "K3"),
    ("K3", 1, "TF1"),
    ("K1", 1, "TF2"),
    ("K2", 1, "K4"),
    ("K4", 1, "TF4"),
]
MEASUREMENTS = {
    "sample": {"TF1": 1.0, "TF2": -1.0, "TF3": 1.0, "Metab__HMDB0000067_c": 1.0}
}


def test_toy_carnival_recovers_planted_route_and_unknown_input_sign():
    result = solve_carnival(
        EDGES,
        MEASUREMENTS,
        lambda_reg=0.1,
        time_limit=60,
        threads=1,
    )
    assert result["accepted"]
    assert result["fit"]["sign_agreement_ratio"].iloc[0] == 1.0
    assert {
        (row.source, row.target)
        for row in result["selected_edges"].itertuples()
    } == {
        ("R1", "K1"), ("R1", "K2"), ("K1", "TF1"), ("K2", "TF2"),
        ("R1", "ENZ"), ("ENZ", "Metab__HMDB0000067_c"),
        ("Metab__HMDB0000067_c", "TF3"),
    }


def test_excluded_inputs_are_not_entered_and_missing_measurements_are_reported():
    edges = EDGES + [("HMGCR", 1, "K1"), ("INSIG1", 1, "K1")]
    data, inputs, missing = make_corneto_data(
        edges,
        {"sample": {"TF1": 1.0, "not_in_graph": -1.0}},
    )
    assert "HMGCR" not in inputs
    assert "INSIG1" not in inputs
    assert missing == {"sample": ["not_in_graph"]}
    assert data is not None


def test_lambda_rule_chooses_largest_point_within_ten_percent():
    results = [
        {"lambda_reg": 0.0, "accepted": True, "weighted_fitting_error": 10.0},
        {"lambda_reg": 0.01, "accepted": True, "weighted_fitting_error": 10.5},
        {"lambda_reg": 0.03, "accepted": True, "weighted_fitting_error": 11.0},
        {"lambda_reg": 0.1, "accepted": True, "weighted_fitting_error": 11.1},
        {"lambda_reg": 0.3, "accepted": False, "weighted_fitting_error": 10.0},
    ]
    assert choose_lambda(results) == 0.03
    assert choose_lambda([results[1]]) is None


def test_measured_source_node_is_an_output_not_an_input():
    edges = EDGES + [("TF5", 1, "K1")]
    _, inputs, _ = make_corneto_data(edges, {"sample": {"TF1": 1.0, "TF5": 1.0}})
    assert "TF5" not in inputs
    assert {"R1", "R2"} <= set(inputs)


def test_node_states_come_from_vertex_value_and_unreachable_measurements_are_reported():
    # U is reachable only from HMGCR, which is never an input, so CORNETO
    # prunes it: it must be reported as unreachable and kept out of the ratio.
    edges = EDGES + [("HMGCR", 1, "U")]
    measurements = {"sample": {**MEASUREMENTS["sample"], "U": 2.0}}
    result = solve_carnival(edges, measurements, lambda_reg=0.1, time_limit=60, threads=1)
    assert result["accepted"]
    fit = result["fit"].iloc[0]
    assert fit["n_unreachable"] == 1
    assert fit["n_reachable"] == 4
    assert fit["sign_agreement_ratio"] == 1.0
    assert fit["n_edge_state_conflicts"] == 0
    nodes = result["nodes"].set_index("node")
    assert nodes.loc["U", "reachability"] == "unreachable"
    assert nodes.loc["U", "state"] == 0
    assert nodes.loc["TF1", "state"] == 1 and nodes.loc["TF2", "state"] == -1
    assert nodes.loc["R1", "is_input"] and not nodes.loc["TF1", "is_input"]
    assert "HMGCR" not in result["inputs"]


def test_lambda_rule_ignores_floating_point_noise_at_zero_baseline():
    results = [
        {"lambda_reg": 0.0, "accepted": True, "weighted_fitting_error": 0.0},
        {"lambda_reg": 0.3, "accepted": True, "weighted_fitting_error": 0.0},
        {"lambda_reg": 1.0, "accepted": True, "weighted_fitting_error": 1.437731e-15},
    ]
    assert choose_lambda(results) == 1.0
    results[2]["weighted_fitting_error"] = 0.5
    assert choose_lambda(results) == 0.3


def test_rejected_time_limit_is_not_an_accepted_solution():
    assert solution_accepted("optimal", float("nan"))
    assert solution_accepted("user_limit", 0.05)
    assert not solution_accepted("user_limit", 0.05001)
    assert not solution_accepted("user_limit", float("nan"))


def test_edge_state_conflicts_check_both_endpoints():
    selected = pd.DataFrame([
        {"source": "A", "target": "B", "edge_value": 1, "pkn_sign": 1},
    ])
    assert _edge_state_conflicts(selected, {"A": 1, "B": 1}) == 0
    assert _edge_state_conflicts(selected, {"A": 1, "B": -1}) == 1


def _stage04_script():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "stage04_script", ROOT / "scripts/analysis/05_corneto_inference.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_backfill_log_parser_reads_only_finished_lambdas(tmp_path):
    log = tmp_path / "stage04.log"
    log.write_text(
        "omnipath: lambda=0\nUnreachable vertices for sample: 0\n"
        "  status=optimal accepted=True gap=0 seconds=30 error=0.000\n"
        "cosmos_tf: lambda=0.3\n  status=optimal accepted=True gap=0.0001485 seconds=564 error=-0.000\n"
        "cosmos_multi: lambda=0\nUnreachable vertices for sample: 0\n"
        "scripts/setup/run_stage04.sh: line 11: 23702 Killed\n"
    )
    records = _stage04_script()._parse_backfill_log(log)
    assert set(records) == {("omnipath", 0.0), ("cosmos_tf", 0.3)}
    assert records[("cosmos_tf", 0.3)]["mip_gap"] == 0.0001485
    assert records[("omnipath", 0.0)]["accepted"] is True


def test_resume_refuses_changed_measurements(tmp_path):
    import pandas as pd
    import pytest

    script = _stage04_script()
    (tmp_path / "inputs").mkdir()
    tfs = pd.DataFrame({"pkn": ["omnipath"], "TF": ["SP1"], "score": [7.29]})
    metabolites = pd.DataFrame({"pkn_node": ["Metab__HMDB0000067_c"], "value": [6.41]})
    tfs.to_csv(tmp_path / "inputs/measured_tfs.tsv", sep="\t", index=False)
    metabolites.to_csv(tmp_path / "inputs/measured_metabolites.tsv", sep="\t", index=False)
    script._check_resume_inputs(tmp_path, tfs, metabolites)
    with pytest.raises(ValueError, match="differs"):
        script._check_resume_inputs(tmp_path, tfs.assign(score=[7.0]), metabolites)
