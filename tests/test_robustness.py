"""Checks for the stage-05 worker, explicit inputs and control bookkeeping."""

from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import corneto_worker  # noqa: E402
from activities import one_sample_moderated_t  # noqa: E402
from inference import make_corneto_data, sample_edges, solve_carnival  # noqa: E402

EDGES = [("A", 1, "B"), ("B", 1, "C"), ("X", 1, "C"), ("C", -1, "D")]


def _stage05_script():
    import importlib.util

    spec = importlib.util.spec_from_file_location("stage05_script", ROOT / "scripts/analysis/06_corneto_robustness.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _stage06_script():
    import importlib.util

    spec = importlib.util.spec_from_file_location("stage06_script", ROOT / "scripts/analysis/08_metabolic_context.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_explicit_inputs_replace_sources_and_may_be_internal_or_excluded():
    _, inputs, _ = make_corneto_data(EDGES, {"s": {"D": 1.0}}, inputs={"B": -1.0, "HMGCR": -1.0})
    assert inputs == ["B"]  # HMGCR is not in the PKN and is dropped
    _, default_inputs, _ = make_corneto_data(EDGES, {"s": {"D": 1.0}})
    assert default_inputs == ["A", "X"]
    _, excluded, _ = make_corneto_data([("HMGCR", 1, "D")], {"s": {"D": 1.0}}, inputs={"HMGCR": -1.0})
    assert excluded == ["HMGCR"]


def test_signed_flux_reports_hmr_even_when_it_is_not_in_the_canonical_steps(tmp_path):
    stage06 = _stage06_script()
    stage06.PKN = __import__("pathlib").Path("results/ionescu_corneto/03_pkn")
    pkn_dir = tmp_path / "results/ionescu_corneto/03_pkn/cosmos/primary"
    pkn_dir.mkdir(parents=True)
    pd.DataFrame([("Gene1916__HDLBP", 1, "Metab__cholesterol")],
                 columns=["source", "sign", "target"]).to_csv(pkn_dir / "pkn.tsv", sep="\t", index=False)
    edge_dir = tmp_path / "results/ionescu_corneto/04_inference/cosmos_multi/lambda_0.1"
    edge_dir.mkdir(parents=True)
    pd.DataFrame([("Gene1916__HDLBP", "Metab__cholesterol")], columns=["source", "target"]).to_csv(
        edge_dir / "edges.tsv", sep="\t", index=False)
    fluxes = pd.DataFrame({
        "condition": ["Ctrl", "PMS"], "reaction_number": [1916, 1916],
        "reaction": ["HMR_1916", "HMR_1916"], "flux": [0.01, -0.01],
    })
    canonical = pd.DataFrame(columns=["recon_id", "reaction_name", "reaction_node"])

    result = stage06._signed_vs_flux(tmp_path, fluxes, canonical)

    row = result.loc[result["reaction"].eq("HMR_1916")].iloc[0]
    assert row["present_in_cosmos"] and row["used_in_stage04_joint"]


def test_known_input_sign_is_enforced():
    consistent = solve_carnival(EDGES, {"s": {"D": 2.0}}, 0.01, time_limit=30, threads=1, inputs={"B": -1.0})
    assert consistent["accepted"] and consistent["weighted_fitting_error"] == pytest.approx(0.0)
    assert consistent["states"]["s"]["B"] == -1 and consistent["states"]["s"]["D"] == 1
    # B = -1 cannot produce D down through C -| D: the fit must pay the full error.
    impossible = solve_carnival(EDGES, {"s": {"D": -2.0}}, 0.01, time_limit=30, threads=1, inputs={"B": -1.0})
    assert impossible["weighted_fitting_error"] == pytest.approx(2.0)


def test_sample_edges_uses_every_sample_column():
    edges = [("A", 1, "T1"), ("A", 1, "T2"), ("B", 1, "T1"), ("B", 1, "T2")]
    result = solve_carnival(edges, {"s1": {"T1": 1.0}, "s2": {"T2": 1.0}}, 0.1, time_limit=30, threads=1)
    alternatives = sample_edges(result["problem"], result["method"], time_limit=30, threads=1, max_samples=3)
    assert list(alternatives.columns) == ["alternative", "source", "target", "edge_value", "n_samples_selected"]
    if not alternatives.empty:
        targets = set(alternatives["target"])
        assert {"T1", "T2"} <= targets  # T2 is only used by the second sample


def test_worker_runs_in_subprocess_and_resumes_only_matching_tasks(tmp_path, monkeypatch):
    pkn = tmp_path / "pkn.tsv"
    pd.DataFrame(EDGES, columns=["source", "sign", "target"]).to_csv(pkn, sep="\t", index=False)
    task = corneto_worker.make_task(label="toy", pkn_path=pkn, measurements={"s": {"D": -1.0}},
                                    lambda_reg=0.01, time_limit=30, root=tmp_path)
    result = corneto_worker.solve_in_subprocess(task, tmp_path / "solve", root=tmp_path)
    assert result["accepted"] and json.loads((tmp_path / "solve/result.json").read_text())["task_signature"] == task["signature"]
    assert {(r.source, r.target) for r in result["selected_edges"].itertuples()} >= {("C", "D")}

    def fail(*_args, **_kwargs):
        raise AssertionError("a matching finished task must not be re-solved")

    monkeypatch.setattr(subprocess, "run", fail)
    resumed = corneto_worker.solve_in_subprocess(task, tmp_path / "solve", resume=True, root=tmp_path)
    assert resumed["weighted_fitting_error"] == result["weighted_fitting_error"]
    changed = corneto_worker.make_task(label="toy", pkn_path=pkn, measurements={"s": {"D": 1.0}},
                                       lambda_reg=0.01, time_limit=30, root=tmp_path)
    with pytest.raises(AssertionError):
        corneto_worker.solve_in_subprocess(changed, tmp_path / "solve", resume=True, root=tmp_path)


def test_worker_refuses_a_pkn_changed_after_the_task_was_written(tmp_path):
    pkn = tmp_path / "pkn.tsv"
    pd.DataFrame(EDGES, columns=["source", "sign", "target"]).to_csv(pkn, sep="\t", index=False)
    task = corneto_worker.make_task(label="toy", pkn_path=pkn, measurements={"s": {"D": 1.0}},
                                    lambda_reg=0.01, time_limit=30, root=tmp_path)
    pkn.write_text("source\tsign\ttarget\nA\t1\tD\n")
    result = corneto_worker.solve_in_subprocess(task, tmp_path / "solve", root=tmp_path)
    assert not result["accepted"] and result["status"].startswith("worker_failed")


def test_one_sample_moderated_t_matches_formula_and_sign():
    rng = np.random.default_rng(1)
    differences = pd.DataFrame(rng.normal(0, 1, size=(200, 7)))
    differences.loc[0] = [-4.1, -3.9, -4.0, -4.2, -3.8, -4.0, -4.0]
    result, prior = one_sample_moderated_t(differences)
    row = result.loc[0]
    posterior = (prior["d0"] * prior["s0_squared"] + 6 * row["s2"]) / (prior["d0"] + 6)
    assert row["moderated_t"] == pytest.approx(row["mean_difference"] / np.sqrt(posterior / 7))
    assert row["moderated_t"] < -5
    with pytest.raises(ValueError):
        one_sample_moderated_t(differences.iloc[:, :2])


def test_core_network_uses_one_denominator_for_every_edge(tmp_path):
    stage05 = _stage05_script()
    frequencies = pd.DataFrame([
        {"run": "omnipath", "source": "A", "target": "B", "n_present": 4, "n_networks": 4, "frequency": 1.0},
    ])
    loo = [pd.DataFrame({"source": ["A", "C"], "target": ["B", "D"]}).assign(run="omnipath", loo="C1"),
           pd.DataFrame({"source": ["A"], "target": ["B"]}).assign(run="omnipath", loo="C2")]
    core = stage05._core_network(frequencies, loo, tmp_path).set_index(["source", "target"])
    assert core.loc[("A", "B"), "frequency"] == pytest.approx(6 / 6)
    assert core.loc[("C", "D"), "frequency"] == pytest.approx(1 / 6)  # not 1/1
    assert not core.loc[("C", "D"), "core"]


def test_comparison_of_unsolved_result_is_missing_not_zero():
    stage05 = _stage05_script()
    failed = corneto_worker.failed_result({"lambda_reg": 0.1}, "worker_failed_exit_1")
    assert all(np.isnan(value) for value in stage05._comparison(failed, {("A", "B")}).values())
    columns = stage05._fit_columns(failed)
    assert not columns["accepted"] and np.isnan(columns["n_edges"])


def test_stage05_defaults_leave_cosmos_controls_off():
    stage05 = _stage05_script()
    args = stage05._parse_args([])
    assert not args.cosmos and args.parts == list(stage05.PARTS)
    assert stage05._parse_args(["--parts", "tight_gap", "sampling"]).parts == ["sampling", "tight_gap"]


def test_context_sv_control_records_a_missing_path_instead_of_crashing(tmp_path, monkeypatch):
    """Pruned PKN without a directed HMGCR -> cholesterol path: no solve, status recorded."""
    module = _stage05_script()
    pkn_dir = tmp_path / "results/ionescu_corneto/06_validation/metabolic_context/context_cosmos"
    pkn_dir.mkdir(parents=True)
    pd.DataFrame([("HMGCR", 1, "Gene1__X"), ("Y", 1, "Metab__HMDB0000067_c")],
                 columns=["source", "sign", "target"]).to_csv(pkn_dir / "pkn.tsv", sep="\t", index=False)
    calls = []
    monkeypatch.setattr(module, "_stage04_measurements", lambda root, run: {"PMS_vs_Ctrl": {"Y": 1.0}})
    monkeypatch.setattr(module, "_sv_measurement",
                        lambda root: ("Metab__HMDB0000067_c", pd.DataFrame({"moderated_t": [-11.2]})))
    monkeypatch.setattr(module, "_fit_columns", lambda result: {"status": "optimal"})

    def fake_solve(root, out, name, *args, **kwargs):
        calls.append(name)
        return {"accepted": True}

    monkeypatch.setattr(module, "_solve", fake_solve)
    table = module._context_cosmos(tmp_path, tmp_path / "out", resume=False)
    assert calls == ["context_cosmos/joint"]
    sv = table[table["scope"].eq("context_sv")].iloc[0]
    assert sv["status"] == "no_directed_path_from_HMGCR"
    assert not sv["accepted"]
