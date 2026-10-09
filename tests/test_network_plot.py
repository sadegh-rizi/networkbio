"""Unit tests for src/network_plot.py on synthetic networks (no solver needed)."""

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from network_plot import edge_table, layered_positions, node_states, pkn_signs  # noqa: E402

PKN = [("A", 1, "B"), ("B", -1, "C"), ("A", 1, "D"), ("X", 1, "D")]


def _sel(rows):
    return pd.DataFrame(rows, columns=["source", "target", "edge_value"])


def test_states_follow_edge_value_and_inhibition():
    # A=+1 -> B=+1 -| C=-1 : edge values are +1 and -1 (CarnivalFlow convention)
    st = node_states(_sel([("A", "B", 1), ("B", "C", -1)]), pkn_signs(PKN))
    assert st == {"A": 1, "B": 1, "C": -1}


def test_source_state_recovered_from_inhibitory_edge():
    # B=-1 sends -1 * -1 = +1 to C, so the source state must come out as -1
    st = node_states(_sel([("B", "C", 1)]), pkn_signs(PKN))
    assert st == {"B": -1, "C": 1}


def test_inconsistent_solution_raises():
    with pytest.raises(ValueError, match="inconsistent state for D"):
        node_states(_sel([("A", "D", 1), ("X", "D", -1)]), pkn_signs(PKN))


def test_edge_not_in_pkn_raises():
    with pytest.raises(ValueError, match="not in the PKN"):
        node_states(_sel([("C", "A", 1)]), pkn_signs(PKN))


def test_conflicting_pkn_signs_raise():
    with pytest.raises(ValueError, match="both signs"):
        pkn_signs([("A", 1, "B"), ("A", -1, "B")])


def test_layers_point_forward_in_a_dag():
    pos = layered_positions(["A", "B", "C", "D", "X"], [(s, t) for s, _, t in PKN])
    assert pos["A"][0] == pos["X"][0] == 0
    assert pos["C"][0] == 2  # longest path A -> B -> C
    assert all(pos[t][0] > pos[s][0] for s, _, t in PKN)


def test_layout_terminates_on_a_cycle():
    pos = layered_positions(["A", "B"], [("A", "B"), ("B", "A")])
    assert set(pos) == {"A", "B"}


def test_edge_frequency_across_samples():
    sol = pd.DataFrame(
        [("s1", "A", "B", 1), ("s2", "A", "B", 1), ("s1", "A", "D", 1)],
        columns=["sample", "source", "target", "edge_value"],
    )
    et = edge_table(PKN, sol).set_index(["source", "target"])
    assert et.loc[("A", "B"), "frequency"] == 1.0
    assert et.loc[("A", "D"), "frequency"] == 0.5
    assert et.loc[("X", "D"), "n_samples_selected"] == 0
    assert len(et) == len(PKN)


def test_toy_constants_parse_without_importing_corneto():
    spec = importlib.util.spec_from_file_location("plot_toy", ROOT / "scripts/analysis/00b_plot_toy_network.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    const = mod.load_toy_constants()
    assert len(const["PKN"]) == 12
    assert const["METAB"] in const["MEASURED"]
    assert const["EXPECTED"] <= {(s, t) for s, _, t in const["PKN"]}
