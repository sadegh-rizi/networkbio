"""Unit tests for src/inputs.py on synthetic values (no data files needed)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from inputs import _parse_factors, parse_mw_sample_id  # noqa: E402


def test_parse_lipidomics_sample_with_sv():
    assert parse_mw_sample_id("AMC A SV P1") == {"name_prefix": "AMC", "line_letter": "A", "replicate": "P1"}


def test_parse_pma_prefix_is_kept_as_deposited():
    assert parse_mw_sample_id("PMA B P2")["name_prefix"] == "PMA"


def test_parse_medium_sample_and_blank():
    assert parse_mw_sample_id("PMS D S9")["line_letter"] == "D"
    assert parse_mw_sample_id("NC 7") == {"name_prefix": "NC", "line_letter": None, "replicate": "7"}


def test_run_id_without_line_returns_none():
    assert parse_mw_sample_id("DN47-01")["line_letter"] is None


def test_parse_factors():
    assert _parse_factors("Disease status:AMC | Treatment:SV") == {"Disease status": "AMC", "Treatment": "SV"}
