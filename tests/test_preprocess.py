"""Synthetic checks for stage-01 scientific transformations; no downloads."""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from preprocess import (  # noqa: E402
    align_measurements, average_lines, normalise_total, pca_scores,
    log2_peak_area, pool_isotopologues, preprocess_rnaseq, reference_to_blanks, size_factor_normalise,
)


def test_pool_keys_missing_values_and_fractional_labelling():
    features = pd.DataFrame({"feature_id": ["a", "b", "c", "d"], "analysis_id": ["neg", "neg", "pos", "pos"],
                             "kegg_id": ["C1"] * 4, "is_13c_isotopologue": [False, True, False, True]})
    values = pd.DataFrame({"s1": [10, 5, 7, 1], "s2": [np.nan, 4, np.nan, np.nan],
                           "s3": [2, np.nan, 5, 5]}, index=features["feature_id"])
    pools, annotation, fractional = pool_isotopologues(values, features)
    np.testing.assert_allclose(pools.loc["a"], [15, 4, 2])
    np.testing.assert_allclose(pools.loc["c"], [8, np.nan, 10], equal_nan=True)
    assert fractional.loc["a", "s1"] == pytest.approx(1 / 3)
    assert np.isnan(fractional.loc["a", "s3"])
    assert annotation["n_pool_rows"].tolist() == [2, 2]
    with pytest.raises(ValueError, match="exactly one unlabelled"):
        pool_isotopologues(values, features.iloc[[1, 2, 3]])


def test_per_analysis_totals_and_missingness():
    values = pd.DataFrame({"s1": [1, 3, 10], "s2": [2, np.nan, 100]}, index=["a", "b", "c"])
    features = pd.DataFrame({"feature_id": values.index, "analysis_id": ["neg", "neg", "pos"]})
    result = normalise_total(values, features)
    np.testing.assert_allclose(result.loc[["a", "b"]].sum(), [3, 3])
    np.testing.assert_allclose(result.loc["c"], [55, 55])
    assert np.isnan(result.loc["b", "s2"])
    assert np.isnan(values.loc["b", "s2"]) and values.loc["c", "s2"] == 100
    with pytest.raises(ValueError, match="positive or NaN"):
        normalise_total(values.fillna(0), features)


def test_line_averaging_uses_two_of_three_and_separates_treatment():
    samples = pd.DataFrame({"sample_id": list("abcdefg"), "line_key": ["st:Ctrl:A"] * 6 + [None],
                            "group": ["Ctrl"] * 6 + ["no_cell_blank"], "line_letter": ["A"] * 6 + [None],
                            "treatment": ["untreated"] * 3 + ["SV"] * 3 + ["untreated"]})
    values = pd.DataFrame([[1, 3, np.nan, 4, np.nan, np.nan, 100],
                           [2, 4, 6, 10, 12, 14, 100]], columns=samples["sample_id"], index=["f1", "f2"])
    result, sheet, counts = average_lines(values, samples, True)
    assert result.loc["f1", "Ctrl_A_untreated"] == 2
    assert np.isnan(result.loc["f1", "Ctrl_A_SV"])
    assert result.loc["f2", "Ctrl_A_SV"] == 12
    assert counts.loc["f1", "Ctrl_A_SV"] == 1
    assert len(sheet) == 2 and sheet["n_replicates"].eq(3).all()
    with pytest.raises(ValueError, match="Ambiguous line"):
        average_lines(values, samples, False)


def test_blank_reference_preserves_uptake_and_no_reference_is_missing():
    values = pd.DataFrame([[1, 3, 3, 3], [5, np.nan, np.nan, np.nan]], columns=["cell", "b1", "b2", "b3"])
    samples = pd.DataFrame({"sample_id": values.columns, "group": ["Ctrl"] + ["no_cell_blank"] * 3})
    result = reference_to_blanks(values, samples)
    assert result.loc[0, "cell"] == -2
    assert result.loc[0, "b1"] == 0
    assert result.loc[1].isna().all()


def test_rnaseq_mapping_collapse_filter_and_two_library_denominators():
    # An unmapped high-count gene must still contribute to the filter denominator.
    counts = pd.DataFrame([[10] * 6, [20] * 6, [1000000] * 6, [1] * 6, [4, 4, 0, 0, 0, 0]],
                          index=["e1", "e2", "unmapped", "low", "two"], columns=[f"C{i}" for i in range(1, 7)])
    hgnc = pd.DataFrame({"ensembl_gene_id": ["e1", "e2", "low", "two", "unmapped"],
                         "symbol": ["G", "G", "LOW", "TWO", "OLD"],
                         "status": ["Approved"] * 4 + ["Entry Withdrawn"]})
    result, mapping, counts_log = preprocess_rnaseq(counts, hgnc)
    assert result.index.tolist() == ["G"]
    np.testing.assert_allclose(result.loc["G"], np.log2(1000001))
    assert counts_log["symbols_with_multiple_ensembl_ids"] == 1
    assert counts_log["genes_cpm_before_mapping"] == 3
    assert counts_log["filtered_mapping_fraction"] == pytest.approx(2 / 3)
    assert mapping.set_index("ensembl_gene_id").loc["unmapped", "reason"] == "no_approved_symbol"


def test_rnaseq_mapping_excludes_ambiguous_approved_ids():
    counts = pd.DataFrame(
        [[10] * 6, [20] * 6, [30] * 6],
        index=["stable", "ambiguous", "unmapped"],
        columns=[f"C{i}" for i in range(1, 7)],
    )
    hgnc = pd.DataFrame(
        {
            "ensembl_gene_id": ["stable", "ambiguous", "ambiguous"],
            "symbol": ["STABLE", "SYMBOL_A", "SYMBOL_B"],
            "status": ["Approved"] * 3,
        }
    )
    result, mapping, counts_log = preprocess_rnaseq(counts, hgnc)

    assert result.index.tolist() == ["STABLE"]
    ambiguous = mapping.set_index("ensembl_gene_id").loc["ambiguous"]
    assert ambiguous["status"] == "ambiguous_approved"
    assert ambiguous["reason"] == "ambiguous_approved_mapping"
    assert not ambiguous["kept"]
    assert counts_log["genes_ambiguous_approved_mapping"] == 1
    assert counts_log["ambiguous_approved_ensembl_ids_in_hgnc"] == 1
    assert counts_log["genes_unmapped_no_approved_symbol"] == 1


def test_size_factor_normalisation_tracks_twofold_library():
    counts = pd.DataFrame({"A": [10, 20, 30], "B": [20, 40, 60]}, index=["g1", "g2", "g3"])
    normalised, factors = size_factor_normalise(counts)
    assert factors["B"] / factors["A"] == pytest.approx(2)
    np.testing.assert_allclose(normalised["A"], normalised["B"])


def test_log2_peak_area_variant_and_blank_exchange_are_unscaled():
    peak = pd.DataFrame([[2, 4, 8, 16]], index=["f"], columns=["cell", "b1", "b2", "b3"])
    samples = pd.DataFrame({"sample_id": peak.columns, "group": ["Ctrl"] + ["no_cell_blank"] * 3})
    logged = log2_peak_area(peak)
    exchanged = reference_to_blanks(logged, samples)
    np.testing.assert_allclose(logged.loc["f"], [1, 2, 3, 4])
    assert exchanged.loc["f", "cell"] == pytest.approx(-2)


def test_sample_alignment_does_not_silently_drop_columns():
    values = pd.DataFrame([[1, 2]], index=["f"], columns=["a", "b"])
    features = pd.DataFrame({"feature_id": ["f"]})
    samples = pd.DataFrame({"sample_id": ["b", "a"], "group": ["PMS", "Ctrl"]})
    assert align_measurements(values, features, samples).columns.tolist() == ["b", "a"]
    with pytest.raises(ValueError, match="columns do not match"):
        align_measurements(values, features, samples.iloc[:1])


def test_pca_centres_features_and_uses_complete_cases():
    values = pd.DataFrame([[1, 2, 3], [2, 4, 6], [1, np.nan, 2]], columns=list("abc"))
    scores, explained, n_features = pca_scores(values)
    assert n_features == 2
    assert explained[0] == pytest.approx(1)
    np.testing.assert_allclose(scores.mean(), 0, atol=1e-14)


def test_hgnc_downloader_fallback_and_immutable_cache(tmp_path):
    """Mock HTTP locally: no network and no files in the real data directory."""
    root = Path(__file__).resolve().parents[1]
    script = tmp_path / "scripts/download/07_resources.sh"
    script.parent.mkdir(parents=True)
    script.write_bytes((root / "scripts/download/07_resources.sh").read_bytes())
    python = tmp_path / ".venv/bin/python"
    python.parent.mkdir(parents=True)
    python.write_text(f'#!/usr/bin/env bash\nexec "{sys.executable}" "$@"\n')
    python.chmod(0o755)
    curl = tmp_path / "bin/curl"
    curl.parent.mkdir()
    curl.write_text('''#!/usr/bin/env bash
set -eu
while (($#)); do
    case "$1" in
        -o) dest="$2"; shift 2 ;;
        *) url="$1"; shift ;;
    esac
done
printf '%s\\n' "$url" >> "$CALL_LOG"
    if [[ "$url" == *2026-10-06* ]]; then
    printf 'not found' > "$dest"
    printf '404'
else
    printf 'ensembl_gene_id\\tsymbol\\tstatus\\nENSG1\\tG1\\tApproved\\n' > "$dest"
    printf '200'
fi
''')
    curl.chmod(0o755)
    log = tmp_path / "curl_calls.txt"
    env = {**os.environ, "PATH": f"{curl.parent}:{os.environ['PATH']}", "CALL_LOG": str(log),
           "DATA_ROOT": str(tmp_path / "data")}
    subprocess.run(["bash", str(script)], env=env, check=True, capture_output=True, text=True)
    resource = tmp_path / "data/resources/hgnc/hgnc_complete_set_2026-07-07.txt"
    metadata = json.loads(resource.with_suffix(".provenance.json").read_text())
    assert metadata["release"] == "2026-07-07" and metadata["row_count"] == 1
    before = resource.read_bytes()
    subprocess.run(["bash", str(script)], env=env, check=True, capture_output=True, text=True)
    assert len(log.read_text().splitlines()) == 2
    assert resource.read_bytes() == before
    resource.write_bytes(before + b"changed")
    result = subprocess.run(["bash", str(script)], env=env, capture_output=True, text=True)
    assert result.returncode != 0 and "checksum mismatch" in result.stderr


@pytest.mark.parametrize("failed_step", ["pytest", "02_preprocess.py", "03_tf_activities.py", "none"])
def test_stage_runner_stops_at_first_failure(tmp_path, failed_step):
    root = Path(__file__).resolve().parents[1]
    script = tmp_path / "scripts/setup/run_stage01_02.sh"
    script.parent.mkdir(parents=True)
    script.write_bytes((root / "scripts/setup/run_stage01_02.sh").read_bytes())
    python = tmp_path / ".venv/bin/python"
    python.parent.mkdir(parents=True)
    python.write_text('''#!/usr/bin/env bash
printf '%s\\n' "$*" >> "$CALL_LOG"
if [[ "$*" == *"$FAILED_STEP"* ]]; then exit 17; fi
''')
    python.chmod(0o755)
    log = tmp_path / "calls.txt"
    result = subprocess.run(["bash", str(script)], capture_output=True, text=True,
                            env={**os.environ, "CALL_LOG": str(log), "FAILED_STEP": failed_step})
    expected_calls = {"pytest": 1, "02_preprocess.py": 2, "03_tf_activities.py": 3, "none": 3}
    assert len(log.read_text().splitlines()) == expected_calls[failed_step]
    assert result.returncode == (0 if failed_step == "none" else 17)
    assert ("STAGES 01-02 FINISHED" in result.stdout) == (failed_step == "none")
