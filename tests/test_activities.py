"""Synthetic activity tests, including a local end-to-end stage-01/02 fixture."""

import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import special, stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import activities  # noqa: E402
import plot_config  # noqa: E402
from activities import (  # noqa: E402
    centre_genes, load_regulon, moderated_t_statistics, run_ulm, welch_statistics,
)
from labelling import TRACED_METABOLITES  # noqa: E402
from preprocess import sha256  # noqa: E402


def toy_regulon() -> pd.DataFrame:
    return pd.DataFrame({"source": ["TF_A"] * 10 + ["TF_B"] * 5 + ["TOO_SMALL"] * 4,
                         "target": [f"g{i}" for i in range(10)] + [f"g{i}" for i in range(5)] + ["g0", "g1", "g2", "outside"],
                         "weight": [1] * 5 + [-1] * 5 + [1, -1, 1, -1, 1] + [1] * 4})


def test_ulm_sign_orientation_target_filter_and_raw_vs_adjusted_pvalues():
    up = np.array([2, 3, 4, 3, 2, -1, -2, -3, -4, -3], dtype=float)
    data = pd.DataFrame([up, -up], index=["up", "down"], columns=[f"g{i}" for i in range(10)])
    # A zero-valued supplied gene stays in the background, including for a contrast.
    data["constant"] = 0.0
    scores, raw, adjusted, summary = run_ulm(data, toy_regulon())
    assert scores.loc["TF_A", "up"] > 0 > scores.loc["TF_A", "down"]
    np.testing.assert_allclose(scores["up"], -scores["down"])
    assert list(scores.columns) == ["up", "down"]
    assert "TOO_SMALL" not in scores.index
    assert summary.loc["TOO_SMALL", "n_targets_expression"] == 3
    expected = 2 * stats.t.sf(abs(scores["up"]), df=9)
    np.testing.assert_allclose(raw["up"], expected)
    np.testing.assert_allclose(adjusted["up"], stats.false_discovery_control(expected), rtol=1e-6)
    assert summary.loc[scores.index, "n_targets_expression"].ge(5).all()


def test_gene_centring_uses_all_lines_and_does_not_mutate_expression():
    expression = pd.DataFrame([[1, 2, 3, 5, 7, 9], [10, 10, 10, 10, 10, 10]],
                              index=["g", "constant"], columns=[f"C{i}" for i in range(1, 7)])
    centred = centre_genes(expression)
    np.testing.assert_allclose(centred.mean(), 0, atol=1e-14)
    assert centred.loc["C6", "g"] == pytest.approx(9 - 4.5)
    assert expression.loc["g", "C1"] == 1


def test_welch_direction_matches_independent_scipy_calculation():
    expression = pd.DataFrame([[1, 2, 3, 5, 7, 9], [8, 9, 10, 1, 3, 4]],
                              index=["up", "down"], columns=[f"C{i}" for i in range(1, 7)])
    groups = pd.Series(["Ctrl"] * 3 + ["PMS"] * 3, index=expression.columns)
    result = welch_statistics(expression, groups.sample(frac=1, random_state=20261006))
    expected = stats.ttest_ind(expression.iloc[:, 3:].to_numpy(), expression.iloc[:, :3].to_numpy(), axis=1, equal_var=False)
    np.testing.assert_allclose(result["welch_t"], expected.statistic)
    np.testing.assert_allclose(result["welch_df"], expected.df)
    assert result.loc["up", "mean_log2cpm_difference"] == 5
    assert result.loc["down", "welch_t"] < 0
    with pytest.raises(ValueError, match="Undefined Welch"):
        welch_statistics(expression * 0, groups)


def test_moderated_t_equal_variances_uses_infinite_prior_df():
    rows = []
    for offset in range(10):
        rows.append([-1 + offset / 10, 1 + offset / 10, offset / 10,
                     -1 + offset / 10 + 2, 1 + offset / 10 + 2, offset / 10 + 2])
    expression = pd.DataFrame(rows, index=[f"g{i}" for i in range(10)])
    groups = pd.Series(["Ctrl"] * 3 + ["PMS"] * 3, index=expression.columns)
    result, prior = moderated_t_statistics(expression, groups)
    expected_s0 = np.exp(-special.digamma(2) + np.log(2))
    assert np.isinf(prior["d0"])
    assert prior["s0_squared"] == pytest.approx(expected_s0)
    assert result["moderated_df"].eq(np.inf).all()


def test_moderated_t_shrinks_a_low_variance_gene():
    expression = pd.DataFrame(
        [[0, 0.1, -0.1, 2, 2.1, 1.9], [0, 1, -1, 2, 3, 1], [0, 2, -2, 2, 4, 0]],
        index=["low", "mid", "high"], columns=[f"C{i}" for i in range(1, 7)], dtype=float
    )
    groups = pd.Series(["Ctrl"] * 3 + ["PMS"] * 3, index=expression.columns)
    result, _ = moderated_t_statistics(expression, groups)
    ordinary = result["mean_log2cpm_difference"] / np.sqrt(result["pooled_s2"] * (2 / 3))
    assert abs(result.loc["low", "moderated_t"]) < abs(ordinary.loc["low"])


def test_resource_fetch_then_offline_reuse_and_checksum_failure(tmp_path, monkeypatch):
    calls = []

    def fetch(**kwargs):
        calls.append(kwargs)
        return toy_regulon()

    monkeypatch.setattr(activities.dc.op, "collectri", fetch)
    first, path, metadata = load_regulon(tmp_path, "collectri")
    second, same_path, same_metadata = load_regulon(tmp_path, "collectri")
    assert calls == [{"organism": "human", "remove_complexes": False}]
    assert path == same_path and metadata == same_metadata
    pd.testing.assert_frame_equal(first, second)
    with path.open("ab") as handle:
        handle.write(b"tampered")
    with pytest.raises(ValueError, match="checksum/settings"):
        load_regulon(tmp_path, "collectri")


def test_collectri_complex_variant_is_requested_separately(tmp_path, monkeypatch):
    calls = []

    def fetch(**kwargs):
        calls.append(kwargs)
        return toy_regulon()

    monkeypatch.setattr(activities.dc.op, "collectri", fetch)
    load_regulon(tmp_path, "collectri")
    load_regulon(tmp_path, "collectri_nocomplex")
    assert calls == [
        {"organism": "human", "remove_complexes": False},
        {"organism": "human", "remove_complexes": True},
    ]


def load_script(filename: str):
    spec = importlib.util.spec_from_file_location(filename.removesuffix(".py"), ROOT / "scripts/analysis" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_stages_end_to_end_on_synthetic_tables(tmp_path, monkeypatch):
    """Exercise schema handoff, figures, provenance, cache reuse and all line counts."""
    stage01, stage02 = load_script("02_preprocess.py"), load_script("03_tf_activities.py")
    monkeypatch.setattr(stage01, "EXPECTED_RNA_GENES", 20)
    monkeypatch.setattr(plot_config, "PRESENTATION_DPI", 30)
    monkeypatch.setattr(plot_config, "MANUSCRIPT_DPI", 40)
    rng = np.random.default_rng(20261006)
    source = tmp_path / "results/ionescu_corneto/00_inputs"
    (source / "rnaseq").mkdir(parents=True)
    (source / "00_inputs.provenance.json").write_text(json.dumps({"stage": "00_inputs", "accessions": {"fixture": "synthetic; not study data"}}))
    genes, symbols = [f"e{i}" for i in range(20)], [f"g{i}" for i in range(20)]
    columns = [f"C{i}" for i in range(1, 7)]
    counts = pd.DataFrame(rng.integers(20, 500, size=(20, 6)), index=genes, columns=columns)
    counts.index.name = "ensembl_gene_id"
    counts.to_csv(source / "rnaseq/counts_raw_baseline.tsv.gz", sep="\t")
    pd.DataFrame({"sample_title": columns, "cell_line": columns, "group": ["Ctrl"] * 3 + ["PMS"] * 3,
                  "use_baseline": True}).to_csv(source / "rnaseq/samples.tsv", sep="\t", index=False)
    pd.DataFrame({"sample": columns, "library_size": counts.sum().to_numpy()}).to_csv(source / "rnaseq/library_sizes.tsv", sep="\t", index=False)
    resource = tmp_path / "data/resources/hgnc/hgnc_complete_set_2026-10-06.txt"
    resource.parent.mkdir(parents=True)
    pd.DataFrame({"ensembl_gene_id": genes, "symbol": symbols, "status": "Approved"}).to_csv(resource, sep="\t", index=False)
    resource.with_suffix(".provenance.json").write_text(json.dumps({"sha256": sha256(resource), "downloaded_utc": "synthetic fixture"}))

    for study in stage01.STUDIES:
        directory = source / "metabolomics" / study
        directory.mkdir(parents=True)
        rows = []
        for group, letters in (("Ctrl", "ABC"), ("PMS", "DEFG")):
            for letter in letters:
                for treatment in (["SV", "untreated"] if study == "ST003328" else ["untreated"]):
                    for rep in range(3):
                        rows.append({"study": study, "sample_id": f"{group}_{letter}_{treatment}_{rep}",
                                     "line_key": f"{study}:{group}:{letter}", "group": group,
                                     "line_letter": letter, "treatment": treatment})
        if study == "ST003332":
            rows.extend({"study": study, "sample_id": f"blank{i}", "group": "no_cell_blank", "treatment": "untreated"} for i in range(3))
        samples = pd.DataFrame(rows)
        feature_ids = [f"f{i}" for i in range(20)]
        kegg_ids = list(TRACED_METABOLITES) + [f"K{i}" for i in range(11)]
        metabolite_names = [TRACED_METABOLITES[kegg][0] for kegg in list(TRACED_METABOLITES)] + [f"fixture_{i}" for i in range(11)]
        features = pd.DataFrame({"feature_id": feature_ids, "analysis_id": ["neg"] * 10 + ["pos"] * 10,
                                 "kegg_id": kegg_ids, "is_13c_isotopologue": False,
                                 "metabolite_name": metabolite_names, "lipid_group": "fixture"})
        if study == "ST003331":
            features.loc[len(features)] = ["labelled", "neg", "C00026", True, "2-Oxoglutarate 13C2", "fixture"]
        values = pd.DataFrame(rng.uniform(1, 100, size=(len(features), len(samples))),
                              index=features["feature_id"], columns=samples["sample_id"])
        samples.to_csv(directory / "samples.tsv", sep="\t", index=False)
        features.to_csv(directory / "features.tsv", sep="\t", index=False)
        values.to_csv(directory / "values.tsv.gz", sep="\t")

    # Exact source snapshots in the fixture permit the real code-hash checks.
    for name in ("scripts/analysis/02_preprocess.py", "scripts/analysis/03_tf_activities.py",
                 "scripts/download/07_resources.sh", "src/preprocess.py", "src/activities.py", "src/labelling.py",
                 "src/plot_config.py", "requirements.lock.txt", stage01.PLAN, stage02.PLAN):
        destination = tmp_path / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, destination)
    monkeypatch.setattr(activities.dc.op, "collectri", lambda **kwargs: toy_regulon())
    monkeypatch.setattr(activities.dc.op, "dorothea", lambda **kwargs: toy_regulon())
    assert stage01.main(tmp_path) == 0
    assert stage02.main(tmp_path) == 0
    first = tmp_path / "results/ionescu_corneto/01_preprocessing"
    second = tmp_path / "results/ionescu_corneto/02_activities"
    assert len(pd.read_csv(first / "rnaseq/log2cpm.tsv", sep="\t", index_col=0)) == 20
    assert (first / "rnaseq/variants/sizefactor/log2norm.tsv").is_file()
    assert (first / "summary/labelling_convention.md").read_text().startswith("# Isotopologue convention")
    for study, expected in stage01.STUDIES.items():
        sheet = pd.read_csv(first / "metabolomics" / study / "line_samples.tsv", sep="\t")
        assert len(sheet) == expected and sheet["n_replicates"].eq(3).all()
        assert "no_cell_blank" not in set(sheet["group"])
    assert (first / "metabolomics/ST003331/fractional_labelling_uncorrected.tsv").is_file()
    assert (first / "metabolomics/ST003331/variants/none/line_log2.tsv").is_file()
    assert (first / "metabolomics/ST003332/variants/pertotal_authors/line_log2_vs_blank.tsv").is_file()
    assert (first / "metabolomics/ST003328/variants/pertotal_all42/line_log2.tsv").is_file()
    for name in activities.RESOURCE_SETTINGS:
        result = pd.read_csv(second / name / "tf_activity_per_line.tsv", sep="\t", index_col=0)
        assert result.shape == (2, 6)
        assert (second / name / "tf_activity_contrast.tsv").is_file()
        assert (second / name / "tf_activity_contrast_logfc.tsv").is_file()
        assert (second / name / "tf_activity_contrast_welch.tsv").is_file()
        assert len(list((second / name / "loo").glob("without_C*.tsv"))) == 6
    assert len(list(first.glob("figures/*.pdf"))) == 16
    assert len(list(second.glob("figures/*.pdf"))) == 2
    for stage in (first, second):
        provenance_path = stage / f"{stage.name}.provenance.json"
        metadata = json.loads(provenance_path.read_text())
        assert all(sha256(stage / name) == digest for name, digest in metadata["outputs"].items())
    # Fully cached reruns retain the original provenance and output bytes.
    before = sha256(second / "02_activities.provenance.json")
    assert stage01.main(tmp_path) == stage02.main(tmp_path) == 0
    assert sha256(second / "02_activities.provenance.json") == before
    with (first / "rnaseq/log2cpm.tsv").open("a") as handle:
        handle.write("changed\n")
    with pytest.raises(ValueError, match="no longer matches provenance"):
        stage02.main(tmp_path)
