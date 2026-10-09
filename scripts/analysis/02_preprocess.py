"""Stage 01: approved RNA-seq/metabolomics preprocessing (run by the user in WSL)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from labelling import annotate_features, whitelist_report  # noqa: E402
from preprocess import (  # noqa: E402
    align_measurements, average_lines, filter_rnaseq_counts, log2_cpm,
    log2_peak_area, normalise_total, plot_layer_qc, pool_isotopologues,
    reference_to_blanks, reuse_stage, sha256, size_factor_normalise,
    stage_provenance, write_stage_provenance,
)

PLAN = "doc/decisions/2026-10-06-stage01-02-preprocessing-tf-activities-plan.md"
EXPECTED_RNA_GENES = 14490
STUDIES = {"ST003331": 7, "ST003332": 7, "ST003328": 14}
BASELINE = {"C1": "Ctrl", "C2": "Ctrl", "C3": "Ctrl", "C4": "PMS", "C5": "PMS", "C6": "PMS"}
PARAMETERS = {
    "unit": "cell line; study-specific identities, no cross-study pairing",
    "rnaseq": {"mapping": "HGNC Approved Ensembl IDs to symbols; exclude IDs with multiple distinct Approved symbols; sum counts for shared symbols",
               "cpm_min": 1, "min_libraries": 3, "filter_denominator": "original unfiltered library sizes",
               "filter_order": "map/collapse then filter; independently check pre-mapping gene count",
                "normalisation": "log2(CPM+1), library sizes recomputed after symbol filtering",
                "sensitivity": "DESeq2 median-of-ratios size factors on filtered symbol counts",
               "expected_pre_mapping_genes": EXPECTED_RNA_GENES, "minimum_mapping_fraction": 0.9,
               "batch_correction": None, "baseline_groups": BASELINE},
    "metabolomics": {"pool_keys": ["analysis_id", "kegg_id"], "pool_sum_min_count": 1,
                      "fractional_labelling": "labelled sum / deposited isotopologue sum, uncorrected",
                      "ST003331_primary": "per-analysis sample total scaled to median total, then log2",
                      "ST003331_variant": "log2 of deposited isotopologue sum without per-total scaling",
                      "ST003332_primary": "log2 peak area, then subtract available blank log2 mean",
                      "ST003332_variant": "per-total sample scaling including blanks, then log2 and blank reference",
                      "ST003328_primary": "log2 peak area without per-total scaling",
                      "ST003328_variant": "per-total sample scaling across all 42 samples, then log2",
                      "ST003328_conditions": ["SV", "untreated"], "imputation": None,
                      "line_aggregation": "mean log2; minimum 2 of 3 non-missing replicates",
                      "max_missing_line_fraction": 0.1, "excluded_layer": "ST003330"},
    "labelling": {"plain_rows": "registered M+0 interpretation; R0 found no total-pool statement",
                   "classes": ["isotopologue_sum_partial", "m0_unconfounded", "m0_confounded"],
                   "natural_abundance_correction": False},
    "qc": {"pca": "complete features, centred not scaled; numpy SVD",
           "ST003332_qc": "normalised log2 peak areas including blanks, before blank subtraction",
           "correlation": "Pearson, pairwise complete log2 features, minimum 2",
           "interpretation": "exploratory only; no statistical tests"},
    "seed": None,
}


def _write_commented_tsv(path: Path, comment: str, table: pd.DataFrame) -> None:
    with path.open("w") as handle:
        handle.write(f"# {comment}\n")
        table.to_csv(handle, sep="\t")


def _write_metabolomics_variant(
    values: pd.DataFrame,
    features: pd.DataFrame,
    samples: pd.DataFrame,
    destination: Path,
    study: str,
    suffix: str,
    include_treatment: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Write one replicate/line variant and return values, sheet and counts."""
    replicate = reference_to_blanks(values, samples) if study == "ST003332" else values
    line_values, line_samples, replicate_counts = average_lines(replicate, samples, include_treatment)
    destination.mkdir(parents=True, exist_ok=True)
    features.to_csv(destination / "features.tsv", sep="\t", index=False)
    replicate.to_csv(destination / f"replicate_{suffix}.tsv", sep="\t")
    line_values.to_csv(destination / f"line_{suffix}.tsv", sep="\t")
    line_samples.to_csv(destination / "line_samples.tsv", sep="\t", index=False)
    return line_values, replicate, line_samples, replicate_counts


LABELLING_CONVENTION = """# Isotopologue convention

The deposited ST003331 and ST003332 mwTab files were checked before stage-01
implementation. Their ST, SP, TR, AN and MS sections document 24-hour
[U-13C6]glucose exposure and the extraction volumes, but contain no explicit
statement that plain metabolite rows are total pools. The data use plain names
such as `D-Glucose` alongside labelled names such as `Glucose 13C6`.

R0 finding: no statement found that plain rows are total pools. D1 is therefore
registered as an interpretation, not a deposit fact: plain rows are treated as
M+0, labelled sums are called `isotopologue_sum_deposited`, and incomplete
labelled coverage is reported. For glucose-reachable metabolites, exchange
values describe changes in the unlabelled fraction rather than abundance.
"""


def main(root: Path = ROOT) -> int:
    source = root / "results/ionescu_corneto/00_inputs"
    out = root / "results/ionescu_corneto/01_preprocessing"
    upstream = source / "00_inputs.provenance.json"
    source_provenance = json.loads(upstream.read_text())  # read before any stage-00 tables
    if source_provenance["stage"] != "00_inputs":
        raise ValueError("Expected stage-00 input provenance")
    hgnc_dir = root / "data/resources/hgnc"
    candidates = [hgnc_dir / f"hgnc_complete_set_{release}.txt" for release in ("2026-10-06", "2026-07-07")]
    hgnc_path = next((p for p in candidates if p.exists()), None)
    if hgnc_path is None:
        raise FileNotFoundError("HGNC resource missing; run bash scripts/download/07_resources.sh first")
    hgnc_sidecar = hgnc_path.with_suffix(".provenance.json")
    resource = json.loads(hgnc_sidecar.read_text())
    if resource["sha256"] != sha256(hgnc_path):
        raise ValueError("HGNC checksum does not match its download provenance")
    resource["path"] = str(hgnc_path.relative_to(root))
    inputs = [upstream] + [source / "rnaseq" / f for f in ("counts_raw_baseline.tsv.gz", "samples.tsv", "library_sizes.tsv")]
    inputs += [source / "metabolomics" / st / f for st in STUDIES for f in ("features.tsv", "samples.tsv", "values.tsv.gz")]
    code = [root / p for p in ("scripts/analysis/02_preprocess.py", "src/preprocess.py", "src/plot_config.py",
                               "src/labelling.py", "scripts/download/07_resources.sh", "requirements.lock.txt", PLAN)]
    provenance = stage_provenance(root, "01_preprocessing", inputs, {"hgnc": resource}, PARAMETERS, code)
    provenance["accessions"] = source_provenance["accessions"]
    if reuse_stage(out, provenance):
        return 0
    for sub in ("rnaseq", "metabolomics", "summary", "figures"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    (out / "summary/labelling_convention.md").write_text(LABELLING_CONVENTION)
    counts_log, exclusions, errors, qc_inputs = [], [], [], []

    # RNA-seq: 6 independent libraries, one per line; group comes from stage 00.
    counts = pd.read_csv(source / "rnaseq/counts_raw_baseline.tsv.gz", sep="\t", index_col=0)
    samples = pd.read_csv(source / "rnaseq/samples.tsv", sep="\t")
    samples = samples.loc[samples["use_baseline"].eq(True)].copy()
    if samples["sample_title"].duplicated().any() or set(samples["sample_title"]) != set(BASELINE):
        raise ValueError("Expected exactly the six annotated baseline RNA-seq lines C1-C6")
    samples = samples.set_index("sample_title").loc[list(BASELINE)].reset_index()
    if samples["cell_line"].tolist() != list(BASELINE):
        raise ValueError("RNA-seq sample-to-line mapping differs from C1-C6; do not treat libraries as independent lines")
    if set(counts.columns) != set(BASELINE) or samples["group"].tolist() != list(BASELINE.values()):
        raise ValueError("RNA-seq columns/group assignments differ from the confirmed baseline")
    counts = counts.loc[:, list(BASELINE)]
    libraries = pd.read_csv(source / "rnaseq/library_sizes.tsv", sep="\t").set_index("sample")
    if not np.array_equal(counts.sum().to_numpy(), libraries.loc[counts.columns, "library_size"].to_numpy()):
        raise ValueError("Stage-00 library sizes do not match the count matrix")
    hgnc = pd.read_csv(hgnc_path, sep="\t", usecols=["ensembl_gene_id", "symbol", "status"])
    filtered_counts, gene_map, rna_summary = filter_rnaseq_counts(counts, hgnc)
    log2cpm = log2_cpm(filtered_counts)
    sizefactor_log2, size_factors = size_factor_normalise(filtered_counts)
    for metric, value in rna_summary.items():
        counts_log.append({"layer": "GSE297192", "metric": metric, "value": value})
    counts_log.append({"layer": "GSE297192", "metric": "sizefactor_genes_all_counts_positive",
                       "value": int(filtered_counts.gt(0).all(axis=1).sum())})
    counts_log.append({"layer": "GSE297192", "metric": "n_lines", "value": len(samples), "expected": 6})
    if rna_summary["genes_cpm_before_mapping"] != EXPECTED_RNA_GENES:
        errors.append(f"RNA-seq pre-mapping CPM count {rna_summary['genes_cpm_before_mapping']} != {EXPECTED_RNA_GENES}")
    if rna_summary["filtered_mapping_fraction"] < 0.9:
        errors.append(f"Approved-symbol coverage {rna_summary['filtered_mapping_fraction']:.1%} is below 90%")
    gene_map.to_csv(out / "rnaseq/gene_map.tsv", sep="\t", index=False)
    log2cpm.to_csv(out / "rnaseq/log2cpm.tsv", sep="\t")
    (out / "rnaseq/variants/sizefactor").mkdir(parents=True, exist_ok=True)
    sizefactor_log2.to_csv(out / "rnaseq/variants/sizefactor/log2norm.tsv", sep="\t")
    size_factors.to_csv(out / "rnaseq/variants/sizefactor/size_factors.tsv", sep="\t", header=True)
    samples["column"] = samples["sample_title"]
    samples["line_key"] = "GSE297192:" + samples["group"] + ":" + samples["cell_line"]
    samples["treatment"] = "untreated"
    samples["n_replicates"] = 1
    samples["library_size"] = samples["column"].map(libraries["library_size"])
    samples.to_csv(out / "rnaseq/samples.tsv", sep="\t", index=False)
    for row in gene_map.loc[~gene_map["kept"]].itertuples(index=False):
        exclusions.append({"layer": "GSE297192", "feature": row.ensembl_gene_id, "line": "", "reason": row.reason})
    qc_inputs.append((log2cpm, samples.rename(columns={"sample_title": "sample_id", "cell_line": "line_letter"}),
                      "GSE297192", "Park et al. 2025, GSE297192; n = 3 Ctrl + 3 PMS lines (one library each). log2(CPM+1)."))
    qc_inputs.append((sizefactor_log2, samples.rename(columns={"sample_title": "sample_id", "cell_line": "line_letter"}),
                      "GSE297192_sizefactor", "Park et al. 2025, GSE297192; n = 3 Ctrl + 3 PMS lines. Size-factor sensitivity."))

    # Metabolomics/lipidomics: primary and registered normalisation variants.
    labelled_features = {}
    for study, expected_lines in STUDIES.items():
        directory = source / "metabolomics" / study
        destination = out / "metabolomics" / study
        destination.mkdir(parents=True, exist_ok=True)
        features = pd.read_csv(directory / "features.tsv", sep="\t")
        samples = pd.read_csv(directory / "samples.tsv", sep="\t")
        values = pd.read_csv(directory / "values.tsv.gz", sep="\t", index_col=0)
        values = align_measurements(values, features, samples)
        counts_log.append({"layer": study, "metric": "features_input", "value": len(features)})
        if not samples["study"].eq(study).all():
            raise ValueError(f"{study}: samples from a different study")
        cells = samples.loc[samples["group"].isin(["Ctrl", "PMS"])]
        expected_treatments = {"SV", "untreated"} if study == "ST003328" else {"untreated"}
        if set(cells["treatment"]) != expected_treatments:
            raise ValueError(f"{study}: unexpected treatment coverage")
        for treatment in expected_treatments:
            line_counts = cells.loc[cells["treatment"].eq(treatment)].groupby("group")["line_key"].nunique().to_dict()
            if line_counts != {"Ctrl": 3, "PMS": 4}:
                raise ValueError(f"{study}/{treatment}: expected 3 Ctrl + 4 PMS lines, found {line_counts}")
        if study == "ST003328":
            keys = cells.groupby("treatment")["line_key"].agg(set)
            if keys["SV"] != keys["untreated"]:
                raise ValueError("Lipidomics treatment arms do not contain the same line keys")
        if study in {"ST003331", "ST003332"}:
            features = annotate_features(features, study)
        if study == "ST003331":
            values, features, fractional = pool_isotopologues(values, features)
            _write_commented_tsv(
                destination / "fractional_labelling_uncorrected.tsv",
                "No natural-abundance correction; missing isotopologues bias fractional labelling downwards.",
                fractional,
            )
            counts_log.append({"layer": study, "metric": "labelled_pools", "value": len(fractional)})
        labelled_features[study] = features

        if study == "ST003331":
            primary_values = log2_peak_area(normalise_total(values, features))
            variants = {"none": log2_peak_area(values)}
        elif study == "ST003332":
            primary_values = log2_peak_area(values)
            variants = {"pertotal_authors": log2_peak_area(normalise_total(values, features))}
        else:
            primary_values = log2_peak_area(values)
            variants = {"pertotal_all42": log2_peak_area(normalise_total(values, features))}

        suffix = "log2_vs_blank" if study == "ST003332" else "log2"
        line_values, primary_replicate, line_samples, replicate_counts = _write_metabolomics_variant(
            primary_values, features, samples, destination, study, suffix, study == "ST003328"
        )
        if study == "ST003332":
            blank_ids = samples.loc[samples["group"].eq("no_cell_blank"), "sample_id"]
            no_reference = primary_replicate.loc[:, blank_ids].isna().all(axis=1)
            counts_log.append({"layer": study, "metric": "features_without_blank_reference", "value": int(no_reference.sum())})
            counts_log.append({"layer": study, "metric": "exchange_interpretation",
                               "value": "m0_unconfounded: relative exchange proxy; m0_confounded: change in unlabelled fraction"})
            for feature in no_reference.index[no_reference]:
                exclusions.append({"layer": study, "feature": feature, "line": "all",
                                   "reason": "all_blank_measurements_missing_reference_undefined"})
        counts_log.extend([
            {"layer": study, "metric": "features_after_pooling", "value": len(features)},
            {"layer": study, "metric": "features_output", "value": len(line_values)},
            {"layer": study, "metric": "n_line_conditions", "value": len(line_samples), "expected": expected_lines},
            {"layer": study, "metric": "missing_replicate_values", "value": int(primary_replicate.isna().sum().sum())},
        ])
        if len(line_samples) != expected_lines:
            errors.append(f"{study}: found {len(line_samples)} line-conditions, expected {expected_lines}")
        if "labelling_class" in features:
            for labelling_class, n_features in features["labelling_class"].value_counts().sort_index().items():
                counts_log.append({"layer": study, "metric": f"labelling_class:{labelling_class}", "value": int(n_features)})
        for line in line_values.columns:
            low = replicate_counts[line] < 2
            counts_log.append({"layer": study, "metric": f"insufficient_replicate_fraction:{line}", "value": float(low.mean())})
            if low.mean() > 0.1:
                errors.append(f"{study}/{line}: {low.mean():.1%} of features have fewer than 2 replicates")
            for feature in low.index[low]:
                exclusions.append({"layer": study, "feature": feature, "line": line,
                                   "reason": "fewer_than_2_nonmissing_replicates",
                                   "n_observed": int(replicate_counts.loc[feature, line])})
        for sid in samples.loc[samples["group"].eq("no_cell_blank"), "sample_id"]:
            exclusions.append({"layer": study, "feature": "", "line": sid, "reason": "blank_reference_only_not_a_line"})
        caption = f"Ionescu et al. 2024, {study}; n = 3 Ctrl + 4 PMS lines, 3 replicates/line."
        if study == "ST003328":
            caption += " SV/untreated paired; SV markers open."
        elif study == "ST003332":
            caption += " Includes 3 no-cell blanks."
        qc_inputs.append((primary_values, samples, study, caption + " Primary log2 peak area before blank subtraction."))

        for variant, variant_values in variants.items():
            variant_suffix = "log2_vs_blank" if study == "ST003332" else "log2"
            variant_destination = destination / "variants" / variant
            _write_metabolomics_variant(
                variant_values, features, samples, variant_destination, study, variant_suffix, study == "ST003328"
            )
            qc_inputs.append((variant_values, samples, f"{study}_{variant}",
                              caption + f" Sensitivity variant {variant}; before blank subtraction."))

    report = whitelist_report(labelled_features)
    report_lines = ["\n## Whitelist coverage\n", "The whitelist is fixed by the revised plan; no IDs were added during implementation.\n"]
    for study in ("ST003331", "ST003332"):
        report_lines.append(f"\n### {study}\n")
        present = report.loc[(report["study"] == study) & report["present"], "kegg_id"].tolist()
        absent = report.loc[(report["study"] == study) & ~report["present"], "kegg_id"].tolist()
        report_lines.append(f"- Found: {', '.join(present) if present else 'none'}\n")
        report_lines.append(f"- Absent: {', '.join(absent) if absent else 'none'}\n")
    (out / "summary/labelling_convention.md").write_text(
        LABELLING_CONVENTION + "".join(report_lines)
    )

    pd.DataFrame(counts_log).to_csv(out / "summary/counts.tsv", sep="\t", index=False)
    pd.DataFrame(exclusions, columns=["layer", "feature", "line", "reason", "n_observed"]).to_csv(
        out / "summary/exclusions.tsv", sep="\t", index=False)
    print(pd.DataFrame(counts_log).to_string(index=False))
    if errors:
        raise ValueError("Stage-01 validation failed (see summary tables):\n" + "\n".join(errors))
    provenance["qc"] = {layer: plot_layer_qc(matrix, metadata, layer, caption, out / "figures")
                        for matrix, metadata, layer, caption in qc_inputs}
    for layer, metrics in provenance["qc"].items():
        counts_log.append({"layer": layer, "metric": "pca_complete_features", "value": metrics["pca_complete_features"]})
    pd.DataFrame(counts_log).to_csv(out / "summary/counts.tsv", sep="\t", index=False)
    write_stage_provenance(out, provenance)
    print(f"Stage 01 complete: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
