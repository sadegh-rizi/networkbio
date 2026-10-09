"""Stage 02: per-line and contrast ULM, with CollecTRI/DoRothEA sensitivity."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from activities import (  # noqa: E402
    RESOURCE_SETTINGS, activity_diagnostics, centre_genes, compare_regulons,
    load_regulon, moderated_t_statistics, plot_activity_heatmap,
    plot_regulon_comparison, run_ulm, welch_statistics,
)
from preprocess import reuse_stage, sha256, stage_provenance, write_stage_provenance  # noqa: E402

PLAN = "doc/decisions/2026-10-09-stage01-02-revision-plan.md"
PARAMETERS = {
    "primary": "ULM on per-gene centred log2(CPM+1), lines x genes",
    "primary_sign": "positive: signed targets higher in a line than the six-line gene mean",
    "contrast_primary": "limma-style moderated t on log2(CPM+1), PMS minus Ctrl",
    "contrast_sensitivities": ["mean log2 difference", "Welch t"],
    "tmin": 5, "decoupler_empty": False,
    "zero_variance_welch": "stop; no invented zero scores or silent gene removal",
    "pvalues": "two-sided ULM slope t-test, df = supplied genes - 2",
    "multiple_testing": "BH over all retained TFs within each resource and each line/contrast separately",
    "inference_selection": None, "resources": RESOURCE_SETTINGS,
    "heatmap": {"selection": "top 30 CollecTRI TFs by within-data variance",
                "scaling": None, "distance": "euclidean", "linkage": "average", "columns": "Ctrl then PMS"},
    "diagnostics": "mean absolute activity vs library size and C1-C6 order: Pearson; shared TFs per line: Spearman",
    "seed": 20261006,
    "decoupler_internal_shuffle": "fixed seed 0 in decoupler 2.2.0; permutation of features only, no stochastic ULM estimator",
    "limitations": "exploratory model outputs; 3 vs 3 has little power; disease confounded with age, sex and genetic background",
}


def _single_vector(values: pd.Series, observation: str) -> pd.DataFrame:
    return pd.DataFrame([values.to_numpy()], index=[observation], columns=values.index)


def _write_contrast(
    directory: Path,
    filename: str,
    scores: pd.DataFrame,
    pvalues: pd.DataFrame,
    padj: pd.DataFrame,
    observation: str,
) -> None:
    pd.DataFrame({"score": scores[observation], "pvalue": pvalues[observation],
                  "padj": padj[observation]}).to_csv(directory / filename, sep="\t")


def main(root: Path = ROOT) -> int:
    source = root / "results/ionescu_corneto/01_preprocessing"
    out = root / "results/ionescu_corneto/02_activities"
    upstream = source / "01_preprocessing.provenance.json"
    previous = json.loads(upstream.read_text())
    if previous["stage"] != "01_preprocessing":
        raise ValueError("Stage-01 success provenance is required before TF activity estimation")
    inputs = [
        upstream,
        source / "rnaseq/log2cpm.tsv",
        source / "rnaseq/samples.tsv",
        source / "rnaseq/variants/sizefactor/log2norm.tsv",
        source / "rnaseq/variants/sizefactor/size_factors.tsv",
    ]
    for path in inputs[1:]:
        if sha256(path) != previous["outputs"][str(path.relative_to(source))]:
            raise ValueError(f"Stage-01 output no longer matches provenance: {path}")
    log2cpm = pd.read_csv(inputs[1], sep="\t", index_col=0)
    sizefactor_log2 = pd.read_csv(inputs[3], sep="\t", index_col=0)
    samples = pd.read_csv(inputs[2], sep="\t").set_index("column")
    expected = [f"C{i}" for i in range(1, 7)]
    if samples.index.tolist() != expected or log2cpm.columns.tolist() != expected or sizefactor_log2.columns.tolist() != expected:
        raise ValueError("Stage 02 expects ordered RNA-seq lines C1-C6")
    groups = samples["group"]
    if groups.tolist() != ["Ctrl"] * 3 + ["PMS"] * 3:
        raise ValueError("RNA-seq groups differ from the confirmed baseline split")

    centred = centre_genes(log2cpm)
    sizefactor_centred = centre_genes(sizefactor_log2)
    welch = welch_statistics(log2cpm, groups)
    moderated, prior = moderated_t_statistics(log2cpm, groups)
    gene_contrast = moderated.join(welch.drop(columns="mean_log2cpm_difference"))

    resources, networks = {}, {}
    for name in RESOURCE_SETTINGS:
        net, path, metadata = load_regulon(root / "data/resources/regulons", name)
        networks[name] = net
        resources[name] = {**metadata, "path": str(path.relative_to(root))}
    serial_prior = {key: ("inf" if np.isinf(value) else value) for key, value in prior.items()}
    parameters = {**PARAMETERS, "moderated_prior": serial_prior}
    code = [root / p for p in ("scripts/analysis/03_tf_activities.py", "src/activities.py",
                               "src/preprocess.py", "src/plot_config.py", "requirements.lock.txt", PLAN)]
    provenance = stage_provenance(root, "02_activities", inputs, resources, parameters, code)
    provenance["accessions"] = {"GSE297192": "Park et al. 2025, Neuron; DOI 10.1016/j.neuron.2025.09.022"}
    if reuse_stage(out, provenance):
        return 0
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "summary").mkdir(parents=True, exist_ok=True)
    gene_contrast.to_csv(out / "gene_contrast.tsv", sep="\t", index_label="symbol")
    counts_log = [
        {"resource": "contrast", "metric": "moderated_d0", "value": serial_prior["d0"]},
        {"resource": "contrast", "metric": "moderated_s0_squared", "value": prior["s0_squared"]},
        {"resource": "contrast", "metric": "moderated_genes_prior", "value": prior["n_genes_prior"]},
    ]
    scores_by_resource, diagnostics = {}, {}
    for name, net in networks.items():
        directory = out / name
        directory.mkdir(parents=True, exist_ok=True)
        scores, pvalue, padj, summary = run_ulm(centred, net, tmin=5)
        scores.to_csv(directory / "tf_activity_per_line.tsv", sep="\t")
        pvalue.to_csv(directory / "tf_pvalue_per_line.tsv", sep="\t")
        padj.to_csv(directory / "tf_padj_per_line.tsv", sep="\t")
        summary.to_csv(directory / "regulon_summary.tsv", sep="\t")

        contrast_vectors = {
            "tf_activity_contrast.tsv": (moderated["moderated_t"], "PMS_vs_Ctrl"),
            "tf_activity_contrast_logfc.tsv": (moderated["mean_log2cpm_difference"], "PMS_vs_Ctrl"),
            "tf_activity_contrast_welch.tsv": (welch["welch_t"], "PMS_vs_Ctrl"),
        }
        contrast_summary = None
        for filename, (vector, observation) in contrast_vectors.items():
            contrast_scores, contrast_p, contrast_q, this_summary = run_ulm(
                _single_vector(vector, observation), net, tmin=5
            )
            if contrast_summary is None:
                contrast_summary = this_summary
            elif not contrast_summary.equals(this_summary):
                raise ValueError("Contrast regulon universes unexpectedly differ")
            _write_contrast(directory, filename, contrast_scores, contrast_p, contrast_q, observation)
        if not summary.equals(contrast_summary):
            raise ValueError("Per-line and contrast regulon universes unexpectedly differ")

        loo_directory = directory / "loo"
        loo_directory.mkdir(parents=True, exist_ok=True)
        for left_out in expected:
            remaining = [line for line in expected if line != left_out]
            loo_centred = centre_genes(log2cpm, remaining)
            loo_scores, _, _, _ = run_ulm(loo_centred, net, tmin=5)
            loo_scores.to_csv(loo_directory / f"without_{left_out}.tsv", sep="\t")

        sf_scores, _, _, _ = run_ulm(sizefactor_centred, net, tmin=5)
        diagnostics[name] = {
            "primary": activity_diagnostics(scores, samples["library_size"], groups),
            "sizefactor_variant": activity_diagnostics(sf_scores, samples["library_size"], groups),
        }
        counts_log.extend({"resource": name, "metric": metric, "value": int(value)} for metric, value in {
            "tfs_resource": len(summary), "edges_resource": len(net),
            "tfs_after_gene_filter": (summary["n_targets_expression"] > 0).sum(),
            "edges_after_gene_filter": summary["n_targets_expression"].sum(),
            "tfs_after_tmin": summary["kept"].sum(), "edges_after_tmin": summary["n_edges_after_tmin"].sum(),
            "minimum_retained_targets": summary.loc[summary["kept"], "n_targets_expression"].min(),
            "complex_sources": summary["is_complex"].sum(),
            "AP1_targets": summary.loc[summary.index == "AP1", "n_targets_expression"].sum(),
            "NFKB_targets": summary.loc[summary.index == "NFKB", "n_targets_expression"].sum(),
            "genes_supplied": len(log2cpm), "n_lines": len(samples),
        }.items())
        scores_by_resource[name] = scores

    diagnostics["library_size_group_means"] = samples["library_size"].groupby(groups).mean().to_dict()
    diagnostics["library_size_warning"] = (
        "Library size is partly aligned with group; at n=6 this correlation cannot separate an artefact from biology."
    )
    diagnostics["comparison"] = compare_regulons(scores_by_resource["collectri"], scores_by_resource["dorothea_ABC"])
    diagnostics["heatmap_row_order"] = plot_activity_heatmap(scores_by_resource["collectri"], groups, out / "figures")
    plot_regulon_comparison(scores_by_resource["collectri"], scores_by_resource["dorothea_ABC"], groups, out / "figures")
    pd.DataFrame(counts_log).to_csv(out / "summary/counts.tsv", sep="\t", index=False)
    (out / "summary/diagnostics.json").write_text(json.dumps(diagnostics, indent=2, allow_nan=False) + "\n")
    print(pd.DataFrame(counts_log).to_string(index=False))
    print(json.dumps(diagnostics, indent=2, allow_nan=False))
    print(PARAMETERS["limitations"])
    provenance["diagnostics"] = diagnostics
    write_stage_provenance(out, provenance)
    print(f"Stage 02 complete: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
