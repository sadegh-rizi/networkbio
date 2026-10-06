"""Stage 00_inputs: tidy sample sheets and feature tables for every layer.

One-time conversions of the deposited files; no filtering, normalisation,
imputation or replicate averaging (those are separate, later decisions).

Inputs (under data/, read only):
  Park/GEO/GSE297192/GSE297192_raw_expression.csv.gz, ..._series_matrix.txt.gz
  Ionescu/MetabolomicsWorkbench/ST003328, ST003331, ST003332

Outputs: results/ionescu_corneto/00_inputs/
  rnaseq/samples.tsv                 all 24 GEO samples, group/use flags
  rnaseq/counts_raw_baseline.tsv.gz  Ensembl x C1-C6 raw counts (as deposited)
  rnaseq/library_sizes.tsv
  metabolomics/<ST>/samples.tsv, features.tsv, values.tsv.gz
  summary/replicates_per_line.tsv    samples per layer x group x treatment x line
  summary/feature_id_coverage.tsv    how many features carry a KEGG ID (PKN mapping)
  summary/lipid_groups.tsv           ST003328 lipid classes
  summary/exclusions.tsv             samples flagged as not usable, with reasons
  00_inputs.provenance.json

Run from the repository root (WSL):
    .venv/bin/python scripts/analysis/01_build_inputs.py
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from inputs import read_geo_series_samples, read_mw_study, read_rnaseq_counts  # noqa: E402

DATA = Path("data")
OUT = Path("results/ionescu_corneto/00_inputs")
GSE = DATA / "Park/GEO/GSE297192"
MW = DATA / "Ionescu/MetabolomicsWorkbench"
MW_STUDIES = {
    "ST003328": "intracellular lipidomics (+/- simvastatin)",
    "ST003331": "intracellular metabolomics, 24 h 13C-glucose",
    "ST003332": "extracellular metabolomics (medium)",
}
# Baseline Ctrl vs PMS split of C1-C6: from the authors' code
# (DARG_PMS src/py/bulk_scripts/09_generate_figures.R), not from GEO.
BASELINE_GROUP = {"C1": "Ctrl", "C2": "Ctrl", "C3": "Ctrl", "C4": "PMS", "C5": "PMS", "C6": "PMS"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    for sub in ("rnaseq", "metabolomics", "summary"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    inputs_used: list[Path] = []
    exclusions: list[dict] = []
    rep_rows: list[pd.DataFrame] = []
    coverage: list[dict] = []

    # ---------------------------------------------------------------- RNA-seq
    series = GSE / "GSE297192_series_matrix.txt.gz"
    counts_path = GSE / "GSE297192_raw_expression.csv.gz"
    inputs_used += [series, counts_path]
    geo = read_geo_series_samples(series)
    counts = read_rnaseq_counts(counts_path)

    geo["in_count_matrix"] = geo["sample_title"].isin(counts.columns)
    geo["group"] = geo["sample_title"].map(BASELINE_GROUP)
    geo["group_source"] = geo["group"].notna().map({True: "authors code (09_generate_figures.R)", False: None})
    geo["use_baseline"] = geo["group"].notna() & geo["in_count_matrix"]
    geo.to_csv(OUT / "rnaseq/samples.tsv", sep="\t", index=False)
    for _, r in geo[~geo["in_count_matrix"]].iterrows():
        exclusions.append({"layer": "rnaseq", "sample": r["sample_title"], "reason": "listed in GEO but not a column of the raw count matrix"})
    for _, r in geo[geo["in_count_matrix"] & geo["group"].isna()].iterrows():
        exclusions.append({"layer": "rnaseq", "sample": r["sample_title"], "reason": "conditioned-medium experiment, not baseline Ctrl vs PMS"})

    baseline = counts[list(BASELINE_GROUP)]
    baseline.to_csv(OUT / "rnaseq/counts_raw_baseline.tsv.gz", sep="\t")
    lib = pd.DataFrame(
        {
            "sample": baseline.columns,
            "group": [BASELINE_GROUP[c] for c in baseline.columns],
            "library_size": baseline.sum().values,
            "genes_detected": (baseline > 0).sum().values,
        }
    )
    lib.to_csv(OUT / "rnaseq/library_sizes.tsv", sep="\t", index=False)
    rep_rows.append(
        lib.assign(layer="rnaseq_GSE297192", treatment="untreated", line_key=lib["sample"])
        .groupby(["layer", "group", "treatment", "line_key"]).size().rename("n_samples").reset_index()
    )
    coverage.append(
        {
            "layer": "rnaseq_GSE297192",
            "n_features": len(counts),
            "id_type": "Ensembl gene ID (symbol mapping needed for regulons/PKN)",
            "n_with_kegg_id": None,
            "n_unique_kegg_ids": None,
            "n_13c_isotopologue_rows": None,
        }
    )

    # ---------------------------------------------------------- Metabolomics
    for st, desc in MW_STUDIES.items():
        sdir = MW / st
        inputs_used += sorted(p for p in sdir.iterdir() if p.is_file())
        samples, features, values = read_mw_study(sdir)
        d = OUT / "metabolomics" / st
        d.mkdir(parents=True, exist_ok=True)
        samples.to_csv(d / "samples.tsv", sep="\t", index=False)
        features.to_csv(d / "features.tsv", sep="\t", index=False)
        values.to_csv(d / "values.tsv.gz", sep="\t")

        for _, r in samples[samples["group"] == "no_cell_blank"].iterrows():
            exclusions.append({"layer": st, "sample": r["sample_id"], "reason": "no-cell medium blank (background reference, not a line)"})
        rep_rows.append(
            samples[samples["group"] != "no_cell_blank"]
            .assign(layer=f"{st} {desc}")
            .groupby(["layer", "group", "treatment", "line_key"]).size().rename("n_samples").reset_index()
        )
        unlabelled = features[~features["is_13c_isotopologue"]]
        coverage.append(
            {
                "layer": f"{st} {desc}",
                "n_features": len(features),
                "id_type": "KEGG compound ID" if features["kegg_id"].notna().any() else "names only (RefMet / lipid shorthand)",
                "n_with_kegg_id": int(unlabelled["kegg_id"].notna().sum()),
                "n_unique_kegg_ids": int(unlabelled["kegg_id"].dropna().nunique()),
                "n_13c_isotopologue_rows": int(features["is_13c_isotopologue"].sum()),
                "n_missing_values": int(values.isna().sum().sum()),
                "n_zero_values": int((values == 0).sum().sum()),
            }
        )
        if st == "ST003328":
            features.groupby("lipid_group", dropna=False).size().rename("n_species").reset_index().sort_values(
                "n_species", ascending=False
            ).to_csv(OUT / "summary/lipid_groups.tsv", sep="\t", index=False)

    reps = pd.concat(rep_rows, ignore_index=True)
    reps.to_csv(OUT / "summary/replicates_per_line.tsv", sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(OUT / "summary/feature_id_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(exclusions).to_csv(OUT / "summary/exclusions.tsv", sep="\t", index=False)

    lines = reps.groupby(["layer", "group", "treatment"])["line_key"].nunique().rename("n_lines").reset_index()
    prov = {
        "stage": "00_inputs",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/analysis/01_build_inputs.py",
        "inputs": {str(p): sha256(p) for p in inputs_used},
        "accessions": {"GSE297192": "Park et al. 2025, Neuron", **{k: "Ionescu et al. 2024, Cell Stem Cell" for k in MW_STUDIES}},
        "decisions": {
            "rnaseq_baseline_groups": "C1-C3 Ctrl, C4-C6 PMS, from the authors' code (not in GEO)",
            "mw_group_labels": "AMC -> Ctrl; deposited 'PMA' sample prefix is PMS by the factors field",
            "line_keys": "study:group:letter; letters are reused across studies and are NOT matched across studies",
            "transformations": "none (values as deposited)",
        },
        "lines_per_layer": lines.to_dict(orient="records"),
        "versions": {p: version(p) for p in ("pandas", "numpy")},
        "python": platform.python_version(),
        "platform": platform.platform(),
    }
    (OUT / "00_inputs.provenance.json").write_text(json.dumps(prov, indent=2, default=str))

    print(lines.to_string(index=False))
    print()
    print(pd.DataFrame(coverage).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
