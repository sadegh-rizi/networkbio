"""Stage-01 transformations and QC for the confirmed 2026-10-06 plan.

Matrices have features in rows and samples in columns. No function imputes
missing measurements or matches line identities between studies.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_config import GROUP_COLORS, OKABE_ITO, apply_plot_style, save_figure


def filter_rnaseq_counts(
    counts: pd.DataFrame, hgnc: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Map approved symbols, sum counts and filter CPM.

    The expression-filter denominator is the original, unfiltered count
    matrix, including unmapped genes. The pre-mapping count is an independent
    QC check; the analysis filter follows symbol collapse. Ensembl IDs with
    multiple distinct approved symbols are excluded rather than assigned
    arbitrarily or duplicated across symbols.
    """
    if not counts.index.is_unique or not counts.columns.is_unique:
        raise ValueError("RNA-seq identifiers must be unique before symbol mapping")
    if not np.isfinite(counts.to_numpy()).all() or (counts < 0).any().any():
        raise ValueError("RNA-seq counts must be finite and non-negative")
    library = counts.sum(axis=0)
    if (library <= 0).any():
        raise ValueError("RNA-seq contains an empty library")
    before_mapping = (counts.div(library, axis=1) * 1e6 >= 1).sum(axis=1) >= 3

    approved = hgnc.loc[hgnc["status"].eq("Approved"), ["ensembl_gene_id", "symbol"]]
    approved = approved.dropna().drop_duplicates()
    symbol_counts = approved.groupby("ensembl_gene_id")["symbol"].nunique()
    ambiguous_ids = set(symbol_counts.index[symbol_counts > 1])
    ambiguous = counts.index.to_series().isin(ambiguous_ids)
    approved = approved.loc[~approved["ensembl_gene_id"].isin(ambiguous_ids)]
    symbol = counts.index.to_series().map(approved.set_index("ensembl_gene_id")["symbol"])
    mapped = counts.loc[symbol.notna()].copy()
    mapped.index = symbol.loc[symbol.notna()].to_numpy()
    collapsed = mapped.groupby(level=0, sort=True).sum()
    keep = (collapsed.div(library, axis=1) * 1e6 >= 1).sum(axis=1) >= 3
    filtered = collapsed.loc[keep]
    filtered_library = filtered.sum(axis=0)
    if filtered.empty or (filtered_library <= 0).any():
        raise ValueError("No usable RNA-seq expression remains after mapping/filtering")
    no_approved_symbol = symbol.isna() & ~ambiguous
    gene_map = pd.DataFrame({"ensembl_gene_id": counts.index, "symbol": symbol.to_numpy()})
    gene_map["status"] = np.select(
        [ambiguous.to_numpy(), symbol.notna().to_numpy()],
        ["ambiguous_approved", "Approved"],
        default="unmapped",
    )
    gene_map["passes_cpm_before_mapping"] = before_mapping.to_numpy()
    gene_map["kept"] = gene_map["symbol"].isin(filtered.index)
    gene_map["reason"] = np.select(
        [ambiguous.to_numpy(), gene_map["kept"].to_numpy(), symbol.notna().to_numpy()],
        ["ambiguous_approved_mapping", "kept", "symbol_below_cpm_filter"],
        default="no_approved_symbol",
    )
    n_filtered = int(before_mapping.sum())
    summary = {
        "genes_input": len(counts),
        "genes_cpm_before_mapping": n_filtered,
        "genes_mapped": int(symbol.notna().sum()),
        "genes_unmapped": int(symbol.isna().sum()),
        "genes_unmapped_no_approved_symbol": int(no_approved_symbol.sum()),
        "genes_ambiguous_approved_mapping": int(ambiguous.sum()),
        "ambiguous_approved_ensembl_ids_in_hgnc": len(ambiguous_ids),
        "filtered_genes_mapped": int((before_mapping & symbol.notna()).sum()),
        "filtered_mapping_fraction": float((before_mapping & symbol.notna()).sum() / n_filtered)
        if n_filtered else 0.0,
        "symbols_before_filter": len(collapsed),
        "symbols_with_multiple_ensembl_ids": int((symbol.value_counts() > 1).sum()),
        "symbols_removed_by_filter": int((~keep).sum()),
        "symbols_output": len(filtered),
    }
    filtered.index.name = "symbol"
    return filtered, gene_map, summary


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Compute log2(CPM + 1), recomputing library sizes on supplied genes."""
    library = counts.sum(axis=0)
    if (library <= 0).any():
        raise ValueError("Cannot compute CPM for an empty filtered library")
    result = np.log2(counts.div(library, axis=1) * 1e6 + 1)
    result.index.name = "symbol"
    return result


def preprocess_rnaseq(
    counts: pd.DataFrame, hgnc: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Map/filter RNA-seq counts and return primary log2(CPM + 1) values."""
    filtered, gene_map, summary = filter_rnaseq_counts(counts, hgnc)
    return log2_cpm(filtered), gene_map, summary


def size_factor_normalise(counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """Apply DESeq2 median-of-ratios normalisation to filtered symbol counts."""
    positive = counts.gt(0).all(axis=1)
    if not positive.any():
        raise ValueError("No gene has positive counts in every library for size factors")
    geometric = np.exp(np.log(counts.loc[positive]).mean(axis=1))
    ratios = counts.loc[positive].div(geometric, axis=0)
    factors = ratios.median(axis=0)
    if not np.isfinite(factors).all() or (factors <= 0).any():
        raise ValueError("Size-factor estimation returned invalid factors")
    normalised = np.log2(counts.div(factors, axis=1) + 1)
    normalised.index.name = "symbol"
    factors.index.name = "sample"
    factors.name = "size_factor"
    return normalised, factors


def log2_peak_area(values: pd.DataFrame) -> pd.DataFrame:
    """Log-transform positive peak areas while preserving missing values."""
    observed = values.to_numpy(dtype=float)
    if not np.isfinite(observed[~np.isnan(observed)]).all() or (observed[~np.isnan(observed)] <= 0).any():
        raise ValueError("Peak areas must be positive or NaN; no pseudocount was approved")
    return np.log2(values.astype(float))


def align_measurements(
    values: pd.DataFrame, features: pd.DataFrame, samples: pd.DataFrame
) -> pd.DataFrame:
    """Require a one-to-one annotation for every measurement row and column."""
    if not values.index.is_unique or not values.columns.is_unique:
        raise ValueError("Duplicate measurement identifiers")
    if features["feature_id"].duplicated().any() or samples["sample_id"].duplicated().any():
        raise ValueError("Duplicate feature/sample annotations")
    if set(values.index) != set(features["feature_id"]):
        raise ValueError("Measurement rows do not match feature annotations")
    if set(values.columns) != set(samples["sample_id"]):
        raise ValueError("Measurement columns do not match sample annotations")
    if not samples["group"].isin(["Ctrl", "PMS", "no_cell_blank"]).all():
        raise ValueError("Unrecognised or missing group in sample sheet")
    return values.loc[features["feature_id"], samples["sample_id"]].copy()


def pool_isotopologues(
    values: pd.DataFrame, features: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Pool by (analysis_id, kegg_id); keep the unlabelled feature ID as the key.

Fractional labelling is defined only for pools with annotated labelled rows.
An all-missing labelled sum stays missing, even if the unlabelled row exists.
The member IDs are retained in the feature table for traceability.
"""
    if features[["analysis_id", "kegg_id"]].isna().any().any():
        raise ValueError("Cannot form an ST003331 pool without analysis_id and KEGG ID")
    if not pd.api.types.is_bool_dtype(features["is_13c_isotopologue"]):
        raise ValueError("is_13c_isotopologue must contain explicit boolean annotations")
    pooled, fractions, annotations = {}, {}, []
    for key, rows in features.groupby(["analysis_id", "kegg_id"], sort=False):
        unlabelled = rows.loc[~rows["is_13c_isotopologue"]]
        if len(unlabelled) != 1:
            raise ValueError(f"Cannot form pool {key}: expected exactly one unlabelled row")
        row = unlabelled.iloc[0].to_dict()
        fid = row["feature_id"]
        pool = values.loc[rows["feature_id"]].sum(axis=0, min_count=1)
        pooled[fid] = pool
        labelled_ids = rows.loc[rows["is_13c_isotopologue"], "feature_id"]
        if len(labelled_ids):
            labelled = values.loc[labelled_ids].sum(axis=0, min_count=1)
            fractions[fid] = labelled / pool.where(pool > 0)
        row["pool_member_ids"] = json.dumps(rows["feature_id"].tolist())
        row["n_pool_rows"] = len(rows)
        row["n_labelled_rows"] = len(labelled_ids)
        annotations.append(row)
    pools = pd.DataFrame.from_dict(pooled, orient="index", columns=values.columns)
    fractional = pd.DataFrame.from_dict(fractions, orient="index", columns=values.columns)
    pools.index.name = fractional.index.name = "feature_id"
    return pools, pd.DataFrame(annotations), fractional


def normalise_total(values: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    """Scale non-missing totals to their median separately within each ion mode."""
    array = values.to_numpy(dtype=float)
    observed = array[~np.isnan(array)]
    if not np.isfinite(observed).all() or (observed <= 0).any():
        raise ValueError("Peak areas must be positive or NaN; no pseudocount was approved")
    analysis = features.set_index("feature_id")["analysis_id"].reindex(values.index)
    if analysis.isna().any():
        raise ValueError("Missing analysis_id for a measurement")
    normalised = values.astype(float).copy()
    for _, ids in analysis.groupby(analysis, sort=False).groups.items():
        block = values.loc[ids]
        total = block.sum(axis=0, min_count=1)
        if total.isna().any() or (total <= 0).any():
            raise ValueError("An analysis/sample has no positive total for normalisation")
        normalised.loc[ids] = block.div(total, axis=1) * total.median()
    return normalised


def reference_to_blanks(log2values: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """Subtract the non-missing mean log2 blank value, without filtering features."""
    blanks = samples.loc[samples["group"].eq("no_cell_blank"), "sample_id"]
    if len(blanks) != 3:
        raise ValueError(f"ST003332 requires 3 no-cell blanks, found {len(blanks)}")
    return log2values.sub(log2values.loc[:, blanks].mean(axis=1), axis=0)


def average_lines(
    log2values: pd.DataFrame, samples: pd.DataFrame, include_treatment: bool
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Mean log2 replicates, requiring 2 of exactly 3 samples per line-condition.

Return values, line annotations and observed replicate counts per feature.
Counts use sample IDs, not the optional deposited replicate-label column.
"""
    if set(log2values.columns) != set(samples["sample_id"]) or samples["sample_id"].duplicated().any():
        raise ValueError("Every replicate needs exactly one sample annotation")
    cells = samples.loc[~samples["group"].eq("no_cell_blank")].copy()
    required = ["line_key", "group", "line_letter", "treatment"]
    if cells[required].isna().any().any():
        raise ValueError("A cell-conditioned sample is missing line/condition metadata")
    if not cells["group"].isin(["Ctrl", "PMS"]).all():
        raise ValueError("Unknown cell-line group")
    means, counts, sheet = {}, {}, []
    for (line_key, treatment), rows in cells.groupby(["line_key", "treatment"], sort=True):
        if len(rows) != 3 or len(rows[["group", "line_letter"]].drop_duplicates()) != 1:
            raise ValueError(f"{line_key}/{treatment}: expected 3 consistently annotated replicates")
        first = rows.iloc[0]
        name = f"{first['group']}_{first['line_letter']}"
        if include_treatment:
            name += f"_{treatment}"
        if name in means:
            raise ValueError(f"Ambiguous line column {name}; studies/conditions cannot be merged")
        block = log2values.loc[:, rows["sample_id"]]
        n = block.notna().sum(axis=1)
        means[name] = block.mean(axis=1).where(n >= 2)
        counts[name] = n
        sheet.append({"column": name, "line_key": line_key, "group": first["group"],
                      "line_letter": first["line_letter"], "treatment": treatment,
                      "n_replicates": len(rows), "sample_ids": json.dumps(rows["sample_id"].tolist())})
    if not means:
        raise ValueError("No cell lines found")
    return pd.DataFrame(means), pd.DataFrame(sheet), pd.DataFrame(counts)


def pca_scores(log2values: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray, int]:
    """PCA by SVD on complete-case, gene/feature-centred data; no scaling."""
    complete = log2values.dropna(axis=0)
    if len(complete) < 2 or complete.shape[1] < 2:
        raise ValueError("PCA needs at least 2 complete features and 2 samples")
    matrix = complete.T.to_numpy(dtype=float, copy=True)
    matrix -= matrix.mean(axis=0)
    u, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    if np.square(singular).sum() == 0:
        raise ValueError("PCA is undefined for a constant matrix")
    scores = pd.DataFrame(u[:, :2] * singular[:2], index=complete.columns, columns=["PC1", "PC2"])
    return scores, np.square(singular[:2]) / np.square(singular).sum(), len(complete)


def plot_layer_qc(
    log2values: pd.DataFrame, samples: pd.DataFrame, layer: str, caption: str, out: Path,
    figure_label: str | None = None,
) -> dict:
    """Exploratory replicate PCA and Pearson correlation (pairwise complete data)."""
    apply_plot_style()
    figure_label = figure_label or layer
    metadata = samples.set_index("sample_id").loc[log2values.columns]
    scores, explained, n_complete = pca_scores(log2values)
    fig, ax = plt.subplots(figsize=(8, 6))
    for (group, treatment), rows in metadata.groupby(["group", "treatment"], dropna=False, sort=True):
        if group not in {"Ctrl", "PMS", "no_cell_blank"}:
            raise ValueError(f"Unknown plot group {group}")
        colour = GROUP_COLORS.get(group, OKABE_ITO[7])
        marker = {"Ctrl": "o", "PMS": "^", "no_cell_blank": "s"}[group]
        points = scores.loc[rows.index]
        ax.scatter(points["PC1"], points["PC2"], marker=marker, edgecolors=colour,
                   facecolors="none" if treatment == "SV" else colour,
                   label=f"{group} / {treatment}", s=45)
        for sid, row in rows.iterrows():
            label = row["line_letter"] if pd.notna(row["line_letter"]) else "blank"
            ax.annotate(str(label), scores.loc[sid], xytext=(3, 3), textcoords="offset points", fontsize=7)
    ax.set(xlabel=f"PC1 ({explained[0]:.1%} variance)", ylabel=f"PC2 ({explained[1]:.1%} variance)",
           title=f"Exploratory PCA — {layer}")
    ax.legend(fontsize=7)
    fig.text(0.02, 0.02, f"{caption}\n{n_complete} complete features; feature-centred, no variance scaling.", fontsize=7)
    fig.tight_layout(rect=(0, 0.11, 1, 1))
    save_figure(fig, out / f"pca_{figure_label}")
    plt.close(fig)

    corr = log2values.corr(method="pearson", min_periods=2)
    labels = [f"{sid} ({metadata.loc[sid, 'group']}, {metadata.loc[sid, 'treatment']})" for sid in corr.index]
    fig, ax = plt.subplots(figsize=(max(7, len(corr) * 0.25), max(6, len(corr) * 0.25)))
    image = ax.imshow(np.ma.masked_invalid(corr), vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(corr)), labels, rotation=90, fontsize=6)
    ax.set_yticks(range(len(corr)), labels, fontsize=6)
    ax.set_title(f"Exploratory replicate correlation — {layer}")
    fig.colorbar(image, ax=ax, label="Pearson r (log2 values; pairwise complete features)")
    fig.text(0.02, 0.01, caption, fontsize=7)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save_figure(fig, out / f"replicate_correlation_{figure_label}")
    plt.close(fig)
    return {"pca_complete_features": n_complete, "pca_variance_ratio": explained.tolist(),
            "correlation_missing_cells": int(corr.isna().sum().sum())}


def sha256(path: Path) -> str:
    """Hash in chunks so provenance does not duplicate large tables in memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stage_provenance(
    root: Path, stage: str, inputs: list[Path], resources: dict, parameters: dict, code: list[Path]
) -> dict:
    """Bind inputs and settings to the exact code, including uncommitted changes."""
    def hashes(paths: list[Path]) -> dict:
        return {str(path.relative_to(root)): sha256(path) for path in paths}

    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, stderr=subprocess.DEVNULL, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).strip())
    except (FileNotFoundError, subprocess.CalledProcessError):
        commit, dirty = "uncommitted", True
    return {"stage": stage, "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "inputs": hashes(inputs), "resources": resources, "parameters": parameters,
            "code": hashes(code), "git_commit": commit, "git_dirty": dirty,
            "versions": {p: version(p) for p in ("numpy", "pandas", "scipy", "matplotlib", "decoupler", "omnipath")},
            "python": platform.python_version(), "platform": platform.platform()}


def write_stage_provenance(out: Path, provenance: dict) -> None:
    """A success sidecar is written last, after tables and figures are complete."""
    sidecar = out / f"{provenance['stage']}.provenance.json"
    provenance["outputs"] = {str(p.relative_to(out)): sha256(p)
                             for p in sorted(out.rglob("*")) if p.is_file() and p != sidecar}
    sidecar.write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")


def reuse_stage(out: Path, provenance: dict) -> bool:
    """Reuse only a complete cache with identical inputs, settings, code and versions."""
    sidecar = out / f"{provenance['stage']}.provenance.json"
    if not sidecar.exists():
        if out.exists() and any(out.iterdir()):
            raise ValueError(f"Incomplete stage directory {out}; inspect it and move it aside before rerunning")
        return False
    previous = json.loads(sidecar.read_text())
    for key in ("inputs", "resources", "parameters", "code", "versions"):
        if previous[key] != provenance[key]:
            raise ValueError(f"{out}: cached {key} differ; archive this stage directory before rerunning")
    for name, digest in previous["outputs"].items():
        path = out / name
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Cached output missing/changed: {path}; inspect before rerunning")
    print(f"Reusing complete stage: {out}")
    return True
