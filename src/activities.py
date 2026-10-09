"""Stage-02 ULM activities, resource caching and exploratory diagnostics.

Activities are model outputs relative to a six-line reference, not measured
TF abundance or evidence of causality. No TF is selected for network inference.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize, special, stats
from scipy.cluster.hierarchy import leaves_list, linkage

from plot_config import GROUP_COLORS, apply_plot_style, save_figure
from preprocess import sha256

RESOURCE_SETTINGS = {
    "collectri": {"function": "decoupler.op.collectri", "organism": "human",
                   "remove_complexes": False, "license": "academic",
                   "source_url": "https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1"},
    "collectri_nocomplex": {"function": "decoupler.op.collectri", "organism": "human",
                             "remove_complexes": True, "license": "academic",
                             "source_url": "https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1"},
    "dorothea_ABC": {"function": "decoupler.op.dorothea", "organism": "human",
                     "levels": ["A", "B", "C"], "license": "academic",
                     "confidence_divisors": {"A": 1, "B": 2, "C": 3},
                     "source_url": "https://omnipathdb.org/interactions/?genesymbols=1&datasets=dorothea&dorothea_levels=A,B,C&fields=dorothea_level&license=academic"},
}
COMPLEX_SOURCES = {"AP1", "NFKB"}


def validate_regulon(net: pd.DataFrame) -> None:
    """Reject ambiguous edges instead of silently inventing signs or aggregating."""
    if not {"source", "target", "weight"}.issubset(net.columns):
        raise ValueError("A regulon must have source, target and weight columns")
    if net.empty or net[["source", "target", "weight"]].isna().any().any():
        raise ValueError("Regulon is empty or contains missing identifiers/weights")
    if net.duplicated(["source", "target"]).any():
        raise ValueError("Duplicate source-target regulon edges need explicit resolution")
    if not np.isfinite(net["weight"].to_numpy(dtype=float)).all() or net["weight"].eq(0).any():
        raise ValueError("Regulon weights must be finite and nonzero")


def load_regulon(cache_dir: Path, name: str) -> tuple[pd.DataFrame, Path, dict]:
    """Fetch once; reuse the dated, checksum-verified resource without network access."""
    settings = RESOURCE_SETTINGS[name]
    cached = sorted(cache_dir.glob(f"{name}_human_*.tsv.gz"))
    if len(cached) > 1:
        raise ValueError(f"Multiple {name} caches found; explicitly retain the intended version")
    if cached:
        path = cached[0]
        metadata = json.loads(path.with_suffix(".provenance.json").read_text())
        if metadata["sha256"] != sha256(path) or metadata["settings"] != settings:
            raise ValueError(f"Resource checksum/settings mismatch: {path}")
        net = pd.read_csv(path, sep="\t")
        if len(net) != metadata["row_count"]:
            raise ValueError(f"Resource row-count mismatch: {path}")
    else:
        if name.startswith("collectri"):
            net = dc.op.collectri(organism="human", remove_complexes=settings["remove_complexes"])
        else:
            net = dc.op.dorothea(organism="human", levels=["A", "B", "C"])
        validate_regulon(net)
        now = datetime.now(timezone.utc)
        cache_dir.mkdir(parents=True, exist_ok=True)
        path = cache_dir / f"{name}_human_{now.date().isoformat()}.tsv.gz"
        net.to_csv(path, sep="\t", index=False, mode="x", compression={"method": "gzip", "mtime": 0})
        metadata = {"resource": name, "downloaded_utc": now.isoformat(), "settings": settings,
                    "row_count": len(net), "sha256": sha256(path),
                    "versions": {p: version(p) for p in ("decoupler", "omnipath")}}
        with path.with_suffix(".provenance.json").open("x") as handle:
            json.dump(metadata, handle, indent=2)
    validate_regulon(net)
    return net, path, metadata


def centre_genes(log2cpm: pd.DataFrame, lines: list[str] | None = None) -> pd.DataFrame:
    """Return selected lines x genes, centred on those lines only."""
    if not log2cpm.index.is_unique or not log2cpm.columns.is_unique:
        raise ValueError("Expression must have unique gene and line identifiers")
    if not np.isfinite(log2cpm.to_numpy()).all():
        raise ValueError("Expression must be finite; no RNA-seq imputation is approved")
    selected = log2cpm if lines is None else log2cpm.loc[:, lines]
    if selected.shape[1] < 3:
        raise ValueError("Gene centring needs at least three lines")
    return selected.sub(selected.mean(axis=1), axis=0).T


def welch_statistics(log2cpm: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """PMS minus Ctrl Welch t input and log2-expression effect/95% interval.

Three independent lines per group; no covariate adjustment. Undefined
statistics stop the run rather than becoming invented zero observations.
"""
    if not groups.index.is_unique or set(groups.index) != set(log2cpm.columns):
        raise ValueError("Every RNA-seq line needs exactly one group annotation")
    groups = groups.reindex(log2cpm.columns)
    if groups.value_counts().to_dict() != {"Ctrl": 3, "PMS": 3}:
        raise ValueError("The confirmed Welch contrast requires 3 Ctrl and 3 PMS lines")
    ctrl, pms = log2cpm.loc[:, groups.eq("Ctrl")], log2cpm.loc[:, groups.eq("PMS")]
    difference = pms.mean(axis=1) - ctrl.mean(axis=1)
    v_ctrl, v_pms = ctrl.var(axis=1, ddof=1) / 3, pms.var(axis=1, ddof=1) / 3
    variance = v_ctrl + v_pms
    se = np.sqrt(variance)
    tvalue = difference / se
    df = variance**2 / (v_ctrl**2 / 2 + v_pms**2 / 2)
    invalid = ~np.isfinite(tvalue) | ~np.isfinite(df)
    if invalid.any():
        raise ValueError(f"Undefined Welch statistic for {int(invalid.sum())} genes (e.g. {list(log2cpm.index[invalid][:5])}); resolve before ULM")
    interval = stats.t.ppf(0.975, df) * se
    return pd.DataFrame({"mean_log2cpm_difference": difference, "standard_error": se,
                         "welch_t": tvalue, "welch_df": df,
                         "ci95_low": difference - interval, "ci95_high": difference + interval})


def moderated_t_statistics(
    log2cpm: pd.DataFrame, groups: pd.Series
) -> tuple[pd.DataFrame, dict[str, float | int]]:
    """Compute a limma-style empirical-Bayes moderated two-group t statistic."""
    if not groups.index.is_unique or set(groups.index) != set(log2cpm.columns):
        raise ValueError("Every RNA-seq line needs exactly one group annotation")
    groups = groups.reindex(log2cpm.columns)
    if groups.value_counts().to_dict() != {"Ctrl": 3, "PMS": 3}:
        raise ValueError("The moderated contrast requires 3 Ctrl and 3 PMS lines")
    ctrl = log2cpm.loc[:, groups.eq("Ctrl")]
    pms = log2cpm.loc[:, groups.eq("PMS")]
    n_ctrl, n_pms = ctrl.shape[1], pms.shape[1]
    df_g = n_ctrl + n_pms - 2
    difference = pms.mean(axis=1) - ctrl.mean(axis=1)
    pooled_s2 = (
        ((ctrl.sub(ctrl.mean(axis=1), axis=0)) ** 2).sum(axis=1)
        + ((pms.sub(pms.mean(axis=1), axis=0)) ** 2).sum(axis=1)
    ) / df_g
    usable = pooled_s2.gt(0) & np.isfinite(pooled_s2)
    if not usable.any():
        raise ValueError("No positive gene variances are available for the moderated-t prior")

    z = np.log(pooled_s2.loc[usable].to_numpy())
    e = z - special.digamma(df_g / 2) + np.log(df_g / 2)
    e_bar = float(np.mean(e))
    v = float(np.var(e, ddof=1) - special.polygamma(1, df_g / 2)) if len(e) > 1 else float("-inf")
    if v > 0:
        root = optimize.brentq(
            lambda x: float(special.polygamma(1, x) - v),
            1e-8,
            1e8,
        )
        d0 = float(2 * root)
        s0_squared = float(np.exp(e_bar + special.digamma(root) - np.log(root)))
    else:
        d0 = float("inf")
        s0_squared = float(np.exp(e_bar))

    if np.isinf(d0):
        posterior = pd.Series(s0_squared, index=pooled_s2.index)
        moderated_df = float("inf")
    else:
        posterior = (d0 * s0_squared + df_g * pooled_s2) / (d0 + df_g)
        moderated_df = d0 + df_g
    moderated_t = difference / np.sqrt(posterior * (1 / n_ctrl + 1 / n_pms))
    result = pd.DataFrame({
        "mean_log2cpm_difference": difference,
        "pooled_s2": pooled_s2,
        "posterior_variance": posterior,
        "moderated_t": moderated_t,
        "moderated_df": moderated_df,
    })
    if not np.isfinite(result.drop(columns="moderated_df").to_numpy()).all():
        raise ValueError("Moderated-t calculation returned undefined values")
    return result, {"df_g": df_g, "d0": d0, "s0_squared": s0_squared,
                    "n_genes_prior": int(usable.sum())}


def run_ulm(
    data: pd.DataFrame, net: pd.DataFrame, tmin: int = 5
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """ULM on observations x genes; return TF x observation scores, raw p, BH p.

Keep zero-valued features (empty=False): a centred constant gene or a zero
Welch t is a valid observation, not an absent gene. This fixes the tmin
universe and the ULM degrees of freedom to the filtered expression genes.
decoupler 2.2.0 returns BH-adjusted p-values; reconstruct raw p from its
documented t distribution, df = number of supplied genes - 2.
"""
    validate_regulon(net)
    if not data.index.is_unique or not data.columns.is_unique or data.shape[1] < 3:
        raise ValueError("ULM requires unique observations/genes and at least 3 genes")
    if not np.isfinite(data.to_numpy()).all() or data.std(axis=1).eq(0).any():
        raise ValueError("ULM input contains nonfinite values or a constant observation")
    shared = net.loc[net["target"].isin(data.columns)]
    n_resource = net.groupby("source")["target"].nunique()
    n_shared = shared.groupby("source")["target"].nunique().reindex(n_resource.index, fill_value=0)
    summary = pd.DataFrame({"n_targets_resource": n_resource, "n_targets_expression": n_shared})
    summary["is_complex"] = summary.index.isin(COMPLEX_SOURCES)
    summary["kept"] = summary["n_targets_expression"] >= tmin
    summary["n_edges_after_tmin"] = summary["n_targets_expression"].where(summary["kept"], 0)
    summary.index.name = "TF"
    filtered = shared.loc[shared["source"].isin(summary.index[summary["kept"]])]
    if filtered.empty:
        raise ValueError(f"No TF has at least {tmin} targets among the filtered genes")
    scores, adjusted = dc.mt.ulm(data, filtered, tmin=tmin, empty=False, verbose=False)
    scores, adjusted = scores.loc[data.index], adjusted.loc[data.index]
    if not np.isfinite(scores.to_numpy()).all() or not np.isfinite(adjusted.to_numpy()).all():
        raise ValueError("ULM returned undefined scores/p-values; inspect expression and regulon")
    raw = pd.DataFrame(2 * stats.t.sf(np.abs(scores), df=data.shape[1] - 2), index=scores.index, columns=scores.columns)
    if set(scores.columns) != set(summary.index[summary["kept"]]):
        raise ValueError("ULM target filtering disagrees with the recorded regulon summary")
    outputs = [frame.T.rename_axis("TF") for frame in (scores, raw, adjusted)]
    return *outputs, summary


def activity_diagnostics(
    scores: pd.DataFrame, library_sizes: pd.Series, groups: pd.Series | None = None
) -> dict:
    """Descriptive correlations only (six observations); no p-value claims."""
    library_sizes = library_sizes.reindex(scores.columns)
    if library_sizes.isna().any():
        raise ValueError("Missing library-size metadata for activity diagnostics")
    mean_absolute = scores.abs().mean(axis=0)
    order = pd.Series(np.arange(1, len(mean_absolute) + 1), index=mean_absolute.index)
    result = {"mean_absolute_score": mean_absolute.to_dict(), "library_sizes": library_sizes.to_dict()}
    if groups is not None:
        groups = groups.reindex(mean_absolute.index)
        result["library_size_group_means"] = library_sizes.groupby(groups).mean().to_dict()
        result["groups"] = groups.to_dict()
    for name, covariate in (("library_size", library_sizes), ("column_order", order)):
        correlation = mean_absolute.corr(covariate) if mean_absolute.std() > 0 and covariate.std() > 0 else np.nan
        result[f"mean_absolute_score_vs_{name}_pearson"] = float(correlation) if np.isfinite(correlation) else None
    return result


def compare_regulons(collectri: pd.DataFrame, dorothea: pd.DataFrame) -> dict:
    """Compare shared TF scores per line with Spearman correlations."""
    shared = collectri.index.intersection(dorothea.index).sort_values()
    if set(collectri.columns) != set(dorothea.columns):
        raise ValueError("Resources must have the same line columns")
    correlations = {}
    for line in collectri.columns:
        a, b = collectri.loc[shared, line], dorothea.loc[shared, line]
        r = a.corr(b, method="spearman") if len(shared) >= 2 and a.nunique() > 1 and b.nunique() > 1 else np.nan
        correlations[line] = float(r) if np.isfinite(r) else None
    finite = [r for r in correlations.values() if r is not None]
    return {"shared_tfs": len(shared), "spearman_per_line": correlations,
            "median_spearman": float(np.median(finite)) if finite else None,
            "n_defined_correlations": len(finite)}


def plot_activity_heatmap(scores: pd.DataFrame, groups: pd.Series, out: Path) -> list[str]:
    """Top-30 variance TFs; raw scores, Euclidean/average-linkage row clustering."""
    apply_plot_style()
    ordered_lines = groups.sort_values(kind="stable").index
    selected = scores.var(axis=1).sort_values(ascending=False, kind="stable").index[:30]
    shown = scores.loc[selected, ordered_lines]
    if len(shown) > 1:
        shown = shown.iloc[leaves_list(linkage(shown.to_numpy(), method="average", metric="euclidean"))]
    limit = max(float(np.abs(shown.to_numpy()).max()), np.finfo(float).eps)
    fig, ax = plt.subplots(figsize=(7, max(5, len(shown) * 0.22)))
    image = ax.imshow(shown, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(range(len(shown.columns)), [f"{line} ({groups[line]})" for line in shown.columns], rotation=45, ha="right")
    ax.set_yticks(range(len(shown)), shown.index)
    ax.set(xlabel="Cell line", ylabel="TF (model output)", title="Exploratory CollecTRI activities")
    fig.colorbar(image, ax=ax, label="ULM t-score (gene-centred expression)")
    fig.text(0.02, 0.01, "Park et al. 2025, GSE297192; n = 3 Ctrl + 3 PMS lines.\n"
             "Top 30 TFs by variance in these data; no row scaling. Euclidean/average linkage.\n"
             "Positive = signed targets higher than six-line mean. No inference-node selection.", fontsize=7)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    save_figure(fig, out / "tf_activity_heatmap_collectri")
    plt.close(fig)
    return shown.index.tolist()


def plot_regulon_comparison(
    collectri: pd.DataFrame, dorothea: pd.DataFrame, groups: pd.Series, out: Path
) -> None:
    """One panel per cell line; points are shared TFs, not independent donors."""
    shared = collectri.index.intersection(dorothea.index).sort_values()
    fig, axes = plt.subplots(2, 3, figsize=(10, 7), sharex=True, sharey=True)
    limit = max(float(collectri.loc[shared].abs().max().max()), float(dorothea.loc[shared].abs().max().max()), 1) if len(shared) else 1
    for line, ax in zip(groups.index, axes.flat):
        group = groups[line]
        ax.scatter(collectri.loc[shared, line], dorothea.loc[shared, line], s=9, alpha=0.5,
                   color=GROUP_COLORS[group], marker="o" if group == "Ctrl" else "^", rasterized=True)
        ax.plot([-limit, limit], [-limit, limit], color="0.5", linewidth=0.7, linestyle="--")
        ax.set(title=f"{line} ({group})", xlim=(-limit, limit), ylim=(-limit, limit))
    fig.supxlabel("CollecTRI ULM t-score")
    fig.supylabel("DoRothEA A-C ULM t-score")
    fig.suptitle("Exploratory regulon sensitivity — shared TFs")
    fig.text(0.02, 0.01, f"Park et al. 2025, GSE297192; n = 3 Ctrl + 3 PMS lines; {len(shared)} shared TFs.\n"
             "Same gene-centred expression; gene/TF observations are not independent lines.", fontsize=7)
    fig.tight_layout(rect=(0.02, 0.09, 1, 0.96))
    save_figure(fig, out / "collectri_vs_dorothea")
    plt.close(fig)
