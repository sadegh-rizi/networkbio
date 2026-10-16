"""Report figures 1-5 (doc/report/main.tex) from stage-01 to stage-05 outputs.

    python scripts/report/make_report_figures.py [--sv-dir <solve directory>]

Figure 3 needs the simvastatin
positive control on the repaired metabolism-only PKN; by default it is read
from stage 05 (05_robustness/sv_control/cosmos_repair_metabolic).  Until that
run exists, ``--sv-dir`` can point at a trial solve directory, and the caption
must say so.  A summary of the solve used is written to fig3_sv_solve.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from plot_config import apply_plot_style, save_figure  # noqa: E402

OUT = ROOT / "doc/report/figures"
STAGE03 = ROOT / "results/ionescu_corneto/03_pkn"
STAGE04 = ROOT / "results/ionescu_corneto/04_inference"
STAGE05 = ROOT / "results/ionescu_corneto/05_robustness"
LAMBDAS = ("0.01", "0.03", "0.1", "0.3", "1")
RUN_LABELS = {"omnipath": "OmniPath, TFs", "cosmos_tf": "COSMOS, TFs",
              "cosmos_multi": "COSMOS, TFs + metab.", "cosmos_repair": "COSMOS repaired"}
BLUE, ORANGE, VERM = "#0072B2", "#E69F00", "#D55E00"
LIGHT_BLUE = "#D6E9F5"
CHOL = "Metab__HMDB0000067_c"
SIZES = (8, 7, 6)  # base, annotation, tick (three-size ladder)


def _style() -> None:
    apply_plot_style()
    plt.rcParams.update({"font.size": SIZES[0], "axes.titlesize": SIZES[0], "axes.labelsize": SIZES[0],
                         "xtick.labelsize": SIZES[2], "ytick.labelsize": SIZES[2], "legend.fontsize": SIZES[1]})


def _letter(ax: plt.Axes, letter: str, x: float = -0.02, y: float = 1.02) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right", va="bottom")


def _load(run: str, lam: str) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    directory = STAGE04 / run / f"lambda_{lam}"
    edges = pd.read_csv(directory / "edges.tsv", sep="\t", dtype={"source": str, "target": str})
    nodes = pd.read_csv(directory / "nodes.tsv", sep="\t", dtype={"node": str})
    return edges, nodes, json.loads((directory / "result.json").read_text())


# --------------------------------------------------------------------------
# Figure 1: workflow, contextualisation and integration levels
# --------------------------------------------------------------------------

def _toy_arrow(ax, p, q, colour, lw, sign=1, z=1):
    style = "-|>" if sign > 0 else "-["
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=8 if sign > 0 else 3,
                                 shrinkA=9, shrinkB=9, color=colour, lw=lw, zorder=z))


SIG_FILL, MET_FILL, DATA_FILL = "#FBE3C0", "#CDEBDD", "#D6E9F5"
MET_EDGE = "#009E73"


def _box(ax, x, y, w, h, text, fc, bold=None, ec="0.5", ls="-", fs=None):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc=fc, ec=ec, lw=0.7, ls=ls, zorder=2))
    if bold:
        top = y + min(h / 2, 0.37)
        ax.text(x, top - 0.07, bold, ha="center", va="top", fontsize=fs or SIZES[1], fontweight="bold", zorder=3)
        ax.text(x, top - 0.27, text, ha="center", va="top", fontsize=SIZES[2], color="0.2", zorder=3,
                linespacing=1.05)
    else:
        ax.text(x, y, text, ha="center", va="center", fontsize=SIZES[2], zorder=3, linespacing=1.1)


def _flow(ax, p, q, colour="0.45", ls="-", lw=0.8, rad=0.0):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=7, color=colour, lw=lw, ls=ls,
                                 shrinkA=0, shrinkB=0, connectionstyle=f"arc3,rad={rad}", zorder=1))


def _workflow(ax) -> None:
    ax.set_xlim(0, 10); ax.set_ylim(0, 3.35); ax.set_axis_off()
    xs = [0.95, 2.95, 4.95, 6.95, 8.95]
    w, h = 1.8, 0.74
    ys = [2.45, 1.5, 0.55]
    for x, head in zip(xs, ["Data (cell lines)", "Measured nodes", "Prior knowledge", "CORNETO method",
                            "Readouts"]):
        ax.text(x, 3.2, head, ha="center", va="center", fontsize=SIZES[1], fontweight="bold")
    _box(ax, xs[0], ys[0], w, h, "cells, medium, lipids\n3 Ctrl vs 4 PMS", DATA_FILL,
         bold="Metabolomics")
    _box(ax, xs[0], (ys[1] + ys[2]) / 2, w, ys[1] - ys[2] + h, "GSE297192\n3 Ctrl vs 3 PMS", DATA_FILL,
         bold="RNA-seq")
    nodes = [("Metabolites", "36 contrasts,\nincl. cholesterol", SIG_FILL),
             ("TF activities", "CollecTRI, top 50", SIG_FILL),
             ("Gene scores", "high / low\nexpression", MET_FILL)]
    pkns = [("COSMOS", "signalling +\nsigned metabolism", SIG_FILL),
            ("OmniPath", "signed\nsignalling", SIG_FILL),
            ("Recon3D", "stoichiometric\nmetabolism", MET_FILL)]
    for (bt, t, c), y in zip(nodes, ys):
        _box(ax, xs[1], y, w, h, t, c, bold=bt)
    for (bt, t, c), y in zip(pkns, ys):
        _box(ax, xs[2], y, w, h, t, c, bold=bt)
    _box(ax, xs[3], (ys[0] + ys[1]) / 2, w, ys[0] - ys[1] + h, "signed paths from\nstarting nodes to\nmeasured nodes,\nfewest edges",
         SIG_FILL, bold="CarnivalFlow")
    _box(ax, xs[3], ys[2], w, h, "flux with\nmass balance", MET_FILL, bold="iMAT / FBA")
    _box(ax, xs[4], ys[1], w, ys[0] - ys[2] + h,
         "cholesterol route:\nin the network?\nselected?\n\nnetwork properties:\nhubs, random\nnetworks, robustness\n\nintegration level:\nlate vs joint", "white")
    R = w / 2
    _flow(ax, (xs[0] + R, ys[0]), (xs[1] - R, ys[0]))
    _flow(ax, (xs[0] + R, ys[1]), (xs[1] - R, ys[1]))
    _flow(ax, (xs[0] + R, ys[2]), (xs[1] - R, ys[2]), colour=MET_EDGE)
    _flow(ax, (xs[1] + R, ys[0]), (xs[2] - R, ys[0]))
    _flow(ax, (xs[1] + R, ys[1]), (xs[2] - R, ys[1]))
    _flow(ax, (xs[1] + R, ys[1] + 0.2), (xs[2] - R, ys[0] - 0.2))
    _flow(ax, (xs[1] + R, ys[2]), (xs[2] - R, ys[2]), colour=MET_EDGE)
    _flow(ax, (xs[2] + R, ys[0]), (xs[3] - R, ys[0]))
    _flow(ax, (xs[2] + R, ys[1]), (xs[3] - R, ys[1]))
    _flow(ax, (xs[2] + R, ys[2]), (xs[3] - R, ys[2]), colour=MET_EDGE)
    _flow(ax, (xs[3] + R, (ys[0] + ys[1]) / 2), (xs[4] - R, (ys[0] + ys[1]) / 2))
    _flow(ax, (xs[3] + R, ys[2]), (xs[4] - R, ys[2]), colour=MET_EDGE)


def _toy(ax) -> None:
    pos = {"S1": (0, 3), "S2": (0, 1.5), "S3": (0, 0), "K1": (1.2, 3.2), "K2": (1.2, 1.9), "K3": (1.2, 0.6),
           "K4": (2.4, 2.6), "K5": (2.4, 0.9), "T1": (3.6, 3.3), "T2": (3.6, 2.2), "T3": (3.6, 1.1), "T4": (3.6, 0.0)}
    pkn = [("S1", "K1", 1), ("S1", "K2", 1), ("S2", "K2", 1), ("S2", "K3", -1), ("S3", "K3", 1), ("K1", "K4", 1),
           ("K1", "T1", 1), ("K2", "K4", 1), ("K2", "K5", -1), ("K3", "K5", 1), ("K4", "T1", 1), ("K4", "T2", 1),
           ("K5", "T3", 1), ("K5", "T4", -1), ("K3", "T4", 1), ("K2", "T3", -1)]
    selected = {("S1", "K1"), ("K1", "K4"), ("K4", "T1"), ("K4", "T2"), ("S1", "K2"), ("K2", "K5"), ("K5", "T3"),
                ("K5", "T4")}
    measured = {"T1": 1, "T2": 1, "T3": -1, "T4": 1}
    state = {"S1": 1, "K1": 1, "K2": 1, "K4": 1, "K5": -1}
    for u, v, sign in pkn:
        if (u, v) in selected:
            _toy_arrow(ax, pos[u], pos[v], "0.1", 1.3, sign, z=2)
        else:
            _toy_arrow(ax, pos[u], pos[v], "0.8", 0.6, sign)
    for n, (x, y) in pos.items():
        if n in measured:
            fc, ec = (VERM if measured[n] > 0 else BLUE), "0.1"
        elif n in state:
            fc, ec = ("#F6C9B0" if state[n] > 0 else "#BFD9EC"), "0.3"
        else:
            fc, ec = "white", "0.45"
        if n.startswith("S"):
            ec = "#CC79A7"
        ax.scatter(x, y, s=150, facecolor=fc, edgecolor=ec, linewidth=1.5 if n.startswith("S") else 0.9, zorder=3)
    ax.set_xlim(-0.5, 4.1); ax.set_ylim(-1.9, 4.2); ax.set_axis_off()
    ax.text(0, 3.5, "starting\nnodes", ha="center", va="bottom", fontsize=SIZES[2], color="#CC79A7")
    ax.text(3.6, 3.75, "measured", ha="center", va="bottom", fontsize=SIZES[2])
    handles = [Line2D([], [], marker="o", ls="", mfc=VERM, mec="0.1", ms=5, label="up in PMS"),
               Line2D([], [], marker="o", ls="", mfc=BLUE, mec="0.1", ms=5, label="down in PMS"),
               Line2D([], [], color="0.1", lw=1.3, label="selected edge"),
               Line2D([], [], color="0.8", lw=1.0, label="unused PKN edge")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(-0.02, -0.02), ncol=2, frameon=False,
              handletextpad=0.3, columnspacing=0.7, handlelength=1.3)


def _mini_graph(ax, origin, kind, scale=1.0):
    """A four-node glyph: signalling (orange, signed arrows) or metabolism (green, reaction square)."""
    x0, y0 = origin
    if kind == "sig":
        pts = [(0, 0.5), (0.55, 0.85), (0.55, 0.15), (1.1, 0.5)]
        for i, j in [(0, 1), (0, 2), (1, 3), (2, 3)]:
            ax.add_patch(FancyArrowPatch((x0 + pts[i][0] * scale, y0 + pts[i][1] * scale),
                                         (x0 + pts[j][0] * scale, y0 + pts[j][1] * scale), arrowstyle="-|>",
                                         mutation_scale=5, shrinkA=4, shrinkB=4, color="0.35", lw=0.7, zorder=3))
        for px, py in pts:
            ax.scatter(x0 + px * scale, y0 + py * scale, s=28, fc=ORANGE, ec="0.2", lw=0.5, zorder=4)
    else:
        mets = [(0, 0.8), (0, 0.2), (1.1, 0.5)]
        rx, ry = 0.55, 0.5
        for mx, my in mets[:2]:
            ax.plot([x0 + mx * scale, x0 + rx * scale], [y0 + my * scale, y0 + ry * scale], color=MET_EDGE,
                    lw=0.8, zorder=3)
        ax.add_patch(FancyArrowPatch((x0 + rx * scale, y0 + ry * scale), (x0 + 1.1 * scale, y0 + 0.5 * scale),
                                     arrowstyle="-|>", mutation_scale=5, shrinkA=2, shrinkB=4, color=MET_EDGE,
                                     lw=0.8, zorder=3))
        for mx, my in mets:
            ax.scatter(x0 + mx * scale, y0 + my * scale, s=28, fc=MET_FILL, ec=MET_EDGE, lw=0.7, zorder=4)
        ax.scatter(x0 + rx * scale, y0 + ry * scale, s=26, marker="s", fc=MET_EDGE, ec=MET_EDGE, zorder=4)


def _integration(ax) -> None:
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.2); ax.set_axis_off()
    ax.text(0.1, 5.85, "Late (this study): separate problems", fontsize=SIZES[1],
            fontweight="bold", va="center")
    _box(ax, 1.6, 4.35, 2.9, 1.8, "", SIG_FILL)
    _box(ax, 8.3, 4.35, 2.9, 1.8, "", MET_FILL)
    _mini_graph(ax, (0.9, 4.3), "sig", scale=1.3)
    _mini_graph(ax, (7.6, 4.3), "met", scale=1.3)
    ax.text(1.6, 3.65, "CarnivalFlow", ha="center", fontsize=SIZES[2])
    ax.text(8.3, 3.65, "iMAT", ha="center", fontsize=SIZES[2])
    _box(ax, 4.95, 4.35, 2.5, 1.3, "compare, or\nprune one with\nthe other", "white")
    _flow(ax, (3.05, 4.35), (3.7, 4.35))
    _flow(ax, (6.85, 4.35), (6.2, 4.35), colour=MET_EDGE)
    ax.text(0.1, 2.75, "Coupled (future): one problem with shared variables", fontsize=SIZES[1],
            fontweight="bold", va="center")
    _box(ax, 4.95, 1.25, 9.7, 2.3, "", "white", ls="--", ec="0.35")
    _box(ax, 1.75, 1.25, 2.9, 1.8, "", SIG_FILL)
    _box(ax, 8.15, 1.25, 2.9, 1.8, "", MET_FILL)
    _mini_graph(ax, (1.05, 1.2), "sig", scale=1.3)
    _mini_graph(ax, (7.45, 1.2), "met", scale=1.3)
    ax.text(1.75, 0.55, "enzyme activity", ha="center", fontsize=SIZES[2])
    ax.text(8.15, 0.55, "reaction flux", ha="center", fontsize=SIZES[2])
    ax.annotate("", xy=(8.16, 1.85), xytext=(2.48, 1.85),
                arrowprops=dict(arrowstyle="<|-|>", ls="--", color="0.15", lw=0.9, mutation_scale=7,
                                shrinkA=4, shrinkB=4), zorder=5)
    ax.text(4.95, 1.05, "shared variable:\ninactive enzyme\n= no flux", ha="center", va="center",
            fontsize=SIZES[2], color="0.15", zorder=6)


def figure1() -> None:
    """Workflow (a), contextualisation toy (b) and integration levels (c)."""
    fig = plt.figure(figsize=(7.2, 5.0))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], width_ratios=[0.8, 1.2], hspace=0.28, wspace=0.1)
    ax_a = fig.add_subplot(grid[0, :]); ax_b = fig.add_subplot(grid[1, 0]); ax_c = fig.add_subplot(grid[1, 1])
    _workflow(ax_a)
    _toy(ax_b)
    _integration(ax_c)
    ax_b.set_title("Contextualisation: fewest edges\nthat fit all measured signs", loc="left")
    ax_c.set_title("How CORNETO methods are combined", loc="left")
    _letter(ax_a, "a", x=0.0, y=0.98); _letter(ax_b, "b", x=-0.02, y=1.0); _letter(ax_c, "c", x=-0.03, y=1.0)
    save_figure(fig, OUT / "fig1_concept")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 2: TF activity and cholesterol genes
# --------------------------------------------------------------------------

FOCAL_TFS = ["SP1", "JUN", "AP1", "NFKB", "JUNB", "EGR1", "SREBF2", "E2F1", "SREBF1", "WT1"]
HYPOTHESIS_TFS = {"E2F1", "EGR1", "WT1", "SP1", "JUN", "JUNB"}
PANEL_GENES = ["HMGCR", "SQLE", "INSIG1", "LDLR", "SREBF2", "SERPINE1", "LMNB1"]


def figure2() -> None:
    activities = ROOT / "results/ionescu_corneto/02_activities"
    contrast = pd.read_csv(activities / "collectri/tf_activity_contrast.tsv", sep="\t", index_col=0)["score"]
    rank = contrast.abs().rank(ascending=False, method="min").astype(int)
    loo = pd.read_csv(STAGE04 / "inputs/tf_selection_loo.tsv", sep="\t").query("pkn == 'cosmos_tf'").set_index("TF")
    loo = loo[[c for c in loo.columns if c.startswith("without_")]]
    expression = pd.read_csv(ROOT / "results/ionescu_corneto/01_preprocessing/rnaseq/log2cpm.tsv", sep="\t", index_col=0)
    samples = pd.read_csv(ROOT / "results/ionescu_corneto/01_preprocessing/rnaseq/samples.tsv", sep="\t").set_index("column")
    groups = samples["group"].reindex(expression.columns)
    genes = pd.read_csv(activities / "gene_contrast.tsv", sep="\t", index_col=0)
    fig = plt.figure(figsize=(7.2, 2.8))
    grid = fig.add_gridspec(1, 2, width_ratios=[1, 1.25], wspace=0.35)
    axa, axb = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    for y, tf in enumerate(FOCAL_TFS):
        colour = ORANGE if tf in HYPOTHESIS_TFS else "0.45"
        if tf in loo.index:
            axa.scatter(loo.loc[tf].to_numpy(dtype=float), [y] * loo.shape[1], s=9, color=colour, alpha=0.45, lw=0, zorder=2)
        axa.scatter(contrast[tf], y, s=36, color=colour, edgecolor="0.1", lw=0.5, zorder=3)
        axa.text(8.3, y, f"{rank[tf]}", va="center", ha="right", fontsize=SIZES[2], color="0.3")
    axa.axvline(0, color="0.6", lw=0.6)
    axa.set_yticks(range(len(FOCAL_TFS)), FOCAL_TFS, fontstyle="italic"); axa.set_ylim(len(FOCAL_TFS) - 0.5, -1.4)
    axa.tick_params(axis="y", length=0)
    axa.set_xlim(-3, 8.6); axa.set_xlabel("TF activity score, PMS vs Ctrl\n(ULM on moderated t)")
    axa.text(8.3, -0.95, f"rank / {len(contrast)}", ha="right", va="center", fontsize=SIZES[2], color="0.3")
    axa.set_title("SP1 and JUN are among the most active TFs", loc="left")
    _letter(axa, "a", x=-0.18)
    handles = [Line2D([], [], marker="o", ls="", mfc=ORANGE, mec="0.1", ms=6, label="hypothesis TF"),
               Line2D([], [], marker="o", ls="", mfc="0.45", mec="0.1", ms=6, label="other TF"),
               Line2D([], [], marker="o", ls="", color="0.45", alpha=0.45, ms=3.5, label="leave one line out")]
    axa.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.05, -0.3), ncol=3, frameon=False,
               handletextpad=0.2, columnspacing=0.8)
    centred = expression.loc[PANEL_GENES].sub(expression.loc[PANEL_GENES, groups.eq("Ctrl")].mean(axis=1), axis=0)
    for x, gene in enumerate(PANEL_GENES):
        for k, group in enumerate(("Ctrl", "PMS")):
            values = centred.loc[gene, groups.eq(group)].to_numpy(dtype=float)
            offsets = np.linspace(-0.08, 0.08, len(values))
            axb.scatter(x + (k - 0.5) * 0.36 + offsets, values, s=16, marker="o" if group == "Ctrl" else "^",
                        color=BLUE if group == "Ctrl" else VERM, zorder=3, label=group if x == 0 else None)
        axb.text(x, centred.loc[gene].max() + 0.45, f"t = {genes.loc[gene, 'moderated_t']:.1f}", ha="center",
                 fontsize=SIZES[2], color="0.3")
    axb.axhline(0, color="0.6", lw=0.6)
    axb.set_xticks(range(len(PANEL_GENES)), PANEL_GENES, fontstyle="italic", rotation=35, ha="right")
    axb.set_ylabel("log2 CPM, relative to\ncontrol mean")
    axb.set_ylim(-5.3, 7.2)
    axb.legend(loc="lower left", frameon=False, ncol=2, handletextpad=0.2)
    axb.set_title("Senescence markers change; sterol genes do not rise", loc="left")
    _letter(axb, "b", x=-0.12)
    save_figure(fig, OUT / "fig2_tf_activity")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 3: Aim 5
# --------------------------------------------------------------------------

STEP_GENES = ["HMGCR", "MVK", "PMVK", "MVD", "IDI1/2", "FDPS", "FDFT1", "SQLE", "LSS", "CYP51A1",
              "TM7SF2/LBR", "MSMO1", "NSDHL", "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24"]


def _reaction_genes(path: Path) -> dict[str, set[str]]:
    table = pd.read_csv(path, sep="\t", dtype=str)
    table = table[table["type"].eq("enzyme_or_reaction")]
    return {str(r.node): {s for s in str(r.approved_symbol).split(";") if s and s != "nan"} for r in table.itertuples()}


def _genes_used(nodes: pd.DataFrame, reaction_genes: dict[str, set[str]]) -> set[str]:
    active = set(nodes.loc[nodes["state"].ne(0), "node"])
    return {g for n in active for g in reaction_genes.get(n, set())}


def _has(entry: str, used: set[str]) -> bool:
    return bool(set(entry.replace("IDI1/2", "IDI1/IDI2").split("/")) & used)


def figure3(sv_dir: Path) -> None:
    """``sv_dir``: repaired metabolism-only SV control; the primary one is read next to it."""
    genes = _reaction_genes(STAGE03 / "cosmos/variants/canonical_repair/nodes.tsv")
    gap = pd.read_csv(STAGE03 / "aim5_route/gap_audit.tsv", sep="\t")
    absent = set(gap.loc[gap["status"].eq("absent_from_cosmos"), "genes"])
    fits = {run: _genes_used(_load(run, "0.1")[1], genes) for run in ("cosmos_multi", "cosmos_repair")}
    sv_used = _genes_used(pd.read_csv(sv_dir / "nodes.tsv", sep="\t", dtype={"node": str}), genes)
    sv_primary_dir = sv_dir.with_name("cosmos_metabolic")
    sv_primary = (_genes_used(pd.read_csv(sv_primary_dir / "nodes.tsv", sep="\t", dtype={"node": str}), genes)
                  if (sv_primary_dir / "nodes.tsv").exists() else None)
    sv_result = json.loads((sv_dir / "result.json").read_text())
    (OUT / "fig3_sv_solve.json").write_text(json.dumps({
        "solve_directory": str(sv_dir), "status": sv_result["status"], "mip_gap": sv_result["mip_gap"],
        "weighted_fitting_error": sv_result["weighted_fitting_error"], "solve_seconds": sv_result["solve_seconds"],
        "canonical_genes_used": sorted(g for g in sv_used if any(g in e.replace("IDI1/2", "IDI1/IDI2").split("/") for e in STEP_GENES)),
    }, indent=2))
    rows = [("Reaction in primary COSMOS", ["partial" if g in absent else "yes" for g in STEP_GENES])]
    if sv_primary is not None:
        rows.append(("SV control, primary", ["yes" if _has(g, sv_primary) else "no" for g in STEP_GENES]))
    rows += [
        ("SV control, repaired", ["yes" if _has(g, sv_used) else "no" for g in STEP_GENES]),
        ("PMS vs Ctrl, primary", ["yes" if _has(g, fits["cosmos_multi"]) else "no" for g in STEP_GENES]),
        ("PMS vs Ctrl, repaired", ["yes" if _has(g, fits["cosmos_repair"]) else "no" for g in STEP_GENES]),
    ]
    fig = plt.figure(figsize=(7.2, 2.9))
    grid = fig.add_gridspec(1, 2, width_ratios=[3.1, 1], wspace=0.12)
    ax = fig.add_subplot(grid[0]); axb = fig.add_subplot(grid[1])
    n_rows = len(rows)
    for j, (_, values) in enumerate(rows):
        y = n_rows - 1 - j
        ax.plot([-0.5, len(STEP_GENES) - 0.5], [y, y], color="0.88", lw=0.6, zorder=0)
        for i, value in enumerate(values):
            if value == "yes":
                ax.scatter(i, y, s=46, color=BLUE if j else "0.35", zorder=3)
            elif value == "partial":
                ax.scatter(i, y, s=46, facecolor="white", edgecolor=VERM, linewidth=1.5, zorder=3)
            else:
                ax.scatter(i, y, s=28, facecolor="white", edgecolor="0.6", linewidth=0.8, zorder=3)
    ax.set_yticks(range(n_rows), [label for label, _ in rows][::-1])
    ax.set_xticks(range(len(STEP_GENES)), STEP_GENES, rotation=55, ha="right", fontstyle="italic")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlim(-0.6, len(STEP_GENES) - 0.4); ax.set_ylim(-0.6, n_rows - 0.25)
    ax.annotate("2 of 4 squalene-synthase half-reactions absent", xy=(6.15, n_rows - 1), xytext=(7.0, n_rows - 0.55),
                fontsize=SIZES[2], color=VERM, va="center", arrowprops=dict(arrowstyle="-", color=VERM, lw=0.6))
    ax.set_title("Recovered only after repair; unused in the PMS fit", loc="left")
    _letter(ax, "a", x=-0.02, y=1.04)
    handles = [Line2D([], [], marker="o", ls="", color="0.35", ms=6, label="reaction present"),
               Line2D([], [], marker="o", ls="", mfc="white", mec=VERM, mew=1.5, ms=6, label="gap (repaired)"),
               Line2D([], [], marker="o", ls="", color=BLUE, ms=6, label="used by CORNETO"),
               Line2D([], [], marker="o", ls="", mfc="white", mec="0.6", ms=5, label="not used")]
    ax.legend(handles=handles, ncol=4, loc="upper left", bbox_to_anchor=(-0.02, -0.33), frameon=False,
              handletextpad=0.2, columnspacing=1.0)
    axb.set_xlim(0, 10); axb.set_ylim(0, 10); axb.set_axis_off()
    nodes = [("HDLBP\n(free input)", 8.4, "o", "white", VERM),
             ("HDLBP reaction\n(Recon3D 4449, reverse)", 5.2, "D", "white", "0.35"),
             ("cholesterol\n(measured, t = 6.4)", 2.0, "s", LIGHT_BLUE, BLUE)]
    for label, y, marker, fc, ec in nodes:
        axb.scatter(1.5, y, s=170, marker=marker, facecolor=fc, edgecolor=ec, linewidth=1.4, zorder=3)
        axb.text(2.6, y, label, va="center", ha="left", fontsize=SIZES[1])
    for y0, y1 in ((8.4, 5.2), (5.2, 2.0)):
        axb.add_patch(FancyArrowPatch((1.5, y0), (1.5, y1), arrowstyle="-|>", mutation_scale=8,
                                      shrinkA=8, shrinkB=8, color="0.2"))
    axb.text(0.2, 0.0, "no selected edge leaves\ncholesterol (both runs)", fontsize=SIZES[2], color="0.3", va="top")
    axb.set_title("PMS fit: transport shortcut", loc="left")
    _letter(axb, "b", x=-0.06, y=1.04)
    save_figure(fig, OUT / "fig3_aim5_route")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 4: cosmos_multi focused subnetwork
# --------------------------------------------------------------------------

def figure4() -> None:
    edges, nodes, _ = _load("cosmos_multi", "0.1")
    graph = nx.DiGraph(list(zip(edges["source"], edges["target"])))
    observed = set(nodes.loc[nodes["is_measured"], "node"])
    tam, rxn = "Metab__5376_c", "Gene4449__HDLBP_reverse"
    targets = ["SP1", "JUNB", "JUN", "EGR1"]
    pos = {tam: (0, 5.4), "MAPK1": (1.2, 6.9), "MAPK14": (1.2, 3.6), "SP1": (2.4, 7.6), "JUNB": (2.4, 4.7),
           "SMAD3": (2.4, 3.6), "TCF3": (2.4, 2.4), "JUN": (3.6, 3.6), "EGR1": (4.8, 3.6), "ID1": (3.6, 2.4),
           "AKT1": (4.8, 2.4), "PRKACA": (6.0, 2.4), "FYN": (7.2, 2.4), "HDLBP": (0, 0.0), rxn: (1.2, 0.0), CHOL: (2.4, 0.0)}
    labels = {tam: "tamoxifen\n(drug node,\nPubChem 5376)", rxn: "HDLBP\nreaction", CHOL: "cholesterol"}
    groups = {
        "MAPK1": ("SP3, JUND, HIF1A, CEBPB, MYC,\nSTAT3, SMAD4, NCOA1, RARA", (4.1, 6.5), (2.85, 6.55)),
        "MAPK14": ("STAT1, ETS2,\nMEF2A, USF1", (1.2, 2.0), (1.2, 2.45)),
        "AKT1": ("CHUK, NRF1", (4.8, 1.15), (4.8, 1.4)),
        "PRKACA": ("SRF, NFKB1, CREB1,\nATF2, ETV1", (6.9, 1.15), (6.55, 1.5)),
        "FYN": ("TGFB1I1", (7.2, 3.5), (7.2, 3.25)),
    }
    for parent, (text, _, _) in groups.items():
        children = [c.strip() for c in text.replace("\n", " ").split(",")]
        if not all(graph.has_edge(parent, c) for c in children):
            raise ValueError(f"Figure 4 grouping no longer matches the network: {parent}")
    on_route = set()
    for t in targets + [CHOL]:
        ancestors = nx.ancestors(graph, t)
        on_route |= {(a, b) for a in ancestors for b in graph.successors(a) if b == t or b in ancestors}
    fig = plt.figure(figsize=(7.2, 4.0))
    grid = fig.add_gridspec(1, 2, width_ratios=[3.2, 1], wspace=0.32)
    ax = fig.add_subplot(grid[0]); axb = fig.add_subplot(grid[1])
    for u, v in graph.edges():
        if u in pos and v in pos:
            on = (u, v) in on_route
            ax.add_patch(FancyArrowPatch(pos[u], pos[v], arrowstyle="-|>", mutation_scale=7, shrinkA=11, shrinkB=11,
                                         color="0.15" if on else "0.6", lw=1.1 if on else 0.8, zorder=2 if on else 1))
    for parent, (text, xy, tip) in groups.items():
        ax.add_patch(FancyArrowPatch(pos[parent], tip, arrowstyle="-|>", mutation_scale=7, color="0.6", lw=0.8,
                                     shrinkA=10, shrinkB=1, zorder=1))
        ax.text(*xy, text, ha="center", va="center", fontsize=SIZES[2], zorder=3,
                bbox=dict(boxstyle="round,pad=0.3", fc=LIGHT_BLUE, ec=BLUE, lw=0.6))
    for n, (x, y) in pos.items():
        marker = "s" if n.startswith("Metab__") else "D" if n.startswith("Gene") else "o"
        fc, ec, lw = (ORANGE, "0.1", 1.0) if n in targets else (LIGHT_BLUE, BLUE, 1.0) if n in observed else ("white", "0.35", 0.8)
        if n in (tam, "HDLBP"):
            ec, lw = VERM, 1.8
        ax.scatter([x], [y], s=240 if marker != "D" else 170, marker=marker, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=4)
        label = labels.get(n, n)
        style = "normal" if n in labels else "italic"
        if n in (tam, "HDLBP"):
            ax.text(x - 0.32, y, label, ha="right", va="center", fontsize=SIZES[1], fontstyle=style, linespacing=0.95)
            ax.text(x - 0.32, y - (0.55 if n == tam else 0.35), "free input", ha="right", va="center",
                    fontsize=SIZES[2], color=VERM)
        elif n == "MAPK14":
            ax.text(x - 0.3, y, label, ha="right", va="center", fontsize=SIZES[1], fontstyle=style)
        elif y == 0.0 or n == "FYN":
            ax.text(x, y - 0.32, label, ha="center", va="top", fontsize=SIZES[1], fontstyle=style)
        else:
            ax.text(x, y + 0.3, label, ha="center", va="bottom", fontsize=SIZES[1], fontstyle=style,
                    fontweight="bold" if n in targets else "normal")
    ax.set_xlim(-1.5, 8.0); ax.set_ylim(-0.9, 8.2); ax.set_axis_off()
    ax.set_title("One drug node feeds every route into the cholesterol-linked TFs", loc="left")
    _letter(ax, "a", x=0.0)
    handles = [Line2D([], [], marker="o", ls="", mfc=ORANGE, mec="0.1", ms=7, label="focal TF (measured)"),
               Line2D([], [], marker="o", ls="", mfc=LIGHT_BLUE, mec=BLUE, ms=7, label="measured node"),
               Line2D([], [], marker="o", ls="", mfc="white", mec="0.35", ms=7, label="inferred, unmeasured"),
               Line2D([], [], marker="o", ls="", mfc="white", mec=VERM, mew=1.8, ms=7, label="free-sign input"),
               Line2D([], [], marker="s", ls="", mfc="white", mec="0.35", ms=6, label="metabolite"),
               Line2D([], [], marker="D", ls="", mfc="white", mec="0.35", ms=5, label="reaction")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.51, -0.02), ncol=2, frameon=False,
              handletextpad=0.3, columnspacing=0.8, labelspacing=0.5)
    reach = set(nodes.loc[nodes["is_measured"] & nodes["reachability"].eq("reachable"), "node"])
    inputs = nodes.loc[nodes["is_input"] & nodes["state"].ne(0), "node"]
    counts = pd.Series({i: len(nx.descendants(graph, i) & reach) for i in inputs if i in graph}).sort_values(ascending=False)
    names = {tam: "tamoxifen (drug)", "Metab__M03167_c": "5-phosphonooxy-\nL-lysine", "Metab__304_c": "PubChem 304",
             "Metab__HMDB0001877_c": "HMDB0001877", "Metab__HMDB0014646_c": "HMDB0014646",
             "Metab__hislysthr_e": "His-Lys-Thr"}
    top, rest = counts.iloc[:6], counts.iloc[6:]
    y = np.arange(len(top))[::-1] + 1
    axb.barh(y, top.values, color=[VERM] + ["0.7"] * (len(top) - 1), height=0.6)
    for yy, v in zip(y, top.values):
        axb.text(v + 0.8, yy, str(v), va="center", fontsize=SIZES[2])
    axb.text(1, 0.1, f"{len(rest)} other inputs: {rest.min()}–{rest.max()} each", fontsize=SIZES[2], color="0.25", va="center")
    axb.set_yticks(y, [names.get(i, i) for i in top.index]); axb.tick_params(axis="y", length=0)
    axb.set_ylim(-0.5, len(top) + 0.6); axb.set_xlim(0, 46)
    axb.set_xlabel(f"Measured nodes downstream\n(of {len(reach)} reachable)")
    axb.set_title("Half the fit hangs off one input", loc="left")
    _letter(axb, "b", x=-0.02)
    save_figure(fig, OUT / "fig4_cosmos_subnetwork")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 5: robustness from the stage-04 lambda grid
# --------------------------------------------------------------------------

def _dominant_inputs() -> pd.DataFrame:
    rows = []
    for run in RUN_LABELS:
        for lam in LAMBDAS:
            edges, nodes, _ = _load(run, lam)
            graph = nx.DiGraph(list(zip(edges["source"], edges["target"])))
            reach = set(nodes.loc[nodes["is_measured"] & nodes["reachability"].eq("reachable"), "node"])
            inputs = [i for i in nodes.loc[nodes["is_input"] & nodes["state"].ne(0), "node"] if i in graph]
            counts = {i: len(nx.descendants(graph, i) & reach) for i in inputs}
            top = max(counts, key=counts.get)
            rows.append({"run": run, "lambda": lam, "top_input": top, "n_explained": counts[top],
                         "n_reachable": len(reach), "share": counts[top] / len(reach),
                         "edges": set(zip(edges["source"], edges["target"]))})
    return pd.DataFrame(rows)


def _rescoring() -> pd.DataFrame:
    """Score every saved solution of a run under every lambda's objective."""
    rows = []
    for run in ("cosmos_multi", "cosmos_repair"):
        solutions = []
        for lam in ("0",) + LAMBDAS:
            edges, nodes, result = _load(run, lam)
            error = sum(v for k, v in result["objectives"].items() if k.startswith("error"))
            penalty = result["objectives"].get("regularization_edge_has_signal")
            if penalty is None:  # lambda = 0 reports no penalty: selected edges plus active input edges
                penalty = len(edges) + int((nodes["is_input"] & nodes["state"].ne(0)).sum())
            solutions.append((float(lam), error, penalty))
        for lam, error, penalty in solutions[1:]:
            own = error + lam * penalty
            best = min(e + lam * p for _, e, p in solutions)
            rows.append({"run": run, "lambda": lam, "own_objective": own, "best_known_objective": best,
                         "excess": (own - best) / best})
    return pd.DataFrame(rows)


def _input_name(node: str) -> str:
    names = {"Metab__5376_c": "tamoxifen", "Metab__HMDB0003441_c": "cinnam-\naldehyde", "Metab__HMDB0011620_c": "HMDB\n0011620"}
    return names.get(node, node.replace("Metab__", ""))


def figure5(properties: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    dominant, rescoring = _dominant_inputs(), _rescoring()
    fig = plt.figure(figsize=(7.2, 5.4))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.7, 1], height_ratios=[1.45, 0.95], wspace=0.45, hspace=0.75)
    axa, axb, axc = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :])
    runs = list(RUN_LABELS)
    for r, run in enumerate(runs):
        sub = dominant[dominant["run"].eq(run)]
        previous = None
        for c, row in enumerate(sub.itertuples()):
            axa.scatter(c, r, s=20 + 200 * row.share, color=ORANGE, zorder=3)
            label = "same" if row.top_input == previous else _input_name(row.top_input)
            axa.text(c, r + 0.33, label, ha="center", va="top", fontsize=SIZES[2],
                     color="0.55" if label == "same" else "0.1")
            previous = row.top_input
    axa.set_xticks(range(len(LAMBDAS)), LAMBDAS); axa.set_yticks(range(len(runs)), [RUN_LABELS[r] for r in runs])
    axa.set_xlabel("λ"); axa.set_ylim(len(runs) - 0.15, -0.6); axa.set_xlim(-0.6, len(LAMBDAS) - 0.4)
    axa.tick_params(axis="y", length=0)
    axa.set_title("The dominant input changes with λ", loc="left")
    _letter(axa, "a", x=-0.02)
    for share in (0.4, 0.9):
        axa.scatter([], [], s=20 + 200 * share, color=ORANGE, label=f"{share:.0%}")
    axa.legend(title="measurements downstream", loc="upper left", bbox_to_anchor=(0.0, -0.22), ncol=2,
               frameon=False, handletextpad=0.2, title_fontsize=SIZES[1])
    for r, run in enumerate(runs):
        sets = list(dominant.loc[dominant["run"].eq(run), "edges"])
        union = set.union(*sets)
        n_all = sum(all(e in s for s in sets) for e in union)
        axb.barh(r, n_all / len(union), color=BLUE, height=0.6)
        axb.text(n_all / len(union) + 0.01, r, f"{n_all}/{len(union)}", va="center", fontsize=SIZES[2])
    axb.set_yticks(range(len(runs)), [""] * len(runs)); axb.set_ylim(len(runs) - 0.3, -0.6); axb.set_xlim(0, 0.42)
    axb.tick_params(axis="y", length=0)
    axb.set_xlabel("Edges in all five\nλ networks (fraction)")
    axb.set_title("Few edges persist", loc="left")
    _letter(axb, "b", x=-0.02)
    # The solver-rescoring panel was removed on 2026-10-16 (user decision); fig5_rescoring.tsv is
    # still written by main() as a record.
    # panel c: degree bias of node selection (07_network_properties.py)
    bias = pd.read_csv(properties / "degree_bias.tsv", sep="\t")
    groups = [("omnipath", "OmniPath"), ("cosmos_tf", "COSMOS, TFs"), ("cosmos_multi", "COSMOS, TFs + metab."),
              ("cosmos_repair", "COSMOS repaired"), ("omnipath_rewired", "random networks")]
    for k, (name, label) in enumerate(groups):
        values = bias.loc[bias["network"].eq(name), "odds_ratio_per_degree_doubling"].dropna().to_numpy()
        jitter = np.linspace(-0.12, 0.12, len(values)) if len(values) > 1 else np.zeros(len(values))
        colour = "0.55" if name == "omnipath_rewired" else BLUE
        axc.scatter(values, k + jitter, s=14, color=colour, zorder=3)
        axc.plot([np.median(values)] * 2, [k - 0.25, k + 0.25], color="0.1", lw=1.2, zorder=4)
    axc.axvline(1, color="0.6", lw=0.6, ls="--")
    axc.set_yticks(range(len(groups)), [label for _, label in groups]); axc.set_ylim(len(groups) - 0.5, -0.6)
    axc.tick_params(axis="y", length=0)
    axc.set_xlim(0.8, 3.3)
    axc.set_xlabel("Odds of selecting a node\nper doubling of its PKN degree")
    axc.set_title("Hubs are favoured, even in random networks", loc="left")
    _letter(axc, "c", x=-0.02)
    save_figure(fig, OUT / "fig5_robustness")
    plt.close(fig)
    return dominant.drop(columns="edges"), rescoring


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sv-dir", type=Path, default=None,
                        help="SV control solve; default: stage-05 output, else the bundled 600 s trial")
    parser.add_argument("--only", nargs="+", choices=["1", "2", "3", "4", "5"], default=["1", "2", "3", "4", "5"])
    parser.add_argument("--network-properties", type=Path,
                        default=ROOT / "results/ionescu_corneto/07_downstream/network_properties",
                        help="output folder of scripts/analysis/07_network_properties.py (Figure 5c)")
    args = parser.parse_args(argv)
    if args.sv_dir is None:
        registered = STAGE05 / "sv_control/cosmos_repair_metabolic"
        args.sv_dir = registered if (registered / "result.json").exists() else OUT / "data/sv_trial_600s"
        if args.sv_dir != registered:
            print("Figure 3: stage-05 SV control not found; using the 600 s trial solve (say so in the caption)")
    # Rows "SV control" are drawn from <sv_dir> (repaired) and its sibling cosmos_metabolic (primary).
    OUT.mkdir(parents=True, exist_ok=True)
    _style()
    if "1" in args.only:
        figure1()
    if "2" in args.only:
        figure2()
    if "3" in args.only:
        figure3(args.sv_dir)
    if "4" in args.only:
        figure4()
    if "5" in args.only:
        dominant, rescoring = figure5(args.network_properties)
        dominant.to_csv(OUT / "fig5_dominant_inputs.tsv", sep="\t", index=False)
        rescoring.to_csv(OUT / "fig5_rescoring.tsv", sep="\t", index=False)
    print(f"Figures written to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
