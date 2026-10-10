"""Run the pre-registered Aim-5 COSMOS canonical-route audit."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from metabolic_route import (  # noqa: E402
    CANONICAL_GENES,
    CHOLESTEROL_NODE_PREFIX,
    canonical_step_table,
    cofactor_nodes,
    fetch_kegg_gene_lists,
    gap_audit,
    genome_absent_by_subsystem,
    metabolism_graph,
    path_variant_tables,
    read_pkn,
    repaired_pkn,
    typed_metabolite_tf_null,
    verify_kegg_canonical_genes,
    waypoint_table,
)
from pkn import TARGET_TFS, fetch_recon3d, load_hgnc_maps  # noqa: E402
from plot_config import apply_plot_style, save_figure  # noqa: E402
from preprocess import stage_provenance, write_stage_provenance  # noqa: E402

ANALYSIS = ROOT / "results/ionescu_corneto/03_pkn/aim5_route"
PKN_ROOT = ROOT / "results/ionescu_corneto/03_pkn/cosmos"
RESOURCE_ROOT = ROOT / "data/resources/pkn"
HGNC_ROOT = ROOT / "data/resources/hgnc"
SEED_BASE = 20261010


def _load_hgnc() -> tuple[dict, Path]:
    paths = sorted(HGNC_ROOT.glob("hgnc_complete_set_*.txt"))
    if len(paths) != 1:
        raise ValueError(f"Expected one HGNC complete set, found {paths}")
    return load_hgnc_maps(paths[0]), paths[0]


def _route_figure(gaps: pd.DataFrame, stem: Path) -> None:
    """Draw the ordered route as a diagnostic, marking missing steps."""
    apply_plot_style()
    rows = gaps.sort_values("step_order").drop_duplicates(
        ["step_order", "previous_intermediate", "next_intermediate"], keep="first"
    )
    rows = rows[rows["branch"].eq("core") | rows["branch"].eq("bloch")]
    if rows.empty:
        raise ValueError("No canonical route rows available for the figure")
    labels = [str(rows.iloc[0]["previous_intermediate"])] + rows["next_intermediate"].astype(str).tolist()
    labels = list(dict.fromkeys(labels))
    positions = {label: (index, 0.0) for index, label in enumerate(labels)}
    fig, ax = plt.subplots(figsize=(max(12, len(labels) * 0.8), 3.5))
    status_colors = {
        "present_forward": "#0072B2",
        "present_reverse_only": "#E69F00",
        "pruned_unexpressed": "#CC79A7",
        "absent_from_cosmos": "#D55E00",
        "present_but_unconnected": "#999999",
    }
    for index, row in enumerate(rows.itertuples(index=False)):
        left = positions.get(str(row.previous_intermediate))
        right = positions.get(str(row.next_intermediate))
        if left is None or right is None:
            continue
        color = status_colors.get(str(row.status), "#999999")
        ax.plot([left[0] + 0.18, right[0] - 0.18], [0, 0], color=color, linewidth=3)
        ax.text((left[0] + right[0]) / 2, 0.12, str(row.status).replace("_", " "), ha="center", va="bottom", fontsize=6, rotation=35)
    for label, (x, y) in positions.items():
        ax.text(x, y, label, ha="center", va="center", fontsize=7,
                bbox={"boxstyle": "round,pad=0.35", "facecolor": "#F0E442", "edgecolor": "#333333"})
    ax.set_xlim(-0.7, max(positions.values(), key=lambda value: value[0])[0] + 0.7)
    ax.set_ylim(-0.8, 0.95)
    ax.set_title("Aim-5 canonical cholesterol route: primary COSMOS gap audit", loc="left")
    ax.legend(
        handles=[Patch(facecolor=color, edgecolor=color, label=label.replace("_", " ")) for label, color in status_colors.items()],
        loc="upper left", bbox_to_anchor=(0, -0.08), ncol=3, fontsize=7,
    )
    ax.set_axis_off()
    fig.tight_layout()
    save_figure(fig, stem)
    plt.close(fig)


def main() -> int:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    hgnc_maps, hgnc_path = _load_hgnc()
    recon_model, recon_path, recon_metadata = fetch_recon3d(RESOURCE_ROOT)
    kegg_ids, kegg_metadata = fetch_kegg_gene_lists(RESOURCE_ROOT)
    kegg_table = verify_kegg_canonical_genes(kegg_ids, hgnc_maps)
    kegg_table.to_csv(ANALYSIS / "canonical_gene_kegg.tsv", sep="\t", index=False)

    primary_path = PKN_ROOT / "primary/pkn.tsv"
    unpruned_path = PKN_ROOT / "variants/unpruned/pkn.tsv"
    primary_edges = read_pkn(primary_path)
    unpruned_edges = read_pkn(unpruned_path)
    all_nodes = set(primary_edges["source"]) | set(primary_edges["target"]) | set(unpruned_edges["source"]) | set(unpruned_edges["target"])
    steps = canonical_step_table(recon_model, hgnc_maps, pkn_nodes=all_nodes)
    gaps = gap_audit(steps, primary_edges, unpruned_edges)
    absent = genome_absent_by_subsystem(recon_model, set(unpruned_edges["source"]) | set(unpruned_edges["target"]))

    primary_metabolism = metabolism_graph(primary_edges)
    cofactors = cofactor_nodes(
        recon_model,
        all_nodes,
        dict(primary_metabolism.degree()),
        set(steps["substrate_model_id"]) | set(steps["product_model_id"]),
    )
    repaired_edges, repaired_nodes, repair_log = repaired_pkn(
        primary_edges,
        pd.read_csv(PKN_ROOT / "primary/nodes.tsv", sep="\t", dtype=str),
        steps,
        gaps,
    )
    repair_dir = PKN_ROOT / "variants/canonical_repair"
    repair_dir.mkdir(parents=True, exist_ok=True)
    repaired_edges.to_csv(repair_dir / "pkn.tsv", sep="\t", index=False)
    repaired_nodes.to_csv(repair_dir / "nodes.tsv", sep="\t", index=False)
    repair_log.to_csv(repair_dir / "filter_log.tsv", sep="\t", index=False)

    variants = {"primary": primary_edges, "unpruned": unpruned_edges, "repaired": repaired_edges}
    waypoints = waypoint_table(variants, steps)
    path_variants, path_edges = path_variant_tables(variants, steps, cofactors)
    # The registered stage-03 null is defined on the primary COSMOS PKN. The
    # repaired PKN is a route diagnostic, not a second null family.
    cholesterol = sorted(node for node in set(primary_edges["source"]) | set(primary_edges["target"]) if node.startswith(CHOLESTEROL_NODE_PREFIX))
    typed_null = typed_metabolite_tf_null(
        primary_edges,
        variant="primary",
        cholesterol_nodes=cholesterol,
        target_tfs=TARGET_TFS,
        n_replicates=100,
        seed_base=SEED_BASE,
    )

    steps.to_csv(ANALYSIS / "canonical_steps.tsv", sep="\t", index=False)
    gaps.to_csv(ANALYSIS / "gap_audit.tsv", sep="\t", index=False)
    absent.to_csv(ANALYSIS / "recon3d_absent_by_subsystem.tsv", sep="\t", index=False)
    cofactors.to_csv(ANALYSIS / "cofactor_nodes.tsv", sep="\t", index=False)
    waypoints.to_csv(ANALYSIS / "waypoints.tsv", sep="\t", index=False)
    path_variants.to_csv(ANALYSIS / "path_variants.tsv", sep="\t", index=False)
    path_edges.to_csv(ANALYSIS / "path_variant_edges.tsv", sep="\t", index=False)
    typed_null.to_csv(ANALYSIS / "null_metabolite_to_tf_typed.tsv", sep="\t", index=False)
    _route_figure(gaps, ANALYSIS / "figures/canonical_route")

    input_paths = [
        recon_path, hgnc_path, primary_path, unpruned_path,
        PKN_ROOT / "primary/nodes.tsv", RESOURCE_ROOT / "cosmos_meta_network.sif",
        RESOURCE_ROOT / "cosmos_metapkn__20200122.txt",
    ] + [RESOURCE_ROOT / f"kegg_{pathway}_gene_links.tsv" for pathway in ("hsa00900", "hsa00100")]
    input_paths = [path for path in input_paths if path.exists()]
    code_paths = [ROOT / "src/metabolic_route.py", ROOT / "src/pkn.py", ROOT / "src/plot_config.py", ROOT / "src/preprocess.py", Path(__file__), ROOT / "doc/decisions/2026-10-10-aim5-metabolic-route-plan.md"]
    parameters = {
        "canonical_genes": list(CANONICAL_GENES),
        "repair_scope": "canonical_only",
        "primary_path_readout": "metabolism_no_cofactors_weighted",
        "path_variants": ["full_unweighted", "metabolism_only", "metabolism_no_cofactors", "metabolism_no_cofactors_weighted", "metabolism_no_cofactors_weighted_top1pct"],
        "typed_null_pkn": "cosmos_primary",
        "typed_null_replicates": 100,
        "seed_base": SEED_BASE,
        "cofactor_rule": "registered BiGG IDs mapped through Recon3D; canonical intermediates protected",
    }
    resources = {"recon3d": recon_metadata, "kegg": kegg_metadata}
    provenance = stage_provenance(ROOT, "aim5_route", input_paths, resources, parameters, code_paths)
    provenance["diagnostics"] = {
        "n_canonical_steps": len(steps),
        "n_absent_canonical_steps": int(gaps["status"].eq("absent_from_cosmos").sum()),
        "n_gene_associated_recon3d": int(absent.loc[absent["subsystem"].eq("__TOTAL__"), "n_gene_associated_reactions"].iloc[0]),
        "n_gene_associated_absent": int(absent.loc[absent["subsystem"].eq("__TOTAL__"), "n_absent_from_cosmos"].iloc[0]),
    }
    write_stage_provenance(ANALYSIS, provenance)
    print(f"Aim-5 route audit complete: {ANALYSIS}")
    print(gaps["status"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
