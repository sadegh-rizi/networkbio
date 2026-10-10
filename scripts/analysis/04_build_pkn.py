"""Stage 03: build the PKNs, identifier maps and Aim-5 structure tables.

Run from the repository root with the project Python.  This script reads the
stage-00/01/02 outputs but does not read TF activity values or run CORNETO.
"""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from network_plot import layered_positions  # noqa: E402
from pkn import (  # noqa: E402
    CHOLESTEROL_HMDB,
    COMPLEX_MEMBERS,
    MEVALONATE_HMDB,
    SOURCE_GENES,
    TARGET_TFS,
    cosmos_node_grammar,
    degree_preserving_null_summary,
    fetch_cosmos,
    fetch_omnipath,
    fetch_recon3d,
    fetch_workbench_lookups,
    filter_cosmos,
    filter_omnipath,
    load_hgnc_maps,
    map_metabolite_features,
    recon3d_hmdb_by_kegg,
    sha256,
    shortest_path_edges,
    signed_shortest_paths,
)
from plot_config import NODE_ROLE_COLORS, apply_plot_style, save_figure  # noqa: E402
from preprocess import stage_provenance, write_stage_provenance  # noqa: E402


ANALYSIS = ROOT / "results/ionescu_corneto/03_pkn"
INPUTS = ROOT / "results/ionescu_corneto/00_inputs"
PREPROCESSED = ROOT / "results/ionescu_corneto/01_preprocessing"
ACTIVITIES = ROOT / "results/ionescu_corneto/02_activities"
PKN_RESOURCES = ROOT / "data/resources/pkn"
METABOLITE_RESOURCES = ROOT / "data/resources/metabolite_ids"
SEED_BASE = 20261009
MAX_PATH_FIGURE_WIDTH_INCHES = 24


def _write_pkn(folder: Path, edges: pd.DataFrame, nodes: pd.DataFrame, log: pd.DataFrame) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    edges[["source", "sign", "target"]].to_csv(folder / "pkn.tsv", sep="\t", index=False)
    nodes.to_csv(folder / "nodes.tsv", sep="\t", index=False)
    log.to_csv(folder / "filter_log.tsv", sep="\t", index=False)


def _load_expressed() -> set[str]:
    expression = pd.read_csv(PREPROCESSED / "rnaseq/log2cpm.tsv", sep="\t", usecols=["symbol"])
    return set(expression["symbol"].astype(str))


def _load_hgnc() -> tuple[dict, Path]:
    candidates = sorted((ROOT / "data/resources/hgnc").glob("hgnc_complete_set_*.txt"))
    if len(candidates) != 1:
        raise ValueError(f"Expected exactly one HGNC complete set, found {candidates}")
    return load_hgnc_maps(candidates[0]), candidates[0]


def _load_features() -> dict[str, pd.DataFrame]:
    paths = {
        "ST003328": INPUTS / "metabolomics/ST003328/features.tsv",
        "ST003331": PREPROCESSED / "metabolomics/ST003331/features.tsv",
        "ST003332": PREPROCESSED / "metabolomics/ST003332/features.tsv",
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Stage-01 feature tables are missing: {missing}")
    return {study: pd.read_csv(path, sep="\t", dtype=str) for study, path in paths.items()}


def _tf_coverage(node_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict] = []
    for resource in ("collectri", "collectri_nocomplex", "dorothea_ABC"):
        summary_path = ACTIVITIES / resource / "regulon_summary.tsv"
        if not summary_path.exists():
            raise FileNotFoundError(f"Stage-02 regulon summary is missing: {summary_path}")
        summary = pd.read_csv(summary_path, sep="\t")
        kept = summary[summary["kept"].map(lambda value: str(value).lower() == "true")]
        for pkn_name, nodes in node_tables.items():
            node_set = set(nodes["node"])
            for row in kept.itertuples(index=False):
                rows.append(
                    {
                        "resource": resource,
                        "pkn": pkn_name,
                        "record_type": "tf",
                        "tf": row.TF,
                        "member_gene": None,
                        "is_complex": bool(row.is_complex),
                        "node_present": row.TF in node_set,
                        "n_targets_expression": int(row.n_targets_expression),
                    }
                )
            for complex_name, members in COMPLEX_MEMBERS.items():
                for member in members:
                    rows.append(
                        {
                            "resource": resource,
                            "pkn": pkn_name,
                            "record_type": "complex_member",
                            "tf": complex_name,
                            "member_gene": member,
                            "is_complex": True,
                            "node_present": member in node_set,
                            "n_targets_expression": None,
                        }
                    )
    return pd.DataFrame(rows)


def _aim5_pairs(pkn_name: str, node_tables: dict[str, pd.DataFrame]) -> list[tuple[str, list[str], list[str], int]]:
    pairs = [("source_tf_to_tf", list(SOURCE_GENES), list(TARGET_TFS), 8)]
    if pkn_name == "cosmos":
        cosmos_nodes = set(node_tables[pkn_name]["node"])
        metabolite_nodes = sorted(
            node
            for node in cosmos_nodes
            if node.startswith(f"Metab__{MEVALONATE_HMDB}_") or node.startswith(f"Metab__{CHOLESTEROL_HMDB}_")
        )
        pairs.insert(0, ("hmgcr_to_metabolite", ["HMGCR"], metabolite_nodes, 80))
        pairs.insert(1, ("metabolite_to_tf", metabolite_nodes, list(TARGET_TFS), 8))
    return pairs


def _aim5_paths(node_tables: dict[str, pd.DataFrame], edge_tables: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path_rows: list[pd.DataFrame] = []
    edge_rows: list[pd.DataFrame] = []
    null_rows: list[pd.DataFrame] = []
    for pkn_name in ("omnipath", "cosmos"):
        n_replicates = 100 if pkn_name == "omnipath" else 20
        for leg, pair_sources, pair_targets, path_cap in _aim5_pairs(pkn_name, node_tables):
            paths = signed_shortest_paths(
                edge_tables[pkn_name],
                pair_sources,
                pair_targets,
                max_length=path_cap,
            )
            paths.insert(0, "path_cap", path_cap)
            paths.insert(0, "leg", leg)
            paths.insert(0, "pkn", pkn_name)
            path_rows.append(paths)
            pkn_paths = paths.drop(columns=["pkn", "leg", "path_cap"])
            path_edges = shortest_path_edges(edge_tables[pkn_name], pkn_paths, max_paths_per_pair_sign=100)
            if not path_edges.empty:
                path_edges.insert(0, "path_cap", path_cap)
                path_edges.insert(0, "leg", leg)
                path_edges.insert(0, "pkn", pkn_name)
                edge_rows.append(path_edges)
            null = degree_preserving_null_summary(
                edge_tables[pkn_name],
                pkn_paths,
                n_replicates=n_replicates,
                seed_base=SEED_BASE,
            )
            null.insert(0, "path_cap", path_cap)
            null.insert(0, "leg", leg)
            null.insert(0, "pkn", pkn_name)
            null_rows.append(null)
    path_table = pd.concat(path_rows, ignore_index=True)
    path_edge_columns = [
        "pkn", "leg", "path_cap", "source", "target", "path_sign", "path_length", "path_index", "edge_index",
        "edge_source", "edge_sign", "edge_target",
    ]
    path_edge_table = pd.concat(edge_rows, ignore_index=True) if edge_rows else pd.DataFrame(columns=path_edge_columns)
    null_table = pd.concat(null_rows, ignore_index=True)
    return path_table, path_edge_table, null_table


def _path_figure(pkn_name: str, edges: pd.DataFrame, paths: pd.DataFrame, path_edges: pd.DataFrame, stem: Path) -> None:
    """Draw the P9 path union without showing stage-02 values."""
    apply_plot_style()
    pkn_paths = path_edges[path_edges["pkn"].eq(pkn_name)]
    pair_paths = paths[paths["pkn"].eq(pkn_name)]
    highlighted = set(zip(pkn_paths["edge_source"], pkn_paths["edge_target"]))
    source_nodes = set(pair_paths["source"]) | set(pair_paths["target"])
    selected_nodes = set(pkn_paths["edge_source"]) | set(pkn_paths["edge_target"]) | source_nodes
    selected_edges = edges[
        edges["source"].isin(selected_nodes) & edges["target"].isin(selected_nodes)
    ].drop_duplicates(["source", "target"])
    nodes = sorted(selected_nodes)
    figure_width = min(max(10, 1.4 * max(len(nodes), 1)), MAX_PATH_FIGURE_WIDTH_INCHES)
    fig, ax = plt.subplots(figsize=(figure_width, 7))
    if not nodes:
        ax.text(0.5, 0.5, "No path within the registered path caps", ha="center", va="center")
        ax.set_axis_off()
        save_figure(fig, stem)
        plt.close(fig)
        return
    positions = layered_positions(nodes, list(zip(selected_edges["source"], selected_edges["target"])))
    for row in selected_edges.itertuples(index=False):
        is_path = (row.source, row.target) in highlighted
        arrow = FancyArrowPatch(
            positions[row.source],
            positions[row.target],
            arrowstyle="-|>" if row.sign > 0 else "-[, widthB=0.45, lengthB=0",
            mutation_scale=12 if is_path else 9,
            shrinkA=17,
            shrinkB=21 if row.sign < 0 else 17,
            color="#333333" if is_path else "#BBBBBB",
            linewidth=2.0 if is_path else 0.8,
            linestyle="-" if is_path else (0, (3, 2)),
            zorder=1,
        )
        ax.add_patch(arrow)
    for node in nodes:
        if node in SOURCE_GENES:
            color = NODE_ROLE_COLORS["input"]
            edge_color = NODE_ROLE_COLORS["input"]
        elif node in TARGET_TFS:
            color = NODE_ROLE_COLORS["measured_up"]
            edge_color = "#000000"
        elif node.startswith("Metab__"):
            color = "#F0E442"
            edge_color = "#555555"
        else:
            color = NODE_ROLE_COLORS["intermediate"]
            edge_color = "#555555"
        ax.scatter(
            *positions[node],
            s=900,
            marker="D" if node.startswith("Metab__") else "o",
            facecolor=color,
            edgecolor=edge_color,
            linewidth=2.0 if node in SOURCE_GENES or node in TARGET_TFS else 0.9,
            zorder=2,
        )
        ax.text(*positions[node], node, ha="center", va="center", fontsize=6.5, zorder=3)
    ax.set_title(
        f"{pkn_name}: shortest signed paths, no activity values shown\n"
        "metabolic leg cap 80; signalling legs cap 8; highlighted edges are shortest paths",
        loc="left",
    )
    ax.legend(
        handles=[
            Patch(facecolor=NODE_ROLE_COLORS["input"], edgecolor=NODE_ROLE_COLORS["input"], label="Aim-5 source"),
            Patch(facecolor=NODE_ROLE_COLORS["measured_up"], edgecolor="#000000", label="Aim-5 target"),
            Patch(facecolor="#F0E442", edgecolor="#555555", label="metabolite"),
            Line2D([], [], color="#333333", linewidth=2, label="shortest-path edge"),
            Line2D([], [], color="#BBBBBB", linewidth=0.8, linestyle=(0, (3, 2)), label="other PKN edge among path nodes"),
        ],
        loc="upper left",
        bbox_to_anchor=(0, -0.02),
        ncol=2,
    )
    ax.set_xlim(min(x for x, _ in positions.values()) - 0.5, max(x for x, _ in positions.values()) + 0.5)
    ax.set_ylim(min(y for _, y in positions.values()) - 0.8, max(y for _, y in positions.values()) + 0.6)
    ax.set_axis_off()
    save_figure(fig, stem)
    plt.close(fig)


def _summary_counts(filter_logs: dict[str, pd.DataFrame], mapping: pd.DataFrame, tf_coverage: pd.DataFrame, paths: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for pkn_name, log in filter_logs.items():
        for record in log.to_dict(orient="records"):
            for metric in ("n_rows", "n_edges", "n_nodes", "n_removed_rows", "n_removed_nodes"):
                rows.append({"section": "filter", "pkn": pkn_name, "metric": f"{record['step']}:{metric}", "value": record[metric]})
    if not mapping.empty:
        grouped = mapping.groupby(["study", "labelling_class"], dropna=False).agg(
            n_features=("feature_id", "nunique"),
            n_mapping_rows=("feature_id", "size"),
            n_rows_mapped=("mapped_to_pkn", "sum"),
        ).reset_index()
        for record in grouped.to_dict(orient="records"):
            for metric in ("n_features", "n_mapping_rows", "n_rows_mapped"):
                rows.append({"section": "mapping", "pkn": record["study"], "metric": f"{record['labelling_class']}:{metric}", "value": record[metric]})
    for pkn_name, group in tf_coverage.groupby("pkn"):
        rows.append({"section": "tf_coverage", "pkn": pkn_name, "metric": "rows_present", "value": int(group["node_present"].sum())})
        rows.append({"section": "tf_coverage", "pkn": pkn_name, "metric": "rows_total", "value": len(group)})
    if not paths.empty:
        for (pkn_name, leg), group in paths.groupby(["pkn", "leg"]):
            rows.append({"section": "aim5", "pkn": pkn_name, "metric": f"{leg}:positive_routes_present", "value": int(group["positive_route_present"].sum())})
            rows.append({"section": "aim5", "pkn": pkn_name, "metric": f"{leg}:pairs", "value": len(group)})
    return pd.DataFrame(rows)


def main() -> int:
    expressed = _load_expressed()
    hgnc_maps, hgnc_path = _load_hgnc()
    hgnc_maps["expressed"] = expressed
    features = _load_features()
    PKN_RESOURCES.mkdir(parents=True, exist_ok=True)
    METABOLITE_RESOURCES.mkdir(parents=True, exist_ok=True)

    omni_raw, omni_resource_path, omni_metadata = fetch_omnipath(PKN_RESOURCES)
    cosmos_raw, cosmos_metadata = fetch_cosmos(PKN_RESOURCES)
    recon_model, recon_path, recon_metadata = fetch_recon3d(PKN_RESOURCES)
    kegg_ids = set()
    for table in features.values():
        kegg_ids.update(table["kegg_id"].dropna().astype(str))
    workbench = fetch_workbench_lookups(METABOLITE_RESOURCES, kegg_ids)
    recon_by_kegg = recon3d_hmdb_by_kegg(recon_model)

    grammar_path = ANALYSIS / "summary/cosmos_node_grammar.md"
    grammar_path.parent.mkdir(parents=True, exist_ok=True)
    grammar_path.write_text(cosmos_node_grammar(cosmos_raw))

    edge_tables: dict[str, pd.DataFrame] = {}
    node_tables: dict[str, pd.DataFrame] = {}
    filter_logs: dict[str, pd.DataFrame] = {}
    for variant, folder in {
        "primary": ANALYSIS / "omnipath/primary",
        "signor_only": ANALYSIS / "omnipath/variants/signor_only",
        "nc_rule": ANALYSIS / "omnipath/variants/nc_rule",
        "unpruned": ANALYSIS / "omnipath/variants/unpruned",
    }.items():
        edges, log, nodes = filter_omnipath(
            omni_raw,
            expressed=expressed,
            hgnc_maps=hgnc_maps,
            variant=variant,
        )
        _write_pkn(folder, edges, nodes, log)
        filter_logs[f"omnipath_{variant}"] = log
        if variant == "primary":
            edge_tables["omnipath"] = edges
            node_tables["omnipath"] = nodes

    for variant, folder in {
        "primary": ANALYSIS / "cosmos/primary",
        "unpruned": ANALYSIS / "cosmos/variants/unpruned",
    }.items():
        edges, log, nodes = filter_cosmos(
            cosmos_raw,
            hgnc_maps=hgnc_maps,
            expressed=expressed,
            variant=variant,
        )
        _write_pkn(folder, edges, nodes, log)
        filter_logs[f"cosmos_{variant}"] = log
        if variant == "primary":
            edge_tables["cosmos"] = edges
            node_tables["cosmos"] = nodes

    mapping = map_metabolite_features(
        features,
        pkn_nodes=set(node_tables["cosmos"]["node"]),
        recon_by_kegg=recon_by_kegg,
        workbench=workbench,
    )
    mapping_dir = ANALYSIS / "mapping"
    mapping_dir.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(mapping_dir / "metabolite_mapping.tsv", sep="\t", index=False)
    mapping[mapping["study"].eq("ST003328")].to_csv(mapping_dir / "lipid_mapping.tsv", sep="\t", index=False)
    tf_coverage = _tf_coverage(node_tables)
    tf_coverage.to_csv(mapping_dir / "tf_node_coverage.tsv", sep="\t", index=False)

    paths, path_edges, null_summary = _aim5_paths(node_tables, edge_tables)
    aim5_dir = ANALYSIS / "aim5"
    aim5_dir.mkdir(parents=True, exist_ok=True)
    paths.to_csv(aim5_dir / "paths.tsv", sep="\t", index=False)
    path_edges.to_csv(aim5_dir / "path_edges.tsv", sep="\t", index=False)
    null_summary.to_csv(aim5_dir / "null_summary.tsv", sep="\t", index=False)
    _path_figure("omnipath", edge_tables["omnipath"], paths, path_edges, ANALYSIS / "figures/aim5_paths_omnipath")
    _path_figure("cosmos", edge_tables["cosmos"], paths, path_edges, ANALYSIS / "figures/aim5_paths_cosmos")

    counts = _summary_counts(filter_logs, mapping, tf_coverage, paths)
    counts.to_csv(ANALYSIS / "summary/counts.tsv", sep="\t", index=False)
    provenance_inputs = [hgnc_path, omni_resource_path, recon_path]
    provenance_inputs.extend(
        path for path in [
            PKN_RESOURCES / "cosmos_meta_network.sif",
            PKN_RESOURCES / "cosmos_metapkn__20200122.txt",
            METABOLITE_RESOURCES / "metabolomics_workbench_kegg_hmdb.tsv",
            INPUTS / "metabolomics/ST003328/features.tsv",
            PREPROCESSED / "metabolomics/ST003331/features.tsv",
            PREPROCESSED / "metabolomics/ST003332/features.tsv",
            PREPROCESSED / "rnaseq/log2cpm.tsv",
            ACTIVITIES / "collectri/regulon_summary.tsv",
            ACTIVITIES / "collectri_nocomplex/regulon_summary.tsv",
            ACTIVITIES / "dorothea_ABC/regulon_summary.tsv",
        ] if path.exists()
    )
    code_paths = [
        ROOT / "src/pkn.py",
        ROOT / "src/network_plot.py",
        ROOT / "src/plot_config.py",
        ROOT / "src/preprocess.py",
        ROOT / "scripts/analysis/04_build_pkn.py",
    ]
    settings = {
        "omnipath_variants": ["primary", "signor_only", "nc_rule", "unpruned"],
        "cosmos_variants": ["primary", "unpruned"],
        "curation_effort_minimum": 2,
        "expression_pruning": "stage-01 log2cpm symbol index, CPM >= 1 in at least 3 libraries",
        "cosmos_parser": {
            "enzyme_gene_tokens": "recognized HGNC suffix tokens",
            "multi_gene_rule": "retain if any recognized gene is expressed",
            "pseudo_reaction_rule": "retain nodes without recognized gene tokens",
            "numeric_metabolite_rule": "preserve numeric model identifiers; do not relabel as HMDB",
        },
        "aim5_path_caps": {"hmgcr_to_metabolite": 80, "metabolite_to_tf": 8, "source_tf_to_tf": 8},
        "rewiring": "3 x |E| accepted directed double-edge swaps, sign-preserving",
        "rewiring_seed_base": SEED_BASE,
        "rewiring_replicates": {"omnipath": 100, "cosmos": 20},
        "source_genes": list(SOURCE_GENES),
        "target_tfs": list(TARGET_TFS),
    }
    provenance = stage_provenance(
        ROOT,
        "03_pkn",
        provenance_inputs,
        {"omnipath": omni_metadata, "cosmos": cosmos_metadata, "recon3d": recon_metadata, "workbench": {"row_count": len(workbench), "sha256": sha256(METABOLITE_RESOURCES / "metabolomics_workbench_kegg_hmdb.tsv")}},
        settings,
        code_paths,
    )
    provenance["diagnostics"] = {"counts_rows": len(counts)}
    write_stage_provenance(ANALYSIS, provenance)

    print("Stage 03 PKN build finished")
    print(counts.to_string(index=False))
    print()
    print(paths[["pkn", "leg", "path_cap", "source", "target", "shortest_length", "shortest_signs", "positive_route_present"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
