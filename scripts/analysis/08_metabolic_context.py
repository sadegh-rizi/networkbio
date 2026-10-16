"""Stage 06: flux-consistent Recon3D context models and COSMOS comparison.

The real-data solve is intentionally explicit and resumable.  Expression is
converted to iMAT scores at the cell-line group level, while all flux outputs
are labelled as capacities because no measured fluxes are available.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from gem import (  # noqa: E402
    BIOMASS_REACTIONS,
    BOUNDARY_PREFIXES,
    DEMAND_REACTION,
    cap_synthesis_direction,
    add_cholesterol_demand,
    context_reactions,
    directional_reaction_nodes,
    gene_scores,
    gpr_token,
    is_active,
    load_recon3d,
    medium_bounds,
    mevalonate_synthesis_sign,
    minimum_extra_exchanges,
    prune_cosmos_to_context,
    reaction_table,
)
from gem_worker import (  # noqa: E402
    ACCEPTANCE_GAP,
    SOLVER,
    SOLVER_RELATIVE_GAP,
    SOLVER_SEED,
    SOLVER_THREADS,
    make_task,
    solve_in_process,
    solve_in_subprocess,
)
from preprocess import sha256, stage_provenance, write_stage_provenance  # noqa: E402


MODEL = ROOT / "data/resources/pkn/Recon3D.json"
HGNC = ROOT / "data/resources/hgnc/hgnc_complete_set_2026-10-06.txt"
INPUTS = ROOT / "results/ionescu_corneto/00_inputs"
PREPROCESSED = ROOT / "results/ionescu_corneto/01_preprocessing"
PKN = ROOT / "results/ionescu_corneto/03_pkn"
OUT = ROOT / "results/ionescu_corneto/06_validation/metabolic_context"
PLAN = ROOT / "doc/decisions/2026-10-13-imat-fba-plan.md"
EPS = 1e-2
TIME_LIMIT = 3600
MEDIUMS = ("open", "defined_chol_plus", "defined_chol_minus")
# FBA readouts run in this process with a cached model (addendum D9); iMAT always
# runs in a subprocess. --fba-subprocess restores one process per FBA solve.
_FBA_IN_PROCESS = True


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stage 06 Recon3D FBA/iMAT context models")
    parser.add_argument("--resume", action="store_true", help="reuse matching worker tasks")
    parser.add_argument("--mediums", nargs="+", choices=MEDIUMS, default=list(MEDIUMS))
    parser.add_argument("--time-limit", type=int, default=TIME_LIMIT)
    parser.add_argument("--fba-subprocess", action="store_true",
                        help="run every FBA solve in its own process (slow; default: in process)")
    return parser.parse_args(argv)


def _blank_kegg(root: Path) -> set[str]:
    features = pd.read_csv(root / INPUTS / "metabolomics/ST003332/features.tsv", sep="\t", dtype=str)
    samples = pd.read_csv(root / INPUTS / "metabolomics/ST003332/samples.tsv", sep="\t", dtype=str)
    values = pd.read_csv(root / INPUTS / "metabolomics/ST003332/values.tsv.gz", sep="\t", index_col=0)
    blanks = samples.loc[samples["group"].eq("no_cell_blank"), "sample_id"]
    detected = values.loc[:, blanks].notna().any(axis=1) & values.loc[:, blanks].gt(0).any(axis=1)
    return set(features.loc[features["feature_id"].isin(values.index[detected]), "kegg_id"].dropna().astype(str))


def _groups(root: Path, expression: pd.DataFrame) -> pd.Series:
    samples = pd.read_csv(root / PREPROCESSED / "rnaseq/samples.tsv", sep="\t").set_index("column")
    groups = samples["group"].reindex(expression.columns)
    if groups.isna().any() or groups.value_counts().to_dict() != {"Ctrl": 3, "PMS": 3}:
        raise ValueError("Expected the six baseline RNA-seq lines: 3 Ctrl and 3 PMS")
    return groups


def _write_medium_table(model: Any, bounds: dict[str, tuple[float, float]], path: Path, detected_kegg: set[str]) -> None:
    table = reaction_table(model)
    table = table[table["reaction"].str.startswith(("EX_", "SK_", "DM_"))].copy()
    table["medium_lower_bound"] = table["reaction"].map(lambda reaction: bounds.get(reaction, (np.nan, np.nan))[0])
    table["medium_upper_bound"] = table["reaction"].map(lambda reaction: bounds.get(reaction, (np.nan, np.nan))[1])
    table["uptake_allowed"] = table["medium_lower_bound"].lt(0)
    table["detected_blank_kegg"] = ",".join(sorted(detected_kegg))
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)


def _task(
    root: Path,
    *,
    label: str,
    method: str,
    gene_expression: dict[str, dict[str, float]] | None = None,
    objectives: dict[str, dict[str, float]] | None = None,
    bounds: dict[str, tuple[float, float]] | None = None,
    output: Path,
    conditions: tuple[str, ...] = ("Ctrl", "PMS"),
    resume: bool,
    time_limit: int,
    add_demand: bool = False,
    gap: float = SOLVER_RELATIVE_GAP,
    in_process: bool = False,
) -> dict[str, Any]:
    reaction_bounds = {condition: dict(bounds or {}) for condition in conditions}
    if objectives is None:
        objectives = {condition: {} for condition in conditions}
    task = make_task(
        label=label, model_path=MODEL, hgnc_path=HGNC, method=method,
        gene_expression=gene_expression, objectives=objectives,
        reaction_bounds=reaction_bounds, time_limit=time_limit,
        solver_relative_gap=gap, eps=EPS, add_demand=add_demand, root=root,
    )
    if in_process and method == "fba":
        return solve_in_process(task, output, resume=resume, root=root)
    return solve_in_subprocess(task, output, resume=resume, root=root)


def _fba(
    root: Path,
    *,
    label: str,
    output: Path,
    reaction: str,
    bounds: dict[str, tuple[float, float]],
    condition: str,
    resume: bool,
    time_limit: int,
    add_demand: bool = False,
    coefficient: float = -1.0,
) -> dict[str, Any]:
    result = _task(
        root, label=label, method="fba", objectives={condition: {reaction: coefficient}},
        bounds=bounds, conditions=(condition,), output=output, resume=resume,
        time_limit=time_limit, add_demand=add_demand, in_process=_FBA_IN_PROCESS,
    )
    return result


def _flux_value(result: dict[str, Any], condition: str, reaction: str) -> float:
    fluxes = result.get("fluxes", pd.DataFrame())
    if fluxes.empty:
        return float("nan")
    rows = fluxes[(fluxes["condition"].eq(condition)) & fluxes["reaction"].eq(reaction)]
    return float(rows["flux"].iloc[0]) if len(rows) else float("nan")


def _model_checks(root: Path, model: Any, output: Path, resume: bool, time_limit: int,
                  chol_minus: dict[str, tuple[float, float]], chol_plus: dict[str, tuple[float, float]]) -> pd.DataFrame:
    """Model sanity checks.

    The sterol checks use the defined media, not the open medium: with every
    uptake open, cholesterol sulfate, cholesterol esters and 7-dehydro-
    cholesterol supply cholesterol even when free-cholesterol uptake and
    HMG-CoA reductase are blocked (first Stage-06 run, 2026-10-14).
    """
    statin = {"HMGCOAR": (0.0, 0.0), "HMGCOARx": (0.0, 0.0)}
    rows = [
        {"check": "reaction_count", "value": model.num_edges, "expected": 10600, "passed": model.num_edges == 10600},
        {"check": "metabolite_count", "value": model.num_vertices, "expected": 5835, "passed": model.num_vertices == 5835},
    ]
    open_biomass = _fba(root, label="open biomass maintenance", output=output / "checks/open_biomass",
                        reaction="BIOMASS_maintenance", bounds={}, condition="open", resume=resume, time_limit=time_limit)
    biomass = _flux_value(open_biomass, "open", "BIOMASS_maintenance")
    rows.append({"check": "open_biomass_maintenance", "value": biomass, "expected": "> 0", "passed": bool(open_biomass["accepted"] and biomass > 0)})
    blocked = {**chol_minus, **statin}
    hmg_blocked = _fba(root, label="sterol sanity blocked", output=output / "checks/sterol_blocked",
                       reaction=DEMAND_REACTION, bounds=blocked, condition="blocked", resume=resume,
                       time_limit=time_limit, add_demand=True)
    hmg_allowed = _fba(root, label="sterol sanity synthesis", output=output / "checks/sterol_synthesis",
                        reaction=DEMAND_REACTION, bounds=dict(chol_minus), condition="synthesis",
                        resume=resume, time_limit=time_limit, add_demand=True)
    uptake_allowed = _fba(root, label="sterol sanity uptake", output=output / "checks/sterol_uptake",
                          reaction=DEMAND_REACTION, bounds={**chol_plus, **statin},
                          condition="uptake", resume=resume, time_limit=time_limit, add_demand=True)
    blocked_value = _flux_value(hmg_blocked, "blocked", DEMAND_REACTION)
    synthesis_value = _flux_value(hmg_allowed, "synthesis", DEMAND_REACTION)
    uptake_value = _flux_value(uptake_allowed, "uptake", DEMAND_REACTION)
    rows.extend([
        {"check": "cholesterol_demand_defined_chol_minus_hmgcr_blocked", "value": blocked_value, "expected": "0", "passed": bool(hmg_blocked["accepted"] and abs(blocked_value) <= EPS)},
        {"check": "cholesterol_demand_defined_chol_minus", "value": synthesis_value, "expected": "> 0", "passed": bool(hmg_allowed["accepted"] and synthesis_value > EPS)},
        {"check": "cholesterol_demand_defined_chol_plus_hmgcr_blocked", "value": uptake_value, "expected": "> 0", "passed": bool(uptake_allowed["accepted"] and uptake_value > EPS)},
    ])
    return pd.DataFrame(rows)


def _medium_support(root: Path, model: Any, output: Path, mediums: dict[str, dict[str, tuple[float, float]]], open_max: float, resume: bool, time_limit: int) -> pd.DataFrame:
    rows = []
    target = 0.1 * open_max
    for variant, bounds in mediums.items():
        result = _fba(root, label=f"{variant} biomass maintenance", output=output / f"medium_support/{variant}",
                      reaction="BIOMASS_maintenance", bounds=bounds, condition=variant, resume=resume, time_limit=time_limit)
        value = _flux_value(result, variant, "BIOMASS_maintenance")
        extra = {"status": "not_needed", "selected_exchanges": ""}
        if not (result["accepted"] and value >= target):
            extra = minimum_extra_exchanges(model, bounds, "BIOMASS_maintenance", target, time_limit=time_limit)
            extra = {"status": extra["status"], "selected_exchanges": ",".join(extra["selected_exchanges"])}
        rows.append({"medium": variant, "open_biomass_maximum": open_max, "target_10_percent": target,
                     "defined_biomass_value": value, "solve_status": result["status"], **extra})
    return pd.DataFrame(rows)


def _active_context_bounds(model: Any, base_bounds: dict[str, tuple[float, float]], fluxes: pd.DataFrame, condition: str) -> dict[str, tuple[float, float]]:
    active = set(context_reactions(fluxes[fluxes["condition"].eq(condition)], EPS).loc[lambda table: table["active"], "reaction"])
    table = reaction_table(model)
    bounds = {}
    for row in table.itertuples(index=False):
        # Boundaries (exchange, sink, demand) keep the medium bounds whether or
        # not iMAT used them: resetting an active sink to its Recon3D default
        # would reopen the supply the defined medium closes.
        if row.reaction.startswith(BOUNDARY_PREFIXES) or row.reaction == DEMAND_REACTION:
            bounds[row.reaction] = base_bounds.get(row.reaction, (float(row.lower_bound), float(row.upper_bound)))
        elif row.reaction not in active:
            bounds[row.reaction] = (0.0, 0.0)
        else:
            bounds[row.reaction] = (float(row.lower_bound), float(row.upper_bound))
    return bounds


def _sterol_flux_ranges(root: Path, model: Any, output: Path, fluxes: pd.DataFrame, medium_bounds_primary: dict[str, tuple[float, float]], canonical: pd.DataFrame, resume: bool, time_limit: int) -> pd.DataFrame:
    reactions = ["HMGCOAR", "HMGCOARx"] + sorted(set(canonical["recon_id"].dropna().astype(str)))
    rows = []
    for condition in ("Ctrl", "PMS"):
        bounds = _active_context_bounds(model, medium_bounds_primary, fluxes, condition)
        for reaction in reactions:
            if reaction not in reaction_table(model)["reaction"].values:
                continue
            maximum = _fba(root, label=f"{condition} max {reaction}", output=output / f"flux_range/{condition}/{reaction}/max",
                           reaction=reaction, bounds=bounds, condition=condition, resume=resume, time_limit=time_limit)
            minimum = _fba(root, label=f"{condition} min {reaction}", output=output / f"flux_range/{condition}/{reaction}/min",
                           reaction=reaction, bounds=bounds, condition=condition, resume=resume, time_limit=time_limit, coefficient=1.0)
            minimum_flux = _flux_value(minimum, condition, reaction)
            maximum_flux = _flux_value(maximum, condition, reaction)
            synthesis = np.nan
            if reaction in ("HMGCOAR", "HMGCOARx"):
                sign = mevalonate_synthesis_sign(model, reaction)
                synthesis = max(0.0, maximum_flux if sign > 0 else -minimum_flux)
            rows.append({"condition": condition, "reaction": reaction,
                         "minimum_flux": minimum_flux,
                         "maximum_flux": maximum_flux,
                         "mevalonate_synthesis_capacity": synthesis,
                         "minimum_status": minimum["status"], "maximum_status": maximum["status"],
                         "accepted": bool(minimum["accepted"] and maximum["accepted"])})
    ranges = pd.DataFrame(rows)
    step_rows = []
    for condition in ("Ctrl", "PMS"):
        for step, group in canonical.groupby("step_order"):
            known = ranges[(ranges["condition"].eq(condition)) & ranges["reaction"].isin(group["recon_id"].astype(str))]
            step_rows.append({"condition": condition, "step_order": step, "step_has_flux": bool(is_active(known["maximum_flux"], EPS).any() or is_active(known["minimum_flux"], EPS).any())})
    return pd.concat([ranges, pd.DataFrame(step_rows).assign(row_type="ordered_step")], ignore_index=True, sort=False)


def _cholesterol_tasks(root: Path, model: Any, output: Path, fluxes: pd.DataFrame, primary_bounds: dict[str, tuple[float, float]], resume: bool, time_limit: int, statin: bool = False, range_table: pd.DataFrame | None = None) -> pd.DataFrame:
    rows = []
    for condition in ("Ctrl", "PMS"):
        bounds = _active_context_bounds(model, primary_bounds, fluxes, condition)
        variants = {"defined": dict(bounds), "cholesterol_blocked": dict(bounds)}
        variants["cholesterol_blocked"]["EX_chsterol_e"] = (0.0, variants["cholesterol_blocked"].get("EX_chsterol_e", (0.0, 1000.0))[1])
        if statin:
            ranges = range_table if range_table is not None else pd.DataFrame()
            for reaction in ("HMGCOAR", "HMGCOARx"):
                maximums = ranges[(ranges.get("condition", pd.Series(dtype=str)).eq(condition)) & ranges.get("reaction", pd.Series(dtype=str)).eq(reaction)]
                if len(maximums):
                    # Cap HMG-CoA reductase in its mevalonate-forming direction,
                    # which is negative flux for Recon3D's HMGCOAR.
                    sign = mevalonate_synthesis_sign(model, reaction)
                    capacity = float(maximums["mevalonate_synthesis_capacity"].iloc[0])
                    for name in ("defined", "cholesterol_blocked"):
                        current = variants[name].get(reaction, (0.0, 1000.0))
                        variants[name][reaction] = cap_synthesis_direction(current, sign, capacity, 0.1)
        for variant, variant_bounds in variants.items():
            result = _fba(root, label=f"{condition} cholesterol {variant}", output=output / f"{condition}/{variant}",
                          reaction=DEMAND_REACTION, bounds=variant_bounds, condition=condition, resume=resume,
                          time_limit=time_limit, add_demand=True)
            rows.append({"condition": condition, "task": variant, "simvastatin_cap": statin,
                         "HMGCOAR_bounds": str(variant_bounds.get("HMGCOAR", "")),
                         "HMGCOARx_bounds": str(variant_bounds.get("HMGCOARx", "")),
                         "demand_flux": _flux_value(result, condition, DEMAND_REACTION), "status": result["status"],
                         "accepted": result["accepted"]})
    return pd.DataFrame(rows)


def _differential_activity(fluxes: pd.DataFrame, model: Any) -> tuple[pd.DataFrame, pd.DataFrame]:
    annotated = context_reactions(fluxes, EPS)
    active = {condition: set(table.loc[table["active"], "reaction"]) for condition, table in annotated.groupby("condition")}
    differential = []
    table = reaction_table(model).set_index("reaction")
    for reaction in sorted(set(active.get("Ctrl", set())) | set(active.get("PMS", set()))):
        differential.append({"reaction": reaction, "active_ctrl": reaction in active.get("Ctrl", set()),
                             "active_pms": reaction in active.get("PMS", set()),
                             "differential": "PMS_only" if reaction in active.get("PMS", set()) - active.get("Ctrl", set()) else "Ctrl_only" if reaction in active.get("Ctrl", set()) - active.get("PMS", set()) else "shared",
                             "subsystem": table.loc[reaction, "subsystem"] if reaction in table.index else ""})
    differential_table = pd.DataFrame(differential)
    universe = table[table["GPR"].ne("")]
    rows = []
    for direction in ("PMS_only", "Ctrl_only"):
        selected = set(differential_table.loc[differential_table["differential"].eq(direction), "reaction"]) & set(universe.index)
        for subsystem, group in universe.groupby("subsystem"):
            M, n, N, k = len(universe), len(group), len(selected), len(set(group.index) & selected)
            rows.append({"direction": direction, "subsystem": subsystem, "universe": M, "subsystem_reactions": n,
                         "selected_reactions": N, "overlap": k, "p_value": hypergeom.sf(k - 1, M, n, N) if N else 1.0})
    enrichment = pd.DataFrame(rows)
    if not enrichment.empty:
        enrichment["padj_bh"] = multipletests(enrichment["p_value"], method="fdr_bh")[1]
    return differential_table, enrichment


def _signed_vs_flux(root: Path, fluxes: pd.DataFrame, canonical: pd.DataFrame, primary_lambda: float = 0.1) -> pd.DataFrame:
    pkn = pd.read_csv(root / PKN / "cosmos/primary/pkn.tsv", sep="\t", dtype=str)
    stage04_dir = root / "results/ionescu_corneto/04_inference/cosmos_multi" / f"lambda_{primary_lambda:g}"
    edges = pd.read_csv(stage04_dir / "edges.tsv", sep="\t", dtype=str) if (stage04_dir / "edges.tsv").exists() else pd.DataFrame(columns=["source", "target"])
    used = set(edges["source"]) | set(edges["target"])
    reaction_nodes = {node for node in set(pkn["source"]) | set(pkn["target"]) if str(node).startswith("Gene")}
    rows = []
    for reaction, group in canonical.groupby("recon_id"):
        nodes = {node for node in reaction_nodes if re.match(rf"^Gene\d+__.*", node) and str(node).split("__", 1)[0] == str(group["reaction_node"].iloc[0]).split("__", 1)[0]}
        subset = fluxes[fluxes["reaction"].eq(str(reaction))]
        rows.append({"reaction": reaction, "reaction_name": group["reaction_name"].iloc[0], "present_in_cosmos": bool(nodes),
                     "used_in_stage04_joint": bool(nodes & used), "active_ctrl": bool(is_active(subset.loc[subset["condition"].eq("Ctrl"), "flux"], EPS).any()),
                     "active_pms": bool(is_active(subset.loc[subset["condition"].eq("PMS"), "flux"], EPS).any()),
                     "minimum_flux": float(subset["flux"].min()) if len(subset) else np.nan,
                     "maximum_flux": float(subset["flux"].max()) if len(subset) else np.nan})
    if "HMR_1916" in set(fluxes["reaction"]):
        subset = fluxes[fluxes["reaction"].eq("HMR_1916")]
        hmr_numbers = set(pd.to_numeric(subset["reaction_number"], errors="coerce").dropna().astype(int))
        hmr_nodes = {
            node for node in reaction_nodes
            if (match := re.match(r"^Gene(\d+)__", node)) and int(match.group(1)) in hmr_numbers
        }
        rows.append({"reaction": "HMR_1916", "reaction_name": "HDLBP cholesterol transport", "present_in_cosmos": bool(hmr_nodes),
                     "used_in_stage04_joint": bool(hmr_nodes & used), "active_ctrl": bool(is_active(subset.loc[subset["condition"].eq("Ctrl"), "flux"], EPS).any()),
                     "active_pms": bool(is_active(subset.loc[subset["condition"].eq("PMS"), "flux"], EPS).any()),
                     "minimum_flux": float(subset["flux"].min()) if len(subset) else np.nan, "maximum_flux": float(subset["flux"].max()) if len(subset) else np.nan})
    return pd.DataFrame(rows)


def _context_pkn(root: Path, output: Path, fluxes: pd.DataFrame) -> pd.DataFrame:
    pkn_path = root / PKN / "cosmos/primary/pkn.tsv"
    pkn = pd.read_csv(pkn_path, sep="\t", dtype={"source": str, "target": str})
    active = directional_reaction_nodes(fluxes, set(pkn["source"]) | set(pkn["target"]), EPS)
    pruned = prune_cosmos_to_context(pkn, active)
    output.mkdir(parents=True, exist_ok=True)
    pruned.to_csv(output / "pkn.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"step": "input", "n_edges": len(pkn), "n_nodes": len(set(pkn["source"]) | set(pkn["target"]))},
        {"step": "directional_iMAT_context", "n_edges": len(pruned), "n_nodes": len(set(pruned["source"]) | set(pruned["target"])), "n_active_reaction_nodes": len(active)},
    ]).to_csv(output / "filter_log.tsv", sep="\t", index=False)
    return pruned


def main(root: Path = ROOT, argv: list[str] | None = None) -> int:
    global _FBA_IN_PROCESS
    args = _parse_args(argv)
    _FBA_IN_PROCESS = not args.fba_subprocess
    root = Path(root)
    output = root / OUT
    if output.exists() and any(output.iterdir()) and not args.resume:
        raise ValueError(f"{output} is not empty; rerun with --resume or archive it")
    output.mkdir(parents=True, exist_ok=True)
    model, gene_map = load_recon3d(root / MODEL, root / HGNC)
    gene_map.to_csv(output / "gene_map.tsv", sep="\t", index=False)
    expression = pd.read_csv(root / PREPROCESSED / "rnaseq/log2cpm.tsv", sep="\t", index_col=0)
    groups = _groups(root, expression)
    # Expression keys must use the same parser-safe tokens as the GPR rules
    # (e.g. MT-ND3 -> MT_ND3); see gem.gpr_token.
    expression.index = expression.index.map(gpr_token)
    if expression.index.duplicated().any():
        raise ValueError(f"gpr_token collisions in expression: {sorted(set(expression.index[expression.index.duplicated()]))[:10]}")
    model_genes = set(gene_map["gpr_token"].astype(str))
    unscored = set(gene_map.loc[gene_map["mapping_status"].ne("unique_approved_symbol"), "gpr_token"].astype(str))
    scores = gene_scores(expression, groups, quantiles=(0.25, 0.75), model_genes=model_genes, unscored_genes=unscored)
    # MultiSampleIMAT validates gene_expression identifiers against GPR tokens;
    # model genes that occur in no reaction are retained in gene_map but are not
    # sent as features to CORNETO.
    known_genes = set(re.findall(r"\b(?!and\b|or\b)[A-Za-z0-9_]+\b", " ".join(reaction_table(model)["GPR"])))
    scores = {condition: {gene: value for gene, value in values.items() if gene in known_genes}
              for condition, values in scores.items()}
    detected_kegg = _blank_kegg(root)
    medium_maps = {variant: medium_bounds(model, variant, detected_kegg) for variant in args.mediums}
    for variant, bounds in medium_maps.items():
        _write_medium_table(model, bounds, output / f"medium/{variant}.tsv", detected_kegg)

    sterol_media = {variant: medium_maps.get(variant) or medium_bounds(model, variant, detected_kegg)
                    for variant in ("defined_chol_minus", "defined_chol_plus")}
    checks = _model_checks(root, model, output, args.resume, args.time_limit,
                           sterol_media["defined_chol_minus"], sterol_media["defined_chol_plus"])
    checks.to_csv(output / "model_checks.tsv", sep="\t", index=False)
    if not checks["passed"].all():
        failed = ", ".join(checks.loc[~checks["passed"].astype(bool), "check"])
        raise RuntimeError(f"Model checks failed ({failed}); see model_checks.tsv. No context model was interpreted.")
    open_result = _fba(root, label="open biomass for medium target", output=output / "checks/open_biomass_target",
                       reaction="BIOMASS_maintenance", bounds={}, condition="open", resume=args.resume, time_limit=args.time_limit)
    open_max = _flux_value(open_result, "open", "BIOMASS_maintenance")
    support = _medium_support(root, model, output, medium_maps, open_max, args.resume, args.time_limit)
    support.to_csv(output / "medium_support.tsv", sep="\t", index=False)

    imat_results = {}
    for variant, bounds in medium_maps.items():
        imat_results[variant] = _task(root, label=f"iMAT {variant}", method="imat", gene_expression=scores,
                                      bounds=bounds, output=output / f"imat/{variant}", resume=args.resume,
                                      time_limit=args.time_limit)
    primary_variant = "defined_chol_plus" if "defined_chol_plus" in imat_results else args.mediums[0]
    primary = imat_results[primary_variant]
    primary_fluxes = primary["fluxes"]
    if not primary["accepted"]:
        raise RuntimeError("Primary iMAT solve was not accepted; no context or readout was interpreted")
    (output / "summary").mkdir(parents=True, exist_ok=True)
    annotated = context_reactions(primary_fluxes, EPS)
    annotated.to_csv(output / "summary/reaction_activity.tsv", sep="\t", index=False)
    differential, enrichment = _differential_activity(primary_fluxes, model)
    differential.to_csv(output / "summary/differential_reactions.tsv", sep="\t", index=False)
    enrichment.to_csv(output / "summary/subsystem_enrichment.tsv", sep="\t", index=False)
    canonical = pd.read_csv(root / PKN / "aim5_route/canonical_steps.tsv", sep="\t", dtype=str)
    ranges = _sterol_flux_ranges(root, model, output, primary_fluxes, medium_maps[primary_variant], canonical, args.resume, args.time_limit)
    ranges.to_csv(output / "summary/sterol_flux_range.tsv", sep="\t", index=False)
    tasks = _cholesterol_tasks(root, model, output / "summary/cholesterol", primary_fluxes, medium_maps[primary_variant], args.resume, args.time_limit)
    tasks.to_csv(output / "summary/cholesterol_tasks.tsv", sep="\t", index=False)
    statin = _cholesterol_tasks(root, model, output / "summary/simvastatin", primary_fluxes, medium_maps[primary_variant], args.resume, args.time_limit, statin=True, range_table=ranges)
    statin.to_csv(output / "summary/simvastatin_task.tsv", sep="\t", index=False)
    signed = _signed_vs_flux(root, primary_fluxes, canonical)
    signed.to_csv(output / "summary/signed_vs_flux.tsv", sep="\t", index=False)
    _context_pkn(root, output / "context_cosmos", primary_fluxes)

    inputs = [root / MODEL, root / HGNC, root / PREPROCESSED / "rnaseq/log2cpm.tsv", root / PREPROCESSED / "rnaseq/samples.tsv",
              root / INPUTS / "metabolomics/ST003332/features.tsv", root / INPUTS / "metabolomics/ST003332/samples.tsv",
              root / INPUTS / "metabolomics/ST003332/values.tsv.gz", root / PKN / "aim5_route/canonical_steps.tsv",
              root / PKN / "cosmos/primary/pkn.tsv", PLAN]
    code = [root / "src/gem.py", root / "src/gem_worker.py", root / "scripts/analysis/08_metabolic_context.py"]
    provenance = stage_provenance(root, "metabolic_context", [path for path in inputs if path.exists()],
                                  {"recon3d": {"sha256": sha256(root / MODEL)}, "hgnc": {"sha256": sha256(root / HGNC)}},
                                  {"mediums": list(args.mediums), "primary_medium": primary_variant, "quantiles": [0.25, 0.75],
                                   "eps": EPS, "time_limit": args.time_limit, "solver": SOLVER, "threads": SOLVER_THREADS,
                                   "seed": SOLVER_SEED, "relative_gap": SOLVER_RELATIVE_GAP, "acceptance_gap": ACCEPTANCE_GAP,
                                   "activity_threshold": EPS - 1e-6, "fba_in_process": _FBA_IN_PROCESS}, code)
    write_stage_provenance(output, provenance)
    print(f"Stage 06 metabolic context finished: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
