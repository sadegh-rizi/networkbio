"""Recon3D import, medium construction, and context-model helpers.

The project does not depend on COBRApy.  Recon3D is read directly from its
COBRA JSON representation and converted to the CORNETO incidence-graph format
used by ``MultiSampleFBA`` and ``MultiSampleIMAT``.
"""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Collection, Mapping, Sequence
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np
import pandas as pd

try:
    from corneto.io import graph_from_vertex_incidence
except ImportError:  # corneto 1.0.0rc8 keeps this reader private
    from corneto.io._metabolism import graph_from_vertex_incidence


EXCHANGE_PREFIX = "EX_"
# Recon3D boundary reactions: exchanges, sinks and demands.  All have the form
# "met -> " (coefficient -1), so a negative flux is an uptake.
SUPPLY_PREFIXES = ("SK_", "DM_")
# Essential amino acids absent from the ST003332 metabolomics panel, so they
# can never be "detected in the blanks"; every culture medium contains them.
# Without them BIOMASS_maintenance is infeasible in the defined media
# (Stage-06 run, 2026-10-14). The other essential amino acids are detected.
UNMEASURED_ESSENTIAL_UPTAKES = ("EX_phe__L_e", "EX_ile__L_e")
MEVALONATE_PREFIX = "mev__R_"
BOUNDARY_PREFIXES = (EXCHANGE_PREFIX, *SUPPLY_PREFIXES)
DEMAND_REACTION = "DM_chsterol_c"
BIOMASS_REACTIONS = ("BIOMASS_reaction", "BIOMASS_maintenance", "BIOMASS_maintenance_noTrTr")
INORGANIC_KEGG = {
    "O2": "C00007",
    "H2O": "C00001",
    "H+": "C00080",
    "Pi": "C00009",
    "Na+": "C00738",
    "K+": "C00238",
    "Cl-": "C00698",
    "Ca2+": "C00076",
    "Fe2+": "C14818",
    "SO4": "C00059",
    "NH4": "C00014",
    "CO2": "C00011",
}


def _split_values(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value if str(item)]
    return [item.strip() for item in re.split(r"[|;,]", str(value)) if item.strip()]


def _hgnc_entrez_symbols(hgnc: Any) -> dict[str, set[str]]:
    if isinstance(hgnc, (str, Path)):
        table = pd.read_csv(hgnc, sep="\t", dtype=str, keep_default_na=False)
        approved = table[table["status"].eq("Approved")]
        result: dict[str, set[str]] = defaultdict(set)
        for row in approved[["symbol", "entrez_id"]].itertuples(index=False):
            for entrez in _split_values(row.entrez_id):
                result[entrez].add(str(row.symbol))
        return dict(result)
    if isinstance(hgnc, pd.DataFrame):
        return _hgnc_entrez_symbols_from_table(hgnc)
    if isinstance(hgnc, Mapping) and "entrez" in hgnc:
        return {str(key): set(map(str, value)) for key, value in hgnc["entrez"].items()}
    raise TypeError("hgnc must be a path, HGNC DataFrame, or a mapping with an 'entrez' map")


def _hgnc_entrez_symbols_from_table(table: pd.DataFrame) -> dict[str, set[str]]:
    required = {"symbol", "status", "entrez_id"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"HGNC table lacks columns: {sorted(missing)}")
    result: dict[str, set[str]] = defaultdict(set)
    for row in table.loc[table["status"].eq("Approved"), ["symbol", "entrez_id"]].itertuples(index=False):
        for entrez in _split_values(row.entrez_id):
            result[entrez].add(str(row.symbol))
    return dict(result)


def _gene_number(gene_id: str) -> str:
    match = re.match(r"(\d+)(?:_AT\d+)?$", str(gene_id))
    return match.group(1) if match else str(gene_id)


def gpr_token(gene: Any) -> str:
    """Return a gene identifier that CORNETO's GPR parser reads as one name.

    CORNETO parses GPR rules as Python expressions.  A hyphen (``MT-ND3``) is
    read as subtraction and the reaction is left unscored, and a bare number
    (an Entrez ID without an approved symbol) is read as a constant whose
    value becomes the reaction score.  Non-identifier characters are therefore
    replaced by ``_`` and identifiers that do not start with a letter or
    underscore are prefixed with ``ENTREZ_``.
    """
    token = re.sub(r"[^A-Za-z0-9_]", "_", str(gene))
    return token if re.match(r"[A-Za-z_]", token) else f"ENTREZ_{token}"


def _rewrite_gpr(rule: str, token_map: Mapping[str, str]) -> str:
    return re.sub(r"\b\d+_AT\d+\b", lambda match: token_map.get(match.group(0), match.group(0)), str(rule or ""))


def load_recon3d(path: Path, hgnc: Any) -> tuple[Any, pd.DataFrame]:
    """Load Recon3D as a CORNETO graph and return its gene mapping table.

    HGNC mappings are accepted only when an Entrez ID has one approved symbol.
    Otherwise the numeric Entrez ID remains in the GPR, which avoids assigning
    an ambiguous reaction to an arbitrary symbol.
    """
    import json

    model = json.loads(Path(path).read_text())
    entrez_symbols = _hgnc_entrez_symbols(hgnc)
    token_map: dict[str, str] = {}
    gene_rows: list[dict[str, Any]] = []
    for gene in model.get("genes", []):
        model_id = str(gene["id"])
        entrez_ids = _split_values(gene.get("annotation", {}).get("ncbigene"))
        entrez = _gene_number(model_id)
        if entrez_ids:
            entrez = sorted(entrez_ids)[0]
        symbols = sorted({symbol for identifier in entrez_ids for symbol in entrez_symbols.get(identifier, set())})
        if len(symbols) == 1:
            mapped = symbols[0]
            status = "unique_approved_symbol"
        else:
            mapped = entrez
            status = "no_unique_approved_symbol" if symbols else "no_approved_symbol"
        token = gpr_token(mapped)
        token_map[model_id] = token
        token_map.setdefault(f"{entrez}_AT1", token)
        gene_rows.append({
            "model_gene_id": model_id,
            "entrez_id": entrez,
            "model_name": str(gene.get("name", "")),
            "mapped_gene": mapped,
            "gpr_token": token,
            "mapping_status": status,
            "approved_symbols": ";".join(symbols),
        })

    metabolite_ids = [str(item["id"]) for item in model["metabolites"]]
    reaction_ids = [str(item["id"]) for item in model["reactions"]]
    metabolite_index = {identifier: index for index, identifier in enumerate(metabolite_ids)}
    stoichiometry = np.zeros((len(metabolite_ids), len(reaction_ids)), dtype=float)
    for column, reaction in enumerate(model["reactions"]):
        for metabolite, coefficient in reaction.get("metabolites", {}).items():
            stoichiometry[metabolite_index[str(metabolite)], column] = float(coefficient)
    graph = graph_from_vertex_incidence(stoichiometry, metabolite_ids, reaction_ids)

    for metabolite in model["metabolites"]:
        identifier = str(metabolite["id"])
        graph.get_attr_vertex(identifier).update({
            "model_id": identifier,
            "name": metabolite.get("name", ""),
            "compartment": metabolite.get("compartment", ""),
            "kegg_ids": _split_values(metabolite.get("annotation", {}).get("kegg.compound")),
            "hmdb_ids": _split_values(metabolite.get("annotation", {}).get("hmdb")),
            "annotation": metabolite.get("annotation", {}),
        })
    for number, reaction in enumerate(model["reactions"], start=1):
        attributes = graph.get_attr_edge(number - 1)
        attributes.update({
            "reaction_number": number,
            "model_id": str(reaction["id"]),
            "name": reaction.get("name", ""),
            "subsystem": reaction.get("subsystem", ""),
            "default_lb": float(reaction["lower_bound"]),
            "default_ub": float(reaction["upper_bound"]),
            "GPR": _rewrite_gpr(str(reaction.get("gene_reaction_rule", "")), token_map),
            "gene_reaction_rule": str(reaction.get("gene_reaction_rule", "")),
            "stoichiometry": {str(k): float(v) for k, v in reaction.get("metabolites", {}).items()},
            "annotation": reaction.get("annotation", {}),
        })
    graph.get_graph_attributes().update({
        "model_id": model.get("id", "Recon3D"),
        "model_version": model.get("version", ""),
        "metabolite_ids": metabolite_ids,
        "reaction_ids": reaction_ids,
    })
    return graph, pd.DataFrame(gene_rows)


def reaction_table(model: Any) -> pd.DataFrame:
    """Return one row per reaction in graph order."""
    rows = []
    for index, attributes in enumerate(model.get_attr_edges(), start=1):
        rows.append({
            "reaction_number": int(attributes.get("reaction_number", index)),
            "reaction": str(attributes["id"]),
            "name": str(attributes.get("name", "")),
            "subsystem": str(attributes.get("subsystem", "")),
            "lower_bound": float(attributes.get("default_lb", 0.0)),
            "upper_bound": float(attributes.get("default_ub", 0.0)),
            "GPR": str(attributes.get("GPR", "")),
        })
    return pd.DataFrame(rows)


def model_stoichiometry(model: Any) -> tuple[np.ndarray, list[str], list[str]]:
    """Extract the signed incidence matrix stored on a loaded graph."""
    metabolites = [str(identifier) for identifier in model.vertices]
    reactions = [str(attributes["id"]) for attributes in model.get_attr_edges()]
    index = {identifier: row for row, identifier in enumerate(metabolites)}
    matrix = np.zeros((len(metabolites), len(reactions)), dtype=float)
    for column, attributes in enumerate(model.get_attr_edges()):
        for metabolite, coefficient in attributes.get("stoichiometry", {}).items():
            matrix[index[str(metabolite)], column] = float(coefficient)
    return matrix, metabolites, reactions


def add_cholesterol_demand(model: Any, reaction_id: str = DEMAND_REACTION) -> Any:
    """Add a non-reversible cholesterol demand reaction to a model graph."""
    if reaction_id in reaction_table(model)["reaction"].tolist():
        return model
    if "chsterol_c" not in model.vertices:
        raise ValueError("Recon3D has no cytosolic cholesterol metabolite")
    model.add_edge(
        {"chsterol_c"}, set(), id=reaction_id, reaction_number=model.num_edges + 1,
        model_id=reaction_id, name="Cholesterol demand", subsystem="Cholesterol metabolism",
        default_lb=0.0, default_ub=1000.0, GPR="", gene_reaction_rule="",
        stoichiometry={"chsterol_c": -1.0}, annotation={},
    )
    return model


def _exchange_reactions(model: Any) -> dict[str, list[dict[str, Any]]]:
    by_kegg: dict[str, list[dict[str, Any]]] = defaultdict(list)
    vertices = {str(identifier): model.get_attr_vertex(identifier) for identifier in model.vertices}
    for attributes in model.get_attr_edges():
        reaction = str(attributes.get("id", ""))
        stoich = attributes.get("stoichiometry", {})
        if not reaction.startswith(EXCHANGE_PREFIX) or len(stoich) != 1:
            continue
        metabolite = str(next(iter(stoich)))
        vertex = vertices.get(metabolite, {})
        for kegg in vertex.get("kegg_ids", []):
            by_kegg[str(kegg).removeprefix("cpd:")].append({
                "reaction": reaction,
                "metabolite": metabolite,
                "lower_bound": float(attributes.get("default_lb", 0.0)),
                "upper_bound": float(attributes.get("default_ub", 0.0)),
            })
    return dict(by_kegg)


def _choose_exchange(candidates: Sequence[dict[str, Any]]) -> dict[str, Any] | None:
    if not candidates:
        return None
    extracellular = [row for row in candidates if str(row["metabolite"]).endswith("_e")]
    return sorted(extracellular or list(candidates), key=lambda row: str(row["reaction"]))[0]


def medium_bounds(model: Any, variant: str, detected_kegg: Collection[str] = ()) -> dict[str, tuple[float, float]]:
    """Return exchange overrides for an open or defined Recon3D medium.

    The defined medium blocks uptake for every exchange, then reopens blank-
    detected compounds, registered inorganic compounds, and optionally
    extracellular cholesterol. Secretion remains at each reaction's default
    upper bound.

    Sink (``SK_``) and demand (``DM_``) reactions are boundaries as well: 92
    sinks and 3 demands in Recon3D can supply their metabolite (lower bound
    < 0).  The defined medium closes that supply direction too; otherwise, for
    example, ``SK_vitd3_c`` supplies vitamin D3, which the reversible
    ``VD3``/``PVD3`` reactions turn into 7-dehydrocholesterol and cholesterol
    with HMG-CoA reductase and cholesterol uptake blocked (Stage-06 check,
    2026-10-14).
    """
    label = str(variant).lower().replace("+", "_plus").replace("-", "_minus").replace(" ", "_")
    if label in {"open", "default"}:
        return {}
    if label not in {"defined_chol_plus", "defined_chol_minus", "defined", "chol_plus", "chol_minus"}:
        raise ValueError(f"Unknown medium variant: {variant}")
    table = reaction_table(model)
    exchanges = table[table["reaction"].str.startswith(EXCHANGE_PREFIX)]
    bounds = {str(row.reaction): (0.0, float(row.upper_bound)) for row in exchanges.itertuples(index=False)}
    supplies = table[table["reaction"].str.startswith(SUPPLY_PREFIXES) & table["lower_bound"].lt(0)]
    bounds.update({str(row.reaction): (0.0, float(row.upper_bound)) for row in supplies.itertuples(index=False)})
    by_kegg = _exchange_reactions(model)
    allowed_kegg = {str(value).removeprefix("cpd:") for value in detected_kegg}
    allowed_kegg.update(INORGANIC_KEGG.values())
    for kegg in allowed_kegg:
        exchange = _choose_exchange(by_kegg.get(kegg, []))
        if exchange is not None:
            bounds[exchange["reaction"]] = (min(-1000.0, exchange["lower_bound"]), exchange["upper_bound"])
    for reaction in UNMEASURED_ESSENTIAL_UPTAKES:
        if reaction in bounds:
            row = table.set_index("reaction").loc[reaction]
            bounds[reaction] = (min(-1000.0, float(row.lower_bound)), float(row.upper_bound))
    cholesterol = "EX_chsterol_e"
    if label in {"defined_chol_plus", "chol_plus"} and cholesterol in bounds:
        row = reaction_table(model).set_index("reaction").loc[cholesterol]
        bounds[cholesterol] = (min(-1000.0, float(row.lower_bound)), float(row.upper_bound))
    elif cholesterol in bounds:
        bounds[cholesterol] = (0.0, bounds[cholesterol][1])
    return bounds


def mevalonate_synthesis_sign(model: Any, reaction: str) -> int:
    """+1 if positive flux of ``reaction`` produces mevalonate, -1 if negative flux does.

    Recon3D writes cytosolic ``HMGCOAR`` as mevalonate -> HMG-CoA, so HMG-CoA
    reductase activity is negative flux there, while peroxisomal ``HMGCOARx``
    is written in the synthesis direction.
    """
    for attributes in model.get_attr_edges():
        if str(attributes.get("id")) == reaction:
            coefficient = sum(float(value) for key, value in attributes.get("stoichiometry", {}).items()
                              if str(key).startswith(MEVALONATE_PREFIX))
            if coefficient == 0:
                raise ValueError(f"{reaction} does not produce or consume mevalonate")
            return 1 if coefficient > 0 else -1
    raise ValueError(f"Unknown reaction: {reaction}")


def cap_synthesis_direction(bounds: tuple[float, float], sign: int, capacity: float, fraction: float = 0.1) -> tuple[float, float]:
    """Limit flux in the synthesis direction to ``fraction`` of ``capacity`` (>= 0)."""
    lower, upper = float(bounds[0]), float(bounds[1])
    cap = max(0.0, fraction * float(capacity))
    return (lower, min(upper, cap)) if sign > 0 else (max(lower, -cap), upper)


def gene_scores(
    expression: pd.DataFrame,
    groups: Mapping[str, str] | pd.Series,
    quantiles: tuple[float, float] = (0.25, 0.75),
    model_genes: Collection[str] = (),
    unscored_genes: Collection[str] = (),
) -> dict[str, dict[str, float]]:
    """Convert group-mean expression to the registered iMAT -1/0/+1 scores.

    Model genes absent from ``expression`` score -1 (not expressed), except
    ``unscored_genes`` (model genes without a unique approved symbol, whose
    expression cannot be looked up), which score 0.
    """
    if len(quantiles) != 2 or not 0 <= quantiles[0] < quantiles[1] <= 1:
        raise ValueError("quantiles must be an increasing pair between 0 and 1")
    labels = pd.Series(groups, index=list(groups.index) if isinstance(groups, pd.Series) else list(groups))
    labels = labels.reindex(expression.columns)
    if labels.isna().any():
        raise ValueError("Every expression column needs a group label")
    means = expression.T.groupby(labels).mean().T
    result: dict[str, dict[str, float]] = {}
    for group in means.columns:
        values = means[group].astype(float)
        low = float(values.quantile(quantiles[0]))
        high = float(values.quantile(quantiles[1]))
        scores = pd.Series(0.0, index=values.index)
        scores.loc[values.le(low)] = -1.0
        scores.loc[values.ge(high)] = 1.0
        score_map = {str(gene): float(score) for gene, score in scores.items()}
        unscored = {str(gene) for gene in unscored_genes}
        for gene in model_genes:
            if str(gene) not in score_map:
                score_map[str(gene)] = 0.0 if str(gene) in unscored else -1.0
        result[str(group)] = score_map
    return result


ACTIVITY_TOLERANCE = 1e-6


def activity_threshold(eps: float = 1e-2) -> float:
    """Smallest |flux| counted as active.

    iMAT places active reactions exactly at its lower bound ``eps``; the solver
    returns that value as 0.01 or as 0.0099999999999909. Comparing with ``eps``
    itself therefore splits the iMAT "on" state at random (defect D9, 2026-10-15).
    """
    return float(eps) - ACTIVITY_TOLERANCE


def is_active(values: Any, eps: float = 1e-2) -> Any:
    """Vectorised activity test with the D9 tolerance."""
    return np.abs(np.asarray(values, dtype=float)) >= activity_threshold(eps)


def context_reactions(fluxes: pd.DataFrame, eps: float = 1e-2) -> pd.DataFrame:
    """Annotate reaction fluxes with activity and direction.

    Activity and direction use the same tolerant threshold, so every active
    reaction has a direction (forward or reverse) and every inactive one has none.
    """
    required = {"condition", "reaction", "flux"}
    if not required.issubset(fluxes.columns):
        raise ValueError(f"Flux table lacks columns: {sorted(required - set(fluxes.columns))}")
    columns = ["condition", "reaction", "flux"]
    if "reaction_number" in fluxes.columns:
        columns.insert(1, "reaction_number")
    result = fluxes[columns].copy()
    result["flux"] = pd.to_numeric(result["flux"], errors="raise")
    threshold = activity_threshold(eps)
    result["active"] = result["flux"].abs().ge(threshold)
    result["direction"] = np.select(
        [result["flux"].ge(threshold), result["flux"].le(-threshold)],
        ["forward", "reverse"], default="inactive",
    )
    return result


def directional_reaction_nodes(fluxes: pd.DataFrame, pkn_nodes: Collection[str], eps: float = 1e-2) -> set[str]:
    """Map positive/negative Recon fluxes to COSMOS forward/reverse nodes."""
    annotated = context_reactions(fluxes, eps)
    if "reaction_number" not in annotated.columns:
        raise ValueError("Fluxes need reaction_number to map to COSMOS reaction nodes")
    directions: dict[str, set[str]] = defaultdict(set)
    for row in annotated[annotated["active"]].itertuples(index=False):
        directions[str(int(row.reaction_number))].add(str(row.direction))
    active: set[str] = set()
    for node in pkn_nodes:
        match = re.match(r"^Gene(\d+)__(.+)$", str(node))
        if not match:
            continue
        number = int(match.group(1))
        node_direction = "reverse" if match.group(2).endswith("_reverse") else "forward"
        if node_direction in directions.get(str(number), set()):
            active.add(str(node))
    return active


def prune_cosmos_to_context(pkn: pd.DataFrame, active: Collection[str]) -> pd.DataFrame:
    """Remove inactive directional ``Gene<N>__...`` nodes from a COSMOS PKN."""
    required = {"source", "target"}
    if not required.issubset(pkn.columns):
        raise ValueError(f"PKN lacks columns: {sorted(required - set(pkn.columns))}")
    active_nodes = {str(node) for node in active}

    def keep(node: Any) -> bool:
        text = str(node)
        return not text.startswith("Gene") or text in active_nodes

    result = pkn.loc[pkn["source"].map(keep) & pkn["target"].map(keep)].copy()
    return result.reset_index(drop=True)


def minimum_extra_exchanges(
    model: Any,
    bounds: Mapping[str, tuple[float, float]],
    objective_reaction: str,
    target: float,
    *,
    time_limit: int = 3600,
) -> dict[str, Any]:
    """Find the smallest set of currently blocked exchanges restoring a target."""
    S, metabolites, reactions = model_stoichiometry(model)
    table = reaction_table(model).set_index("reaction")
    if objective_reaction not in table.index:
        raise ValueError(f"Unknown objective reaction: {objective_reaction}")
    exchange_ids = [reaction for reaction in reactions if reaction.startswith(BOUNDARY_PREFIXES)]
    # Candidates: boundary reactions whose uptake the medium closed but the
    # model allows.  Their lower bound is relaxed to the model default and the
    # uptake is switched by ``opened``.  (Keeping the medium lower bound of 0
    # here made every candidate unusable and the MILP always infeasible.)
    candidates = [reaction for reaction in exchange_ids
                  if bounds.get(reaction, (float(table.loc[reaction, "lower_bound"]), 0.0))[0] >= 0
                  and float(table.loc[reaction, "lower_bound"]) < 0]
    v = cp.Variable(len(reactions))
    opened = cp.Variable(len(candidates), boolean=True)
    lower = np.array([bounds.get(reaction, (table.loc[reaction, "lower_bound"], table.loc[reaction, "upper_bound"]))[0] for reaction in reactions], dtype=float)
    upper = np.array([bounds.get(reaction, (table.loc[reaction, "lower_bound"], table.loc[reaction, "upper_bound"]))[1] for reaction in reactions], dtype=float)
    position = {reaction: index for index, reaction in enumerate(reactions)}
    for reaction in candidates:
        lower[position[reaction]] = float(table.loc[reaction, "lower_bound"])
    constraints = [S @ v == 0, v >= lower, v <= upper]
    for index, reaction in enumerate(candidates):
        constraints.append(v[position[reaction]] >= float(table.loc[reaction, "lower_bound"]) * opened[index])
    constraints.append(v[position[objective_reaction]] >= float(target))
    problem = cp.Problem(cp.Minimize(cp.sum(opened)), constraints)
    problem.solve(solver="HIGHS", time_limit=time_limit, mip_rel_gap=0.01)
    selected = [reaction for index, reaction in enumerate(candidates) if opened.value is not None and float(opened.value[index]) > 0.5]
    return {"status": str(problem.status), "objective_value": problem.value, "selected_exchanges": selected}
