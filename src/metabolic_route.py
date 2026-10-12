"""Aim-5 canonical cholesterol-route audit and graph variants.

The route audit is structural: it compares Recon3D's reaction definitions
with the two cached COSMOS PKN variants.  It does not use samples or infer
biological activity.  Reaction direction is the registered pathway direction,
not necessarily Recon3D's written stoichiometric direction (HMGCOAR and the
FDFT1 half-reactions are important examples).
"""

from __future__ import annotations

from collections import defaultdict, deque
from pathlib import Path
import hashlib
import json
import re
import random
from typing import Any, Iterable, Mapping
from urllib.request import Request, urlopen

import networkx as nx
import numpy as np
import pandas as pd

from pkn import (
    Edge,
    positive_shortest_lengths,
    rewire_type_preserving_signed_graph,
    sha256,
)

CANONICAL_GENES = (
    "HMGCR", "MVK", "PMVK", "MVD", "IDI1", "IDI2", "FDPS", "FDFT1",
    "SQLE", "LSS", "CYP51A1", "TM7SF2", "LBR", "MSMO1", "NSDHL",
    "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24",
)
COFACTOR_BIGG_IDS = frozenset({
    "h", "h2o", "atp", "adp", "amp", "gtp", "gdp", "utp", "udp", "ctp", "cdp",
    "nad", "nadh", "nadp", "nadph", "fad", "fadh2", "coa", "pi", "ppi", "o2",
    "co2", "hco3", "nh4", "h2o2", "na1", "k", "cl", "q10", "q10h2", "amet", "ahcys",
})
KEGG_PATHWAYS = ("hsa00900", "hsa00100")
KEGG_URL = "https://rest.kegg.jp/link/hsa/{pathway}"
CHOLESTEROL_NODE_PREFIX = "Metab__HMDB0000067_"

# Base IDs are Recon3D/BiGG IDs without the compartment suffix.  The pairs
# are deliberately explicit so that a co-substrate cannot create a fake route.
_STEP_SPECS = (
    (1, "HMG-CoA", "Mevalonate", "hmgcoa", "mev__R", ("HMGCR",), "core"),
    (2, "Mevalonate", "Phosphomevalonate", "mev__R", "5pmev", ("MVK",), "core"),
    (3, "Phosphomevalonate", "Diphosphomevalonate", "5pmev", "5dpmev", ("PMVK",), "core"),
    (4, "Diphosphomevalonate", "IPP", "5dpmev", "ipdp", ("MVD",), "core"),
    (5, "IPP", "DMAPP", "ipdp", "dmpp", ("IDI1", "IDI2"), "core"),
    (6, "DMAPP", "GPP", "dmpp", "grdp", ("FDPS",), "core"),
    (7, "GPP", "FPP", "grdp", "frdp", ("FDPS",), "core"),
    # Recon3D also contains a direct DMAPP -> FPP representation in cytosol.
    (7.5, "DMAPP", "FPP", "dmpp", "frdp", ("FDPS",), "core_direct"),
    (8, "FPP", "Presqualene diphosphate", "frdp", "HC01118", ("FDFT1",), "core"),
    (9, "Presqualene diphosphate", "Squalene", "HC01118", "sql", ("FDFT1",), "core"),
    (10, "Squalene", "2,3-oxidosqualene", "sql", "Ssq23epx", ("SQLE",), "core"),
    (11, "2,3-oxidosqualene", "Lanosterol", "Ssq23epx", "lanost", ("LSS",), "core"),
    (12, "Lanosterol", "M00939", "lanost", "M00939", ("CYP51A1",), "core"),
    (13, "M00939", "M00937", "M00939", "M00937", ("CYP51A1",), "core"),
    (14, "M00937", "44mctr", "M00937", "44mctr", ("CYP51A1",), "core"),
    (15, "44mctr", "44mzym", "44mctr", "44mzym", ("TM7SF2", "LBR"), "core"),
    (16, "44mzym", "M00961", "44mzym", "M00961", ("MSMO1",), "core"),
    (17, "M00961", "M00957", "M00961", "M00957", ("MSMO1",), "core"),
    (18, "M00957", "4mzym_int1", "M00957", "4mzym_int1", ("MSMO1",), "core"),
    (19, "4mzym_int1", "4mzym_int2", "4mzym_int1", "4mzym_int2", ("NSDHL",), "core"),
    (20, "4mzym_int2", "HC02110", "4mzym_int2", "HC02110", ("HSD17B7",), "core"),
    (21, "HC02110", "M00963", "HC02110", "M00963", ("MSMO1",), "core"),
    (22, "M00963", "M00959", "M00963", "M00959", ("MSMO1",), "core"),
    (23, "M00959", "M00955", "M00959", "M00955", ("MSMO1",), "core"),
    (24, "M00955", "M01067", "M00955", "M01067", ("NSDHL",), "core"),
    (25, "M01067", "Zymosterol", "M01067", "zymst", ("HSD17B7",), "core"),
    (26, "Zymosterol", "Cholest-7-enol", "zymst", "chlstol", ("EBP",), "bloch"),
    (27, "Cholest-7-enol", "7-dehydrodesmosterol", "chlstol", "ddsmsterol", ("SC5D",), "bloch"),
    (28, "7-dehydrodesmosterol", "Desmosterol", "ddsmsterol", "dsmsterol", ("DHCR7",), "bloch"),
    (29, "Desmosterol", "Cholesterol", "dsmsterol", "chsterol", ("DHCR24",), "bloch"),
    (26.1, "Zymosterol", "Zymosterol intermediate", "zymst", "zymstnl", ("DHCR24",), "kandutsch_russell"),
    (27.1, "Zymosterol intermediate", "Lathosterol", "zymstnl", "lthstrl", ("EBP",), "kandutsch_russell"),
    (28.1, "Lathosterol", "7-dehydrocholesterol", "lthstrl", "7dhchsterol", ("SC5D",), "kandutsch_russell"),
    (29.1, "7-dehydrocholesterol", "Cholesterol", "7dhchsterol", "chsterol", ("DHCR7",), "kandutsch_russell"),
    (27.2, "Cholest-7-enol", "Lathosterol", "chlstol", "lthstrl", ("DHCR24",), "kandutsch_russell"),
    (29.2, "Lathosterol", "Cholesterol", "lthstrl", "chsterol", ("DHCR24",), "kandutsch_direct"),
)

WAYPOINTS = (
    ("HMGCR", "HMG-CoA"),
    ("HMG-CoA", "Mevalonate"),
    ("Mevalonate", "Phosphomevalonate"),
    ("Phosphomevalonate", "Diphosphomevalonate"),
    ("Diphosphomevalonate", "IPP/DMAPP"),
    ("IPP/DMAPP", "GPP"),
    ("GPP", "FPP"),
    ("FPP", "Presqualene diphosphate"),
    ("Presqualene diphosphate", "Squalene"),
    ("Squalene", "2,3-oxidosqualene"),
    ("2,3-oxidosqualene", "Lanosterol"),
    ("Lanosterol", "Cholesterol"),
)


def _base_id(model_id: str) -> str:
    return str(model_id).rsplit("_", 1)[0]


def _compartment(model_id: str) -> str:
    return str(model_id).rsplit("_", 1)[-1]


def _normalise_hmdb(value: Any) -> str | None:
    match = re.search(r"HMDB\s*0*(\d+)", str(value), flags=re.IGNORECASE)
    return f"HMDB{int(match.group(1)):07d}" if match else None


def _gene_rule_ids(rule: str) -> set[str]:
    """Extract Recon3D gene identifiers without treating the ``AT1`` suffix as a gene."""
    ids = set(re.findall(r"(\d+)_AT\d+", str(rule)))
    if not ids and str(rule).strip():
        ids = set(re.findall(r"\b\d+\b", str(rule)))
    return ids


def reaction_symbols(reaction: Mapping[str, Any], hgnc_maps: Mapping[str, Any]) -> tuple[str, ...]:
    """Return approved HGNC symbols represented in a Recon3D gene rule."""
    symbols: set[str] = set()
    for identifier in _gene_rule_ids(str(reaction.get("gene_reaction_rule", ""))):
        symbols.update(hgnc_maps.get("entrez", {}).get(identifier, set()))
    return tuple(sorted(symbols))


def _reaction_node_candidates(number: int, symbols: Iterable[str]) -> tuple[str, ...]:
    symbols = tuple(sorted(set(symbols)))
    candidates = [f"Gene{number}__{symbol}" for symbol in symbols]
    candidates.extend(f"Gene{number}__{symbol}_reverse" for symbol in symbols)
    if symbols:
        candidates.append(f"Gene{number}__{'_'.join(symbols)}")
        candidates.append(f"Gene{number}__{'_'.join(symbols)}_reverse")
    return tuple(dict.fromkeys(candidates))


def _metabolite_candidates(metabolite: Mapping[str, Any]) -> tuple[str, ...]:
    model_id = str(metabolite["id"])
    compartment = str(metabolite.get("compartment") or _compartment(model_id))
    candidates = []
    for value in metabolite.get("annotation", {}).get("hmdb", []):
        if normalized := _normalise_hmdb(value):
            candidates.append(f"Metab__{normalized}_{compartment}")
    candidates.append(f"Metab__{model_id}")
    return tuple(dict.fromkeys(candidates))


def _model_node_map(model: Mapping[str, Any]) -> dict[str, tuple[str, ...]]:
    return {str(item["id"]): _metabolite_candidates(item) for item in model.get("metabolites", [])}


def _resolve_model_node(model_id: str, model_nodes: Mapping[str, tuple[str, ...]], pkn_nodes: set[str]) -> str | None:
    for candidate in model_nodes.get(model_id, (f"Metab__{model_id}",)):
        if candidate in pkn_nodes:
            return candidate
    return None


def _present(value: Any) -> bool:
    """True for a real node name; False for None, NaN and empty values."""
    if value is None:
        return False
    if isinstance(value, float) and np.isnan(value):
        return False
    return str(value) not in {"", "nan", "None"}


def infer_names_from_cosmos_reactions(
    model: Mapping[str, Any],
    edges: pd.DataFrame,
    reaction_numbers: Iterable[int],
    resolved: Mapping[str, str],
) -> tuple[dict[str, str], dict[str, str]]:
    """Name Recon3D metabolites from the COSMOS neighbours of their reactions.

    COSMOS ``Gene<N>__...`` nodes are Recon3D reaction N, with in-edges from
    the reactants and out-edges to the products of Recon3D's written
    direction.  Some intermediates (e.g. 7-dehydrocholesterol, lathosterol,
    desmosterol) lack HMDB annotations in Recon3D but are HMDB-named in
    COSMOS.  On each side of a reaction, a name is assigned only when exactly
    one non-cofactor Recon3D metabolite is unresolved and exactly one COSMOS
    metabolite neighbour is unclaimed.  The rule is repeated until no new
    name is assigned; ambiguous cases stay unresolved.
    """
    predecessors: dict[str, list[str]] = defaultdict(list)
    successors: dict[str, list[str]] = defaultdict(list)
    forward_nodes: dict[int, list[str]] = defaultdict(list)
    for source, target in zip(edges["source"].astype(str), edges["target"].astype(str)):
        predecessors[target].append(source)
        successors[source].append(target)
    for node in set(predecessors) | set(successors):
        match = re.match(r"Gene(\d+)__", node)
        if match and not node.endswith("_reverse"):
            forward_nodes[int(match.group(1))].append(node)

    reactions = list(model.get("reactions", []))
    assigned = dict(resolved)
    sources: dict[str, str] = {}
    changed = True
    while changed:
        changed = False
        for number in sorted(set(int(value) for value in reaction_numbers)):
            stoichiometry = reactions[number - 1].get("metabolites", {})
            for reaction_node in forward_nodes.get(number, []):
                for side, neighbours in ((-1, predecessors[reaction_node]), (1, successors[reaction_node])):
                    model_ids = [str(m) for m, c in stoichiometry.items() if np.sign(float(c)) == side]
                    claimed = {assigned[m] for m in model_ids if m in assigned}
                    unresolved = [m for m in model_ids if m not in assigned and _base_id(m) not in COFACTOR_BIGG_IDS]
                    unclaimed = [n for n in neighbours if n.startswith("Metab__") and n not in claimed]
                    if len(unresolved) == 1 and len(unclaimed) == 1:
                        assigned[unresolved[0]] = unclaimed[0]
                        sources[unresolved[0]] = f"cosmos_neighbour_of_Gene{number}"
                        changed = True

    by_node: dict[str, set[str]] = defaultdict(set)
    for model_id, node in assigned.items():
        if model_id in sources:
            by_node[node].add(_base_id(model_id) + "_" + _compartment(model_id))
    conflicts = {node: ids for node, ids in by_node.items() if len(ids) > 1}
    if conflicts:
        raise ValueError(f"COSMOS neighbour inference assigned one node to several metabolites: {conflicts}")
    inferred = {model_id: node for model_id, node in assigned.items() if model_id in sources}
    return inferred, sources


def canonical_step_table(
    model: Mapping[str, Any],
    hgnc_maps: Mapping[str, Any],
    *,
    pkn_nodes: set[str] | None = None,
    canonical_genes: Iterable[str] = CANONICAL_GENES,
    edges: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build one row per Recon3D reaction matching a registered route step.

    Substrate and product names come from Recon3D annotations first. When
    ``edges`` (a COSMOS PKN) is given, names still missing are inferred from
    the COSMOS neighbours of the canonical reactions
    (:func:`infer_names_from_cosmos_reactions`); the source of each name is
    recorded in ``substrate_source`` / ``product_source``.
    """
    canonical = set(canonical_genes)
    model_nodes = _model_node_map(model)
    rows: list[dict[str, Any]] = []
    for number, reaction in enumerate(model.get("reactions", []), start=1):
        symbols = reaction_symbols(reaction, hgnc_maps)
        matched_genes = tuple(sorted(set(symbols) & canonical))
        if not matched_genes:
            continue
        metabolite_ids = {str(item) for item in reaction.get("metabolites", {})}
        bases = {_base_id(item) for item in metabolite_ids}
        for order, previous, following, substrate_base, product_base, genes, branch in _STEP_SPECS:
            if not (set(genes) & set(matched_genes)):
                continue
            if substrate_base not in bases or product_base not in bases:
                continue
            substrate_ids = sorted(item for item in metabolite_ids if _base_id(item) == substrate_base)
            product_ids = sorted(item for item in metabolite_ids if _base_id(item) == product_base)
            stoichiometry = reaction.get("metabolites", {})
            for substrate_id in substrate_ids:
                for product_id in product_ids:
                    # A conversion step needs substrate and product on opposite
                    # sides (excludes e.g. HMR_1532, GGPP + PPi <=> DMAPP + FPP).
                    if np.sign(float(stoichiometry[substrate_id])) == np.sign(float(stoichiometry[product_id])):
                        continue
                    expected_nodes = _reaction_node_candidates(number, matched_genes)
                    rows.append({
                        "step_order": order,
                        "previous_intermediate": previous,
                        "next_intermediate": following,
                        "recon_number": number,
                        "recon_id": str(reaction["id"]),
                        "reaction_name": str(reaction.get("name", "")),
                        "subsystem": str(reaction.get("subsystem", "")),
                        "substrate_model_id": substrate_id,
                        "product_model_id": product_id,
                        "substrate": _resolve_model_node(substrate_id, model_nodes, pkn_nodes or set()),
                        "product": _resolve_model_node(product_id, model_nodes, pkn_nodes or set()),
                        "compartment": _compartment(substrate_id),
                        "genes": ";".join(matched_genes),
                        "branch": branch,
                        "reaction_node": expected_nodes[0] if expected_nodes else f"Gene{number}",
                        "reaction_node_candidates": ";".join(expected_nodes),
                    })
    table = pd.DataFrame(rows)
    if table.empty:
        raise ValueError("No canonical Recon3D route reactions matched the registered step table")
    table["substrate_source"] = np.where(table["substrate"].map(_present), "recon3d_annotation", "unresolved")
    table["product_source"] = np.where(table["product"].map(_present), "recon3d_annotation", "unresolved")
    if edges is not None:
        reactions = list(model.get("reactions", []))
        numbers = sorted(set(table["recon_number"].astype(int)))
        resolved: dict[str, str] = {}
        for number in numbers:
            for model_id in reactions[number - 1].get("metabolites", {}):
                node = _resolve_model_node(str(model_id), model_nodes, pkn_nodes or set())
                if node is not None:
                    resolved[str(model_id)] = node
        inferred, sources = infer_names_from_cosmos_reactions(model, edges, numbers, resolved)
        for column, id_column, source_column in (
            ("substrate", "substrate_model_id", "substrate_source"),
            ("product", "product_model_id", "product_source"),
        ):
            missing = ~table[column].map(_present)
            table.loc[missing, column] = table.loc[missing, id_column].map(inferred)
            filled = missing & table[column].map(_present)
            table.loc[filled, source_column] = table.loc[filled, id_column].map(sources)
    return table.sort_values(["step_order", "compartment", "recon_number", "recon_id"]).reset_index(drop=True)


def read_pkn(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype={"source": str, "target": str})
    required = {"source", "sign", "target"}
    if not required.issubset(table.columns):
        raise ValueError(f"PKN lacks columns: {sorted(required - set(table.columns))}")
    return table[["source", "sign", "target"]].astype({"source": str, "target": str, "sign": int})


def _nodes(edges: pd.DataFrame) -> set[str]:
    return set(edges["source"]) | set(edges["target"])


def _edge_pairs(edges: pd.DataFrame) -> set[tuple[str, str]]:
    return set(edges[["source", "target"]].itertuples(index=False, name=None))


def gap_audit(
    steps: pd.DataFrame,
    primary_edges: pd.DataFrame,
    unpruned_edges: pd.DataFrame,
) -> pd.DataFrame:
    """Classify each canonical reaction against primary and unpruned COSMOS."""
    primary_nodes = _nodes(primary_edges)
    unpruned_nodes = _nodes(unpruned_edges)
    primary_pairs = _edge_pairs(primary_edges)
    unpruned_pairs = _edge_pairs(unpruned_edges)
    rows = []
    for row in steps.itertuples(index=False):
        candidates = str(row.reaction_node_candidates).split(";") if row.reaction_node_candidates else []
        forward_node = next((
            node for node in candidates
            if node in primary_nodes and _present(row.substrate) and _present(row.product)
            and (row.substrate, node) in primary_pairs and (node, row.product) in primary_pairs
        ), None)
        reverse_node = next((
            node for node in candidates
            if node in primary_nodes and _present(row.substrate) and _present(row.product)
            and (row.product, node) in primary_pairs and (node, row.substrate) in primary_pairs
        ), None)
        primary_reaction = forward_node or reverse_node or next((node for node in candidates if node in primary_nodes), None)
        unpruned_reaction = next((node for node in candidates if node in unpruned_nodes), None)
        forward = forward_node is not None
        reverse = reverse_node is not None
        if forward:
            status = "present_forward"
        elif reverse:
            status = "present_reverse_only"
        elif unpruned_reaction is not None and primary_reaction is None:
            status = "pruned_unexpressed"
        elif unpruned_reaction is None:
            status = "absent_from_cosmos"
        else:
            status = "present_but_unconnected"
        rows.append({
            **row._asdict(),
            "primary_reaction_node": primary_reaction,
            "unpruned_reaction_node": unpruned_reaction,
            "primary_forward": forward,
            "primary_reverse": reverse,
            "status": status,
            "expected_repair": status == "absent_from_cosmos",
            "primary_node_present": primary_reaction is not None,
            "unpruned_node_present": unpruned_reaction is not None,
        })
    return pd.DataFrame(rows)


def genome_absent_by_subsystem(model: Mapping[str, Any], unpruned_nodes: set[str]) -> pd.DataFrame:
    """Count gene-associated Recon3D reactions absent from unpruned COSMOS."""
    rows = []
    for number, reaction in enumerate(model.get("reactions", []), start=1):
        if not str(reaction.get("gene_reaction_rule", "")).strip():
            continue
        prefix = f"Gene{number}__"
        rows.append({
            "recon_number": number,
            "subsystem": str(reaction.get("subsystem") or "Unspecified"),
            "present_in_unpruned": any(node.startswith(prefix) for node in unpruned_nodes),
        })
    table = pd.DataFrame(rows)
    summary = table.groupby("subsystem", as_index=False).agg(
        n_gene_associated_reactions=("recon_number", "size"),
        n_absent_from_cosmos=("present_in_unpruned", lambda values: int((~values).sum())),
    )
    total = pd.DataFrame([{
        "subsystem": "__TOTAL__",
        "n_gene_associated_reactions": len(table),
        "n_absent_from_cosmos": int((~table["present_in_unpruned"]).sum()),
    }])
    return pd.concat([summary, total], ignore_index=True)


def cofactor_nodes(
    model: Mapping[str, Any],
    pkn_nodes: set[str],
    graph_degree: Mapping[str, int],
    protected_model_ids: Iterable[str],
) -> pd.DataFrame:
    """Map registered BiGG cofactors to COSMOS nodes and annotate degree."""
    protected = {_base_id(item) for item in protected_model_ids}
    rows: list[dict[str, Any]] = []
    for metabolite in model.get("metabolites", []):
        model_id = str(metabolite["id"])
        bigg_ids = {str(value) for value in metabolite.get("annotation", {}).get("bigg.metabolite", [])}
        bigg_ids.add(_base_id(model_id))
        matched = sorted(bigg_ids & COFACTOR_BIGG_IDS)
        if not matched:
            continue
        candidates = _metabolite_candidates(metabolite)
        mapped = [node for node in candidates if node in pkn_nodes]
        for bigg_id in matched:
            for node in mapped or [None]:
                rows.append({
                    "model_id": model_id,
                    "bigg_id": bigg_id,
                    "compartment": str(metabolite.get("compartment") or _compartment(model_id)),
                    "pkn_node": node,
                    "degree": int(graph_degree.get(node, 0)) if node else np.nan,
                    "protected": _base_id(model_id) in protected,
                    "mapped": node is not None,
                })
    return pd.DataFrame(rows, columns=["model_id", "bigg_id", "compartment", "pkn_node", "degree", "protected", "mapped"])


def _reaction_node_present(row: Mapping[str, Any], nodes: set[str]) -> str | None:
    candidates = str(row.get("reaction_node_candidates", "")).split(";")
    return next((node for node in candidates if node in nodes), None)


def repaired_pkn(
    primary_edges: pd.DataFrame,
    primary_nodes: pd.DataFrame,
    steps: pd.DataFrame,
    gaps: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Add only canonical steps absent from both COSMOS variants."""
    edges = [tuple(row) for row in primary_edges[["source", "sign", "target"]].itertuples(index=False, name=None)]
    nodes = _nodes(primary_edges)
    additions = []
    for row in gaps[gaps["status"].eq("absent_from_cosmos")].itertuples(index=False):
        reaction_node = f"Gene{int(row.recon_number)}__{row.genes.replace(';', '_')}"
        if not (_present(row.substrate) and _present(row.product)):
            raise ValueError(f"Cannot repair {row.recon_id}: substrate/product is not mapped to COSMOS")
        additions.extend([
            (str(row.substrate), 1, reaction_node),
            (reaction_node, 1, str(row.product)),
        ])
        nodes.update((str(row.substrate), reaction_node, str(row.product)))
    edge_set = {(source, target) for source, _, target in edges}
    for source, sign, target in additions:
        if (source, target) not in edge_set:
            edges.append((source, sign, target))
            edge_set.add((source, target))
    repaired_edges = pd.DataFrame(edges, columns=["source", "sign", "target"]).sort_values(["source", "target", "sign"]).reset_index(drop=True)

    # Stage-03 node tables may be read with Arrow string columns; use object
    # columns because the repaired rows add boolean and numeric fields.
    base_nodes = primary_nodes.set_index("node").copy().astype(object)
    for node in sorted(nodes - set(base_nodes.index)):
        base_nodes.loc[node, "type"] = "enzyme_or_reaction" if node.startswith("Gene") else "metabolite_model_id" if node.startswith("Metab__") else "other"
        base_nodes.loc[node, "approved_symbol"] = node.split("__", 1)[1] if "__" in node else None
        base_nodes.loc[node, "expressed"] = True
    indegree = repaired_edges.groupby("target").size().to_dict()
    outdegree = repaired_edges.groupby("source").size().to_dict()
    base_nodes["in_degree"] = [int(indegree.get(node, 0)) for node in base_nodes.index]
    base_nodes["out_degree"] = [int(outdegree.get(node, 0)) for node in base_nodes.index]
    repaired_nodes = base_nodes.reset_index().rename(columns={"index": "node"})
    filter_log = pd.DataFrame([{
        "variant": "canonical_repair",
        "step": "canonical_repair",
        "description": "add canonical reactions classified absent from primary and unpruned COSMOS",
        "n_rows": len(repaired_edges),
        "n_edges": len(repaired_edges),
        "n_nodes": len(nodes),
        "n_added_reaction_nodes": len(additions) // 2,
        "n_added_edges": len(additions),
    }])
    return repaired_edges, repaired_nodes, filter_log


def _is_metabolite(node: str) -> bool:
    return str(node).startswith("Metab__")


def _is_reaction(node: str) -> bool:
    return str(node).startswith("Gene")


def _reaction_has_gene(node: str, gene: str) -> bool:
    if "__" not in node:
        return False
    return gene in node.split("__", 1)[1].split("_")


def metabolism_graph(edges: pd.DataFrame) -> nx.DiGraph:
    """Keep reaction/metabolite edges plus a gene's own reaction edges."""
    graph = nx.DiGraph()
    for source, _, target in edges[["source", "sign", "target"]].itertuples(index=False, name=None):
        source, target = str(source), str(target)
        if (_is_metabolite(source) and _is_reaction(target)) or (_is_reaction(source) and _is_metabolite(target)):
            graph.add_edge(source, target)
        elif _is_reaction(target) and not _is_reaction(source) and not _is_metabolite(source) and _reaction_has_gene(target, source):
            graph.add_edge(source, target)
        elif _is_reaction(source) and not _is_reaction(target) and not _is_metabolite(target) and _reaction_has_gene(source, target):
            graph.add_edge(source, target)
    return graph


def route_restricted_graph(
    edges: pd.DataFrame,
    steps: pd.DataFrame,
) -> nx.DiGraph:
    """Restrict COSMOS to canonical reactions/intermediates and valid transports."""
    graph = nx.DiGraph()
    pairs = _edge_pairs(edges)
    reaction_nodes = set()
    metabolite_nodes = set()
    for row in steps.itertuples(index=False):
        if _present(row.substrate):
            metabolite_nodes.add(str(row.substrate))
        if _present(row.product):
            metabolite_nodes.add(str(row.product))
        reaction_nodes.update(str(row.reaction_node_candidates).split(";"))
    allowed = {node for node in reaction_nodes if node in _nodes(edges)} | metabolite_nodes | {"HMGCR"}
    incoming_transport = {target for source, target in pairs if source in metabolite_nodes}
    outgoing_transport = {source for source, target in pairs if target in metabolite_nodes}
    transports = {
        node for node in incoming_transport & outgoing_transport
        if _is_reaction(node) and node not in reaction_nodes
    }
    allowed.update(transports)
    # Sorted insertion makes tie-breaking between equal-length paths independent
    # of Python's string-hash randomisation (otherwise reruns report different paths).
    for source, target in sorted(pairs):
        if source in allowed and target in allowed:
            graph.add_edge(source, target)
    return graph


def _label_nodes(steps: pd.DataFrame) -> dict[str, set[str]]:
    labels: dict[str, set[str]] = defaultdict(set)
    for row in steps.itertuples(index=False):
        previous = getattr(row, "previous_intermediate", None)
        following = getattr(row, "next_intermediate", None)
        if _present(row.substrate) and _present(previous):
            labels[str(previous)].add(str(row.substrate))
        if _present(row.product) and _present(following):
            labels[str(following)].add(str(row.product))
    # The waypoint table merges IPP and DMAPP into one checkpoint.
    labels["IPP/DMAPP"] = labels.get("IPP", set()) | labels.get("DMAPP", set())
    labels["Cholesterol"] = {node for nodes in labels.values() for node in nodes if node.startswith(CHOLESTEROL_NODE_PREFIX)}
    return labels


def waypoint_table(variants: Mapping[str, pd.DataFrame], steps: pd.DataFrame) -> pd.DataFrame:
    """Test ordered canonical waypoints on primary, unpruned and repaired graphs."""
    labels = _label_nodes(steps)
    rows = []
    for variant, edges in variants.items():
        graph = route_restricted_graph(edges, steps)
        first_failure = None
        checks = []
        for previous, following in WAYPOINTS:
            sources = {"HMGCR"} if previous == "HMGCR" else labels.get(previous, set())
            targets = labels.get(following, set())
            reachable = any(
                source in graph and target in graph and nx.has_path(graph, source, target)
                for source in sources for target in targets
            ) if sources and targets else False
            if not reachable and first_failure is None:
                first_failure = f"{previous}->{following}"
            checks.append((previous, following, reachable))
        for previous, following, reachable in checks:
            rows.append({
                "variant": variant,
                "previous_intermediate": previous,
                "next_intermediate": following,
                "reachable": reachable,
                "first_failing_waypoint": first_failure or "",
            })
    return pd.DataFrame(rows)


def _shortest_path(graph: nx.DiGraph, source: str, target: str, weight: str | None = None) -> tuple[list[str] | None, float]:
    try:
        path = nx.shortest_path(graph, source, target, weight=weight)
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return None, np.nan
    return path, float(nx.path_weight(graph, path, weight)) if weight else float(len(path) - 1)


def _path_branch(path: list[str], steps: pd.DataFrame) -> str:
    matched = []
    triples = {}
    for row in steps.itertuples(index=False):
        if not (_present(row.substrate) and _present(row.product)):
            continue
        for reaction_node in str(row.reaction_node_candidates).split(";"):
            triples[(str(row.substrate), reaction_node, str(row.product))] = str(row.branch)
    for index in range(len(path) - 2):
        branch = triples.get((path[index], path[index + 1], path[index + 2]))
        if branch:
            matched.append(branch)
    return matched[-1] if matched else "none"


def _path_metrics(
    path: list[str] | None,
    cost: float,
    steps: pd.DataFrame,
    graph: nx.DiGraph,
) -> dict[str, Any]:
    canonical_reactions = set(steps["reaction_node_candidates"].str.split(";").explode())
    canonical_metabolites = {str(node) for column in ("substrate", "product") for node in steps[column].dropna()}
    route_nodes = canonical_reactions | canonical_metabolites | {"HMGCR"}
    if path is None:
        return {
            "path_found": False, "steps": np.nan, "cost": np.nan,
            "n_canonical_enzymes": 0, "n_waypoints_in_order": 0,
            "fraction_waypoints_in_order": 0.0,
            "first_noncanonical_node": "", "branch": "none", "path": "",
        }
    # Ordered intermediates from mevalonate to cholesterol (HMG-CoA is skipped:
    # a path starting at the HMGCR gene enters the reaction node directly).
    labels = _label_nodes(steps)
    ordered = [following for _, following in WAYPOINTS if following != "HMG-CoA"]
    position = {node: index for index, node in reversed(list(enumerate(path)))}
    last, in_order = -1, 0
    for label in ordered:
        hits = sorted(
            index for node in labels.get(label, set())
            for index in [position.get(node)] if index is not None and index > last
        )
        if hits:
            last, in_order = hits[0], in_order + 1
    noncanonical = next((node for node in path[1:] if node not in route_nodes), "")
    return {
        "path_found": True,
        "steps": len(path) - 1,
        "cost": cost,
        "n_canonical_enzymes": sum(node in canonical_reactions for node in path),
        "n_waypoints_in_order": in_order,
        "fraction_waypoints_in_order": in_order / len(ordered),
        "first_noncanonical_node": noncanonical,
        "branch": _path_branch(path, steps),
        "path": " -> ".join(path),
    }


def path_variant_tables(
    variants: Mapping[str, pd.DataFrame],
    steps: pd.DataFrame,
    cofactors: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run full, metabolism-only and weighted path variants."""
    rows = []
    edge_rows = []
    protected_nodes = set(cofactors.loc[cofactors["protected"] & cofactors["mapped"], "pkn_node"].astype(str)) if not cofactors.empty else set()
    intermediate_nodes = {
        str(node) for column in ("substrate", "product") if column in steps.columns
        for node in steps[column] if _present(node)
    }
    for variant, edges in variants.items():
        if variant == "unpruned":
            continue
        full = nx.DiGraph()
        full.add_edges_from((str(source), str(target)) for source, _, target in edges.itertuples(index=False, name=None))
        metabolism = metabolism_graph(edges)
        no_cofactor = metabolism.copy()
        remove = set(cofactors.loc[cofactors["mapped"] & ~cofactors["protected"], "pkn_node"].astype(str)) if not cofactors.empty else set()
        no_cofactor.remove_nodes_from(remove)
        degrees = dict(no_cofactor.degree())
        weighted = no_cofactor.copy()
        for source, target in weighted.edges:
            weighted[source][target]["weight"] = float(degrees.get(target, 1) if _is_metabolite(target) else 0.0)
        # Top 1% of metabolite *nodes* by degree (not of distinct degree
        # values); canonical intermediates and HMGCR are always protected.
        metabolite_degrees = sorted((degree for node, degree in degrees.items() if _is_metabolite(node)), reverse=True)
        cutoff_index = max(0, int(np.ceil(len(metabolite_degrees) * 0.01)) - 1)
        cutoff = metabolite_degrees[cutoff_index] if metabolite_degrees else np.inf
        route_protected = protected_nodes | intermediate_nodes | {"HMGCR"}
        top_degree = {
            node for node, degree in degrees.items()
            if _is_metabolite(node) and degree >= cutoff and node not in route_protected
        }
        weighted_top = weighted.copy()
        weighted_top.remove_nodes_from(top_degree)
        restricted = route_restricted_graph(edges, steps)
        graphs = {
            "restricted_route_unweighted": (restricted, None),
            "full_unweighted": (full, None),
            "metabolism_only": (metabolism, None),
            "metabolism_no_cofactors": (no_cofactor, None),
            "metabolism_no_cofactors_weighted": (weighted, "weight"),
            "metabolism_no_cofactors_weighted_top1pct": (weighted_top, "weight"),
        }
        cholesterol_targets = sorted(node for node in _nodes(edges) if node.startswith(CHOLESTEROL_NODE_PREFIX))
        for path_variant, (graph, weight) in graphs.items():
            for target in cholesterol_targets:
                path, cost = _shortest_path(graph, "HMGCR", target, weight)
                metrics = _path_metrics(path, cost, steps, graph)
                rows.append({"variant": variant, "path_variant": path_variant, "source": "HMGCR", "target": target, **metrics})
                if path:
                    for index, (source, child) in enumerate(zip(path, path[1:]), start=1):
                        edge_rows.append({
                            "variant": variant,
                            "path_variant": path_variant,
                            "source": "HMGCR",
                            "target": target,
                            "edge_index": index,
                            "edge_source": source,
                            "edge_target": child,
                        })
    return pd.DataFrame(rows), pd.DataFrame(edge_rows)


def typed_metabolite_tf_null(
    edges: pd.DataFrame,
    *,
    variant: str,
    cholesterol_nodes: Iterable[str],
    target_tfs: Iterable[str],
    n_replicates: int = 100,
    seed_base: int = 20261010,
    max_length: int = 8,
) -> pd.DataFrame:
    """Compare cholesterol-to-TF positive paths with type-preserving rewires."""
    sources = sorted(set(cholesterol_nodes))
    targets = sorted(set(target_tfs))
    observed = positive_shortest_lengths(edges, sources, targets, max_length=max_length)
    finite = [length for length in observed.values() if length is not None]
    search_cap = max(finite, default=max_length)
    counts = {pair: 0 for pair, length in observed.items() if length is not None}
    valid = {pair: 0 for pair, length in observed.items() if length is not None}
    for replicate in range(n_replicates):
        rewired = rewire_type_preserving_signed_graph(edges, seed=seed_base + replicate)
        lengths = positive_shortest_lengths(rewired, sources, targets, max_length=search_cap)
        for pair, observed_length in observed.items():
            if observed_length is None:
                continue
            valid[pair] += 1
            if lengths[pair] is not None and lengths[pair] <= observed_length:
                counts[pair] += 1
    rows = []
    for (source, target), observed_length in observed.items():
        rows.append({
            "variant": variant,
            "source": source,
            "target": target,
            "observed_positive_shortest_length": observed_length,
            "n_replicates": n_replicates,
            "n_rewired_no_longer_than_observed": counts.get((source, target), 0),
            "fraction_no_longer_than_observed": counts.get((source, target), 0) / valid.get((source, target), 1) if observed_length is not None else np.nan,
            "seed_base": seed_base,
            "max_length": max_length,
        })
    return pd.DataFrame(rows)


def fetch_kegg_gene_lists(cache_dir: Path) -> tuple[dict[str, set[str]], dict[str, dict[str, Any]]]:
    """Fetch or verify KEGG hsa gene-link files and return Entrez IDs."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    gene_ids: dict[str, set[str]] = {}
    metadata: dict[str, dict[str, Any]] = {}
    for pathway in KEGG_PATHWAYS:
        path = cache_dir / f"kegg_{pathway}_gene_links.tsv"
        sidecar = path.with_suffix(".provenance.json")
        if path.exists():
            if not sidecar.exists():
                raise ValueError(f"Cached KEGG file has no provenance sidecar: {path}")
            sidecar_data = json.loads(sidecar.read_text())
        else:
            request = Request(KEGG_URL.format(pathway=pathway), headers={"User-Agent": "networkbio-aim5/1.0"})
            with urlopen(request, timeout=120) as response:
                payload = response.read()
            path.write_bytes(payload)
            sidecar_data = {"url": KEGG_URL.format(pathway=pathway), "sha256": sha256(path), "byte_count": len(payload)}
            sidecar.write_text(json.dumps(sidecar_data, indent=2) + "\n")
        if sidecar_data.get("sha256") != sha256(path):
            raise ValueError(f"KEGG cache checksum mismatch: {path}")
        ids = set()
        for line in path.read_text().splitlines():
            fields = line.split("\t")
            for field in fields:
                match = re.search(r"hsa:(\d+)", field)
                if match:
                    ids.add(match.group(1))
        gene_ids[pathway] = ids
        metadata[pathway] = sidecar_data
    return gene_ids, metadata


def verify_kegg_canonical_genes(
    gene_ids: Mapping[str, set[str]],
    hgnc_maps: Mapping[str, Any],
    *,
    canonical_genes: Iterable[str] = CANONICAL_GENES,
) -> pd.DataFrame:
    """Fail if a registered gene is absent from both KEGG pathway lists."""
    rows = []
    for gene in canonical_genes:
        pathway_hits = {
            pathway: any(gene in hgnc_maps.get("entrez", {}).get(identifier, set()) for identifier in ids)
            for pathway, ids in gene_ids.items()
        }
        if not any(pathway_hits.values()):
            raise ValueError(f"Canonical gene is absent from both KEGG pathway lists: {gene}")
        rows.append({"gene": gene, **pathway_hits})
    return pd.DataFrame(rows)
