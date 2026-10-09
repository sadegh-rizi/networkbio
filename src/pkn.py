"""Build, map and inspect the signed prior-knowledge networks for stage 03.

The functions in this module deliberately keep resource retrieval separate
from filtering.  A cached resource is never silently replaced, and every
filter returns a row in a log that can be written next to the resulting PKN.
The COSMOS resource is distributed with internal numeric identifiers; these
are converted to HGNC symbols and the canonical metabolite names used by the
rest of the project before expression pruning and path analysis.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
import hashlib
import json
from importlib.metadata import version
from pathlib import Path
import random
import re
from typing import Any, Iterable, Iterator, Mapping
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

Edge = tuple[str, int, str]

COSMOS_NETWORKCOMMONS_URL = "https://commons.omnipathdb.org/prior_knowledge/meta_network.sif"
COSMOS_OMNIPATH_URL = "https://metapkn.omnipathdb.org/metapkn__20200122.txt"
RECON3D_URL = "https://bigg.ucsd.edu/static/models/Recon3D.json"
WORKBENCH_URL = "https://www.metabolomicsworkbench.org/rest/compound/kegg_id/{kegg}/all"

SOURCE_GENES = ("HMGCR", "SREBF1", "SREBF2", "SCAP", "INSIG1")
TARGET_TFS = ("E2F1", "EGR1", "WT1", "SP1", "JUN", "JUNB")
MEVALONATE_HMDB = "HMDB0000227"
CHOLESTEROL_HMDB = "HMDB0000067"
COMPLEX_MEMBERS = {
    "AP1": ("ATF3", "FOS", "FOSB", "JUN", "JUNB"),
    "NFKB": ("NFKB1", "NFKB2", "RELA", "RELB", "REL"),
}


def sha256(path: Path) -> str:
    """Return the SHA256 digest of a local resource."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resource_sidecar(path: Path) -> Path:
    return path.with_suffix(".provenance.json")


def _fetch_bytes(path: Path, url: str, *, timeout: int = 300) -> tuple[Path, dict[str, Any]]:
    """Fetch one immutable resource or verify its cached copy."""
    sidecar = _resource_sidecar(path)
    if path.exists():
        if not sidecar.exists():
            raise ValueError(f"Cached resource has no provenance sidecar: {path}")
        metadata = json.loads(sidecar.read_text())
        if metadata.get("sha256") != sha256(path):
            raise ValueError(f"Cached resource checksum mismatch: {path}")
        return path, metadata

    path.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": "networkbio-stage03/1.0"})
    with urlopen(request, timeout=timeout) as response:
        payload = response.read()
    path.write_bytes(payload)
    metadata = {
        "url": url,
        "downloaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256": sha256(path),
        "byte_count": len(payload),
    }
    sidecar.write_text(json.dumps(metadata, indent=2))
    return path, metadata


def _as_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def _normalise_hmdb(value: Any) -> str | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    match = re.search(r"HMDB\s*0*(\d+)", str(value), flags=re.IGNORECASE)
    if not match:
        return None
    return f"HMDB{int(match.group(1)):07d}"


def _split_identifiers(value: Any) -> Iterator[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return
    for item in re.split(r"[|;,]", str(value)):
        item = item.strip()
        if item:
            yield item


def load_hgnc_maps(path: Path) -> dict[str, Any]:
    """Load approved-symbol, alias/previous-symbol and Entrez mappings."""
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"symbol", "status", "alias_symbol", "prev_symbol", "entrez_id"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"HGNC table lacks columns: {sorted(missing)}")
    approved_table = table[table["status"].eq("Approved")].copy()
    approved = set(approved_table["symbol"])

    def collect(column: str) -> dict[str, set[str]]:
        mapping: dict[str, set[str]] = defaultdict(set)
        for row in approved_table[["symbol", column]].itertuples(index=False):
            for identifier in _split_identifiers(row[1]):
                mapping[identifier].add(row[0])
        return dict(mapping)

    entrez: dict[str, set[str]] = defaultdict(set)
    for row in approved_table[["symbol", "entrez_id"]].itertuples(index=False):
        for identifier in _split_identifiers(row[1]):
            entrez[identifier].add(row[0])

    return {
        "approved": approved,
        "prev_symbol": collect("prev_symbol"),
        "alias_symbol": collect("alias_symbol"),
        "entrez": dict(entrez),
        "n_approved": len(approved),
    }


def harmonise_symbol(symbol: str, maps: Mapping[str, Any]) -> tuple[str, str]:
    """Map one symbol through approved, previous and alias HGNC names."""
    symbol = str(symbol).strip()
    if symbol in maps["approved"]:
        return symbol, "approved"
    previous = maps["prev_symbol"].get(symbol, set())
    if len(previous) == 1:
        return next(iter(previous)), "previous_symbol"
    aliases = maps["alias_symbol"].get(symbol, set())
    if len(aliases) == 1:
        return next(iter(aliases)), "alias_symbol"
    if len(previous) > 1 or len(aliases) > 1:
        return symbol, "ambiguous"
    return symbol, "unmapped"


def _node_count(edges: pd.DataFrame) -> int:
    if edges.empty:
        return 0
    return len(set(edges["source"]) | set(edges["target"]))


def _log_filter(
    rows: list[dict[str, Any]],
    variant: str,
    step: str,
    description: str,
    before: pd.DataFrame,
    after: pd.DataFrame,
    **extra: Any,
) -> None:
    rows.append(
        {
            "variant": variant,
            "step": step,
            "description": description,
            "n_rows": len(after),
            "n_edges": len(after.drop_duplicates(["source", "target"])),
            "n_nodes": _node_count(after),
            "n_removed_rows": len(before) - len(after),
            "n_removed_nodes": max(_node_count(before) - _node_count(after), 0),
            **extra,
        }
    )


def _deduplicate_edges(edges: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """Remove conflicting pairs and merge exact duplicate signed edges."""
    if edges.empty:
        return edges.copy(), {"conflicting_pairs": 0, "conflicting_rows": 0, "identical_rows_merged": 0}
    grouped = edges.groupby(["source", "target"], sort=False)["sign"]
    conflicting_pairs = int((grouped.nunique() > 1).sum())
    conflicting_keys = set(grouped.nunique()[grouped.nunique() > 1].index)
    key = list(zip(edges["source"], edges["target"]))
    conflicting_rows = sum(pair in conflicting_keys for pair in key)
    kept = edges.loc[[pair not in conflicting_keys for pair in key]].copy()
    identical_rows_merged = len(kept) - len(kept.drop_duplicates(["source", "target"], keep="first"))
    kept = kept.drop_duplicates(["source", "target"], keep="first").reset_index(drop=True)
    return kept, {
        "conflicting_pairs": conflicting_pairs,
        "conflicting_rows": conflicting_rows,
        "identical_rows_merged": int(identical_rows_merged),
    }


def _prepare_omnipath(raw: pd.DataFrame) -> pd.DataFrame:
    source = "source_genesymbol" if "source_genesymbol" in raw.columns else "source"
    target = "target_genesymbol" if "target_genesymbol" in raw.columns else "target"
    required = {source, target}
    missing = required - set(raw.columns)
    if missing:
        raise ValueError(f"OmniPath table lacks columns: {sorted(missing)}")
    out = raw.copy()
    out["source"] = out[source].astype(str).str.strip()
    out["target"] = out[target].astype(str).str.strip()
    return out


def _source_has_signor(value: Any) -> bool:
    return bool(re.search(r"(?:^|[;,])\s*SIGNOR\s*(?:$|[;,])", str(value), flags=re.IGNORECASE))


def _gene_nodes_table(edges: pd.DataFrame, expressed: set[str]) -> pd.DataFrame:
    nodes = sorted(set(edges["source"]) | set(edges["target"]))
    indegree = edges.groupby("target").size().to_dict()
    outdegree = edges.groupby("source").size().to_dict()
    return pd.DataFrame(
        {
            "node": nodes,
            "type": ["gene"] * len(nodes),
            "approved_symbol": nodes,
            "expressed": [node in expressed for node in nodes],
            "in_degree": [int(indegree.get(node, 0)) for node in nodes],
            "out_degree": [int(outdegree.get(node, 0)) for node in nodes],
        }
    )


def filter_omnipath(
    raw: pd.DataFrame,
    *,
    expressed: set[str],
    hgnc_maps: Mapping[str, Any],
    variant: str = "primary",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Apply P1-P6 to OmniPath and return PKN, filter log and node table."""
    edges = _prepare_omnipath(raw)
    log_rows: list[dict[str, Any]] = []
    _log_filter(log_rows, variant, "raw", "OmniPath response", edges, edges)

    if variant == "signor_only":
        before = edges
        edges = edges[edges["sources"].map(_source_has_signor)].copy()
        _log_filter(log_rows, variant, "P1_signor", "retain rows citing SIGNOR", before, edges)
    elif variant not in {"primary", "nc_rule", "unpruned"}:
        raise ValueError(f"Unknown OmniPath variant: {variant}")
    else:
        _log_filter(log_rows, variant, "P1_source", "all OmniPath resources", edges, edges)

    if not {"consensus_direction", "consensus_stimulation", "consensus_inhibition"}.issubset(edges.columns):
        raise ValueError("OmniPath response does not contain the consensus direction/sign columns")
    direction = edges["consensus_direction"].map(_as_bool)
    stimulation = edges["consensus_stimulation"].map(_as_bool)
    inhibition = edges["consensus_inhibition"].map(_as_bool)
    if variant == "nc_rule":
        keep = direction & (stimulation | inhibition)
        sign = np.where(stimulation, 1, -1)
    else:
        keep = direction & (stimulation ^ inhibition)
        sign = np.where(stimulation, 1, -1)
    before = edges
    edges = edges.loc[keep].copy()
    edges["sign"] = sign[keep.to_numpy()]
    _log_filter(
        log_rows,
        variant,
        "P2_direction_sign",
        "directed and signed consensus rows",
        before,
        edges,
        n_both_sign_rows=int((direction & stimulation & inhibition).sum()),
        n_no_consensus_direction=int((~direction).sum()),
    )

    before = edges
    curation = pd.to_numeric(edges["curation_effort"], errors="coerce") if "curation_effort" in edges else pd.Series(0, index=edges.index)
    edges = edges.loc[curation.ge(2)].copy()
    _log_filter(log_rows, variant, "P3_curation", "curation_effort >= 2", before, edges)

    before = edges
    complex_mask = edges["source"].map(_is_complex_node) | edges["target"].map(_is_complex_node)
    edges = edges.loc[~complex_mask].copy()
    _log_filter(log_rows, variant, "P4_complexes", "remove complex identifiers", before, edges, n_complex_rows=int(complex_mask.sum()))

    before = edges
    loop_mask = edges["source"].eq(edges["target"])
    edges = edges.loc[~loop_mask].copy()
    _log_filter(log_rows, variant, "P4_self_loops", "remove self-loops", before, edges, n_self_loop_rows=int(loop_mask.sum()))

    before = edges
    edges, duplicate_stats = _deduplicate_edges(edges[["source", "target", "sign"]])
    _log_filter(log_rows, variant, "P4_duplicates", "remove sign conflicts and merge identical pairs", before, edges, **duplicate_stats)

    before = edges
    mapped = edges["source"].map(lambda value: harmonise_symbol(value, hgnc_maps))
    mapped_target = edges["target"].map(lambda value: harmonise_symbol(value, hgnc_maps))
    edges["source"] = mapped.map(lambda item: item[0])
    edges["target"] = mapped_target.map(lambda item: item[0])
    mapping_reasons = list(mapped) + list(mapped_target)
    edges, duplicate_stats = _deduplicate_edges(edges)
    _log_filter(
        log_rows,
        variant,
        "P5_hgnc",
        "approved, previous and unique alias symbol harmonisation",
        before,
        edges,
        n_approved_symbol_mappings=sum(reason[1] == "approved" for reason in mapping_reasons),
        n_previous_symbol_mappings=sum(reason[1] == "previous_symbol" for reason in mapping_reasons),
        n_alias_symbol_mappings=sum(reason[1] == "alias_symbol" for reason in mapping_reasons),
        n_ambiguous_symbols=sum(reason[1] == "ambiguous" for reason in mapping_reasons),
        n_unmapped_symbols=sum(reason[1] == "unmapped" for reason in mapping_reasons),
        **duplicate_stats,
    )

    if variant == "unpruned":
        _log_filter(log_rows, variant, "P6_expression", "expression pruning skipped", edges, edges, pruned=False)
    else:
        before = edges
        keep = edges["source"].isin(expressed) & edges["target"].isin(expressed)
        edges = edges.loc[keep].copy()
        _log_filter(log_rows, variant, "P6_expression", "retain genes expressed in stage 01", before, edges, pruned=True)

    edges = edges[["source", "sign", "target"]].sort_values(["source", "target", "sign"]).reset_index(drop=True)
    nodes = _gene_nodes_table(edges.rename(columns={"sign": "sign"}), expressed)
    return edges, pd.DataFrame(log_rows), nodes


def _is_complex_node(value: Any) -> bool:
    text = str(value)
    return text.startswith("COMPLEX:") or "_" in text


def read_cosmos_edges(path: Path) -> pd.DataFrame:
    """Read either the NetworkCommons SIF or OmniPath COSMOS text format."""
    table = pd.read_csv(path, sep="\t", dtype=str)
    columns = {str(column).strip().lower(): column for column in table.columns}
    source_column = columns.get("source")
    target_column = columns.get("target")
    sign_column = columns.get("interaction") or columns.get("sign")
    if source_column is None or target_column is None or sign_column is None:
        table = pd.read_csv(path, sep="\t", dtype=str, header=None, names=["source", "sign", "target"])
        source_column, sign_column, target_column = "source", "sign", "target"
    out = pd.DataFrame(
        {
            "source": table[source_column].astype(str).str.strip(),
            "sign": pd.to_numeric(table[sign_column], errors="coerce"),
            "target": table[target_column].astype(str).str.strip(),
        }
    )
    out = out[out["sign"].isin([-1, 1])].copy()
    out["sign"] = out["sign"].astype(int)
    if out.isna().any().any():
        raise ValueError(f"COSMOS resource contains missing fields: {path}")
    return out.reset_index(drop=True)


def clean_cosmos_edges(raw: pd.DataFrame) -> pd.DataFrame:
    """Apply the P7 duplicate-sign mean rule without identifier mapping."""
    edges = raw[["source", "sign", "target"]].copy()
    edges = edges[~edges["source"].eq(edges["target"])]
    mean_sign = edges.groupby(["source", "target"], as_index=False)["sign"].mean()
    mean_sign = mean_sign[mean_sign["sign"].isin([-1, 1])].copy()
    mean_sign["sign"] = mean_sign["sign"].astype(int)
    return mean_sign[["source", "sign", "target"]].sort_values(["source", "target"]).reset_index(drop=True)


def parse_cosmos_node(node: str, hgnc_maps: Mapping[str, Any]) -> dict[str, Any]:
    """Parse and canonicalize one COSMOS node.

    The downloaded 2020 file uses internal prefixes.  Numeric XMetab IDs are
    HMDB numbers, X<Entrez> nodes are proteins, and XGene<Entrez>__<reaction>
    nodes are enzyme/reaction intermediates associated with that gene.
    """
    raw = str(node).strip()
    metabolite = re.fullmatch(r"XMetab__(.+)___([a-z])____", raw)
    canonical_metabolite = re.fullmatch(r"Metab__([^_]+)_([a-z])", raw)
    if metabolite:
        identifier, compartment = metabolite.groups()
        if identifier.isdigit():
            hmdb = f"HMDB{int(identifier):07d}"
            canonical = f"Metab__{hmdb}_{compartment}"
        else:
            hmdb = _normalise_hmdb(identifier)
            canonical = f"Metab__{hmdb or identifier}_{compartment}"
        return {
            "raw_node": raw,
            "node": canonical,
            "type": "metabolite",
            "approved_symbol": None,
            "expressed": True,
            "hmdb_id": hmdb,
            "compartment": compartment,
        }
    if canonical_metabolite:
        identifier, compartment = canonical_metabolite.groups()
        hmdb = _normalise_hmdb(identifier)
        if hmdb is None and identifier.isdigit():
            hmdb = f"HMDB{int(identifier):07d}"
        return {
            "raw_node": raw,
            "node": f"Metab__{hmdb or identifier}_{compartment}",
            "type": "metabolite",
            "approved_symbol": None,
            "expressed": True,
            "hmdb_id": hmdb,
            "compartment": compartment,
        }

    enzyme = re.fullmatch(r"(?:X)?Gene(\d+)__(.+)", raw)
    if enzyme:
        gene_id = enzyme.group(1)
        symbols = hgnc_maps.get("entrez", {}).get(gene_id, set())
        approved = next(iter(symbols)) if len(symbols) == 1 else None
        return {
            "raw_node": raw,
            "node": raw,
            "type": "enzyme_or_reaction",
            "approved_symbol": approved,
            "expressed": approved in hgnc_maps.get("expressed", set()),
            "hmdb_id": None,
            "compartment": None,
        }

    gene = re.fullmatch(r"X(\d+)", raw)
    if gene:
        symbols = hgnc_maps.get("entrez", {}).get(gene.group(1), set())
        approved = next(iter(symbols)) if len(symbols) == 1 else None
        return {
            "raw_node": raw,
            "node": approved or raw,
            "type": "gene",
            "approved_symbol": approved,
            "expressed": approved in hgnc_maps.get("expressed", set()),
            "hmdb_id": None,
            "compartment": None,
        }

    canonical_symbol, mapping_reason = harmonise_symbol(raw, hgnc_maps)
    if mapping_reason in {"approved", "previous_symbol", "alias_symbol"}:
        return {
            "raw_node": raw,
            "node": canonical_symbol,
            "type": "gene",
            "approved_symbol": canonical_symbol,
            "expressed": canonical_symbol in hgnc_maps.get("expressed", set()),
            "hmdb_id": None,
            "compartment": None,
        }
    return {
        "raw_node": raw,
        "node": raw,
        "type": "other",
        "approved_symbol": None,
        "expressed": None,
        "hmdb_id": None,
        "compartment": None,
    }


def cosmos_node_grammar(raw: pd.DataFrame) -> str:
    """Return a compact grammar report with counts and up to five examples."""
    nodes = sorted(set(raw["source"]) | set(raw["target"]))
    rows: dict[str, list[str]] = defaultdict(list)
    for node in nodes:
        if re.fullmatch(r"XMetab__\d+___[a-z]____", node):
            category = "metabolite: numeric HMDB"
        elif re.fullmatch(r"XMetab__[^_]+___[a-z]____", node):
            category = "metabolite: model identifier"
        elif re.fullmatch(r"(?:X)?Gene\d+__.+", node):
            category = "enzyme_or_reaction: XGene Entrez plus reaction"
        elif re.fullmatch(r"X\d+", node):
            category = "gene: X Entrez"
        elif re.fullmatch(r"Metab__.+_[a-z]", node):
            category = "metabolite: canonical HMDB or model identifier"
        elif re.fullmatch(r"[A-Za-z0-9.-]+", node):
            category = "gene: plain symbol"
        elif "_" in node:
            category = "complex: underscore-joined symbols"
        elif node.startswith("X"):
            category = "other: X-prefixed"
        else:
            category = "other"
        if len(rows[category]) < 5:
            rows[category].append(node)
    lines = [
        "# COSMOS meta-PKN node grammar",
        "",
        "The 2020 resource was inspected before filtering. Node identifiers are",
        "kept in `raw_node` provenance; PKN outputs use canonical names where a",
        "gene or numeric HMDB identifier can be resolved.",
        "",
        "| Pattern | Unique nodes | Examples |",
        "|---|---:|---|",
    ]
    counts = Counter()
    for node in nodes:
        if re.fullmatch(r"XMetab__\d+___[a-z]____", node):
            counts["metabolite: numeric HMDB"] += 1
        elif re.fullmatch(r"XMetab__[^_]+___[a-z]____", node):
            counts["metabolite: model identifier"] += 1
        elif re.fullmatch(r"(?:X)?Gene\d+__.+", node):
            counts["enzyme_or_reaction: XGene Entrez plus reaction"] += 1
        elif re.fullmatch(r"X\d+", node):
            counts["gene: X Entrez"] += 1
        elif re.fullmatch(r"Metab__.+_[a-z]", node):
            counts["metabolite: canonical HMDB or model identifier"] += 1
        elif re.fullmatch(r"[A-Za-z0-9.-]+", node):
            counts["gene: plain symbol"] += 1
        elif "_" in node:
            counts["complex: underscore-joined symbols"] += 1
        elif node.startswith("X"):
            counts["other: X-prefixed"] += 1
        else:
            counts["other"] += 1
    for category in sorted(counts):
        lines.append(f"| `{category}` | {counts[category]} | {', '.join(rows[category])} |")
    lines.extend(["", f"Unique nodes: {len(nodes)}", f"Rows inspected: {len(raw)}"])
    return "\n".join(lines) + "\n"


def _cosmos_nodes_table(edges: pd.DataFrame, hgnc_maps: Mapping[str, Any]) -> pd.DataFrame:
    nodes = sorted(set(edges["source"]) | set(edges["target"]))
    info = {node: parse_cosmos_node(node, hgnc_maps) for node in nodes}
    indegree = edges.groupby("target").size().to_dict()
    outdegree = edges.groupby("source").size().to_dict()
    return pd.DataFrame(
        {
            "node": nodes,
            "type": [info[node]["type"] for node in nodes],
            "approved_symbol": [info[node]["approved_symbol"] for node in nodes],
            "expressed": [info[node]["expressed"] for node in nodes],
            "in_degree": [int(indegree.get(node, 0)) for node in nodes],
            "out_degree": [int(outdegree.get(node, 0)) for node in nodes],
        }
    )


def filter_cosmos(
    raw: pd.DataFrame,
    *,
    hgnc_maps: Mapping[str, Any],
    expressed: set[str],
    variant: str = "primary",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Apply P7, P5 and P6 to the COSMOS meta-PKN."""
    if variant not in {"primary", "unpruned"}:
        raise ValueError(f"Unknown COSMOS variant: {variant}")
    local_maps = dict(hgnc_maps)
    local_maps["expressed"] = expressed
    edges = raw[["source", "sign", "target"]].copy()
    log_rows: list[dict[str, Any]] = []
    _log_filter(log_rows, variant, "raw", "COSMOS resource", edges, edges)

    before = edges
    loop_mask = edges["source"].eq(edges["target"])
    edges = edges.loc[~loop_mask].copy()
    _log_filter(log_rows, variant, "P7_self_loops", "remove self-loops", before, edges, n_self_loop_rows=int(loop_mask.sum()))

    before = edges
    grouped = edges.groupby(["source", "target"], as_index=False)["sign"].agg(["mean", "nunique"])
    grouped = grouped.reset_index()
    conflict_pairs = int(grouped["nunique"].gt(1).sum())
    edges = grouped[grouped["mean"].isin([-1, 1])][["source", "target", "mean"]].rename(columns={"mean": "sign"})
    edges["sign"] = edges["sign"].astype(int)
    _log_filter(log_rows, variant, "P7_duplicate_sign_average", "average duplicate signs and retain exactly +/-1", before, edges, n_conflicting_pairs=conflict_pairs)

    before = edges.copy()
    source_info = edges["source"].map(lambda node: parse_cosmos_node(node, local_maps))
    target_info = edges["target"].map(lambda node: parse_cosmos_node(node, local_maps))
    edges["source"] = source_info.map(lambda item: item["node"])
    edges["target"] = target_info.map(lambda item: item["node"])
    edges, duplicate_stats = _deduplicate_edges(edges)
    mapping_reasons = [item["type"] for item in list(source_info) + list(target_info)]
    _log_filter(
        log_rows,
        variant,
        "P5_hgnc",
        "canonicalize COSMOS genes and metabolites",
        before,
        edges,
        n_gene_nodes=sum(kind == "gene" for kind in mapping_reasons),
        n_metabolite_nodes=sum(kind == "metabolite" for kind in mapping_reasons),
        n_enzyme_reaction_nodes=sum(kind == "enzyme_or_reaction" for kind in mapping_reasons),
        **duplicate_stats,
    )

    if variant == "unpruned":
        _log_filter(log_rows, variant, "P6_expression", "expression pruning skipped", edges, edges, pruned=False)
    else:
        before = edges
        info = {node: parse_cosmos_node(node, local_maps) for node in set(edges["source"]) | set(edges["target"])}
        keep = edges["source"].map(lambda node: info[node]["type"] == "metabolite" or info[node]["approved_symbol"] in expressed)
        keep &= edges["target"].map(lambda node: info[node]["type"] == "metabolite" or info[node]["approved_symbol"] in expressed)
        edges = edges.loc[keep].copy()
        _log_filter(log_rows, variant, "P6_expression", "retain expressed genes and all metabolites", before, edges, pruned=True)

    edges = edges[["source", "sign", "target"]].sort_values(["source", "target", "sign"]).reset_index(drop=True)
    nodes = _cosmos_nodes_table(edges, local_maps)
    return edges, pd.DataFrame(log_rows), nodes


def fetch_omnipath(cache_dir: Path) -> tuple[pd.DataFrame, Path, dict[str, Any]]:
    """Fetch and cache the OmniPath interaction table."""
    path = cache_dir / "omnipath_interactions.tsv.gz"
    if path.exists():
        _, metadata = _fetch_bytes(path, metadata_url(path))
        return pd.read_csv(path, sep="\t"), path, metadata

    from omnipath.interactions import OmniPath

    raw = OmniPath.get(
        genesymbols=True,
        datasets="omnipath",
        organisms="human",
        directed=True,
        signed=True,
        license="academic",
    )
    if not isinstance(raw, pd.DataFrame) or raw.empty:
        raise ValueError("OmniPath returned no interactions")
    path.parent.mkdir(parents=True, exist_ok=True)
    raw.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    metadata = {
        "resource": "OmniPath omnipath dataset",
        "url": "omnipath.interactions.OmniPath.get",
        "downloaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256": sha256(path),
        "row_count": len(raw),
        "columns": list(raw.columns),
        "omnipath_version": version("omnipath"),
    }
    _resource_sidecar(path).write_text(json.dumps(metadata, indent=2))
    return raw, path, metadata


def metadata_url(path: Path) -> str:
    """Return the original URL recorded for a cached local resource.

    The URL is only used for the required argument of `_fetch_bytes` when a
    caller asks to verify a cache; an existing sidecar remains authoritative.
    """
    sidecar = _resource_sidecar(path)
    if sidecar.exists():
        return str(json.loads(sidecar.read_text()).get("url", "cached"))
    return "cached"


def fetch_cosmos(cache_dir: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Fetch COSMOS, preferring NetworkCommons and checking both sources."""
    resources = [
        ("networkcommons", "cosmos_meta_network.sif", COSMOS_NETWORKCOMMONS_URL),
        ("omnipath", "cosmos_metapkn__20200122.txt", COSMOS_OMNIPATH_URL),
    ]
    successful: dict[str, tuple[Path, dict[str, Any]]] = {}
    failures: dict[str, str] = {}
    for label, filename, url in resources:
        path = cache_dir / filename
        try:
            successful[label] = _fetch_bytes(path, url)
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            failures[label] = f"{type(exc).__name__}: {exc}"
    if not successful:
        raise RuntimeError(f"Could not fetch either COSMOS resource: {failures}")

    parsed = {label: read_cosmos_edges(item[0]) for label, item in successful.items()}
    cleaned = {label: set(map(tuple, clean_cosmos_edges(edges).itertuples(index=False, name=None))) for label, edges in parsed.items()}
    comparison = None
    if len(cleaned) == 2:
        comparison = cleaned["networkcommons"] == cleaned["omnipath"]
        # The current endpoints expose different namespaces and edge counts.
        # Keep this result in provenance rather than silently calling them
        # interchangeable; NetworkCommons remains the registered primary.
    selected = "networkcommons" if "networkcommons" in parsed else "omnipath"
    metadata = {
        "selected": selected,
        "sources_available": sorted(parsed),
        "source_failures": failures,
        "same_after_cleanup": comparison,
        "clean_edge_counts": {label: len(edges) for label, edges in cleaned.items()},
        "resource_metadata": {label: item[1] for label, item in successful.items()},
    }
    return parsed[selected], metadata


def fetch_recon3d(cache_dir: Path) -> tuple[dict[str, Any], Path, dict[str, Any]]:
    """Fetch and cache Recon3D JSON."""
    path, metadata = _fetch_bytes(cache_dir / "Recon3D.json", RECON3D_URL)
    model = json.loads(path.read_text())
    if not {"metabolites", "genes"}.issubset(model):
        raise ValueError("Recon3D resource lacks metabolites or genes")
    return model, path, metadata


def recon3d_hmdb_by_kegg(model: Mapping[str, Any]) -> dict[str, list[str]]:
    """Return all normalized HMDB IDs annotated for each KEGG compound."""
    out: dict[str, set[str]] = defaultdict(set)
    for metabolite in model.get("metabolites", []):
        annotations = metabolite.get("annotation", {})
        kegg_values = annotations.get("kegg.compound", [])
        hmdb_values = annotations.get("hmdb", [])
        if isinstance(kegg_values, str):
            kegg_values = [kegg_values]
        if isinstance(hmdb_values, str):
            hmdb_values = [hmdb_values]
        hmdb_ids = [normalized for value in hmdb_values if (normalized := _normalise_hmdb(value))]
        for kegg in kegg_values:
            out[str(kegg)].update(hmdb_ids)
    return {key: sorted(value) for key, value in out.items()}


def _extract_hmdb_values(payload: Any) -> set[str]:
    found: set[str] = set()
    if isinstance(payload, Mapping):
        for key, value in payload.items():
            if str(key).lower() in {"hmdb_id", "hmdb", "hmdbid"}:
                if isinstance(value, list):
                    values = value
                else:
                    values = [value]
                for item in values:
                    if normalized := _normalise_hmdb(item):
                        found.add(normalized)
            found.update(_extract_hmdb_values(value))
    elif isinstance(payload, list):
        for item in payload:
            found.update(_extract_hmdb_values(item))
    elif isinstance(payload, str):
        for match in re.findall(r"HMDB\d+", payload, flags=re.IGNORECASE):
            if normalized := _normalise_hmdb(match):
                found.add(normalized)
    return found


def fetch_workbench_lookups(cache_dir: Path, kegg_ids: Iterable[str]) -> pd.DataFrame:
    """Fetch one Metabolomics Workbench lookup per KEGG ID and cache the table."""
    requested = sorted({str(value) for value in kegg_ids if str(value) and str(value) != "nan"})
    path = cache_dir / "metabolomics_workbench_kegg_hmdb.tsv"
    if path.exists():
        table = pd.read_csv(path, sep="\t", dtype=str)
        if set(requested).issubset(set(table["kegg_id"])):
            return table

    rows: list[dict[str, Any]] = []
    for kegg in requested:
        url = WORKBENCH_URL.format(kegg=quote(kegg))
        try:
            with urlopen(Request(url, headers={"User-Agent": "networkbio-stage03/1.0"}), timeout=60) as response:
                payload = json.loads(response.read().decode("utf-8"))
            hmdb_ids = sorted(_extract_hmdb_values(payload))
            for hmdb in hmdb_ids or [None]:
                rows.append({"kegg_id": kegg, "hmdb_id": hmdb, "source": "metabolomics_workbench", "error": None})
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
            rows.append({"kegg_id": kegg, "hmdb_id": None, "source": "metabolomics_workbench", "error": f"{type(exc).__name__}: {exc}"})
    table = pd.DataFrame(rows, columns=["kegg_id", "hmdb_id", "source", "error"])
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)
    _resource_sidecar(path).write_text(
        json.dumps(
            {
                "resource": "Metabolomics Workbench KEGG compound REST lookups",
                "url_template": WORKBENCH_URL,
                "downloaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "sha256": sha256(path),
                "row_count": len(table),
            },
            indent=2,
        )
    )
    return table


def map_metabolite_features(
    feature_tables: Mapping[str, pd.DataFrame],
    *,
    pkn_nodes: set[str],
    recon_by_kegg: Mapping[str, list[str]],
    workbench: pd.DataFrame,
) -> pd.DataFrame:
    """Map stage-01/00 features to canonical COSMOS metabolite nodes."""
    workbench_by_kegg: dict[str, list[str]] = defaultdict(list)
    if not workbench.empty:
        for row in workbench.itertuples(index=False):
            if pd.notna(row.hmdb_id):
                workbench_by_kegg[str(row.kegg_id)].append(str(row.hmdb_id))
    compartments = {"ST003331": ("c", "m"), "ST003332": ("e",), "ST003328": ("c", "m")}
    rows: list[dict[str, Any]] = []
    for study, features in feature_tables.items():
        for feature in features.to_dict(orient="records"):
            raw_kegg = feature.get("kegg_id")
            kegg = "" if raw_kegg is None or pd.isna(raw_kegg) else str(raw_kegg).strip()
            name = str(feature.get("metabolite_name") or "")
            refmet = str(feature.get("refmet_name") or "")
            label = feature.get("labelling_class")
            if pd.isna(label) or label is None or str(label) == "":
                label = "not_labelled" if study == "ST003328" else "unknown"
            mapping_source = None
            hmdb_ids: list[str] = []
            if kegg and kegg in recon_by_kegg:
                hmdb_ids = recon_by_kegg[kegg]
                mapping_source = "recon3d"
            elif kegg and kegg in workbench_by_kegg:
                hmdb_ids = sorted(set(workbench_by_kegg[kegg]))
                mapping_source = "metabolomics_workbench"
            elif study == "ST003328" and (name.lower() == "cholesterol" or refmet.lower() == "cholesterol"):
                hmdb_ids = [CHOLESTEROL_HMDB]
                mapping_source = "lipid_name_hmdb"

            if not hmdb_ids:
                rows.append(
                    {
                        "study": study,
                        "feature_id": feature.get("feature_id"),
                        "feature_name": name,
                        "kegg_id": kegg or None,
                        "hmdb_id": None,
                        "compartment": None,
                        "pkn_node": None,
                        "mapping_source": "unmapped",
                        "labelling_class": label,
                        "mapped_to_pkn": False,
                    }
                )
                continue
            for hmdb in sorted(set(hmdb_ids)):
                for compartment in compartments.get(study, ("c",)):
                    pkn_node = f"Metab__{hmdb}_{compartment}"
                    rows.append(
                        {
                            "study": study,
                            "feature_id": feature.get("feature_id"),
                            "feature_name": name,
                            "kegg_id": kegg or None,
                            "hmdb_id": hmdb,
                            "compartment": compartment,
                            "pkn_node": pkn_node,
                            "mapping_source": mapping_source,
                            "labelling_class": label,
                            "mapped_to_pkn": pkn_node in pkn_nodes,
                        }
                    )
    return pd.DataFrame(rows)


def _is_metabolite(node: str) -> bool:
    return str(node).startswith("Metab__") or str(node).startswith("XMetab__")


def _signed_bfs(
    source: str,
    adjacency: Mapping[str, list[tuple[str, int]]],
    *,
    max_length: int,
    keep_predecessors: bool = False,
) -> tuple[dict[tuple[str, int], int], dict[tuple[str, int], int], dict[tuple[str, int], bool], dict[tuple[str, int], list[tuple[str, int, int]]]]:
    start = (source, 1)
    distance = {start: 0}
    count = {start: 1}
    through_metabolite = {start: _is_metabolite(source)}
    predecessors: dict[tuple[str, int], list[tuple[str, int, int]]] = defaultdict(list)
    queue: deque[tuple[str, int]] = deque([start])
    while queue:
        node, state_sign = queue.popleft()
        current_distance = distance[(node, state_sign)]
        if current_distance >= max_length:
            continue
        for target, edge_sign in adjacency.get(node, []):
            next_state = (target, state_sign * edge_sign)
            next_distance = current_distance + 1
            next_through = through_metabolite[(node, state_sign)] or _is_metabolite(target)
            if next_state not in distance:
                distance[next_state] = next_distance
                count[next_state] = count[(node, state_sign)]
                through_metabolite[next_state] = next_through
                if keep_predecessors:
                    predecessors[next_state].append((node, state_sign, edge_sign))
                queue.append(next_state)
            elif distance[next_state] == next_distance:
                count[next_state] += count[(node, state_sign)]
                through_metabolite[next_state] |= next_through
                if keep_predecessors:
                    predecessors[next_state].append((node, state_sign, edge_sign))
    return distance, count, through_metabolite, predecessors


def signed_shortest_paths(
    edges: Iterable[Edge] | pd.DataFrame,
    sources: Iterable[str],
    targets: Iterable[str],
    *,
    max_length: int = 8,
) -> pd.DataFrame:
    """Summarize shortest paths for every source-target pair and sign."""
    if isinstance(edges, pd.DataFrame):
        edge_list = list(edges[["source", "sign", "target"]].itertuples(index=False, name=None))
    else:
        edge_list = [(str(source), int(sign), str(target)) for source, sign, target in edges]
    adjacency: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for source, sign, target in edge_list:
        adjacency[source].append((target, int(sign)))
    for source in adjacency:
        adjacency[source].sort()
    rows: list[dict[str, Any]] = []
    for source in sources:
        source = str(source)
        distance, count, through_metabolite, _ = _signed_bfs(source, adjacency, max_length=max_length)
        for target in targets:
            target = str(target)
            positive = (target, 1)
            negative = (target, -1)
            positive_length = distance.get(positive)
            negative_length = distance.get(negative)
            lengths = [length for length in (positive_length, negative_length) if length is not None]
            shortest_length = min(lengths) if lengths else None
            shortest_signs = []
            if shortest_length is not None and positive_length == shortest_length:
                shortest_signs.append("+1")
            if shortest_length is not None and negative_length == shortest_length:
                shortest_signs.append("-1")
            rows.append(
                {
                    "source": source,
                    "target": target,
                    "shortest_length": shortest_length,
                    "shortest_signs": ";".join(shortest_signs),
                    "positive_shortest_length": positive_length,
                    "negative_shortest_length": negative_length,
                    "n_shortest_positive": count.get(positive, 0) if shortest_length is not None and positive_length == shortest_length else 0,
                    "n_shortest_negative": count.get(negative, 0) if shortest_length is not None and negative_length == shortest_length else 0,
                    "n_positive_shortest_paths": count.get(positive, 0),
                    "n_negative_shortest_paths": count.get(negative, 0),
                    "positive_shortest_passes_metabolite": bool(through_metabolite.get(positive, False)),
                    "negative_shortest_passes_metabolite": bool(through_metabolite.get(negative, False)),
                    "any_shortest_path_passes_metabolite": any(
                        through_metabolite.get(state, False)
                        for state in (positive, negative)
                        if distance.get(state) == shortest_length
                    ),
                    "positive_route_present": positive_length is not None and positive_length <= max_length,
                }
            )
    return pd.DataFrame(rows)


def _iter_paths(
    state: tuple[str, int],
    source: str,
    predecessors: Mapping[tuple[str, int], list[tuple[str, int, int]]],
    *,
    limit: int,
) -> Iterator[list[tuple[str, int, str]]]:
    if limit <= 0:
        return
    if state == (source, 1):
        yield []
        return
    emitted = 0
    for previous_node, previous_sign, edge_sign in predecessors.get(state, []):
        for prefix in _iter_paths((previous_node, previous_sign), source, predecessors, limit=limit - emitted):
            yield prefix + [(previous_node, edge_sign, state[0])]
            emitted += 1
            if emitted >= limit:
                return


def shortest_path_edges(
    edges: Iterable[Edge] | pd.DataFrame,
    paths: pd.DataFrame,
    *,
    max_paths_per_pair_sign: int = 100,
) -> pd.DataFrame:
    """Return the union of representative shortest positive/negative paths."""
    if isinstance(edges, pd.DataFrame):
        edge_list = list(edges[["source", "sign", "target"]].itertuples(index=False, name=None))
    else:
        edge_list = [(str(source), int(sign), str(target)) for source, sign, target in edges]
    adjacency: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for source, sign, target in edge_list:
        adjacency[source].append((target, int(sign)))
    for source in adjacency:
        adjacency[source].sort()
    rows: list[dict[str, Any]] = []
    for source in paths["source"].drop_duplicates():
        source = str(source)
        source_path_lengths = paths[paths["source"].eq(source)][
            ["positive_shortest_length", "negative_shortest_length"]
        ].stack().dropna()
        max_path_length = int(source_path_lengths.max()) if not source_path_lengths.empty else 0
        distance, _, _, predecessors = _signed_bfs(
            source,
            adjacency,
            max_length=max_path_length,
            keep_predecessors=True,
        )
        source_paths = paths[paths["source"].eq(source)]
        for row in source_paths.itertuples(index=False):
            for path_sign, length in ((1, row.positive_shortest_length), (-1, row.negative_shortest_length)):
                if pd.isna(length):
                    continue
                state = (str(row.target), path_sign)
                if distance.get(state) != int(length):
                    continue
                for path_index, path in enumerate(
                    _iter_paths(state, source, predecessors, limit=max_paths_per_pair_sign), start=1
                ):
                    for edge_index, (edge_source, edge_sign, edge_target) in enumerate(path, start=1):
                        rows.append(
                            {
                                "source": source,
                                "target": str(row.target),
                                "path_sign": path_sign,
                                "path_length": int(length),
                                "path_index": path_index,
                                "edge_index": edge_index,
                                "edge_source": edge_source,
                                "edge_sign": edge_sign,
                                "edge_target": edge_target,
                            }
                        )
    return pd.DataFrame(rows)


def rewire_signed_graph(edges: Iterable[Edge] | pd.DataFrame, *, seed: int, accepted_swaps: int | None = None) -> list[Edge]:
    """Degree- and sign-preserving directed double-edge swaps."""
    if isinstance(edges, pd.DataFrame):
        current = [(str(s), int(sign), str(t)) for s, sign, t in edges[["source", "sign", "target"]].itertuples(index=False, name=None)]
    else:
        current = [(str(s), int(sign), str(t)) for s, sign, t in edges]
    if len({(s, t) for s, _, t in current}) != len(current):
        raise ValueError("Rewiring requires unique source-target pairs")
    target_swaps = accepted_swaps if accepted_swaps is not None else 3 * len(current)
    if target_swaps <= 0:
        return sorted(current)
    rng = random.Random(seed)
    accepted = 0
    attempts = 0
    max_attempts = max(1000, target_swaps * 200)
    edge_keys = {(source, target) for source, _, target in current}
    while accepted < target_swaps and attempts < max_attempts:
        attempts += 1
        first_index, second_index = rng.sample(range(len(current)), 2)
        first = current[first_index]
        second = current[second_index]
        source_a, sign_a, target_b = first
        source_c, sign_c, target_d = second
        proposed = [(source_a, sign_a, target_d), (source_c, sign_c, target_b)]
        proposed_keys = [(source, target) for source, _, target in proposed]
        if source_a == target_d or source_c == target_b or proposed_keys[0] == proposed_keys[1]:
            continue
        old_keys = {(source_a, target_b), (source_c, target_d)}
        if any(key in edge_keys and key not in old_keys for key in proposed_keys):
            continue
        current[first_index], current[second_index] = proposed
        edge_keys.difference_update(old_keys)
        edge_keys.update(proposed_keys)
        accepted += 1
    if accepted < target_swaps:
        raise RuntimeError(f"Could not complete {target_swaps} accepted swaps after {attempts} attempts")
    return sorted(current)


def positive_shortest_lengths(
    edges: Iterable[Edge] | pd.DataFrame,
    sources: Iterable[str],
    targets: Iterable[str],
    *,
    max_length: int,
) -> dict[tuple[str, str], int | None]:
    """Find only positive-sign path lengths for the P10 null efficiently."""
    if isinstance(edges, pd.DataFrame):
        edge_list = list(edges[["source", "sign", "target"]].itertuples(index=False, name=None))
    else:
        edge_list = [(str(source), int(sign), str(target)) for source, sign, target in edges]
    adjacency: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for source, sign, target in edge_list:
        adjacency[source].append((target, int(sign)))
    out: dict[tuple[str, str], int | None] = {}
    for source in sources:
        distance, _, _, _ = _signed_bfs(str(source), adjacency, max_length=max_length)
        for target in targets:
            out[(str(source), str(target))] = distance.get((str(target), 1))
    return out


def degree_preserving_null_summary(
    edges: pd.DataFrame,
    observed_paths: pd.DataFrame,
    *,
    n_replicates: int,
    seed_base: int,
) -> pd.DataFrame:
    """Run the registered P10 rewiring null for the Aim-5 pairs."""
    pairs = list(observed_paths[["source", "target"]].itertuples(index=False, name=None))
    target_nodes = sorted(set(observed_paths["target"]))
    source_nodes = sorted(set(observed_paths["source"]))
    counts = {pair: 0 for pair in pairs}
    valid = {pair: 0 for pair in pairs}
    for replicate in range(n_replicates):
        rewired = rewire_signed_graph(edges, seed=seed_base + replicate)
        observed_by_pair = observed_paths.set_index(["source", "target"])
        observed_lengths = {
            pair: observed_by_pair.loc[pair, "positive_shortest_length"]
            for pair in pairs
        }
        finite = [int(length) for length in observed_lengths.values() if pd.notna(length)]
        max_length = max(finite, default=0)
        rewired_lengths = positive_shortest_lengths(rewired, source_nodes, target_nodes, max_length=max_length)
        for pair in pairs:
            observed_length = observed_lengths[pair]
            if pd.isna(observed_length):
                continue
            valid[pair] += 1
            rewired_length = rewired_lengths[pair]
            if rewired_length is not None and rewired_length <= int(observed_length):
                counts[pair] += 1
    rows = []
    for source, target in pairs:
        rows.append(
            {
                "source": source,
                "target": target,
                "n_replicates": n_replicates,
                "n_observed_positive": valid[(source, target)],
                "n_rewired_no_longer_than_observed": counts[(source, target)],
                "fraction_no_longer_than_observed": (
                    counts[(source, target)] / valid[(source, target)] if valid[(source, target)] else np.nan
                ),
                "seed_base": seed_base,
            }
        )
    return pd.DataFrame(rows)
