"""Readers for the deposited input tables (stage 00_inputs).

These functions convert the deposited files into tidy tables. They do not
filter, normalise, impute or average anything: those are later, separate
decisions (see doc/agent-rules.md, "Preprocessing and quality control").
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import pandas as pd

# Group labels as deposited -> project labels.
GROUP_MAP = {"AMC": "Ctrl", "PMS": "PMS", "no cell ctrl": "no_cell_blank"}


# ---------------------------------------------------------------- RNA-seq (GEO)


def read_geo_series_samples(series_matrix: Path) -> pd.DataFrame:
    """Sample table from a GEO series matrix: title, GSM, characteristics."""
    fields: dict[str, list[list[str]]] = {}
    with gzip.open(series_matrix, "rt") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            fields.setdefault(key, []).append([v.strip('"') for v in vals])
    out = pd.DataFrame(
        {
            "sample_title": fields["!Sample_title"][0],
            "gsm": fields["!Sample_geo_accession"][0],
        }
    )
    for row in fields.get("!Sample_characteristics_ch1", []):
        key = row[0].split(":", 1)[0].strip().replace(" ", "_")
        out[key] = [v.split(":", 1)[1].strip() if ":" in v else v for v in row]
    return out


def read_rnaseq_counts(path: Path) -> pd.DataFrame:
    """Raw count matrix: rows = Ensembl gene IDs, columns = sample titles."""
    counts = pd.read_csv(path, index_col=0)
    counts.index.name = "ensembl_gene_id"
    return counts


# ------------------------------------------------------ Metabolomics Workbench


def _mwtab_blocks(path: Path) -> tuple[list[list[str]], list[list[str]]]:
    """Return (SUBJECT_SAMPLE_FACTORS rows, METABOLITES block rows) of one mwTab file."""
    ssf, met, in_met = [], [], False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("SUBJECT_SAMPLE_FACTORS"):
            ssf.append(line.split("\t")[1:])
        elif line.startswith("METABOLITES_START"):
            in_met = True
        elif line.startswith("METABOLITES_END"):
            in_met = False
        elif in_met:
            met.append(line.split("\t"))
    return ssf, met


def _parse_factors(text: str) -> dict[str, str]:
    out = {}
    for part in text.split("|"):
        if ":" in part:
            k, v = part.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def parse_mw_sample_id(sample_id: str) -> dict[str, str | None]:
    """Line letter and replicate from a deposited sample name.

    Patterns seen: 'AMC A SV P1', 'PMA B P2', 'PMS D S9', 'NC 7'. ST003331 uses
    run IDs ('DN47-01') with no line in the name; those return None here and
    the line comes from the mwTab subject column instead.
    """
    m = re.fullmatch(r"(AMC|PMA|PMS)\s+([A-Z])\s+(?:SV\s+)?([PS]\d+)", sample_id)
    if m:
        return {"name_prefix": m.group(1), "line_letter": m.group(2), "replicate": m.group(3)}
    m = re.fullmatch(r"NC\s+(\d+)", sample_id)
    if m:
        return {"name_prefix": "NC", "line_letter": None, "replicate": m.group(1)}
    return {"name_prefix": None, "line_letter": None, "replicate": None}


def read_mw_study(study_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read one Metabolomics Workbench study.

    Returns:
        samples: one row per sample with group, treatment, line and replicate.
        features: one row per (analysis, metabolite) with the deposited IDs.
        values: features x samples, peak areas as deposited (float, NaN = missing).
    """
    study = study_dir.name
    data = json.loads((study_dir / f"{study}_data.json").read_text())
    factors = json.loads((study_dir / f"{study}_factors.json").read_text())

    # Subject IDs and the per-analysis metabolite annotation live in the mwTab files.
    subject_of: dict[str, str] = {}
    annot: dict[tuple[str, str], dict[str, str]] = {}
    for mwtab in sorted(study_dir.glob(f"{study}_AN*.mwtab.txt")):
        analysis = mwtab.name.split("_")[1].split(".")[0]
        ssf, met = _mwtab_blocks(mwtab)
        for row in ssf:
            if len(row) >= 2:
                subject_of[row[1]] = row[0]
        if met:
            header = met[0]
            for row in met[1:]:
                rec = dict(zip(header, row))
                annot[(analysis, rec.get("metabolite_name", ""))] = rec

    sample_rows = []
    for rec in factors.values():
        sid = rec["local_sample_id"]
        fac = _parse_factors(rec["factors"])
        disease = next((v for k, v in fac.items() if k.lower().startswith("disease")), None)
        parsed = parse_mw_sample_id(sid)
        subject = subject_of.get(sid)
        line_letter = parsed["line_letter"] or (subject if subject not in (None, "NA", "-") else None)
        group = GROUP_MAP.get(disease, disease)
        sample_rows.append(
            {
                "study": study,
                "sample_id": sid,
                "mb_sample_id": rec.get("mb_sample_id"),
                "sample_source": rec.get("sample_source"),
                "disease_as_deposited": disease,
                "group": group,
                "treatment": fac.get("Treatment", "untreated"),
                "subject_as_deposited": subject,
                "name_prefix": parsed["name_prefix"],
                "line_letter": line_letter,
                # Line key is only meaningful within one study (letters are reused).
                "line_key": f"{study}:{group}:{line_letter}" if line_letter else None,
                "replicate": parsed["replicate"],
            }
        )
    samples = pd.DataFrame(sample_rows).sort_values(["group", "treatment", "line_letter", "sample_id"])

    feat_rows, value_rows = [], {}
    for key, rec in data.items():
        analysis = rec["analysis_id"]
        name = rec["metabolite_name"]
        fid = f"{analysis}|{name}"
        a = annot.get((analysis, name), {})
        feat_rows.append(
            {
                "feature_id": fid,
                "study": study,
                "analysis_id": analysis,
                "analysis_summary": rec.get("analysis_summary"),
                "metabolite_name": name,
                "refmet_name": rec.get("refmet_name") or None,
                "kegg_id": (a.get("KEGG ID") or "").strip() or None,
                "lipid_group": (a.get("LipidGroup") or "").strip() or None,
                "is_13c_isotopologue": bool(re.search(r"13C\d*", name)),
                "units": rec.get("units"),
            }
        )
        value_rows[fid] = {s: _to_float(v) for s, v in rec["DATA"].items()}
    features = pd.DataFrame(feat_rows)
    values = pd.DataFrame.from_dict(value_rows, orient="index")
    values = values[samples["sample_id"].tolist()]
    values.index.name = "feature_id"
    return samples, features, values


def _to_float(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")
