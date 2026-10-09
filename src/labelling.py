"""Labelling annotations for the 24-hour [U-13C6]glucose measurements.

The deposited files do not state whether a plain metabolite row is a total
pool or M+0.  The revised stage plan therefore keeps deposited isotope sums
explicitly labelled as partial and treats other glucose-reachable plain rows
as unlabelled-fraction measurements, not abundance measurements.
"""

from __future__ import annotations

import re

import pandas as pd


TRACED_METABOLITES = {
    "C00026": ("2-Oxoglutarate", 5),
    "C00417": ("cis-Aconitate", 6),
    "C00158": ("Citrate", 6),
    "C00122": ("Fumarate", 4),
    "C00031": ("Glucose", 6),
    "C00490": ("Itaconate", 5),
    "C00149": ("Malate", 4),
    "C01107": ("Phosphomevalonate", 6),
    "C00042": ("Succinate", 4),
}

WHITELIST_REASONS = {
    "C00135": "essential amino acid; glucose carbon cannot be incorporated by human synthesis",
    "C00123": "essential amino acid; glucose carbon cannot be incorporated by human synthesis",
    "C00047": "essential amino acid; glucose carbon cannot be incorporated by human synthesis",
    "C00073": "essential amino acid; glucose carbon cannot be incorporated by human synthesis",
    "C00079": "essential amino acid; glucose carbon cannot be incorporated by human synthesis",
    "C00188": "essential amino acid; glucose carbon cannot be incorporated by human synthesis",
    "C00078": "essential amino acid; glucose carbon cannot be incorporated by human synthesis",
    "C00183": "essential amino acid; glucose carbon cannot be incorporated by human synthesis",
    "C00082": "formed from phenylalanine; no direct glucose-carbon synthesis route expected",
    "C02989": "methionine oxidation product; no direct glucose-carbon synthesis route expected",
    "C05335": "selenometionine; no direct glucose-carbon synthesis route expected",
    "C00318": "carnitine derives from lysine and methionine; no direct glucose-carbon synthesis route expected",
    "C00072": "ascorbate is not synthesised in human cells",
    "C05422": "dehydroascorbate is not synthesised in human cells",
    "C00864": "vitamin B5 derivative; no direct glucose-carbon synthesis route expected",
    "C03492": "vitamin B5 derivative; no direct glucose-carbon synthesis route expected",
    "C01595": "essential fatty acid; no human de novo synthesis",
    "C02700": "tryptophan catabolite without added carbon in the listed route",
    "C00108": "tryptophan catabolite without added carbon in the listed route",
    "C10164": "tryptophan catabolite without added carbon in the listed route",
    "C03722": "tryptophan catabolite without added carbon in the listed route",
    "C00322": "tryptophan catabolite without added carbon in the listed route",
    "C00463": "tryptophan catabolite without added carbon in the listed route",
    "C05658": "tryptophan catabolite without added carbon in the listed route",
    "C00331": "tryptophan catabolite without added carbon in the listed route",
    "C00637": "tryptophan catabolite without added carbon in the listed route",
    "C05635": "tryptophan catabolite without added carbon in the listed route",
    "C02796": "tryptophan or tyrosine catabolite without added carbon in the listed route",
    "C00788": "tyrosine catabolite without added carbon in the listed route",
}


def _isotope_number(name: str) -> int:
    match = re.search(r"13C(\d+)", str(name))
    if match is None:
        raise ValueError(f"Could not read a 13C number from isotopologue name: {name}")
    return int(match.group(1))


def validate_traced_isotopologues(features: pd.DataFrame) -> dict[str, dict[str, object]]:
    """Check deposited isotope forms and return their present/missing forms."""
    required = {"kegg_id", "is_13c_isotopologue", "metabolite_name"}
    if not required.issubset(features.columns):
        raise ValueError(f"Missing labelling columns: {sorted(required - set(features.columns))}")
    report: dict[str, dict[str, object]] = {}
    for kegg_id, (name, carbon_count) in TRACED_METABOLITES.items():
        rows = features.loc[features["kegg_id"].eq(kegg_id)]
        if rows.empty:
            raise ValueError(f"Traced metabolite {kegg_id} ({name}) is absent")
        numbers = sorted({_isotope_number(value) for value in rows.loc[
            rows["is_13c_isotopologue"], "metabolite_name"
        ]})
        expected = set(range(1, carbon_count + 1))
        if not set(numbers).issubset(expected):
            raise ValueError(f"Unexpected 13C form for {kegg_id}: {numbers}")
        report[kegg_id] = {
            "name": name,
            "carbon_count": carbon_count,
            "deposited_13c": numbers,
            "missing_13c": sorted(expected - set(numbers)),
        }
    return report


def annotate_features(features: pd.DataFrame, study: str) -> pd.DataFrame:
    """Assign exactly one labelling class and traceability fields per row."""
    if study not in {"ST003331", "ST003332"}:
        raise ValueError(f"Labelling classes are only defined for ST003331/ST003332, got {study}")
    required = {"kegg_id", "is_13c_isotopologue", "metabolite_name"}
    if not required.issubset(features.columns):
        raise ValueError(f"Missing labelling columns: {sorted(required - set(features.columns))}")

    result = features.copy()
    result["labelling_class"] = "m0_confounded"
    result["labelling_reason"] = "plain row may contain glucose-derived carbon after 24 hours"
    result["deposited_13c"] = ""
    result["missing_13c"] = ""

    traced = validate_traced_isotopologues(features) if study == "ST003331" else {}
    for kegg_id, details in traced.items():
        mask = result["kegg_id"].eq(kegg_id)
        result.loc[mask, "labelling_class"] = "isotopologue_sum_partial"
        result.loc[mask, "labelling_reason"] = "deposited labelled forms are incomplete; sum is not a total pool"
        result.loc[mask, "deposited_13c"] = ",".join(map(str, details["deposited_13c"]))
        result.loc[mask, "missing_13c"] = ",".join(map(str, details["missing_13c"]))

    for kegg_id, reason in WHITELIST_REASONS.items():
        mask = result["kegg_id"].eq(kegg_id) & result["labelling_class"].eq("m0_confounded")
        result.loc[mask, "labelling_class"] = "m0_unconfounded"
        result.loc[mask, "labelling_reason"] = reason

    if study == "ST003332":
        mask = result["kegg_id"].eq("C00031")
        result.loc[mask, "labelling_reason"] = "13C6-glucose medium; uptake not measured by the M+0 row"

    if result["labelling_class"].isna().any() or not result["labelling_class"].isin(
        {"isotopologue_sum_partial", "m0_unconfounded", "m0_confounded"}
    ).all():
        raise ValueError("Every metabolomics feature must receive exactly one labelling class")
    return result


def whitelist_report(feature_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """List every whitelist ID as present or absent in each labelled study."""
    rows = []
    for study, features in feature_tables.items():
        present = set(features["kegg_id"].dropna())
        for kegg_id, reason in WHITELIST_REASONS.items():
            rows.append({"study": study, "kegg_id": kegg_id, "present": kegg_id in present,
                         "reason": reason})
    return pd.DataFrame(rows)
