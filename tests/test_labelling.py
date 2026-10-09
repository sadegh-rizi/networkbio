"""Synthetic checks for the registered 13C labelling interpretation."""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from labelling import (  # noqa: E402
    TRACED_METABOLITES,
    WHITELIST_REASONS,
    annotate_features,
    validate_traced_isotopologues,
    whitelist_report,
)


def traced_fixture() -> pd.DataFrame:
    deposited = {
        "C00026": [2, 3, 4, 5], "C00417": [2, 3, 4], "C00158": [2, 3, 4, 5],
        "C00122": [2, 3, 4], "C00031": [6], "C00490": [2], "C00149": [2, 3, 4],
        "C01107": [4], "C00042": [2, 3, 4],
    }
    rows = []
    for kegg_id, (name, carbon_count) in TRACED_METABOLITES.items():
        rows.append({"feature_id": f"{kegg_id}_m0", "kegg_id": kegg_id,
                     "is_13c_isotopologue": False, "metabolite_name": name,
                     "analysis_id": "neg"})
        for number in deposited[kegg_id]:
            rows.append({"feature_id": f"{kegg_id}_13c{number}", "kegg_id": kegg_id,
                         "is_13c_isotopologue": True, "metabolite_name": f"{name} 13C{number}",
                         "analysis_id": "neg"})
    rows.append({"feature_id": "whitelist", "kegg_id": "C00135",
                 "is_13c_isotopologue": False, "metabolite_name": "L-histidine",
                 "analysis_id": "neg"})
    rows.append({"feature_id": "unknown", "kegg_id": "C00031",
                 "is_13c_isotopologue": False, "metabolite_name": "D-Glucose",
                 "analysis_id": "neg"})
    return pd.DataFrame(rows)


def test_carbon_counts_and_missing_forms_are_explicit():
    features = traced_fixture()
    report = validate_traced_isotopologues(features)
    assert report["C00026"]["deposited_13c"] == [2, 3, 4, 5]
    assert report["C00026"]["missing_13c"] == [1]
    assert report["C00417"]["missing_13c"] == [1, 5, 6]


def test_each_labelled_feature_gets_one_class_and_whitelist_absences_are_listed():
    features = traced_fixture()
    annotated = annotate_features(features, "ST003331")
    assert annotated["labelling_class"].notna().all()
    assert annotated["labelling_class"].isin(
        {"isotopologue_sum_partial", "m0_unconfounded", "m0_confounded"}
    ).all()
    assert annotated.loc[annotated["kegg_id"].eq("C00135"), "labelling_class"].eq("m0_unconfounded").all()

    medium = annotate_features(features.loc[features["is_13c_isotopologue"].eq(False)], "ST003332")
    glucose = medium.loc[medium["kegg_id"].eq("C00031")].iloc[0]
    assert glucose["labelling_class"] == "m0_confounded"
    assert "M+0" in glucose["labelling_reason"]

    report = whitelist_report({"ST003331": features, "ST003332": medium})
    absent = report.loc[~report["present"]]
    assert set(absent["kegg_id"]) == set(WHITELIST_REASONS) - {"C00135"}
