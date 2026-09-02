#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "data/planning/chapter3_nationwide_species_tree_sampling_v4.json"
BUDGET = ROOT / "data/planning/chapter3_nationwide_species_tree_budget_v4.csv"
DOC = ROOT / "docs/CHAPTER3_NATIONWIDE_SPECIES_TREE_SAMPLING_V4.md"
README = ROOT / "README.md"


def validate() -> dict:
    d = json.loads(CONTRACT.read_text(encoding="utf-8"))
    assert d["contract_version"] == "chapter3_nationwide_species_tree_sampling_v4"
    assert d["authoritative_for_new_sampling"] is True
    assert d["primary_product"] == "JAPAN_WIDE_NUCLEAR_SPECIES_TREE_AND_NETWORK"
    c = d["taxonomic_census"]
    assert c["nmns_authority_records_snapshot"] == 161
    assert c["planning_species_baseline"] == 120
    assert c["planning_species_upper_envelope"] == 125
    assert c["count_status"] == "PLANNING_BASELINE_NOT_FINAL_OPERATIONAL_CENSUS"
    n = d["nationwide_sampling"]
    assert n["base_individuals_per_species"] == 2
    assert n["nominal_base_total_at_120_species"] == 240
    assert n["nominal_full_bank_at_120_species"] == 290
    assert n["full_bank_envelope_if_125_species"] == 300
    assert [x["nominal_additional_individuals"] for x in n["enrichment"]] == [30,20]
    tc = d["target_capture"]
    assert tc["role"] == "PRIMARY_CROSS_SPECIES_SCAFFOLD"
    assert tc["wave_1"]["nominal_individuals_at_120_species"] == 240
    assert tc["wave_2"]["nominal_max_additional_individuals"] == 50
    assert tc["nominal_full_target_capture_individuals"] == 290
    r = d["population_rad"]
    assert r["target_focal_primary_rad_total"] == 168
    assert r["nominal_tree_rad_overlap_individuals"] == 20
    assert r["nominal_additional_unique_focal_plants_beyond_290_bank"] == 148
    assert r["nominal_total_unique_physical_plants_if_full_v4_is_completed"] == 438
    assert sum(x["primary_rad"] for x in r["target_focal_design"]) == 168
    assert d["analysis_levels"] == [
        "LEVEL_0_TAXONOMIC_CENSUS",
        "LEVEL_1_NATIONWIDE_TARGET_CAPTURE_SPECIES_TREE_NETWORK",
        "LEVEL_2_TRAIT_HISTORY_MAPPING_ON_NATIONWIDE_ENSEMBLE",
        "LEVEL_3_NESTED_POPULATION_RAD_FOR_FOCAL_TRANSITIONS",
        "LEVEL_4_M01_MECHANISTIC_DECOMPOSITION",
    ]
    s = d["current_state"]
    assert s["physical_samples"] == 0
    assert all(v is False for k,v in s.items() if k != "physical_samples")

    rows = list(csv.DictReader(BUDGET.open(encoding="utf-8", newline="")))
    by = {x["component"]: x for x in rows}
    for comp, expected in {
        "T1_BASE_ALL_SPECIES":"240",
        "NATIONWIDE_TARGET_CAPTURE_FULL":"290",
        "NATIONWIDE_UPPER_ENVELOPE":"300",
        "RAD_FOCAL_TOTAL":"168",
        "TREE_RAD_OVERLAP":"20",
        "FULL_UNIQUE_PHYSICAL_NOMINAL":"438",
    }.items():
        assert by[comp]["nominal_individuals"] == expected

    doc = DOC.read_text(encoding="utf-8")
    for term in ("Japan-wide Cirsium species tree/network", "120 species × 2 individuals", "290 individuals", "300 individuals", "168", "438", "target capture, not RAD"):
        assert term in doc, term
    readme = README.read_text(encoding="utf-8")
    assert "CHAPTER3_NATIONWIDE_SPECIES_TREE_SAMPLING_V4.md" in readme
    assert "nationwide species tree" in readme.lower()
    return d


def main() -> int:
    for p in (CONTRACT,BUDGET,DOC,README):
        assert p.exists() and p.stat().st_size > 0, p
    d = validate()
    print("chapter3_nationwide_species_tree_sampling_v4_valid=true")
    print(f"planning_species_baseline={d['taxonomic_census']['planning_species_baseline']}")
    print(f"target_capture_nominal={d['nationwide_sampling']['nominal_full_bank_at_120_species']}")
    print(f"focal_rad_primary={d['population_rad']['target_focal_primary_rad_total']}")
    print("sampling_authorized=false")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
