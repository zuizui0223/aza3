#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "data/planning/chapter3_field_priority_v5.json"
LEDGER = ROOT / "data/planning/chapter3_field_campaign_priority_v5.csv"
DOC = ROOT / "docs/CHAPTER3_FIELD_PRIORITY_CALENDAR_V5.md"
README = ROOT / "README.md"


def main() -> int:
    d = json.loads(CONTRACT.read_text(encoding="utf-8"))
    assert d["contract_version"] == "chapter3_field_priority_v5"
    assert d["planning_species_baseline"] == 120
    assert d["nationwide_tree"]["wave1_nominal_individuals"] == 240
    assert d["nationwide_tree"]["full_nominal_target_capture"] == 290
    starter = d["focal_rad_staging"]["starter"]
    assert starter["populations_total"] == 10
    assert starter["primary_rad_total"] == 120
    assert starter["physical_total"] == 150
    completion = d["focal_rad_staging"]["conditional_completion"]
    assert completion["additional_primary_rad"] == 48
    assert completion["full_primary_rad"] == 168
    assert "Level-1 tree admission" in completion["opening_rule"]
    tech = d["target_capture_technical_batch"]
    assert tech == {
        "species": 12,
        "individuals_per_species": 2,
        "libraries": 24,
        "counts_toward_wave1": True,
        "purpose": "qualify extraction/capture/QC across geography, DNA quality and cytotype/genome-size classes without creating a separate biological quota"
    }
    assert [x["id"] for x in d["campaigns"]] == ["C1","C2","C3","C4","C5","C6"]
    assert d["current_state"]["physical_samples"] == 0
    assert all(v is False for k,v in d["current_state"].items() if k != "physical_samples")

    rows = list(csv.DictReader(LEDGER.open(encoding="utf-8", newline="")))
    assert len(rows) == 6
    assert rows[0]["campaign_id"] == "C1" and rows[0]["priority_class"] == "A"
    assert "missing Level-1 species always outranks extra RAD depth" in rows[-1]["stop_rule"]

    doc = DOC.read_text(encoding="utf-8")
    for term in (
        "tree breadth first",
        "Priority A",
        "September 2026",
        "late February to late March 2027",
        "Starter total = **120 primary RAD**",
        "Wave 1 — breadth",
        "collect the missing species-tree taxon",
    ):
        assert term in doc, term
    readme = README.read_text(encoding="utf-8")
    assert "CHAPTER3_FIELD_PRIORITY_CALENDAR_V5.md" in readme
    print("chapter3_field_priority_v5_valid=true")
    print("wave1_target_capture=240")
    print("starter_rad=120")
    print("field_sampling_authorized=false")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
