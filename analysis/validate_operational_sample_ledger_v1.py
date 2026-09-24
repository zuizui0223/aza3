#!/usr/bin/env python3
"""Validate Chapter 3 operational biological-sample ledgers."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "data" / "planning" / "chapter3_core_operational_sample_ledger_v1.csv"
M01_PATH = ROOT / "data" / "planning" / "m01_operational_population_ledger_v1.csv"
DOC_PATH = ROOT / "docs" / "CHAPTER3_OPERATIONAL_SAMPLE_LEDGER_V1.md"
README_PATH = ROOT / "README.md"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_core() -> list[dict[str, str]]:
    data = rows(CORE_PATH)
    ids = [row["member_id"] for row in data]
    if ids != [f"JPN_{i:02d}" for i in range(1, 39)]:
        raise AssertionError("core operational ledger must contain ordered JPN_01-JPN_38 exactly once")

    minimum = sum(int(row["core_min_primary"]) for row in data)
    recommended = sum(int(row["core_recommended_primary"]) for row in data)
    if (minimum, recommended) != (167, 193):
        raise AssertionError(f"operational core totals drift: {(minimum, recommended)}")

    lookup = {row["member_id"]: row for row in data}
    expected_focal = {
        "JPN_06": (16, 24, "P02_FOCAL_HIGH"),
        "JPN_15": (16, 24, "P02_FOCAL_HIGH"),
        "JPN_36": (30, 40, "P01_FOCAL_VERY_HIGH"),
    }
    for member, (minimum_n, recommended_n, status) in expected_focal.items():
        row = lookup[member]
        if (int(row["core_min_primary"]), int(row["core_recommended_primary"]), row["priority_status"]) != (
            minimum_n,
            recommended_n,
            status,
        ):
            raise AssertionError(f"focal operational target drift: {member}")

    blocked = {
        "JPN_29": "BLOCKED_TAXON_CONCEPT_IDENTITY_GATE",
        "JPN_31": "BLOCKED_IDENTITY_GEOGRAPHY_RESOLUTION",
        "JPN_33": "BLOCKED_NAME_VOUCHER_CONCEPT_RESOLUTION",
    }
    for member, gate in blocked.items():
        row = lookup[member]
        if row["identity_gate"] != gate or "NO_COLLECTION_UNTIL_IDENTITY_GATE" not in row["acquisition_route"]:
            raise AssertionError(f"identity-blocked concept was opened: {member}")
        if row["replacement_rule"] != "NO_CONVENIENCE_REPLACEMENT":
            raise AssertionError(f"blocked concept allows convenience replacement: {member}")

    wild_repair = {"JPN_32", "JPN_34", "JPN_35", "JPN_36", "JPN_37", "JPN_38"}
    for member in wild_repair:
        if "WILD_JAPAN_RESAMPLE" not in lookup[member]["identity_gate"]:
            raise AssertionError(f"wild-Japan provenance repair missing: {member}")

    forbidden_columns = {"latitude", "longitude", "exact_locality", "permit_document"}
    if forbidden_columns.intersection(data[0].keys()):
        raise AssertionError("sensitive locality fields entered public core operational ledger")

    if "two Moreyra tree codes" not in lookup["JPN_20"]["population_rule"]:
        raise AssertionError("JPN20 paper-concept collapse rule drift")
    jpn38 = lookup["JPN_38"]
    if jpn38["existing_reference_class"] != "Japanese_distributed_taxon_sampled_outside_Japan":
        raise AssertionError("JPN38 outside-Japan reference provenance drift")
    if "wild Japanese individuals" not in jpn38["population_rule"]:
        raise AssertionError("JPN38 wild-Japan replacement target drift")
    return data


def validate_m01() -> list[dict[str, str]]:
    data = rows(M01_PATH)
    lookup = {row["sample_unit_id"]: row for row in data}

    discovery = [
        "M01_BREV_OKI",
        "M01_BREV_AMAMI",
        "M01_IRUM_MIYAKO",
        "M01_IRUM_ISHIGAKI",
    ]
    for sample_id in discovery:
        row = lookup.get(sample_id)
        if row is None:
            raise AssertionError(f"missing discovery population: {sample_id}")
        counts = tuple(int(row[field]) for field in (
            "primary_individuals",
            "reserve_individuals",
            "pigment_primary_n",
            "rna_collect_n",
            "rna_sequence_primary_n",
        ))
        if counts != (15, 3, 3, 6, 5):
            raise AssertionError(f"discovery nested counts drift: {sample_id} {counts}")

    e3_focal = [row for row in data if row["taxon_or_role"] in {"Cirsium brevicaule", "Cirsium irumtiense"}]
    if len(e3_focal) != 8 or sum(int(row["primary_individuals"]) for row in e3_focal) != 120:
        raise AssertionError("M01 final focal E3 panel must remain 8 populations x 15 = 120")
    if sum(int(lookup[s]["primary_individuals"]) for s in discovery) != 60:
        raise AssertionError("M01 discovery primary bank must remain 60")

    required_outgroups = {
        "Cirsium morii",
        "Cirsium pengii",
        "Cirsium tatakaense",
        "Cirsium kawakamii",
        "Cirsium japonicum var. albescens",
        "Cirsium japonicum var. fukienense",
    }
    out_rows = [row for row in data if row["population_role"] in {"E1_SISTER_CLADE_BRACKET", "E1_STATE_BALANCED_DEEP_BRACKET"}]
    if {row["taxon_or_role"] for row in out_rows} != required_outgroups:
        raise AssertionError("M01 six-concept ancestral-state bracket drift")
    if sum(int(row["primary_individuals"]) for row in out_rows) != 18:
        raise AssertionError("M01 minimum nonfocal E1 bracket must remain 18 individuals")

    controls = {row["taxon_or_role"] for row in data if row["population_role"] == "OPTIONAL_MOLECULAR_POSITIVE_CONTROL"}
    expected_controls = {
        "Cirsium japonicum var. takaoense white morph",
        "Cirsium japonicum var. takaoense bluish-purple morph",
    }
    if controls != expected_controls:
        raise AssertionError("takaoense optional molecular-control pair drift")

    if lookup["M01_IRUM_ISHIGAKI"]["population_role"] != "DISCOVERY_AND_E3":
        raise AssertionError("Ishigaki must remain preregistered Yaeyama discovery population")
    if lookup["M01_IRUM_IRIOMOTE"]["population_role"] != "E3_EXPANSION":
        raise AssertionError("Iriomote was silently promoted into discovery after outcomes")

    forbidden_columns = {"latitude", "longitude", "exact_locality", "permit_document"}
    if forbidden_columns.intersection(data[0].keys()):
        raise AssertionError("sensitive locality fields entered public M01 operational ledger")
    return data


def validate_narrative() -> None:
    text = DOC_PATH.read_text(encoding="utf-8")
    required = [
        "target count is not admission",
        "JPN_29 Cirsium verutum",
        "JPN_31 Cirsium yuki-uenoanum",
        "JPN_33 Cirsium effusum",
        "Ishigaki is frozen as the discovery Yaeyama population",
        "C. morii",
        "C. pengii",
        "C. tatakaense",
        "C. kawakamii",
        "C. japonicum var. albescens",
        "C. japonicum var. fukienense",
        "60 of the final 120",
        "does not invent one",
    ]
    missing = [item for item in required if item not in text]
    if missing:
        raise AssertionError(f"operational ledger narrative missing: {missing}")

    readme = README_PATH.read_text(encoding="utf-8")
    if "CHAPTER3_OPERATIONAL_SAMPLE_LEDGER_V1.md" not in readme:
        raise AssertionError("README does not route to operational sample ledger")


def main() -> int:
    for path in (CORE_PATH, M01_PATH, DOC_PATH, README_PATH):
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing or empty operational ledger file: {path.relative_to(ROOT)}")
    core = validate_core()
    m01 = validate_m01()
    validate_narrative()
    print("chapter3_operational_sample_ledger_valid=true")
    print(f"core_concepts={len(core)}")
    print(f"core_minimum_primary={sum(int(r['core_min_primary']) for r in core)}")
    print(f"core_recommended_primary={sum(int(r['core_recommended_primary']) for r in core)}")
    print(f"m01_rows={len(m01)}")
    print("blocked_identity_concepts=JPN_29,JPN_31,JPN_33")
    print("sampling_authorized=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
