#!/usr/bin/env python3
"""Fail-closed validation of the EAzami Chapter 2 to aza3 Chapter 3 handoff."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "data" / "contracts" / "chapter3_eazami_handoff_contract_v1.json"
PRIORITY_PATH = ROOT / "data" / "planning" / "chapter3_sampling_priorities_v1.csv"
PRIOR_PATH = ROOT / "data" / "planning" / "chapter3_bounded_prior_registry_v1.csv"
PROTOCOL_PATH = ROOT / "data" / "planning" / "chapter3_protocol_registry_v1.csv"
INTAKE_PATH = ROOT / "data" / "intake" / "chapter3_individual_intake_v1.csv"
SCOPE_PATH = ROOT / "docs" / "CHAPTER3_SCOPE_AND_HANDOFF_V1.md"
README_PATH = ROOT / "README.md"

SOURCE_MERGE_SHA = "4fc03f128a7ec05ce9e16e1daedef23b61104b89"
SOURCE_PRIORITY_SHA256 = "48066b17bb3767f8aa23e4945184125d83d3f71e09ecdfd4aa94352e15fec252"

INTAKE_FIELDS = [
    "individual_id", "taxon_concept", "population_id", "deidentified_locality_key",
    "voucher_id", "rad_tissue_id", "rad_library_id", "orientation_state",
    "phyllary_state", "stickiness_state", "cytotype", "access_authorization_id",
    "collection_authorization_id", "conservation_review_id", "admission_status",
    "exclusion_reason",
]


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_contract() -> dict:
    contract = load_json(CONTRACT_PATH)
    if contract.get("contract_version") != "chapter3_eazami_handoff_contract_v1":
        raise AssertionError("contract version drift")
    if contract.get("source_repository") != "https://github.com/zuizui0223/EAzami":
        raise AssertionError("EAzami source repository drift")
    if contract.get("source_merge_sha") != SOURCE_MERGE_SHA:
        raise AssertionError("handoff is not pinned to the merged EAzami Chapter 2 commit")
    if contract.get("frozen_text_hash_semantics") != (
        "SHA-256 of UTF-8 text after CRLF-to-LF normalization"
    ):
        raise AssertionError("source hash semantics drift")

    sources = contract.get("source_files", [])
    if len(sources) != 6:
        raise AssertionError("source manifest must contain six audited EAzami inputs")
    paths = [row["path"] for row in sources]
    if len(paths) != len(set(paths)):
        raise AssertionError("duplicate source path in provenance manifest")
    for row in sources:
        if not re.fullmatch(r"[0-9a-f]{64}", row.get("canonical_text_sha256", "")):
            raise AssertionError(f"invalid source hash for {row.get('path')}")
        if not row.get("role", "").strip():
            raise AssertionError(f"empty source role for {row.get('path')}")
    source_priority = next(
        row for row in sources
        if row["path"] == "data/evidence/chapter2_to_chapter3_sampling_priorities_v1.csv"
    )
    if source_priority["canonical_text_sha256"] != SOURCE_PRIORITY_SHA256:
        raise AssertionError("source sampling-priority provenance drift")

    locked = contract["chapter2_locked_results"]
    expected = {
        "orientation_minimum_changes": (4, [4, 6]),
        "phyllary_posture_minimum_changes": (3, [3, 3]),
        "stickiness_minimum_changes": (5, [5, 5]),
    }
    for key, (lower, range_) in expected.items():
        if locked[key] != {"lower_bound": lower, "ensemble_range": range_}:
            raise AssertionError(f"locked Chapter 2 result drift: {key}")
    if locked["orientation_forced_ml_edges"] != 0:
        raise AssertionError("orientation event-resolution result drift")
    if abs(locked["jpn36_phyllary_terminal_fraction"] - 0.754) > 1e-12:
        raise AssertionError("JPN36 event-resolution result drift")
    if (locked["species_tip_compression_systems"], locked["morph_linked_testable_systems"]) != (4, 1):
        raise AssertionError("cross-scale identifiability result drift")

    state = contract["current_state"]
    required_false = (
        "field_execution_authorized", "tissue_collection_authorized",
        "definitive_japan_wide_species_tree_claim_authorized",
        "sensitive_coordinates_allowed_in_repository",
    )
    if state["own_biological_records_admitted"] != 0:
        raise AssertionError("own biological data were admitted without a new contract")
    for field in required_false:
        if state[field] is not False:
            raise AssertionError(f"fail-closed current state was opened: {field}")
    if state["jpn36_protocol_state"] != "TECHNICALLY_READY_NOT_FIELD_AUTHORIZED":
        raise AssertionError("JPN36 technical readiness was confused with authorization")
    if "not authorized" not in contract["failure_actions"]["field_protocol_gate"]:
        raise AssertionError("field-protocol failure action drift")
    return contract


def validate_sampling_priorities() -> list[dict[str, str]]:
    rows = read_rows(PRIORITY_PATH)
    if [row["priority_id"] for row in rows] != ["P01", "P02", "P03", "P04", "P05"]:
        raise AssertionError("sampling priorities must contain ordered P01-P05")
    required = {
        "priority_id", "rank", "trait_module", "focal_concepts", "chapter2_locked_result",
        "chapter3_discriminator", "minimum_own_data", "predeclared_falsifier",
        "authorization_gate", "claim_boundary", "source_priority_id",
    }
    if not rows or set(rows[0]) != required:
        raise AssertionError("sampling-priority schema drift")
    for rank, row in enumerate(rows, start=1):
        if int(row["rank"]) != rank or row["source_priority_id"] != row["priority_id"]:
            raise AssertionError(f"sampling-priority rank/source drift for {row['priority_id']}")
        for field in required:
            if not row[field].strip():
                raise AssertionError(f"{row['priority_id']} has empty {field}")
    if "0.754" not in rows[0]["chapter2_locked_result"] or "JPN_36" not in rows[0]["focal_concepts"]:
        raise AssertionError("JPN36 phyllary discriminator is not priority 1")
    if "100/100" not in rows[1]["chapter2_locked_result"]:
        raise AssertionError("JPN06-JPN15 public sister contrast drift")
    if "at least 4" not in rows[2]["chapter2_locked_result"]:
        raise AssertionError("orientation lower-bound language drift")
    if "same voucher-linked RAD individuals" not in rows[3]["minimum_own_data"]:
        raise AssertionError("same-individual cross-module linkage gate drift")
    return rows


def validate_bounded_priors() -> list[dict[str, str]]:
    rows = read_rows(PRIOR_PATH)
    if [row["prior_id"] for row in rows] != [f"B{i:02d}" for i in range(1, 11)]:
        raise AssertionError("bounded-prior registry must contain ordered B01-B10")
    allowed = {
        "BOUNDED_PRIOR", "MEASUREMENT_REQUIREMENT", "DESIGN_CONSTRAINT",
        "FEASIBILITY_PRIOR", "STOP_RULE", "SAMPLING_REQUIREMENT", "REFERENCE_ONLY",
    }
    for row in rows:
        if row["admission_status"] not in allowed:
            raise AssertionError(f"unknown prior admission status: {row['prior_id']}")
        for field, value in row.items():
            if not value.strip():
                raise AssertionError(f"{row['prior_id']} has empty {field}")
    lookup = {row["prior_id"]: row for row in rows}
    if "RR=2.67364" not in lookup["B01"]["exact_result"]:
        raise AssertionError("herbivory quantitative prior drift")
    if "antagonist 2 mixed 2 pollinator 1" not in lookup["B02"]["exact_result"]:
        raise AssertionError("selection-mosaic prior drift")
    if "6/6" not in lookup["B08"]["exact_result"] or lookup["B08"]["admission_status"] != "STOP_RULE":
        raise AssertionError("generic-meta stop rule drift")
    if "0/64" not in lookup["B10"]["exact_result"] or lookup["B10"]["admission_status"] != "REFERENCE_ONLY":
        raise AssertionError("simulation boundary was promoted to a biological prior")
    return rows


def validate_protocol_registry() -> list[dict[str, str]]:
    rows = read_rows(PROTOCOL_PATH)
    if [row["protocol_id"] for row in rows] != ["F01", "F02"]:
        raise AssertionError("protocol registry membership drift")
    for row in rows:
        if row["field_execution_authorized"] != "false" or row["tissue_collection_authorized"] != "false":
            raise AssertionError(f"unauthorized protocol was opened: {row['protocol_id']}")
        if not row["stop_rule"].strip() or not row["claim_ceiling"].strip():
            raise AssertionError(f"protocol boundary missing: {row['protocol_id']}")
    if rows[0]["technical_state"] != "TECHNICALLY_READY_NOT_FIELD_AUTHORIZED":
        raise AssertionError("JPN36 readiness state drift")
    if rows[1]["technical_state"] != "DESIGN_CANDIDATE_NOT_PROTOCOL_READY":
        raise AssertionError("JPN15 design candidate was promoted to a protocol")
    return rows


def validate_empty_intake() -> None:
    with INTAKE_PATH.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != INTAKE_FIELDS:
            raise AssertionError("individual-intake schema drift")
        rows = list(reader)
    if rows:
        raise AssertionError("v1 handoff must remain empty until a separately authorized intake contract")
    forbidden = {"latitude", "longitude", "exact_locality", "permit_document"}
    if forbidden.intersection(INTAKE_FIELDS):
        raise AssertionError("sensitive or protected fields entered the public intake schema")


def validate_narrative() -> None:
    scope = SCOPE_PATH.read_text(encoding="utf-8")
    required_scope = [
        "minimum state-change lower bounds of four", "JPN36", "JPN06-JPN15",
        "RR=2.674", "0/64", "field_execution_authorized=false",
        "no definitive Japan-wide species-tree claim", "no adaptation or convergence claim",
    ]
    missing = [needle for needle in required_scope if needle not in scope]
    if missing:
        raise AssertionError(f"Chapter 3 scope missing: {missing}")
    readme = README_PATH.read_text(encoding="utf-8")
    required_readme = [SOURCE_MERGE_SHA, "at least four", "admitted: **0**", "authorized: **false**"]
    missing = [needle for needle in required_readme if needle not in readme]
    if missing:
        raise AssertionError(f"README fail-closed state missing: {missing}")


def main() -> int:
    paths = [CONTRACT_PATH, PRIORITY_PATH, PRIOR_PATH, PROTOCOL_PATH, INTAKE_PATH, SCOPE_PATH, README_PATH]
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing or empty Chapter 3 handoff file: {path.relative_to(ROOT)}")
    contract = validate_contract()
    priorities = validate_sampling_priorities()
    priors = validate_bounded_priors()
    protocols = validate_protocol_registry()
    validate_empty_intake()
    validate_narrative()
    print("chapter3_eazami_handoff_valid=true")
    print(f"source_merge_sha={contract['source_merge_sha']}")
    print(f"sampling_priorities={len(priorities)}")
    print(f"bounded_prior_rows={len(priors)}")
    print(f"protocol_rows={len(protocols)}")
    print("own_biological_records_admitted=0")
    print("field_execution_authorized=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
