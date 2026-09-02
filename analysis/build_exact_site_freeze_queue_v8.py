#!/usr/bin/env python3
"""Build the Chapter 3 one-slot-per-sample field-site freeze queue.

The public repository stores range sectors and deidentified locality IDs only.
Exact coordinates and sensitive locality descriptions are deliberately prohibited.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLANNING = ROOT / "data" / "planning"
OUT = PLANNING / "chapter3_exact_site_freeze_queue_v8.csv"

SECTOR_FILES = {
    1: PLANNING / "chapter3_holefill_priority1_sample_sectors_v7.csv",
    2: PLANNING / "chapter3_holefill_priority2_sample_sectors_v7.csv",
    3: PLANNING / "chapter3_holefill_priority3_sample_sectors_v7.csv",
    4: PLANNING / "chapter3_holefill_priority4_sample_sectors_v7.csv",
    5: PLANNING / "chapter3_holefill_priority5_sample_sectors_v7.csv",
    6: PLANNING / "chapter3_holefill_priority6_sample_sectors_v7.csv",
    7: PLANNING / "chapter3_holefill_priority7_sample_sectors_v7.csv",
}

FIELDS = [
    "slot_id",
    "priority",
    "species_binomial",
    "japanese_names",
    "sample_slot",
    "required_for_tree",
    "range_sector",
    "nmns_distribution_source",
    "sector_design",
    "public_locality_audit_required",
    "public_locality_audit_status",
    "current_occurrence_status",
    "current_occurrence_evidence_id",
    "current_occurrence_checked_date",
    "conservation_review_status",
    "land_manager_status",
    "collection_permission_status",
    "tissue_permission_status",
    "target_locality_id",
    "private_exact_site_record_status",
    "field_window_status",
    "collaborator_or_collector",
    "freeze_status",
    "exclusion_or_block_reason",
]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def safe_id(species: str) -> str:
    return species.replace(" ", "_").replace(".", "")


def base_row(priority: int, row: dict[str, str], slot: str, required: bool, sector: str, design: str,
             public_audit_required: bool = False) -> dict[str, str]:
    return {
        "slot_id": f"P{priority}_{safe_id(row['species_binomial'])}_{slot}",
        "priority": f"P{priority}",
        "species_binomial": row["species_binomial"],
        "japanese_names": row.get("japanese_names", ""),
        "sample_slot": slot,
        "required_for_tree": "true" if required else "false",
        "range_sector": sector,
        "nmns_distribution_source": row.get("nmns_distribution_source", ""),
        "sector_design": design,
        "public_locality_audit_required": "true" if public_audit_required else "false",
        "public_locality_audit_status": "NOT_REQUIRED" if not public_audit_required else "NOT_CHECKED",
        "current_occurrence_status": "NOT_CHECKED",
        "current_occurrence_evidence_id": "",
        "current_occurrence_checked_date": "",
        "conservation_review_status": "NOT_CHECKED",
        "land_manager_status": "NOT_CHECKED",
        "collection_permission_status": "NOT_CHECKED",
        "tissue_permission_status": "NOT_CHECKED",
        "target_locality_id": "",
        "private_exact_site_record_status": "NOT_CREATED",
        "field_window_status": "NOT_CHECKED",
        "collaborator_or_collector": "",
        "freeze_status": "NOT_FROZEN",
        "exclusion_or_block_reason": "",
    }


def build() -> list[dict[str, str]]:
    out: list[dict[str, str]] = []

    # P1-P5: two required new representatives, A and B.
    for priority in range(1, 6):
        for row in read_rows(SECTOR_FILES[priority]):
            out.append(base_row(priority, row, "A", True, row["sample_A_sector"], row["sector_design"]))
            out.append(base_row(priority, row, "B", True, row["sample_B_sector"], row["sector_design"]))

    # P6: one required own complement; final range sector cannot be frozen until
    # the existing public accession locality is audited.
    for row in read_rows(SECTOR_FILES[6]):
        out.append(base_row(
            6,
            row,
            "OWN",
            True,
            row["own_sample_sector_rule"],
            row["sector_design"],
            public_audit_required=True,
        ))

    # P7: tree is already filled; one own phenotype-linked sample is optional.
    for row in read_rows(SECTOR_FILES[7]):
        out.append(base_row(
            7,
            row,
            "OPTIONAL_OWN",
            False,
            row["optional_own_trait_link_sector_rule"],
            "OPTIONAL_COMPLEMENT_AFTER_PUBLIC_LOCALITY_AUDIT",
            public_audit_required=True,
        ))

    return out


def main() -> int:
    rows = build()
    if len(rows) != 229:
        raise AssertionError(f"expected 229 execution slots (228 required +1 optional), got {len(rows)}")
    required = sum(r["required_for_tree"] == "true" for r in rows)
    if required != 228:
        raise AssertionError(f"expected 228 required tree slots, got {required}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print("chapter3_exact_site_freeze_queue_v8_built=true")
    print(f"queue_rows={len(rows)}")
    print(f"required_tree_slots={required}")
    print("optional_trait_link_slots=1")
    print("exact_coordinates_in_public_repo=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
