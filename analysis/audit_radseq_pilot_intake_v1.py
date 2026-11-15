#!/usr/bin/env python3
"""Aggregate, fail-closed RAD pilot material audit. Never reports authorization."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NATURE = Path("data/intake/aza3_nature_individual_intake_v1.csv")
CHAPTER3 = Path("data/intake/chapter3_individual_intake_v1.csv")
HANDOFF = Path("data/intake/grow_ito_sample_handoff_v1.csv")

REQUIRED = {
    NATURE: {
        "individual_id", "taxon_concept", "population_id", "voucher_id",
        "leaf_dna_tissue_id", "orientation_deg_gravity", "phyllary_posture",
        "stickiness_state", "access_authorization_id",
        "collection_authorization_id", "conservation_review_id",
        "cytotype_x", "flow_cytometry_record_id", "admission_status",
    },
    CHAPTER3: {
        "individual_id", "taxon_concept", "voucher_id", "rad_tissue_id",
        "access_authorization_id", "collection_authorization_id",
        "conservation_review_id", "admission_status",
    },
    HANDOFF: {
        "handoff_id", "individual_id", "voucher_id", "material_type",
        "transfer_status", "material_use_permission", "sequencing_permission",
        "conservation_review", "dna_qc_status", "linked_head_image_id",
        "quantity", "preservation_type",
        "admission_status",
    },
}

POSITIVE = {"yes", "approved", "allowed", "granted", "true"}
RECEIVED = {"received", "confirmed_received"}
DNA_PASS = {"pass", "passed", "qc_pass"}
ADMITTED = {"admitted", "approved"}


def read_csv(root: Path, relpath: Path) -> list[dict[str, str]]:
    path = root / relpath
    if not path.is_file():
        raise ValueError(f"required intake file missing: {relpath}")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        missing = REQUIRED[relpath] - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{relpath}: missing required columns {sorted(missing)}")
        return list(reader)


def val(row: dict[str, str], key: str) -> str:
    return (row.get(key) or "").strip()


def norm(row: dict[str, str], key: str) -> str:
    return val(row, key).lower()


def handoff_qc_pass(row: dict[str, str]) -> bool:
    """Only a documented handoff; not a physical/lab authorization decision."""
    fields = ("individual_id", "voucher_id", "handoff_id")
    return (
        all(val(row, k) for k in fields)
        and norm(row, "material_type") in {
            "dna", "dna_extract", "leaf", "leaf_tissue", "silica_leaf",
            "silica-dried leaf", "frozen_leaf",
        }
        and val(row, "quantity")
        and val(row, "preservation_type")
        and norm(row, "transfer_status") in RECEIVED
        and all(norm(row, k) in POSITIVE for k in (
            "material_use_permission", "sequencing_permission", "conservation_review"
        ))
        and norm(row, "dna_qc_status") in DNA_PASS
        and norm(row, "admission_status") in ADMITTED
    )


def audit(root: Path) -> dict:
    nature = read_csv(root, NATURE)
    chapter3 = read_csv(root, CHAPTER3)
    handoffs = read_csv(root, HANDOFF)

    for name, records in (("nature", nature), ("chapter3", chapter3)):
        ids = [val(row, "individual_id") for row in records if val(row, "individual_id")]
        if len(ids) != len(set(ids)):
            raise ValueError(f"duplicate individual_id within {name} intake")

    tissue_ids = {
        val(row, "individual_id")
        for row in nature
        if val(row, "individual_id") and val(row, "leaf_dna_tissue_id")
    }
    tissue_ids |= {
        val(row, "individual_id")
        for row in chapter3
        if val(row, "individual_id") and val(row, "rad_tissue_id")
    }

    complete_handoffs = [row for row in handoffs if handoff_qc_pass(row)]
    handoff_ids = {
        val(row, "individual_id") for row in complete_handoffs
    }
    phenotype_ids = {
        val(row, "individual_id")
        for row in nature
        if val(row, "individual_id")
        and (
            val(row, "orientation_deg_gravity")
            or val(row, "phyllary_posture")
            or val(row, "stickiness_state")
        )
    }
    linked_ids = handoff_ids & phenotype_ids

    n = len(handoff_ids)
    if n == 0:
        gate = "HOLD_NO_DOCUMENTED_RECEIVED_PERMISSIONED_QC_PASS_HANDOFF"
    elif n < 24:
        gate = "HOLD_BELOW_PROVISIONAL_24_SAMPLE_PILOT_FLOOR"
    else:
        gate = "REVIEW_PILOT_COHORT_AND_PRIVATE_PERMITS_NOT_AUTO_APPROVED"

    return {
        "contract_version": "radseq_pilot_material_preflight_v1",
        "source_row_counts": {
            "aza3_nature_individual_intake": len(nature),
            "aza3_chapter3_individual_intake": len(chapter3),
            "grow_ito_sample_handoff": len(handoffs),
        },
        "individuals_with_repository_recorded_rad_or_dna_tissue_id": len(tissue_ids),
        "received_permissioned_qc_pass_handoff_individuals_in_public_ledger": n,
        "handoff_individuals_with_same_id_focal_trait_record": len(linked_ids),
        "provisional_pilot_sample_window": [24, 48],
        "gate": gate,
        "pilot_sequencing_authorized_by_this_audit": False,
        "physical_inventory_exhaustive": False,
        "sensitive_localities_exposed": False,
        "note": (
            "Only public ledger documentation was checked. Zero records does not "
            "prove zero physical specimens; even positive counts require private "
            "chain-of-custody, permissions, DNA and cytotype review."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument(
        "--require-documentary-pilot-floor", action="store_true",
        help="Exit 2 when fewer than 24 handoff individuals pass the documented gate.",
    )
    args = parser.parse_args()
    report = audit(args.root)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if args.require_documentary_pilot_floor and (
        report["received_permissioned_qc_pass_handoff_individuals_in_public_ledger"] < 24
    ):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
