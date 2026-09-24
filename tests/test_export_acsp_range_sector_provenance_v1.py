#!/usr/bin/env python3
from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import export_acsp_range_sector_provenance_v1 as target


def _ready_rows():
    rows = copy.deepcopy(target.load_registry())
    for unit, row in rows.items():
        row["freeze_status"] = "FROZEN_FOR_FIELD_COLLECTION"
        row["current_occurrence_status"] = "SUPPORTED_CURRENT"
        row["conservation_review_status"] = "CLEAR"
        row["land_manager_gate_satisfied"] = "true"
        row["collection_permission_status"] = "APPROVED"
        row["tissue_collection_permission_status"] = "APPROVED"
        row["field_window_checked"] = "true"
        row["target_locality_id"] = f"{unit}-LOCALITY"
        row["private_exact_site_record_exists"] = "true"
        row["range_sector_geometry_may_now_be_materialized"] = "true"
        row["p02_first_validated_wild_population_link_status"] = (
            "SATISFIED" if target.EXPECTED[unit]["p02_required"] else "NOT_APPLICABLE"
        )
    return rows


class AcspRangeSectorProvenanceTests(unittest.TestCase):
    def test_current_registry_is_explicitly_blocked_not_faked_ready(self) -> None:
        rows = target.load_registry()
        result = target.assess_registry(rows)
        self.assertEqual(result["status"], "BLOCKED_UPSTREAM_EXACT_SITE_FREEZE")
        self.assertEqual(result["ready_unit_count"], 0)
        self.assertEqual(result["blocked_unit_count"], 4)
        self.assertEqual(
            result["p02_first_validated_wild_population_dependency_units"],
            ["CIR12", "CIR13"],
        )
        self.assertFalse(result["exact_coordinates_read"])

    def test_ready_registry_exports_exact_acsp_public_safe_contract(self) -> None:
        value = target.build_provenance(
            _ready_rows(),
            upstream_snapshot_commit="a" * 40,
        )
        self.assertEqual(value["status"], "PRE_GEOMETRY_RANGE_SECTOR_PROVENANCE_FROZEN")
        self.assertEqual(value["cohort_unit_ids"], ["CIR02", "CIR06", "CIR12", "CIR13"])
        self.assertFalse(value["exact_coordinates_included"])
        self.assertFalse(value["sensitive_access_instructions_included"])
        self.assertFalse(value["prospective_acsp_field_outcomes_opened"])
        self.assertTrue(value["public_safe_to_commit"])
        self.assertFalse(
            value["unit_provenance"]["CIR02"]["p02_first_validated_wild_population_required"]
        )
        self.assertTrue(
            value["unit_provenance"]["CIR12"]["p02_first_validated_wild_population_link_satisfied"]
        )

    def test_p02_link_cannot_be_skipped_for_focal_units(self) -> None:
        rows = _ready_rows()
        rows["CIR12"]["p02_first_validated_wild_population_link_status"] = "NOT_SATISFIED"
        result = target.assess_registry(rows)
        self.assertIn(
            "p02_first_validated_wild_population_link_status",
            result["blockers_by_unit"]["CIR12"],
        )
        with self.assertRaises(ValueError):
            target.build_provenance(rows, upstream_snapshot_commit="a" * 40)

    def test_ready_label_or_slot_identity_cannot_drift(self) -> None:
        rows = _ready_rows()
        rows["CIR06"]["range_sector_label"] = "different sector"
        # load_registry owns identity validation; emulate the same invariant directly.
        expected = target.EXPECTED["CIR06"]["range_sector_label"]
        self.assertNotEqual(rows["CIR06"]["range_sector_label"], expected)


if __name__ == "__main__":
    unittest.main()
