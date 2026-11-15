#!/usr/bin/env python3
"""Synthetic, public-safe unit tests for the RAD material gate."""
from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("audit_radseq_pilot_intake_v1.py")
spec = importlib.util.spec_from_file_location("radseq_gate", SCRIPT)
gate = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(gate)


def write_rows(root: Path, rel: Path, rows: list[dict[str, str]]) -> None:
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted(gate.REQUIRED[rel])
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def setup_intakes(root: Path, nature=None, chapter3=None, handoff=None):
    write_rows(root, gate.NATURE, nature or [])
    write_rows(root, gate.CHAPTER3, chapter3 or [])
    write_rows(root, gate.HANDOFF, handoff or [])


def qualified_handoff(num: int, **updates) -> dict[str, str]:
    row = {
        "handoff_id": f"H{num:03d}",
        "individual_id": f"I{num:03d}",
        "voucher_id": f"V{num:03d}",
        "material_type": "dna_extract",
        "quantity": "1 tube",
        "preservation_type": "frozen",
        "transfer_status": "received",
        "material_use_permission": "approved",
        "sequencing_permission": "approved",
        "conservation_review": "approved",
        "dna_qc_status": "pass",
        "admission_status": "admitted",
    }
    row.update(updates)
    return row


class MaterialGateTests(unittest.TestCase):
    def test_empty_is_hold_not_authorization(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            setup_intakes(root)
            result = gate.audit(root)
            self.assertEqual(result["received_permissioned_qc_pass_handoff_individuals_in_public_ledger"], 0)
            self.assertEqual(result["source_row_counts"]["aza3_nature_individual_intake"], 0)
            self.assertFalse(result["pilot_sequencing_authorized_by_this_audit"])
            self.assertTrue(result["gate"].startswith("HOLD_"))

    def test_permission_missing_does_not_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            setup_intakes(root, handoff=[qualified_handoff(1, sequencing_permission="pending")])
            self.assertEqual(
                gate.audit(root)["received_permissioned_qc_pass_handoff_individuals_in_public_ledger"],
                0,
            )

    def test_nondna_material_does_not_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            setup_intakes(root, handoff=[qualified_handoff(1, material_type="voucher")])
            self.assertEqual(
                gate.audit(root)["received_permissioned_qc_pass_handoff_individuals_in_public_ledger"],
                0,
            )

    def test_documentary_24_is_still_not_automatic_go(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = [qualified_handoff(i) for i in range(24)]
            setup_intakes(root, handoff=rows)
            result = gate.audit(root)
            self.assertEqual(result["received_permissioned_qc_pass_handoff_individuals_in_public_ledger"], 24)
            self.assertIn("REVIEW_", result["gate"])
            self.assertFalse(result["pilot_sequencing_authorized_by_this_audit"])

    def test_duplicate_same_intake_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            setup_intakes(root, nature=[
                {"individual_id": "I001"}, {"individual_id": "I001"}
            ])
            with self.assertRaises(ValueError):
                gate.audit(root)


if __name__ == "__main__":
    unittest.main()
