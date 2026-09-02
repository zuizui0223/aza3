#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import validate_radseq_physical_template_readiness_v1 as target


class RadPhysicalTemplateReadinessTests(unittest.TestCase):
    def test_inventory_has_16_fail_closed_templates(self) -> None:
        rows = target.validate_inventory()
        self.assertEqual(len(rows), 16)
        self.assertEqual(sum(r["template_role"] == "TEMPLATE_1" for r in rows), 8)
        self.assertEqual(sum(r["template_role"] == "TEMPLATE_2" for r in rows), 8)
        self.assertTrue(all(r["physical_material_status"] == "UNKNOWN_NOT_VERIFIED" for r in rows))

    def test_no_physical_possession_is_invented(self) -> None:
        rows = target.validate_inventory()
        self.assertTrue(all(not r["physical_sample_id"] for r in rows))
        self.assertTrue(all(r["primary_blocker"] == "PHYSICAL_SAMPLE_NOT_CONFIRMED" for r in rows))
        self.assertTrue(all(r["stage_a_eligible"] == "false" for r in rows))
        self.assertTrue(all(r["stage_b_eligible"] == "false" for r in rows))

    def test_stage_a_requires_one_template_per_anchor(self) -> None:
        rows = target.validate_inventory()
        stage_a = [r for r in rows if r["template_role"] == "TEMPLATE_1"]
        self.assertEqual({r["anchor_id"] for r in stage_a}, {f"RAD_A{i:02d}" for i in range(1, 9)})
        self.assertTrue(all(r["stage_requirement"] == "STAGE_A_REQUIRED" for r in stage_a))

    def test_narrative_preserves_public_vs_physical_boundary(self) -> None:
        target.validate_narrative()


if __name__ == "__main__":
    unittest.main()
