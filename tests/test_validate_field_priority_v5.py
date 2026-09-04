#!/usr/bin/env python3
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

class FieldPriorityV5Tests(unittest.TestCase):
    def setUp(self):
        self.d = json.loads((ROOT / "data/planning/chapter3_field_priority_v5.json").read_text(encoding="utf-8"))

    def test_tree_breadth_precedes_rad_depth(self):
        self.assertEqual(self.d["primary_field_objective"], "MAXIMIZE_NATIONWIDE_SPECIES_TREE_TAXON_BREADTH_BEFORE_POPULATION_DEPTH")

    def test_rad_is_staged(self):
        self.assertEqual(self.d["focal_rad_staging"]["starter"]["primary_rad_total"], 120)
        self.assertEqual(self.d["focal_rad_staging"]["conditional_completion"]["full_primary_rad"], 168)

    def test_campaign_order(self):
        self.assertEqual([x["id"] for x in self.d["campaigns"]], ["C1","C2","C3","C4","C5","C6"])

    def test_fail_closed(self):
        self.assertEqual(self.d["current_state"]["physical_samples"], 0)
        self.assertFalse(self.d["current_state"]["field_sampling_authorized"])

if __name__ == "__main__":
    unittest.main()
