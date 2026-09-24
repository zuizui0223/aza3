#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import validate_inverse_sampling_design_v1 as target


class InverseSamplingDesignTests(unittest.TestCase):
    def test_core_floor_and_focal_topups(self) -> None:
        d = target.validate_design()
        core = d["core_breadth_panel"]
        self.assertEqual(core["all_japan_floor"]["total_individuals"], 114)
        self.assertEqual(core["derived_totals"]["minimum_primary_individuals"], 167)
        self.assertEqual(core["derived_totals"]["recommended_primary_individuals"], 193)

    def test_p05_adds_no_fresh_quota(self) -> None:
        d = target.validate_design()
        self.assertIn("no additional fresh-collection quota", d["core_breadth_panel"]["p05_rule"])

    def test_m01_discovery_uses_two_populations_per_lineage(self) -> None:
        d = target.validate_design()
        roles = d["m01_depth_panel"]["stage_m1_history_and_expression_discovery"]["focal_population_roles"]
        self.assertEqual(len(roles), 2)
        self.assertTrue(all("role_1" in row and "role_2" in row for row in roles))

    def test_rna_and_outgroup_targets(self) -> None:
        d = target.validate_design()
        m1 = d["m01_depth_panel"]["stage_m1_history_and_expression_discovery"]
        self.assertEqual(m1["rna_subset"]["planned_libraries_at_one_frozen_developmental_stage"], 20)
        self.assertEqual(m1["ancestral_state_outgroup_panel"]["minimum_nonfocal_individuals"], 18)
        self.assertIn("per observed colour state", m1["ancestral_state_outgroup_panel"]["polymorphism_rule"])

    def test_e3_is_120_and_gated(self) -> None:
        d = target.validate_design()
        e3 = d["m01_depth_panel"]["stage_m2_population_genomic_selection"]
        self.assertEqual(e3["primary_individuals_total"], 120)
        self.assertIn("RAD-only cross-species FST is not sufficient", e3["sequencing_platform_gate"])
        self.assertFalse(d["current_state"]["e3_population_genomic_sequencing_open"])

    def test_narrative_is_linked(self) -> None:
        target.validate_narrative()


if __name__ == "__main__":
    unittest.main()
