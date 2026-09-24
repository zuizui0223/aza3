#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import validate_zero_baseline_sampling_v2 as target


class ZeroBaselineSamplingV2Tests(unittest.TestCase):
    def test_zero_material_start(self) -> None:
        d = target.validate_contract()
        self.assertEqual(d["starting_material"], "ZERO_PHYSICAL_SAMPLES")

    def test_s1_is_population_derived_not_38x3(self) -> None:
        d = target.validate_contract()
        self.assertEqual(d["S1_focal_launch"]["derived_minimum"]["physical_plants"], 165)
        self.assertEqual(d["S1_focal_launch"]["derived_minimum"]["initial_rad_individuals"], 132)
        self.assertEqual(d["S4_optional_all_japan_rad"]["priority"], "LAST_OPTIONAL")

    def test_focal_pilot_uses_required_systems_only(self) -> None:
        d = target.validate_contract()
        pilot = d["S1_rad_pilot_nested"]
        self.assertEqual(len(pilot["systems"]), 5)
        self.assertEqual(pilot["stage_A"]["libraries"], 15)
        self.assertEqual(pilot["stage_B"]["libraries_total"], 25)
        self.assertIn("does not qualify high-ploidy/polyploid", pilot["polyploid_boundary"])

    def test_s2_is_uncertainty_cover_design(self) -> None:
        d = target.validate_contract()
        self.assertTrue(any(">=80%" in x for x in d["S2_transition_neighbourhood"]["nomination_algorithm"]))
        self.assertFalse(d["S2_transition_neighbourhood"]["concept_count_is_fixed_in_advance"])

    def test_m01_full_is_gated(self) -> None:
        d = target.validate_contract()
        self.assertEqual(d["S3_m01_expansion"]["added_physical"], 60)
        self.assertEqual(d["S3_m01_expansion"]["final_m01_initial_rad"], 96)
        self.assertIn("non-RAD confirmation assay", d["S3_m01_expansion"]["selection_boundary"])

    def test_population_target_rows_match_contract(self) -> None:
        rows = target.validate_targets()
        self.assertGreaterEqual(len(rows), 18)

    def test_readme_and_doc_are_v2(self) -> None:
        target.validate_narrative()


if __name__ == "__main__":
    unittest.main()
