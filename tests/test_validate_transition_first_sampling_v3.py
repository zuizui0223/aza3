#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import validate_transition_first_sampling_v3 as target


class TransitionFirstV3Tests(unittest.TestCase):
    def test_japan38_is_not_sampling_universe(self) -> None:
        d = target.validate_contract()
        u = d["candidate_universe"]
        self.assertEqual(u["source_record_count_observed_2026_09_02"], 161)
        self.assertEqual(u["moreyra_japan_subset"], 38)
        self.assertEqual(u["moreyra_subset_role"], "HYPOTHESIS_ORIGIN_SUBSET_NOT_SAMPLING_UNIVERSE")

    def test_focal_population_replication_is_population_based(self) -> None:
        d = target.validate_contract()
        self.assertEqual(d["A1_fixed_focal_population_layer"]["minimum_launch_total"]["primary_rad"], 132)
        self.assertEqual(d["A1_fixed_focal_population_layer"]["target_total"]["primary_rad"], 168)

    def test_p03_requires_nonoverlapping_neighbourhoods(self) -> None:
        d = target.validate_contract()
        p03 = d["A3_transition_neighbourhoods"]["P03"]
        self.assertEqual(p03["minimum_neighbourhoods"], 4)
        self.assertTrue(p03["must_be_nonoverlapping"])
        self.assertIn("JPN36 can count at most once", p03["requirements"])

    def test_nonbackbone_taxa_require_nuclear_placement_before_rad(self) -> None:
        d = target.validate_contract()
        a2 = d["A2_backbone_augmentation"]
        self.assertFalse(a2["rad_before_stable_placement"])
        self.assertEqual(a2["placement_representatives_per_taxon"], 2)

    def test_focal_assay_pilot_is_nested(self) -> None:
        d = target.validate_contract()
        p = d["A1_rad_assay_pilot"]
        self.assertEqual(p["stage_A"]["libraries"], 15)
        self.assertEqual(p["stage_B"]["total_libraries"], 25)

    def test_narrative_and_builder_are_linked(self) -> None:
        target.validate_narrative()


if __name__ == "__main__":
    unittest.main()
