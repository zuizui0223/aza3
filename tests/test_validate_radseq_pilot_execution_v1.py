#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import validate_radseq_pilot_execution_v1 as target


class RadseqPilotExecutionTests(unittest.TestCase):
    def test_three_candidates_are_unselected(self) -> None:
        d = target.validate_contract()
        self.assertEqual(len(d["candidate_protocols"]), 3)
        self.assertTrue(all(x["status"] == "SCREEN_CANDIDATE_NOT_WINNER" for x in d["candidate_protocols"]))

    def test_stage_counts_are_24_and_40(self) -> None:
        d = target.validate_contract()
        self.assertEqual(d["stage_A_screen"]["planned_libraries"], 24)
        self.assertEqual(d["stage_B_reproducibility"]["planned_libraries_total"], 40)

    def test_reproducibility_gates_are_frozen(self) -> None:
        d = target.validate_contract()
        g = d["stage_B_reproducibility"]["technical_concordance_gate"]
        self.assertEqual(g["minimum_overlap_genotype_concordance"], 0.95)
        self.assertEqual(g["minimum_core_locus_recovery"], 0.90)

    def test_r80_and_cross_concept_occupancy_are_separate(self) -> None:
        d = target.validate_contract()
        s = d["stacks_pilot_rule"]
        self.assertIn("within focal or population-like subsets", s["within_population_parameter_optimization"])
        self.assertIn("6/8, 7/8, 8/8", s["cross_concept_gate"])

    def test_library_allocation_is_balanced(self) -> None:
        rows = target.validate_allocation()
        self.assertEqual(len([r for r in rows if r["pilot_stage"] == "STAGE_A"]), 24)
        self.assertEqual(len([r for r in rows if r["pilot_stage"] == "STAGE_B"]), 32)
        repeats = [r for r in rows if r["pilot_stage"] == "STAGE_B_REPEAT"]
        self.assertEqual(len(repeats), 8)
        self.assertEqual(sum(r["protocol_slot"] == "TOP1" for r in repeats), 4)
        self.assertEqual(sum(r["protocol_slot"] == "TOP2" for r in repeats), 4)

    def test_fail_closed_state(self) -> None:
        d = target.validate_contract()
        self.assertTrue(all(v is False for v in d["current_state"].values()))

    def test_narrative(self) -> None:
        target.validate_narrative()


if __name__ == "__main__":
    unittest.main()
