#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import validate_radseq_end_to_end_v1 as target


class RadseqEndToEndTests(unittest.TestCase):
    def test_anchor_panel_is_16_templates_across_eight_systems(self) -> None:
        rows = target.validate_anchor_ledger()
        self.assertEqual(len(rows), 8)
        self.assertEqual(sum(int(r["biological_templates"]) for r in rows), 16)
        self.assertTrue(all(r["counts_toward_core_biological_n"] == "false" for r in rows))

    def test_library_intake_is_empty_before_production(self) -> None:
        target.validate_empty_library_intake()

    def test_rad_products_keep_bounded_roles(self) -> None:
        d = target.validate_contract()
        products = d["primary_products"]
        self.assertIn("secondary strict-shared-locus", products["R2_all_japan_sensitivity"])
        self.assertIn("no pooled diploid matrix", products["R3_mixed_ploidy_overlay"])
        self.assertIn("not final E3 evidence", products["R4_m01_selection_support"])

    def test_pilot_precedes_production(self) -> None:
        d = target.validate_contract()
        self.assertTrue(d["pilot"]["required"])
        self.assertIn("at least three candidate", d["pilot"]["stage_A_enzyme_screen"]["rule"])
        self.assertEqual(d["pilot"]["stage_B_reproducibility_screen"]["technical_repeat_minimum"], 8)
        self.assertIn("prohibit concatenating different RAD protocols", d["pilot"]["stratification_fallback"])

    def test_stacks_parameters_are_outcome_blind(self) -> None:
        d = target.validate_contract()
        stacks = d["bioinformatics"]["stacks_parameter_freeze"]
        self.assertEqual(stacks["m_grid"], [2, 3, 4, 5, 6])
        self.assertEqual(stacks["M_grid"], [1, 2, 3, 4, 5, 6])
        self.assertIn("r80-style", stacks["optimization_basis"])
        self.assertIn("before inspecting P01-P05", stacks["freeze_timing"])

    def test_mixed_ploidy_is_not_forced_diploid(self) -> None:
        d = target.validate_contract()
        mixed = d["analysis_matrices"]["mixed_ploidy_matrices"]
        self.assertIn("polyRAD", mixed["polyploid"])
        self.assertIn("no single diploid genotype caller", mixed["cross_ploidy"])

    def test_m01_selection_uses_rad_as_background_not_final_proof(self) -> None:
        d = target.validate_contract()
        m01 = d["m01_selection_boundary"]
        self.assertIn("two-species FST outlier", m01["rad_not_sufficient"])
        self.assertIn("non-RAD confirmation route", m01["promotion_rule"])
        self.assertIn("restriction-site dropout", m01["dropout_audit"])

    def test_rad_starts_fail_closed(self) -> None:
        d = target.validate_contract()
        self.assertTrue(all(value is False for value in d["current_state"].values()))

    def test_rad_narrative_is_linked(self) -> None:
        target.validate_narrative()


if __name__ == "__main__":
    unittest.main()
