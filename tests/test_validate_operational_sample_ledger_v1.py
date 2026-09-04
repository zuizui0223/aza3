#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import validate_operational_sample_ledger_v1 as target


class OperationalSampleLedgerTests(unittest.TestCase):
    def test_core_has_38_ordered_concepts_and_frozen_totals(self) -> None:
        rows = target.validate_core()
        self.assertEqual(len(rows), 38)
        self.assertEqual(rows[0]["member_id"], "JPN_01")
        self.assertEqual(rows[-1]["member_id"], "JPN_38")
        self.assertEqual(sum(int(r["core_min_primary"]) for r in rows), 167)
        self.assertEqual(sum(int(r["core_recommended_primary"]) for r in rows), 193)

    def test_identity_blocks_do_not_allow_replacement(self) -> None:
        rows = {row["member_id"]: row for row in target.validate_core()}
        for member in ("JPN_29", "JPN_31", "JPN_33"):
            self.assertIn("BLOCKED", rows[member]["identity_gate"])
            self.assertEqual(rows[member]["replacement_rule"], "NO_CONVENIENCE_REPLACEMENT")

    def test_wild_japan_repairs_are_explicit(self) -> None:
        rows = {row["member_id"]: row for row in target.validate_core()}
        for member in ("JPN_32", "JPN_34", "JPN_35", "JPN_36", "JPN_37", "JPN_38"):
            self.assertIn("WILD_JAPAN_RESAMPLE", rows[member]["identity_gate"])

    def test_m01_discovery_is_four_populations_with_nested_assays(self) -> None:
        rows = {row["sample_unit_id"]: row for row in target.validate_m01()}
        for sample_id in ("M01_BREV_OKI", "M01_BREV_AMAMI", "M01_IRUM_MIYAKO", "M01_IRUM_ISHIGAKI"):
            self.assertEqual(int(rows[sample_id]["primary_individuals"]), 15)
            self.assertEqual(int(rows[sample_id]["pigment_primary_n"]), 3)
            self.assertEqual(int(rows[sample_id]["rna_sequence_primary_n"]), 5)

    def test_m01_e1_bracket_is_six_concepts_and_18_individuals(self) -> None:
        rows = target.validate_m01()
        bracket = [r for r in rows if r["population_role"] in {"E1_SISTER_CLADE_BRACKET", "E1_STATE_BALANCED_DEEP_BRACKET"}]
        self.assertEqual(len(bracket), 6)
        self.assertEqual(sum(int(r["primary_individuals"]) for r in bracket), 18)

    def test_m01_final_e3_is_120_and_discovery_is_nested(self) -> None:
        rows = target.validate_m01()
        focal = [r for r in rows if r["taxon_or_role"] in {"Cirsium brevicaule", "Cirsium irumtiense"}]
        self.assertEqual(sum(int(r["primary_individuals"]) for r in focal), 120)
        discovery = [r for r in focal if r["population_role"] == "DISCOVERY_AND_E3"]
        self.assertEqual(sum(int(r["primary_individuals"]) for r in discovery), 60)

    def test_narrative_routes_operational_constraints(self) -> None:
        target.validate_narrative()


if __name__ == "__main__":
    unittest.main()
