#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import validate_chapter3_handoff_v1 as target


class Chapter3HandoffTests(unittest.TestCase):
    def test_source_and_chapter2_results_are_frozen(self) -> None:
        contract = target.validate_contract()
        self.assertEqual(contract["source_merge_sha"], target.SOURCE_MERGE_SHA)
        self.assertEqual(
            contract["chapter2_locked_results"]["orientation_minimum_changes"]["lower_bound"],
            4,
        )

    def test_five_sampling_priorities_are_ranked(self) -> None:
        rows = target.validate_sampling_priorities()
        self.assertEqual(len(rows), 5)
        self.assertIn("JPN_36", rows[0]["focal_concepts"])
        self.assertIn("100/100", rows[1]["chapter2_locked_result"])

    def test_meta_and_simulation_results_are_bounded(self) -> None:
        rows = target.validate_bounded_priors()
        lookup = {row["prior_id"]: row for row in rows}
        self.assertEqual(len(rows), 14)
        self.assertEqual(lookup["B10"]["admission_status"], "REFERENCE_ONLY")
        self.assertIn("not biological", lookup["B10"]["claim_ceiling"])
        self.assertEqual(lookup["B12"]["admission_status"], "REFERENCE_ONLY")
        self.assertEqual(lookup["B14"]["admission_status"], "ROUTING_ONLY")

    def test_protocols_are_not_authorized(self) -> None:
        rows = target.validate_protocol_registry()
        self.assertTrue(all(row["field_execution_authorized"] == "false" for row in rows))
        self.assertIn("NOT_FIELD_AUTHORIZED", rows[0]["technical_state"])

    def test_intake_is_empty_and_narrative_is_bounded(self) -> None:
        target.validate_empty_intake()
        target.validate_narrative()


if __name__ == "__main__":
    unittest.main()
