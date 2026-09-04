from __future__ import annotations

import unittest

from analysis.validate_acsp_prospective_validation_v1 import validate


class AcspProspectiveValidationMirrorTests(unittest.TestCase):
    def test_frozen_cohort_and_prefill_validate(self) -> None:
        summary = validate()
        self.assertEqual(summary["status"], "OK")
        self.assertEqual(summary["cohort_size"], 13)
        self.assertEqual(summary["regime_counts"], {"LOCAL_CONTINUATION": 9, "SENTINEL": 4})
        self.assertEqual(
            summary["anchor_replication_counts"],
            {
                "MULTIPLE_PRIMARY_ANCHORS": 5,
                "ZERO_PRIMARY_ANCHOR": 4,
                "SINGLE_PRIMARY_ANCHOR": 4,
            },
        )
        self.assertEqual(
            summary["method_arm_counts"],
            {"STRUCTURAL_LOCAL": 8, "STRUCTURAL_SENTINEL": 4, "SPATIAL_BASELINE_ONLY": 1},
        )
        self.assertFalse(summary["candidate_patches_built"])
        self.assertFalse(summary["field_outcomes_opened"])


if __name__ == "__main__":
    unittest.main()
