#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import validate_mechanistic_demonstration_v1 as target


class MechanisticDemonstrationTests(unittest.TestCase):
    def test_core_priorities_remain_p01_to_p05(self) -> None:
        rows = target.validate_core_priority_independence()
        self.assertEqual([row["priority_id"] for row in rows], ["P01", "P02", "P03", "P04", "P05"])

    def test_m01_is_embedded_not_p06(self) -> None:
        m01 = target.validate_m01_contract()
        self.assertEqual(m01["demonstration_id"], "M01")
        self.assertEqual(m01["role"]["architecture"], "BREADTH_TO_DEPTH")
        self.assertEqual(m01["role"]["core_priority_status"], "EMBEDDED_CASE_NOT_P06")
        self.assertIn("same Chapter 3 identifiability problem", m01["role"]["dissertation_role"])
        self.assertIn("P01-P05 provide the breadth test", m01["role"]["general_problem_bridge"])
        self.assertIn("worked example", m01["role"]["not_primary_subject"])

    def test_competing_histories_are_not_preferred(self) -> None:
        m01 = target.validate_m01_contract()
        histories = m01["competing_histories"]
        self.assertEqual([row["history_id"] for row in histories], ["H1", "H2", "H3"])
        self.assertTrue(all(row["prior_status"] == "COMPETING_NOT_PREFERRED" for row in histories))
        self.assertIn("gained or regained", histories[1]["statement"])

    def test_core_connection_runs_from_species_tip_to_causal_ecology(self) -> None:
        m01 = target.validate_m01_contract()
        connection = m01["core_connection"]
        self.assertIn("not yet a historical event", connection["stage_1_species_tip_contrast"])
        self.assertIn("individual/population ancestry", connection["stage_2_linked_history"])
        self.assertIn("pigment chemistry", connection["stage_3_transition_decomposition"])
        self.assertIn("population genomics", connection["stage_4_selection_layer"])
        self.assertIn("separate functional or agent-specific ecological design", connection["stage_5_causal_ecology"])

    def test_evidence_ladder_preserves_claim_ceilings(self) -> None:
        m01 = target.validate_m01_contract()
        ladder = {row["level"]: row for row in m01["evidence_ladder"]}
        self.assertIn("no loss, regain", ladder["E0"]["maximum_claim"])
        self.assertIn("not a causal regulatory variant", ladder["E2"]["maximum_claim"])
        self.assertIn("does not identify pollinators", ladder["E3"]["maximum_claim"])
        self.assertIn("direct agent-specific ecological evidence", ladder["E4"]["maximum_claim"])

    def test_selection_lane_rejects_shortcuts(self) -> None:
        m01 = target.validate_m01_contract()
        self.assertFalse(m01["sampling_contract"]["field_manipulation_required"])
        self.assertFalse(m01["sampling_contract"]["pollinator_observation_required"])
        self.assertIn(
            "two-species FST outlier alone",
            m01["selection_contract"]["prohibited_shortcuts"],
        )
        self.assertIn(
            "genomics alone cannot identify",
            m01["selection_contract"]["selective_agent_boundary"],
        )

    def test_m01_starts_fail_closed(self) -> None:
        m01 = target.validate_m01_contract()
        state = m01["current_state"]
        self.assertEqual(state["own_m01_biological_records_admitted"], 0)
        self.assertFalse(state["tissue_collection_authorized"])
        self.assertFalse(state["field_manipulation_authorized"])
        self.assertFalse(state["regain_claim_authorized"])
        self.assertFalse(state["selection_claim_authorized"])
        self.assertFalse(state["pollinator_agent_claim_authorized"])

    def test_m01_narrative_is_bounded(self) -> None:
        target.validate_narrative()


if __name__ == "__main__":
    unittest.main()
