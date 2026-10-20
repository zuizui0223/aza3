"""Do not collapse source-described variable phyllary direction to a fake fixed tip."""
import unittest
import pathlib
import sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"analysis"))
import audit_jpn19_source_phyllary_range_v1 as v


def source(state="ascending_or_recurved",orientation="upward_or_erect"):
    return {"paper_japan_member_id":"JPN_19",
            "paper_taxon_concept":"Cirsium nippoense",
            "orientation_state":orientation,
            "phyllary_posture":state,
            "stickiness_state":"sticky",
            "source_url":"https://www.kahaku.go.jp/research/db/botany/azami/detail.html?no=110",
            "authority_match_status":"exact_authority_concept_match"}

BASE={"erect_not_outward":3,"erect_outward":0,
      "nodding_not_outward":2,"nodding_outward":1}


class SourceAmbiguity(unittest.TestCase):
    def test_partial_identification_two_possible_extensions(self):
        r=v.outcome_envelope(source(),BASE)
        self.assertEqual(r["forced_new_source_tip_state"],None)
        self.assertEqual(r["if_ascending"]["four_cell_counts"]["erect_not_outward"],4)
        self.assertEqual(r["if_ascending"]["n_observed_combinations"],3)
        self.assertEqual(r["if_recurved"]["four_cell_counts"]["erect_outward"],1)
        self.assertEqual(r["if_recurved"]["n_observed_combinations"],4)
        self.assertEqual(r["n_observed_combinations_bound"],[3,4])
        self.assertEqual(r["n_existing_confirmed_dominant_taxa"],6)
        self.assertFalse(r["within_single_plant_polymorphism_demonstrated"])

    def test_invalid_source_recode_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"ORIGINAL_AMBIGUOUS_SOURCE_CHANGED"):
            v.outcome_envelope(source("spreading"),BASE)

    def test_invalid_orientation_is_not_force_coded(self):
        with self.assertRaisesRegex(ValueError,"UNRESOLVED_ERECT_SOURCE_STATE"):
            v.outcome_envelope(source(orientation="upward_or_ascending"),BASE)

    def test_only_prior_source_main6_matrix_is_valid(self):
        other=dict(BASE);other["erect_outward"]=1
        with self.assertRaisesRegex(ValueError,"PRIOR_SOURCE_COHORT_CHANGED"):
            v.outcome_envelope(source(),other)

    def test_no_biological_reactivation_claim(self):
        x=v.claim_limit()
        self.assertFalse(x["jpn19_fixed_outward_morph_assigned"])
        self.assertFalse(x["within_plant_reversible_phyllary_change_measured"])
        self.assertFalse(x["actual_adaptive_selection_or_pollination_measured"])


if __name__=="__main__":
    unittest.main()
