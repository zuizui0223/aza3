"""Head module mosaic: synthetic tests never count as new field biology."""
import copy
import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]/"analysis"))
import audit_cirsium_two_radiation_head_mosaic_v1 as m


def na(taxon, orientation, outer):
    return {"Species":taxon, "HeadPosition":orientation,
            "OuterPhyllaryPosition":outer}


def jp(taxon, orientation, outer):
    return {"paper_japan_member_id":taxon,
            "paper_taxon_concept":"Cirsium " + taxon,
            "orientation_state":orientation,
            "phyllary_posture":outer,
            "authority_match_status":"exact_authority_concept_match"}


class TestMosaic(unittest.TestCase):
    def test_north_american_binary_geometry_codes(self):
        self.assertEqual(m.na_state(na("x","0","0")), (0,0))
        self.assertEqual(m.na_state(na("x","1","2")), (1,1))
        self.assertEqual(m.na_state(na("x","0","1")), (0,1))
        with self.assertRaisesRegex(ValueError,"INVALID_NA_HEAD_ORIENTATION"):
            m.na_state(na("x","2","0"))
        with self.assertRaisesRegex(ValueError,"INVALID_NA_OUTER_PHYLLARY"):
            m.na_state(na("x","0","unknown"))

    def test_japanese_coarse_categories_are_validated_not_force_coded(self):
        self.assertEqual(m.jp_state(jp("A","upward_or_erect","appressed_or_ascending")), (0,0))
        self.assertEqual(m.jp_state(jp("A","upward_or_erect","spreading_or_recurved")), (0,1))
        self.assertEqual(m.jp_state(jp("B","downward_or_nodding","spreading")), (1,1))
        self.assertEqual(m.jp_state(jp("C","upward_or_erect","ascending_or_recurved")), (0,None))
        self.assertEqual(m.jp_state(jp("D","source_conflict_index_downward_detail_erect","appressed")), (None,0))
        self.assertEqual(m.jp_state(jp("E","upward_or_ascending","unknown")), (None,None))
        with self.assertRaisesRegex(ValueError,"INVALID_JP_PHYLLARY_SOURCE_CODE"):
            m.jp_state(jp("F","upward_or_erect","angry"))
        with self.assertRaisesRegex(ValueError,"INVALID_JP_ORIENTATION_SOURCE_CODE"):
            m.jp_state(jp("F","horizontal","ascending"))

    def test_four_state_matrix_and_leave_one_tip_guards(self):
        records=[(0,0),(0,0),(0,1),(1,0),(1,1)]
        r=m.matrix_summary(records)
        self.assertEqual(r["n_taxa"],5)
        self.assertEqual(r["four_cell_counts"],{
            "erect_not_outward":2,"erect_outward":1,
            "nodding_not_outward":1,"nodding_outward":1})
        self.assertTrue(r["all_four_observed"])
        self.assertEqual(r["leave_one_taxon_out_preserves_four_n"],2)
        self.assertFalse(r["all_single_taxon_deletions_preserve_four"])

    def test_completely_nonsynchronous_structure_with_either_orientation(self):
        records=[(0,0),(0,1),(1,0),(1,1)]
        r=m.matrix_summary(records)
        self.assertFalse(r["deterministic_outer_given_orientation"])
        self.assertFalse(r["deterministic_orientation_given_outer"])

    def test_three_state_subset_must_not_claim_four(self):
        r=m.matrix_summary([(0,0),(0,0),(1,0),(1,1)])
        self.assertEqual(r["n_observed_combinations"],3)
        self.assertFalse(r["all_four_observed"])
        self.assertTrue(r["all_single_taxon_deletions_preserve_four"] is False)

    def test_fisher_manual_cumulative_exact(self):
        self.assertAlmostEqual(m.fisher_exact_two_sided(1,0,0,1),1.0)
        self.assertAlmostEqual(m.fisher_exact_two_sided(3,0,0,3),0.1)
        self.assertAlmostEqual(m.fisher_exact_two_sided(1,2,1,3),1.0)
        self.assertAlmostEqual(m.fisher_exact_two_sided(1,2,0,3),1.0)
        self.assertAlmostEqual(m.fisher_exact_two_sided(6,1,10,46),
                               0.0007009516685550101,places=12)

    def test_refuses_duplicated_source_taxon(self):
        with self.assertRaisesRegex(ValueError,"DUPLICATE_SOURCE_TAXON"):
            m.unique_rows([na("X","0","0"),na("X","1","1")],"Species")

    def test_dominant_japan_does_not_silently_include_secondary_history(self):
        self.assertEqual(m.role("JPN_06"),"secondary")
        self.assertEqual(m.role("JPN_15"),"secondary")
        self.assertEqual(m.role("JPN_13"),"dominant")

    def test_missing_japan_outer_is_not_no_armature_or_ascending(self):
        r=m.jp_state(jp("X","downward_or_nodding","unknown"))
        self.assertEqual(r,(1,None))

    def test_claim_ceiling_fails_closed(self):
        limits=m.claim_ceiling()
        self.assertFalse(limits["observed_morphology_equals_adaptive_module"])
        self.assertFalse(limits["four_tip_combinations_prove_independent_transitions"])
        self.assertFalse(limits["jp_and_na_ontologies_are_identical_quantitative_angles"])
        self.assertFalse(limits["same_head_visitor_access_fitness_verified"])
        self.assertFalse(limits["genetic_independence_or_ancestral_trait_memory_demonstrated"])


if __name__=="__main__":
    unittest.main()
