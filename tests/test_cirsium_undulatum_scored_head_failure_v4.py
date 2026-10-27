"""Source-authenticated Cirsium undulatum real head-recorded alternative routes.

Tests MUST distinguish a scored lateral head, a lateral flower reaching
anthesis, and a lateral head making >0 author-qualified good seed.
Only the first and third are direct from the 323-row head CSV.

Original West & Louda (2021), DOI 10.1002/ajb2.1744;
USDA ADC DOI 10.15482/USDA.ADC/1522649. No new samples.
"""
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from audit_cirsium_undulatum_scored_head_failure_v4 import (
    audit_scored_heads, head_mother_crosswalk, score_damage_as_number,
    source_locked_tables,
)


class CirsiumRealHeadFailureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mothers,cls.heads=source_locked_tables()
        cls.r=audit_scored_heads(cls.mothers,cls.heads,strict=True)

    def test_source_globally_has_106_mothers_323_heads_and_frozen_sha(self):
        s=self.r["source"]
        self.assertEqual(s["original_mother_rows"],106)
        self.assertEqual(s["original_head_rows"],323)
        self.assertEqual(s["original_2008_mothers"],67)
        self.assertEqual(s["original_2008_head_rows"],178)
        self.assertEqual(s["original_2008_lateral_head_rows"],111)
        self.assertEqual(s["original_2008_apical_head_rows"],67)
        self.assertTrue(s["mother_head_reconciliation_all_2008"])
        self.assertEqual(s["mother_source_sha1"],"458800272c3d0a03c7653702027a8f8578d9e88c")
        self.assertEqual(s["head_source_sha1"],"637a4dcd5194d0bf5b26e7f9cf0dde56a50d7616")

    def test_2008_damaged_ambient_joint_failures_have_scored_lateral_heads(self):
        t=self.r["year2008_by_original_treatment"]["dni"]
        self.assertEqual(t["mothers"],12)
        self.assertEqual(t["mothers_both_seed_routes_zero"],5)
        self.assertEqual(t["joint_failure_mothers_with_scored_lateral_head"],5)
        self.assertEqual(t["joint_failure_mothers_with_two_plus_scored_lateral_heads"],1)
        self.assertEqual(t["joint_failure_mothers_with_3plus_source_buds"],4)
        self.assertEqual(t["joint_failure_mothers_with_lateral_scored_good_seed_zero"],5)
        self.assertEqual(t["mothers_with_lateral_head_record"],12)
        self.assertEqual(t["mothers_with_positive_lateral_good_seed"],7)
        self.assertFalse(self.r["scientific_claim_limits"]["scored_lateral_head_means_confirmed_anthesis"])

    def test_intact_ambient_failure_not_exclusively_absent_side_heads(self):
        t=self.r["year2008_by_original_treatment"]["ndni"]
        self.assertEqual(t["mothers"],15)
        self.assertEqual(t["mothers_both_seed_routes_zero"],12)
        self.assertEqual(t["joint_failure_mothers_with_scored_lateral_head"],9)
        self.assertEqual(t["joint_failure_mothers_with_two_plus_scored_lateral_heads"],5)
        self.assertEqual(t["joint_failure_mothers_with_3plus_source_buds"],9)
        self.assertEqual(t["mothers_with_lateral_head_record"],11)
        self.assertEqual(t["mothers_with_positive_lateral_good_seed"],2)

    def test_insecticide_groups_include_mothers_with_apex_but_no_scored_lateral(self):
        t=self.r["year2008_by_original_treatment"]["di"]
        self.assertEqual(t["mothers"],19)
        self.assertEqual(t["mothers_with_lateral_head_record"],15)
        self.assertEqual(t["mothers_with_positive_lateral_good_seed"],15)
        self.assertEqual(t["mothers_both_seed_routes_zero"],0)
        self.assertEqual(t["mothers_without_scored_lateral_head"],4)
        self.assertEqual(t["mothers_with_scored_lateral_but_zero_lateral_seed"],0)
        self.assertEqual(t["apical_seed_positive_mothers"],11)

    def test_highest_recorded_damage_does_not_force_seed_zero(self):
        t=self.r["year2008_by_original_treatment"]["dni"]
        self.assertEqual(t["scored_lateral_head_rows"],23)
        self.assertEqual(t["scored_lateral_good_seed_positive_head_rows"],7)
        self.assertEqual(t["scored_lateral_good_seed_zero_head_rows"],16)
        self.assertEqual(t["lateral_head_rows_at_max_observed_damage_088"],17)
        self.assertEqual(t["lateral_positive_seed_heads_at_max_damage_088"],6)
        self.assertEqual(t["lateral_zero_seed_heads_at_max_damage_088"],11)
        self.assertAlmostEqual(self.r["source"]["max_recorded_head_damage_score"],0.88)
        self.assertFalse(self.r["scientific_claim_limits"]["damage_score_deterministically_classifies_reproductive_outcome"])

    def test_within_same_mother_damage_not_consistently_lower_on_positive_lateral_head(self):
        t=self.r["within_mother_lateral_head_damage_comparison"]
        self.assertEqual(t["n_2008_mothers_with_both_positive_and_zero_lateral_heads"],12)
        self.assertEqual(t["successful_lateral_heads_lower_damage_mothers"],5)
        self.assertEqual(t["successful_lateral_heads_higher_damage_mothers"],5)
        self.assertEqual(t["successful_lateral_heads_same_damage_mothers"],2)
        self.assertEqual(t["by_treatment"]["dni"],{"n":5,"lower":1,"higher":3,"same":1})
        self.assertFalse(self.r["scientific_claim_limits"]["within_mother_damage_association_is_causal"])

    def test_2007_individual_head_archive_does_not_fully_justify_crossgrain_inference(self):
        q=self.r["other_year_archive_quality"]
        self.assertEqual(q["year"],"2007")
        self.assertEqual(q["mother_rows"],39)
        self.assertEqual(q["head_rows"],145)
        self.assertEqual(q["mother_ids_without_scored_head_rows"],["2007|329","2007|367"])
        self.assertEqual(q["mother_ids_with_lateral_seed_sum_disagreement"],["2007|327","2007|329","2007|367"])
        self.assertEqual(q["mother_ids_with_treatment_label_conflicts"],["2007|323","2007|375"])
        self.assertFalse(q["full_2007_mother_head_joint_state_authenticated"])

    def test_source_flowering_variable_semantics_not_taken_as_entire_whole_plant(self):
        c=self.r["source_semantics"]
        self.assertEqual(c["source_mothers_with_zero_no_flw_but_positive_total_seed"],10)
        self.assertEqual(c["source_2008_mothers_zero_no_flw_positive_seed"],6)
        self.assertFalse(c["no_flw_is_complete_verified_plant_head_anthesis_inventory"])
        self.assertFalse(c["score_hd_gsd_zero_guarantees_no_flower_opened"])

    def test_score_damage_missing_is_not_undegraded_or_zero(self):
        self.assertIsNone(score_damage_as_number("NA"))
        self.assertIsNone(score_damage_as_number(""))
        self.assertEqual(score_damage_as_number("0.88"),0.88)
        with self.assertRaisesRegex(ValueError,"DAMAGE_SCORE"):
            score_damage_as_number("-0.15")
        with self.assertRaisesRegex(ValueError,"DAMAGE_SCORE"):
            score_damage_as_number("wounded")

    def test_real_head_seed_total_mutation_fails_closed(self):
        mothers=[dict(r) for r in self.mothers]
        heads=[dict(r) for r in self.heads]
        target=next(r for r in heads if r["year"]=="2008"
                    and r["plt"]=="ES68" and r["hd"]=="300")
        target["hd.gsd"]=str(int(target["hd.gsd"])+2)
        with self.assertRaisesRegex(ValueError,"MOTHER_HEAD_GOOD_SEED_SUM_CONFLICT"):
            head_mother_crosswalk(mothers,heads,strict2008=True)

    def test_2008_source_treatment_conflict_fails_closed(self):
        mothers=[dict(r) for r in self.mothers]
        heads=[dict(r) for r in self.heads]
        target=next(r for r in heads if r["year"]=="2008" and r["plt"]=="ES68")
        target["ins"]="i"
        with self.assertRaisesRegex(ValueError,"SOURCE_2008_TREATMENT_CONFLICT"):
            head_mother_crosswalk(mothers,heads,strict2008=True)

    def test_duplicate_head_number_for_one_plant_rejected(self):
        mothers=[dict(r) for r in self.mothers]
        heads=[dict(r) for r in self.heads]
        first=next(r for r in heads if r["year"]=="2008")
        heads.append(dict(first))
        with self.assertRaisesRegex(ValueError,"DUPLICATE_HEAD_POSITION"):
            head_mother_crosswalk(mothers,heads,strict2008=False)

    def test_no_unmeasured_stamen_genetics_or_exact_slk_payoff_claim(self):
        flags=self.r["scientific_claim_limits"]
        for key,val in flags.items():
            self.assertFalse(val,key)
        self.assertEqual(self.r["conclusion"],
            "REAL_LATERAL_HEADS_RECORDED_DURING_JOINT_SEED_FAILURE_CAUSE_UNIDENTIFIED")


if __name__=="__main__":
    unittest.main()
