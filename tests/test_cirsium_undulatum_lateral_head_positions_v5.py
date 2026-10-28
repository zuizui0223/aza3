"""Real 2008 Cirsium lateral-head positional rerouting; source head codes are NOT timestamps.

Source: West & Louda (2021) USDA-ADC CC0 DOI 10.15482/USDA.ADC/1522649.
One plant is one observed opportunity set; all plant/head original records are
from the SAME original experimental cohort as earlier Cirsium reports.
"""
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from audit_cirsium_undulatum_lateral_head_positions_v5 import (
    analyse_lateral_positions, classify_lateral_position_routes,
    source_prefix_group
)
from audit_cirsium_undulatum_scored_head_failure_v4 import source_locked_tables


class TestOriginalLateralHeadPosition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        moms,heads=source_locked_tables()
        cls.rows=moms
        cls.heads=heads
        cls.x=analyse_lateral_positions(moms,heads,strict=True)

    def test_2008_actual_source_cohort_and_not_independent_heads(self):
        s=self.x["source"]
        self.assertEqual(s["original_mothers"],106)
        self.assertEqual(s["original_heads"],323)
        self.assertEqual(s["jointly_validated_2008_mothers"],67)
        self.assertEqual(s["original_2008_scored_lateral_heads"],111)
        self.assertEqual(s["mothers_with_2plus_scored_lateral_heads"],33)
        self.assertTrue(s["all_2008_head_good_seed_sums_match_maternal_source"])
        self.assertFalse(self.x["claim_limits"]["head_position_codes_are_proven_chronological_order"])

    def test_actual_four_treatment_position_contingency(self):
        expected={
            "ndni":(6,1,0,0,5),
            "dni": (6,2,0,3,1),
            "ndi": (8,1,6,1,0),
            "di":  (13,0,12,1,0),
        }
        for arm,(n,first_no_later_yes,first_yes_later_yes,
                 first_yes_later_no,first_no_later_no) in expected.items():
            a=self.x["by_2008_treatment"][arm]["lateral_rank_presence_2plus"]
            self.assertEqual(a["mothers"],n)
            self.assertEqual(a["first_zero_later_positive"],first_no_later_yes)
            self.assertEqual(a["first_positive_later_positive"],first_yes_later_yes)
            self.assertEqual(a["first_positive_later_zero"],first_yes_later_no)
            self.assertEqual(a["first_zero_later_zero"],first_no_later_no)
            self.assertEqual(sum(v for k,v in a.items() if k!="mothers"),n)

    def test_apical_and_first_lateral_failure_with_higher_code_success(self):
        a=self.x["by_2008_treatment"]["dni"]
        self.assertEqual(a["first_zero_later_positive_mothers"],["ES68","SD05"])
        self.assertEqual(a["first_zero_later_positive_mothers_with_apical_good_seed_zero"],2)
        self.assertEqual(a["lateral_seed_positive_mothers"],7)
        self.assertEqual(a["first_zero_later_positive_in_positive_mothers"],2)
        self.assertFalse(self.x["claim_limits"]["higher_numbered_head_is_proven_later_flowering"])

    def test_source_id_prefixes_protect_against_one_prefix_only_story(self):
        p=self.x["original_id_prefix_strata"]
        self.assertEqual(sorted(p),["ES","SD","WG","WR"])
        self.assertEqual([p[k]["total_mothers"] for k in ("ES","SD","WG","WR")],
                         [17,17,12,21])
        zero={k:p[k]["original_arms"]["dni"]["both_seed_routes_zero"]
              for k in p}
        self.assertEqual(zero,{"ES":0,"SD":2,"WG":2,"WR":1})
        intact_zero={k:p[k]["original_arms"]["ndni"]["both_seed_routes_zero"]
                    for k in p}
        self.assertEqual(intact_zero,{"ES":2,"SD":4,"WG":3,"WR":3})
        self.assertFalse(self.x["claim_limits"]["source_id_prefix_verified_as_geographic_site"])

    def test_nonapical_position_is_author_number_not_arbitrary_string_sort(self):
        self.assertEqual(classify_lateral_position_routes(
            [{"hd":"1300","hd.gsd":"4"},{"hd":"200","hd.gsd":"0"},
             {"hd":"310","hd.gsd":"0"}])["class"],"first_zero_later_positive")
        self.assertEqual(classify_lateral_position_routes(
            [{"hd":"200","hd.gsd":"0"}])["class"],"one_scored_lateral")
        self.assertEqual(classify_lateral_position_routes([])["class"],"none_scored_lateral")
        self.assertEqual(source_prefix_group("WG32.1"),"WG")
        with self.assertRaisesRegex(ValueError,"INVALID_HEAD_POSITION"):
            classify_lateral_position_routes(
                [{"hd":"odd","hd.gsd":"4"},{"hd":"200","hd.gsd":"0"}])

    def test_no_claim_of_adaptive_architecture_or_genetic_afterlife(self):
        self.assertTrue(self.x["biology"]["original_nonapical_heads_do_have_multiple_source_position_codes"])
        for name,val in self.x["claim_limits"].items():
            self.assertFalse(val,name)
        self.assertEqual(self.x["conclusion"],
          "ACTUAL_NONAPICAL_HEAD_CODE_ALTERNATIVES_WITH_POSITIONAL_FALLBACK_NOT_EVOLUTIONARY_REGENERATION")


if __name__=="__main__":
    unittest.main()
