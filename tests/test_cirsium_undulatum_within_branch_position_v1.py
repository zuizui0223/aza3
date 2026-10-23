"""Original USDA Cirsium head location vs injury and seed output, source guards."""
from pathlib import Path
import copy
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_undulatum_within_branch_position_v1 as v

DATA=ROOT/"data/external/usda_wavyleaf_2007_2008"

class TestRealWithinBranchPositions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mothers=v.read_locked(DATA/"mother_plants.csv","mother_plants.csv",
                                   ("plt","year","dam","ins","trt"))
        cls.heads=v.read_locked(DATA/"flower_heads.csv","flower_heads.csv",
                                 ("plt","year","dam","ins","hd","av.dam","hd.gsd"))
        cls.dictionary=v.read_locked(DATA/"source_dictionary.txt","source_dictionary.txt")

    def test_original_position_code_is_anatomical_order_not_angle(self):
        self.assertEqual(v.position({"hd":"100"}),"apical")
        self.assertEqual(v.position({"hd":"200"}),"branch_terminal")
        self.assertEqual(v.position({"hd":"1300"}),"branch_terminal")
        self.assertEqual(v.position({"hd":"210"}),"intermediate")
        self.assertEqual(v.position({"hd":"410"}),"intermediate")
        with self.assertRaisesRegex(ValueError,"INVALID_HEAD_POSITION_CODE"):
            v.position({"hd":"311"})

    def test_uses_297_source_aligned_real_head_records(self):
        ix,clean,audit=v.raw_audit(self.mothers,self.heads)
        self.assertEqual(audit["n_original_mothers"],106)
        self.assertEqual(audit["n_original_head_rows"],323)
        self.assertEqual(audit["n_clean_head_rows"],297)
        self.assertEqual(audit["source_position_counts"],
                         {"apical":103,"branch_terminal":173,"intermediate":21})
        self.assertEqual(audit["n_original_rows_with_conflicting_treatment"],12)
        self.assertEqual(audit["source_conflict_mother_ids"],["2007|323","2007|375"])

    def test_pair_conditioning_not_194_independent_laterals(self):
        _,clean,_=v.raw_audit(self.mothers,self.heads)
        records=v.pair_mothers(self.mothers,clean)
        self.assertEqual(len(records),14)
        self.assertEqual([sum(x["year"]==year for x in records)
                          for year in v.YEARS],[7,7])
        self.assertTrue(all(x["n_branch_terminal"]>=1 and
                            x["n_intermediate"]>=1 for x in records))

    def test_2008_intermediate_seed_shortfall_is_all_seven_but_exploratory(self):
        d=v.analyze(self.mothers,self.heads,self.dictionary,reps=100)
        r=d["paired_within_mother_year_results"]["2008"]["position_differences"]
        x=r["good_seed_per_head"]
        self.assertEqual((x["n_increase"],x["n_decrease"],x["n_equal"]),(0,7,0))
        self.assertAlmostEqual(x["posthoc_exact_two_sided_sign_p"],.015625)
        self.assertAlmostEqual(x["BH_q_six_exploratory_sign_tests"],.09375)
        self.assertAlmostEqual(x["mean_intermediate_minus_branch_terminal"],-27.476190476190474)
        self.assertTrue(d["limits"]["posthoc_six_tests_do_not_confirm_architectural_adaptation"])

    def test_strongest_mother_outlier_does_not_drive_sign_pattern(self):
        _,clean,_=v.raw_audit(self.mothers,self.heads)
        rr=v.pair_mothers(self.mothers,clean)
        z=[x["delta"]["good_seed_per_head"] for x in rr
           if x["year"]=="2008" and x["id"]!="2008|WG25"]
        self.assertEqual(len(z),6)
        self.assertTrue(all(x<0 for x in z))
        self.assertAlmostEqual(sum(z)/len(z),-13.97222222222222,places=8)

    def test_original_damage_differences_do_not_follow_seed_differences(self):
        d=v.analyze(self.mothers,self.heads,self.dictionary,reps=100)
        x=d["paired_within_mother_year_results"]["2008"]["position_differences"]["damage_score"]
        self.assertEqual((x["n_increase"],x["n_decrease"],x["n_equal"]),(4,3,0))
        self.assertEqual(x["posthoc_exact_two_sided_sign_p"],1.0)
        self.assertEqual(
          d["paired_within_mother_year_results"]["2008"]["position_differences"]["seed_positive_fraction"]["n_decrease"],5)

    def test_2007_mixed_signs_not_copy_of_2008(self):
        d=v.analyze(self.mothers,self.heads,self.dictionary,reps=100)
        r=d["paired_within_mother_year_results"]["2007"]["position_differences"]["good_seed_per_head"]
        self.assertEqual((r["n_increase"],r["n_decrease"]),(3,4))
        self.assertAlmostEqual(r["mean_intermediate_minus_branch_terminal"],-0.5214285714285715)

    def test_exact_sign_on_rare_paired_clusters(self):
        self.assertAlmostEqual(v.exact_sign_p(0,7),.015625)
        self.assertEqual(v.exact_sign_p(4,3),1.0)
        self.assertEqual(v.exact_sign_p(0,0),1.0)
        self.assertEqual(v.bh({"a":.015625,"b":1.0})["a"],.03125)

    def test_source_missing_damage_cannot_be_imputed_as_zero(self):
        x=v.head_stats([{"hd.gsd":"0","av.dam":"NA"},
                        {"hd.gsd":"4","av.dam":"0.88"}])
        self.assertEqual(x["n_scored_damage"],1)
        self.assertAlmostEqual(x["mean_damage"],.88)
        self.assertEqual(x["n_positive_seed_heads"],1)

    def test_invalid_source_dictionary_fails(self):
        with self.assertRaisesRegex(ValueError,"SOURCE_HEAD_POSITION_ONTOLOGY_UNVERIFIED"):
            v.analyze(self.mothers,self.heads,"unverified",reps=100)

    def test_original_data_are_not_experiment_on_head_placement(self):
        d=v.analyze(self.mothers,self.heads,self.dictionary,reps=100)
        self.assertTrue(d["limits"]["not_original_head_position_randomization"])
        self.assertTrue(d["limits"]["very_small_seven_paired_mothers_in_2008"])
        self.assertEqual(d["limits"]["new_individual_head_observations_collected"],0)

if __name__=="__main__":
    unittest.main()
