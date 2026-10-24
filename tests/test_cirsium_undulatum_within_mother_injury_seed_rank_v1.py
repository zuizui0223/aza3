"""Only original USDA Cirsium maternal and head records, not synthetic study findings."""
import copy
import math
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_undulatum_within_mother_injury_seed_rank_v1 as v

SRC=ROOT/"data/external/usda_wavyleaf_2007_2008"

class OriginalHeadInjurySeedRank(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m=v.parse_locked(SRC/"mother_plants.csv","mother_plants.csv")
        cls.h=v.parse_locked(SRC/"flower_heads.csv","flower_heads.csv")
        cls.d=v.parse_locked(SRC/"source_dictionary.txt","source_dictionary.txt")
        cls.report=v.analyze(cls.m,cls.h,cls.d)

    def test_real_usda_mother_head_and_label_counts(self):
        s=self.report["source_audit"]
        self.assertEqual(s["original_mother_rows"],106)
        self.assertEqual(s["original_head_rows"],323)
        self.assertEqual(s["head_with_source_mother"],309)
        self.assertEqual(s["head_without_source_mother"],14)
        self.assertEqual(s["treatment_mismatch_head_rows"],12)
        self.assertEqual(s["consistent_heads"],297)
        self.assertEqual(s["consistent_lateral_head_rows"],194)
        self.assertEqual(s["scored_lateral_head_rows"],192)
        self.assertEqual(s["scored_terminal_head_rows"],171)
        self.assertEqual(s["scored_intermediate_head_rows"],21)

    def test_author_damage_is_a_scored_fraction_of_quarters_not_germination(self):
        o=self.report["damage_ontology"]
        self.assertEqual(o["maximum_recorded_score"],.88)
        self.assertTrue(o["max_score_is_coded_not_exact_fraction_of_tissue_lost"])
        self.assertTrue(o["author_good_seeds_are_not_verified_germination"])

    def test_actual_48_mothers_all_lateral_head_extrema(self):
        z=self.report["analyses"]["lateral_all_classes"]["summary"]
        self.assertEqual(z["n_original_mothers"],48)
        self.assertEqual(z["directions"],{
            "lower_damage_has_more_good_seeds":26,
            "higher_damage_has_more_good_seeds":17,
            "equal_good_seed_counts":5})
        self.assertAlmostEqual(z["mean_lower_minus_higher_good_seed"],
                               6.840277777777778)
        self.assertAlmostEqual(z["exact_two_sided_p_for_symmetry_exploratory"],
                               .2220528201560228)

    def test_actual_terminal_class_alone_defeats_monotone_predictive_rule(self):
        z=self.report["analyses"]["terminal_only"]["summary"]
        self.assertEqual(z["n_original_mothers"],46)
        self.assertEqual(z["directions"],{
            "lower_damage_has_more_good_seeds":23,
            "higher_damage_has_more_good_seeds":18,
            "equal_good_seed_counts":5})
        self.assertAlmostEqual(z["exact_two_sided_p_for_symmetry_exploratory"],
                               .5327092552361137)
        self.assertAlmostEqual(z["mean_lower_minus_higher_good_seed"],
                               4.380434782608695)
        self.assertEqual(z["n_higher_damage_head_produces_any_good_seed"],36)
        self.assertEqual(z["n_lower_damage_head_fails_while_higher_damage_succeeds"],5)

    def test_each_original_mother_once_even_with_many_head_records(self):
        for key in ("terminal_only","lateral_all_classes"):
            rs=self.report["analyses"][key]["per_mother"]
            self.assertEqual(len(rs),len({r["mother"] for r in rs}))
            for row in rs:
                self.assertGreaterEqual(row["eligible_scored_heads"],2)
                self.assertGreater(row["damage_score_gap"],0)
                self.assertTrue(all(int(h["head_code"])>100
                                    for h in row["low_heads"]+row["high_heads"]))
                self.assertAlmostEqual(
                    row["least_damaged_mean_good_seed"]-
                    row["most_damaged_mean_good_seed"],
                    row["lower_damage_minus_higher_damage_good_seed"])

    def test_true_source_year_split_without_pooled_pseudoreplication(self):
        years=self.report["analyses"]["terminal_only"]["by_year_arm"]
        self.assertEqual(years["2007"]["total"]["n_original_mothers"],22)
        self.assertEqual(years["2007"]["total"]["directions"],{
            "lower_damage_has_more_good_seeds":12,
            "higher_damage_has_more_good_seeds":9,
            "equal_good_seed_counts":1})
        self.assertEqual(years["2008"]["total"]["n_original_mothers"],24)
        self.assertEqual(years["2008"]["total"]["directions"],{
            "lower_damage_has_more_good_seeds":11,
            "higher_damage_has_more_good_seeds":9,
            "equal_good_seed_counts":4})
        self.assertAlmostEqual(
            years["2007"]["total"]["exact_two_sided_p_for_symmetry_exploratory"],
            .6636238098144531)
        self.assertAlmostEqual(
            years["2008"]["total"]["exact_two_sided_p_for_symmetry_exploratory"],
            .8238029479980469)

    def test_individual_true_mother_injury_reversal_not_synthetic(self):
        mother={z["mother"]:z
                 for z in self.report["analyses"]["terminal_only"]["per_mother"]}
        a=mother["2007|351"]
        self.assertEqual(a["least_damaged_score"],.205)
        self.assertEqual(a["most_damaged_score"],.755)
        self.assertEqual(a["least_damaged_mean_good_seed"],0)
        self.assertEqual(a["most_damaged_mean_good_seed"],81)
        self.assertEqual(a["low_heads"][0]["head_code"],"400")
        self.assertEqual(a["high_heads"][0]["head_code"],"300")
        b=mother["2008|ES69"]
        self.assertEqual(b["least_damaged_score"],.205)
        self.assertEqual(b["most_damaged_score"],.88)
        self.assertEqual(b["least_damaged_mean_good_seed"],0)
        self.assertEqual(b["most_damaged_mean_good_seed"],42)
        self.assertTrue(a["severe_success_vs_low_failure"])
        self.assertTrue(b["severe_success_vs_low_failure"])

    def test_identical_scored_damage_different_seed_counts_on_same_mother(self):
        rr=self.report["equal_scored_damage_distinct_seed_outputs"]
        s=rr["summary"]
        self.assertEqual(s["n_same_source_damage_score_terminal_groups"],13)
        self.assertEqual(s["n_distinct_source_mothers"],13)
        self.assertEqual(s["n_same_damage_but_different_good_seed_counts"],9)
        self.assertEqual(s["n_same_damage_but_seed_positive_vs_zero"],8)
        self.assertEqual(s["by_year"]["2007"],{
            "eligible_groups":3,"variable_good_seed_groups":3})
        self.assertEqual(s["by_year"]["2008"],{
            "eligible_groups":10,"variable_good_seed_groups":6})
        same={x["mother"]:x for x in rr["groups"]}
        self.assertEqual(same["2008|SD04"]["damage_score"],.5175)
        self.assertEqual((same["2008|SD04"]["min_good_seed"],
                          same["2008|SD04"]["max_good_seed"]),(134,148))

    def test_prior_original_same_branch_code_contrast_is_sensitivity_not_independent_study(self):
        z=self.report["prior_source_nested_branch_germination_not_measured"]
        self.assertEqual(z["n_parent_source_code_matched_pairs"],20)
        self.assertEqual(z["by_year"]["2007"]["n"],11)
        self.assertEqual(z["by_year"]["2008"]["n"],9)
        self.assertEqual(
            z["by_year"]["2007"]["damage_and_seed_direction_cross_tab"]
             ["positive_damage__positive_good_seed"],1)
        self.assertEqual(
            z["by_year"]["2008"]["damage_and_seed_direction_cross_tab"]
             ["equal_damage__negative_good_seed"],2)

    def test_scored_injury_not_sufficient_for_seed_success_prediction(self):
        rr=self.report["equal_scored_damage_distinct_seed_outputs"]["groups"]
        pairs=[x for x in rr if x["mixed_zero_and_positive"]]
        self.assertEqual(len(pairs),8)
        self.assertTrue(all(x["min_good_seed"]==0 and x["max_good_seed"]>0
                            for x in pairs))

    def test_negative_duplicate_head_id_rejected(self):
        heads=copy.deepcopy(self.h)
        heads[1]["hd"]=heads[0]["hd"]
        with self.assertRaisesRegex(ValueError,"DUPLICATE_ORIGINAL_HEAD_RECORD"):
            v.analyze(self.m,heads,self.d)

    def test_negative_conflicting_original_maternal_treatment_rejected(self):
        heads=copy.deepcopy(self.h)
        heads[0]["ins"]="i" if heads[0]["ins"]=="ni" else "ni"
        with self.assertRaisesRegex(ValueError,
                                     "WRONG_ORIGINAL_PARENT_HEAD_JOIN_OR_TREATMENT_CONFLICTS"):
            v.analyze(self.m,heads,self.d)

    def test_negative_injury_code_exceeds_original_scale(self):
        heads=copy.deepcopy(self.h)
        heads[1]["av.dam"]="1.01"
        with self.assertRaisesRegex(ValueError,"INVALID_CODED_RECEPTACLE_DAMAGE"):
            v.analyze(self.m,heads,self.d)

    def test_negative_fabricated_good_seed_count_rejected(self):
        heads=copy.deepcopy(self.h)
        heads[1]["hd.gsd"]="-15"
        with self.assertRaisesRegex(ValueError,"INVALID_AUTHOR_GOOD_SEEDS"):
            v.analyze(self.m,heads,self.d)

    def test_negative_dictionary_semantics_changed(self):
        with self.assertRaisesRegex(ValueError,"AUTHOR_DICTIONARY_POSITION_OR_DAMAGE_CHANGED"):
            v.analyze(self.m,self.h,self.d.replace("four quarters of the receptacle",
                                                  "unverified")) 

    def test_source_naive_vs_within_mother_association_not_conflated(self):
        z=self.report["terminal_score_vs_good_seed_pearson_diagnostic"]
        for year,n,m,rep,naive,within in [
            ("2007",69,31,60,-.06400727450417101,-.1393467469601251),
            ("2008",102,55,79,-.4539108689115904,-.26984588634097434),
            ("all_years_descriptive_only",171,86,139,
              -.2938450971688947,-.20954884460100218)]:
            x=z[year]
            self.assertEqual(x["n_original_head_records"],n)
            self.assertEqual(x["n_original_mothers"],m)
            self.assertEqual(x["n_terminal_heads_in_mothers_at_least_two"],rep)
            self.assertAlmostEqual(x["pooled_head_level_r_not_independent"],
                                   naive,places=12)
            self.assertAlmostEqual(x["within_mother_demeaned_r_descriptive"],
                                   within,places=12)
            self.assertLess(x["pooled_head_level_r_not_independent"],0)
            self.assertLess(x["within_mother_demeaned_r_descriptive"],0)
        self.assertEqual(z["all_years_descriptive_only"]
                         ["n_mothers_at_least_two_heads"],54)

    def test_exact_sign_reference_never_measures_functional_randomized_damage(self):
        self.assertAlmostEqual(v.exact_sign(23,18),.5327092552361137)
        self.assertEqual(v.exact_sign(0,0),1)
        with self.assertRaisesRegex(ValueError,"INVALID_EXACT_SIGNS"):
            v.exact_sign(-1,2)
        lim=self.report["limits"]
        self.assertTrue(lim["min_max_selection_on_observed_posttreatment_damage_not_randomized"])
        self.assertTrue(lim["no_causal_insect_seed_fitness_or_genetic_LEGO_support"])
        self.assertEqual(lim["n_new_field_observations"],0)

if __name__=="__main__":
    unittest.main()
