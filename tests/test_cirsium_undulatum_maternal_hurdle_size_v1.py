"""Synthetic correctness guards and original USDA Cirsium hurdle observations."""
from __future__ import annotations
import copy
import pathlib
import sys
import unittest

ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_undulatum_maternal_hurdle_size_v1 as v

BASE=ROOT/"data/external/usda_wavyleaf_2007_2008"


class TestOriginalUndulatumTwoStages(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m=v.read_locked(BASE/"mother_plants.csv","mother_plants.csv",
                             ("plt","year","trt","plt.gsd","no.flw","r.diam"))
        cls.h=v.read_locked(BASE/"flower_heads.csv","flower_heads.csv",
                             ("plt","year","dam","ins","hd","hd.gsd"))
        cls.d=v.read_locked(BASE/"source_dictionary.txt","source_dictionary.txt")
        cls.groups=v.source_groups(cls.m)

    def test_original_2007_2008_four_arms(self):
        self.assertEqual(len(self.m),106)
        self.assertEqual(len(self.h),323)
        self.assertEqual({k:len(x) for k,x in self.groups["2008"].items()},
                         {"ndni":15,"dni":12,"ndi":21,"di":19})

    def test_2008_reduced_herbivory_never_zero_good_seeds(self):
        r={a:v.group_stats(z) for a,z in self.groups["2008"].items()}
        self.assertEqual([r[a]["n_seed_bearing_mothers"] for a in v.ARMS],
                         [3,7,21,19])
        self.assertAlmostEqual(r["ndni"]["p_any_good_seeds"],.2)
        self.assertAlmostEqual(r["dni"]["p_any_good_seeds"],7/12)
        self.assertEqual(r["di"]["p_any_good_seeds"],1)
        self.assertAlmostEqual(r["ndni"]["mean_good_seeds_given_positive"],40/3)

    def test_exact_2008_fisher_and_non_support_of_damage_ambient(self):
        t=v.exact_seed_presence_tests(self.groups)
        self.assertAlmostEqual(t["2008"]["insecticide_intact"]["fisher_exact_two_sided_p"],
                               3.635121085883381e-7,places=12)
        self.assertAlmostEqual(t["2008"]["insecticide_damaged"]["fisher_exact_two_sided_p"],
                               .0046612638381270085,places=12)
        self.assertAlmostEqual(t["2008"]["damage_ambient"]["fisher_exact_two_sided_p"],
                               .05675057208237994,places=12)
        self.assertGreater(t["2008"]["damage_ambient"]["BH_q_across_eight_predefined_contrasts"],.05)
        self.assertLess(t["2008"]["insecticide_damaged"]["BH_q_across_eight_predefined_contrasts"],.05)
        self.assertEqual(t["2007"]["insecticide_intact"]["fisher_exact_two_sided_p"],1.0)

    def test_two_part_identity_without_causal_mediation(self):
        d=v.hurdle_identity(self.groups)["2008"]
        self.assertAlmostEqual(d["insecticide_intact"]["probability_any_seed_ratio"],5.)
        self.assertAlmostEqual(d["insecticide_intact"]["positive_mother_conditional_mean_seed_ratio"],
                               (2906/21)/(40/3))
        self.assertAlmostEqual(d["insecticide_intact"]["observed_mean_good_seed_ratio"],
                               (2906/21)/(40/15))
        self.assertFalse(d["insecticide_intact"]["reproductive_processes_are_causally_decomposed"])

    def test_baseline_size_complete_case_n_is_36_and_56(self):
        z=v.analyze(self.m,self.h,self.d)["rosette_complete_case_log1p_HC3_sensitivity"]
        self.assertEqual(z["2007"]["complete_cases_with_rosette_covariate"]["n_mothers"],36)
        self.assertEqual(z["2008"]["complete_cases_with_rosette_covariate"]["n_mothers"],56)
        self.assertEqual(z["2008"]["n_missing_or_invalid_rosette"],11)

    def test_size_adjusted_model_no_doubled_plant_head_sampling(self):
        z=v.analyze(self.m,self.h,self.d)
        q=z["rosette_complete_case_log1p_HC3_sensitivity"]["2008"]
        self.assertAlmostEqual(
            q["complete_cases_with_rosette_covariate"]["HC3_log1p_contrasts"]["insecticide_intact"]["difference_log1p_good_seeds"],
            4.0330906837633105,places=9)
        self.assertAlmostEqual(
            q["complete_cases_with_rosette_covariate"]["HC3_log1p_contrasts"]["insecticide_damaged"]["difference_log1p_good_seeds"],
            3.0669436961496075,places=9)
        self.assertGreater(
            q["complete_cases_with_rosette_covariate"]["HC3_log1p_contrasts"]["insecticide_damaged"]["normal_approximate_95_interval"][0],
            0)
        self.assertLess(
            q["complete_cases_with_rosette_covariate"]["HC3_log1p_contrasts"]["damage_ambient"]["normal_approximate_95_interval"][0],
            0)

    def test_no_flowered_head_when_apical_positive_is_semantic_not_zero_seed(self):
        q=v.audit_head_semantics(self.m,self.h)
        self.assertEqual(q["n_positive_mothers_no_reported_flowered_heads"],10)
        self.assertEqual(q["n_positive_mothers_zero_reported_seed_bearing_heads"],17)
        self.assertTrue(q["all_such_positive_seed_cases_are_apical_only"])
        self.assertFalse(q["no_flw_and_n_fl_sds_count_full_flowering_capitula_confirmed"])
        self.assertFalse(q["these_rows_are_proven_data_errors"])

    def test_12_head_treatment_conflicts_are_not_mutated(self):
        q=v.audit_head_semantics(self.m,self.h)
        self.assertEqual(q["n_head_rows_with_treatment_disagreement"],12)
        self.assertEqual(q["n_head_rows_matching_original_mothers"],309)
        self.assertTrue(q["original_head_treatment_discrepancies_are_not_autocorrected"])

    def test_rosette_na_is_missing_not_zero_or_excluded_from_whole_population(self):
        self.assertIsNone(v.value({"r.diam":"NA"},"r.diam"))
        self.assertIsNone(v.value({"r.diam":""},"r.diam"))
        with self.assertRaisesRegex(ValueError,"SIZE_MISSING_IN_ADJUSTED_ANALYSIS"):
            v.fit_log1p_hc3([dict(z, **{"r.diam":"NA"}) if i==0 else z
                             for i,z in enumerate(self.m[:20])],True)
        self.assertEqual(len(self.m),106)

    def test_numerical_fisher_corner(self):
        self.assertAlmostEqual(v.fisher_two_sided(1,2,1,3),1.)
        self.assertAlmostEqual(v.fisher_two_sided(3,0,0,3),.1)
        self.assertEqual(v.fisher_two_sided(0,3,0,7),1)
        with self.assertRaisesRegex(ValueError,"INVALID_FISHER"):
            v.fisher_two_sided(-1,2,1,3)

    def test_invalid_treatment_is_not_reinterpreted(self):
        s=copy.deepcopy(self.m)
        s[0]["trt"]="newarm"
        with self.assertRaisesRegex(ValueError,"SOURCE_TREATMENT_NOT_VERIFIED"):
            v.source_groups(s)

    def test_no_life_history_or_causal_trait_promotion(self):
        z=v.analyze(self.m,self.h,self.d)
        self.assertTrue(z["scientific_limits"]["year_group_treatment_populations_are_not_matched_repeated_mothers"])
        self.assertTrue(z["scientific_limits"]["rosette_measurement_timing_is_not_verified_pre_treatment"])
        self.assertTrue(z["scientific_limits"]["independent_C_canescens_original_134_mother_rows_unavailable"])
        self.assertEqual(z["scientific_limits"]["new_field_plants_or_heads_collected"],0)


if __name__=="__main__":
    unittest.main()
