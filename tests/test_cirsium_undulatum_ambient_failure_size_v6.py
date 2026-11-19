"""Exploratory exact source-level check: is maternal seed failure just plant size?

Actual 2008 Cirsium undulatum USDA mother and head source. Does NOT prove
causal maternal height, pollutant/pollinator mechanism or SLK specialization.
All six features and tail test are fixed in the source-specific receipt.
"""
import math
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from audit_cirsium_undulatum_ambient_failure_size_v6 import (
    audit_ambient_size, exact_stratified_binary,
    bh_qvalues_six, source_ambient_mothers,
)


class OriginalAmbientFailureSizeV6(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=audit_ambient_size(strict=True)

    def test_exact_original_source_cohort_and_mother_grain(self):
        s=self.report["source"]
        self.assertEqual(s["original_mothers_all_years"],106)
        self.assertEqual(s["original_heads_all_years"],323)
        self.assertEqual(s["complete_2008_mothers"],67)
        self.assertEqual(s["ambient_insect_2008_mothers"],27)
        self.assertEqual(s["ambient_apical_intact"],15)
        self.assertEqual(s["ambient_apical_damaged"],12)
        self.assertEqual(s["any_seed_positive"],10)
        self.assertEqual(s["both_seed_routes_failed"],17)
        self.assertEqual(s["original_source_id_prefixes"],["ES","SD","WG","WR"])

    def test_plant_height_ambient_mean_difference_is_real_but_not_claimed_causal(self):
        r=self.report["six_exploratory_covariates"]["plant_height_cm"]
        self.assertEqual((r["n_measured_mothers"],r["seed_positive"],
                          r["seed_zero"]),(27,10,17))
        self.assertAlmostEqual(r["mean_in_zero_seed_mothers"],46.294117647058826)
        self.assertAlmostEqual(r["mean_in_seed_positive_mothers"],53.8)
        self.assertAlmostEqual(r["mean_difference_positive_minus_zero"],7.505882352941174)
        self.assertEqual(r["exact_stratified_permutation"]["total_assignments"],360360)
        self.assertAlmostEqual(r["exact_stratified_permutation"]["p_two_sided"],.012218337218337219,places=13)
        self.assertAlmostEqual(r["q_BH_six_exploratory_covariates"],.07331002331002331,places=13)
        self.assertGreater(r["q_BH_six_exploratory_covariates"],.05)

    def test_non_height_variables_fail_5_percent_screen(self):
        six=self.report["six_exploratory_covariates"]
        self.assertEqual(set(six),{
            "plant_height_cm","basal_rosette_diameter_cm",
            "flower_bud_count","branch_count",
            "mean_head_damage_fraction","scored_lateral_head_count"})
        for k,r in six.items():
            if k!="plant_height_cm":
                self.assertGreater(r["exact_stratified_permutation"]["p_two_sided"],.15)
                self.assertGreater(r["q_BH_six_exploratory_covariates"],.05)
        self.assertEqual(six["basal_rosette_diameter_cm"]["n_measured_mothers"],24)
        self.assertEqual(six["plant_height_cm"]["n_measured_mothers"],27)
        self.assertEqual(six["flower_bud_count"]["n_measured_mothers"],27)
        self.assertEqual(six["mean_head_damage_fraction"]["n_measured_mothers"],27)

    def test_source_prefix_leave_one_out_flags_height_instability(self):
        d=self.report["height_sensitivity"]
        self.assertAlmostEqual(d["leave_out_source_prefix"]["WR"]["p_two_sided"],.2945165945165945,places=12)
        self.assertAlmostEqual(d["leave_out_source_prefix"]["SD"]["p_two_sided"],.0659090909090909,places=12)
        self.assertAlmostEqual(d["leave_out_source_prefix"]["ES"]["p_two_sided"],.0007215007215007215,places=12)
        self.assertAlmostEqual(d["leave_out_source_prefix"]["WG"]["p_two_sided"],.010833333333333334,places=12)
        self.assertFalse(self.report["interpretive_ceiling"]["source_prefixes_are_verified_field_sites"])

    def test_highest_successful_height_is_influential(self):
        d=self.report["height_sensitivity"]["leave_out_highest_successful_height"]
        self.assertEqual(d["source_plant_identifier"],"WR51")
        self.assertEqual(d["height_cm"],77)
        self.assertEqual(d["n_mothers"],26)
        self.assertAlmostEqual(d["p_two_sided"],.07473082473082474,places=12)
        self.assertGreater(d["p_two_sided"],.05)

    def test_exact_integer_stratification_does_not_use_asymptotic_fake_precision(self):
        z=[
          {"arm":"a","positive":1,"value":"1"},
          {"arm":"a","positive":0,"value":"0"},
          {"arm":"b","positive":1,"value":"2"},
          {"arm":"b","positive":0,"value":"0"},
        ]
        r=exact_stratified_binary(z,["arm"])
        self.assertEqual(r["total_assignments"],4)
        self.assertEqual(r["p_two_sided"],.5)
        self.assertEqual(r["n_informative_strata"],2)
        self.assertAlmostEqual(r["observed_centered_positive_sum"],1.5)

    def test_exact_prefix_and_treatment_conditioning_is_not_site_proof(self):
        d=self.report["height_sensitivity"]["source_prefix_and_apical_treatment_conditional"]
        self.assertEqual(d["total_assignments"],432)
        self.assertEqual(d["n_informative_strata"],4)
        self.assertAlmostEqual(d["p_two_sided"],.023148148148148147,places=12)
        self.assertEqual(d["uninformative_stratum_mothers"],12)

    def test_source_flowering_field_excluded_due_existing_semantic_conflict(self):
        self.assertEqual(self.report["source"]["no_flw_not_in_six_covariates"],
                         "ORIGINAL_ZERO_FLOWERED_BUT_POSITIVE_APICAL_SEED_CONFLICT")
        self.assertTrue(self.report["source"]["not_seven_posthoc_favorable_covariates"])
        self.assertFalse(self.report["interpretive_ceiling"]["flowering_fields_are_complete_verified_anthesis"])

    def test_synthetic_missingness_stratum_and_zero_variance_fail_closed(self):
        z=[
          {"arm":"a","positive":1,"value":"NA"},
          {"arm":"a","positive":1,"value":"1"},
          {"arm":"a","positive":0,"value":"0"},
          {"arm":"b","positive":1,"value":"3"},
          {"arm":"b","positive":0,"value":"1"},
        ]
        r=exact_stratified_binary(z,["arm"])
        self.assertEqual(r["n_measured_mothers"],4)
        self.assertEqual(r["n_missing"],1)
        self.assertEqual(r["total_assignments"],4)
        with self.assertRaisesRegex(ValueError,"INVALID_STRATUM"):
            exact_stratified_binary([{"arm":"a","positive":1,"value":"2"}],["arm"])
        with self.assertRaisesRegex(ValueError,"INVALID_VALUE"):
            exact_stratified_binary([{"arm":"a","positive":1,"value":"malformed"},
                {"arm":"a","positive":0,"value":"0"}],["arm"])

    def test_multiple_testing_and_biological_boundaries(self):
        q=bh_qvalues_six({
            "x":.01,"y":.02,"z":.4,"q":.5,"r":.7,"s":.9
        })
        self.assertAlmostEqual(q["x"],.06)
        self.assertAlmostEqual(q["y"],.06)
        for claim,value in self.report["interpretive_ceiling"].items():
            self.assertFalse(value,claim)
        self.assertEqual(self.report["scientific_conclusion"],
            "AMBIENT_MATERNAL_SEED_FAILURE_NOT_DETERMINISTICALLY_EXPLAINED_BY_SIZE_OR_SCORED_DAMAGE_NO_CAUSE_IDENTIFIED")


if __name__=="__main__":
    unittest.main()
