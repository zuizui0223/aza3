"""Source-origin and false-causal-selection regression tests (no raw copy)."""
import json
import math
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_original_ipomopsis_two_stage_selection_v1 as v
RECEIPT=ROOT/"data/evidence/ipomopsis_two_stage_trait_selection_independent_original_js_v1.json"

class RealIndependentOutgroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt=json.loads(RECEIPT.read_text(encoding="utf-8"))

    def test_author_archive_source_and_life_stage_definitions(self):
        p=self.receipt["source"]
        self.assertEqual(p["doi"],"10.1086/716740")
        self.assertEqual(p["open_data_doi"],"10.7280/D1KM49")
        self.assertEqual(p["metadata_author_md5"],v.EXPECTED_MD5)
        self.assertTrue(p["published_author_data_not_an_aza3_field_study"])
        self.assertTrue(self.receipt["limits"]["Ipomopsis_not_Cirsium"])
        self.assertTrue(self.receipt["limits"]["flyescape_measures_fraction_of_fruits_avoiding_predation"])

    def test_original_437_rows_and_absent_plant_id_preserved(self):
        s=self.receipt["audit"]
        self.assertEqual(s["n_original_data_records"],437)
        self.assertEqual(s["n_original_fields"],12)
        self.assertEqual(s["source_header"],v.COLS)
        self.assertFalse(s["source_plant_identifier_field_exists"])
        self.assertEqual(s["source_exact_duplicate_full_rows"],0)
        self.assertEqual(s["site_counts"],{"agg":221,"hyb":144,"ten":72})
        self.assertEqual(s["type_counts"],{"H":246,"A":175,"T":16})
        self.assertEqual(sum(s["year_counts"].values()),437)

    def test_original_life_stage_missingness_is_not_zero(self):
        s=self.receipt["audit"]
        m=s["missing_count_by_column"]
        self.assertEqual(
            (m["width"],m["sepal"],m["seeds"],m["seedsinit"],m["flyescape"]),
            (18,329,32,28,51))
        self.assertEqual(s["seeds_greater_than_seedsinit"],41)
        self.assertTrue(self.receipt["limits"][
            "staged_outcomes_related_not_seeds_equals_initiated_times_flyescape"])

    def test_original_complete_cohort_within_real_source_strata(self):
        z=self.receipt["analyses"]
        w=z["width"]
        s=z["sepal"]
        self.assertEqual(w["summary"]["n_source_complete_records"],367)
        self.assertEqual(w["summary"]["n_within_group_comparable_records"],354)
        self.assertEqual(w["summary"]["n_groups"],34)
        self.assertEqual(s["summary"]["n_source_complete_records"],96)
        self.assertEqual(s["summary"]["n_within_group_comparable_records"],90)
        self.assertEqual(s["summary"]["n_groups"],15)
        self.assertEqual(sum(x["n"] for x in w["design"]["source_eligible_strata"]),354)
        self.assertEqual(sum(x["n"] for x in s["design"]["source_eligible_strata"]),90)

    def test_raw_opposite_correlations_at_same_original_width_rows(self):
        w=self.receipt["analyses"]["width"]
        r=w["summary"]["raw_correlations"]
        self.assertAlmostEqual(r["width_or_sepal_vs_log1p_initiated"],
                               .28283346836233314,places=11)
        self.assertAlmostEqual(r["width_or_sepal_vs_fly_escape"],
                               -.14404478690024206,places=11)
        self.assertGreater(r["width_or_sepal_vs_log1p_initiated"],0)
        self.assertLess(r["width_or_sepal_vs_fly_escape"],0)

    def test_conditioning_on_site_year_and_plant_type_shrinks_assoc(self):
        w=self.receipt["analyses"]["width"]
        self.assertAlmostEqual(w["r_init"],.05349503292406693,places=11)
        self.assertAlmostEqual(w["r_escape"],-.028808559820120226,places=11)
        self.assertLess(abs(w["r_init"]),
                        abs(w["summary"]["raw_correlations"]
                            ["width_or_sepal_vs_log1p_initiated"]))
        self.assertLess(abs(w["r_escape"]),
                        abs(w["summary"]["raw_correlations"]
                            ["width_or_sepal_vs_fly_escape"]))
        self.assertEqual(w["design"]["permutations"],10000)
        self.assertEqual(w["design"]["bootstrap_draws"],3000)

    def test_conditional_permutation_and_uncertainty_do_not_force_a_sign(self):
        a=self.receipt["analyses"]
        w=a["width"]
        self.assertAlmostEqual(
            w["permutation_p_two_sided"]["init"],.3347665233476652)
        self.assertAlmostEqual(
            w["permutation_p_two_sided"]["escape"],.6106389361063893)
        self.assertAlmostEqual(
            w["permutation_p_two_sided"]["stage_difference"],.33266673332666735)
        for k in ("init","escape","stage_difference"):
            interval=w["within_stratum_row_resample_ci95"][k]
            self.assertLess(interval[0],0)
            self.assertGreater(interval[1],0)
        self.assertGreater(a["four_analyses_multiplicity_BH_q"]["width_init"],.05)

    def test_sepal_external_flower_proxy_has_different_missingness(self):
        a=self.receipt["analyses"]
        s=a["sepal"]
        self.assertAlmostEqual(s["r_init"],-.005002846883399559,places=11)
        self.assertAlmostEqual(s["r_escape"],-.2079573761878985,places=11)
        self.assertAlmostEqual(
            s["permutation_p_two_sided"]["escape"],.0826917308269173)
        self.assertAlmostEqual(
            a["four_analyses_multiplicity_BH_q"]["sepal_escape"],
            .3307669233076692)
        self.assertTrue(s["within_stratum_row_resample_ci95"]["escape"][1]>0)

    def test_original_upper_site_merging_preserves_weak_associations(self):
        a=self.receipt["analyses"]
        m=a["merged_actual_upper_site_sensitivity"]["width"]
        self.assertEqual((m["n_within_group_comparable_records"],
                          m["n_groups"]),(358,28))
        self.assertAlmostEqual(
            m["within_group_correlations"]["trait_vs_log1p_initiated"],
            .06392148857280598,places=11)
        self.assertAlmostEqual(
            m["within_group_correlations"]["trait_vs_fly_escape"],
            -.030951879296620827,places=11)

    def test_one_year_removed_makes_predation_stage_direction_unstable(self):
        a=self.receipt["analyses"]["leave_one_year_out_width"]
        self.assertEqual(len(a),9)
        self.assertGreater(
            a["2017"]["within_group_correlations"]["trait_vs_fly_escape"],0)
        self.assertLess(
            a["2018"]["within_group_correlations"]["trait_vs_fly_escape"],0)
        self.assertGreater(
            a["2017"]["within_group_correlations"]["trait_vs_log1p_initiated"],0)

    def test_authentic_source_md5_requires_exact_binary_not_html_challenge(self):
        with self.assertRaisesRegex(ValueError,"ORIGINAL_ZENODO_MD5_MISMATCH"):
            v.load(b"Site,Year,Type,width\nfake\n")

    def test_synthetic_degenerate_group_does_not_impute_artificial_correlation(self):
        rr=[
            {"Site":"agg","Year":2015,"Type":"A","width":2.0,
             "seedsinit":5,"flyescape":.8},
            {"Site":"agg","Year":2015,"Type":"A","width":2.0,
             "seedsinit":9,"flyescape":.9},
            {"Site":"ten","Year":2015,"Type":"T","width":4.0,
             "seedsinit":3,"flyescape":.1}
        ]
        source,groups=v.source_groups(rr,"width")
        self.assertEqual(len(source),3)
        self.assertEqual(len(groups),0)
        self.assertIsNone(v.pearson([1,1,1,1],[1,2,3,4]))
        with self.assertRaisesRegex(ValueError,"INVALID_OR_ZERO_XORSHIFT_SEED"):
            next(v.xorshift32(0))

    def test_trial_real_sites_not_independent_cirsium_guild_expansions(self):
        lim=self.receipt["limits"]
        self.assertTrue(lim["original_source_has_no_unique_plant_identifier"])
        self.assertTrue(lim["traits_are_corolla_and_sepal_width_not_Cirsium_spines"])
        self.assertTrue(lim["selection_model_after_original_article_not_reconstructed"])
        self.assertTrue(lim["permutation_exchangeability_not_proven"])
        self.assertTrue(lim["analyses_exploratory_after_source_selection"])
        self.assertEqual(lim["new_own_field_observations"],0)

    def test_science_not_an_ephemeral_single_stage_binary_tradeoff(self):
        z=self.receipt["analyses"]["width"]
        self.assertGreater(z["r_init"]-z["r_escape"],0)
        self.assertGreater(z["permutation_p_two_sided"]["stage_difference"],.05)
        self.assertTrue(self.receipt["limits"][
            "observational_within_source_strata_not_true_trait_randomized_intervention"])

if __name__=="__main__":
    unittest.main()
