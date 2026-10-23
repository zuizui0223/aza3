"""Fail-close exact conditional maternal reproductive-route covariance interpretation."""
import copy
import pathlib
import sys
import unittest

ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_undulatum_route_covariance_v1 as v
SOURCE=ROOT/"data/external/usda_wavyleaf_2007_2008/mother_plants.csv"


class TestMeasuredRouteOverlap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original=v.load_source(SOURCE)
        cls.coherent,cls.bad=v.classify_original(cls.original)

    def test_106_plant_source_remains_105_for_joint_route_only(self):
        self.assertEqual(len(self.original),106)
        self.assertEqual(len(self.coherent),105)
        self.assertEqual(self.bad,["2007|327"])
        self.assertEqual(v.SOURCE_GIT_BLOB,
                         "458800272c3d0a03c7653702027a8f8578d9e88c")

    def test_two_routes_have_four_observed_outcome_states(self):
        a=v.exact_route_analysis(self.coherent)
        self.assertEqual((a["pooled_ignoring_year_and_treatment"]["apical_positive"],
                          a["pooled_ignoring_year_and_treatment"]["lateral_positive"],
                          a["pooled_ignoring_year_and_treatment"]["both_positive"],
                          a["pooled_ignoring_year_and_treatment"]["neither_positive"]),
                         (53,68,35,19))
        self.assertEqual(a["eight_original_year_by_treatment_strata"][5]["route_counts"],
            {"neither":5,"apical_only":0,"lateral_only":7,"both":0})

    def test_hypergeom_exact_probabilities_no_random_simulations(self):
        x=v.distribution(4,2,2)
        self.assertAlmostEqual(x[0],1/6)
        self.assertAlmostEqual(x[1],4/6)
        self.assertAlmostEqual(x[2],1/6)
        self.assertAlmostEqual(sum(x.values()),1)
        with self.assertRaisesRegex(ValueError,"INVALID_HYPERGEOMETRIC_MARGINAL"):
            v.distribution(4,5,2)

    def test_convolution_is_normalized_and_matches_manual(self):
        a=v.convolve({0:.5,1:.5},{0:.5,1:.5})
        self.assertEqual(a,{0:.25,1:.5,2:.25})
        with self.assertRaisesRegex(ValueError,"NULL_CONVOLUTION_MASS_DRIFT"):
            v.convolve({0:.6},{0:.5,1:.5})

    def test_real_sign_reversal_is_descriptive_not_sig_pledge(self):
        a=v.exact_route_analysis(self.coherent)
        unadjusted=a["pooled_ignoring_year_and_treatment"]
        adjusted=a["stratified_conditional_marginal_independence"]
        self.assertGreater(unadjusted["observed_minus_independence_expected_both"],0)
        self.assertLess(adjusted["observed_minus_expected_both"],0)
        self.assertAlmostEqual(unadjusted["pooled_odds_ratio_both_vs_mutual_failure"],
                               (35*19)/(18*33),places=13)
        self.assertAlmostEqual(adjusted["observed_minus_expected_both"],
                              -2.87410951621478,places=12)
        self.assertAlmostEqual(adjusted["exact_two_inclusive_tails_doubled_min"],
                               0.09078087375300997,places=12)
        self.assertTrue(a["scientific_claim_limits"][
            "negative_covariance_direction_is_not_significant_at_doubled_tail_0p05"])

    def test_one_sided_045_is_not_two_sided_045(self):
        a=v.exact_route_analysis(self.coherent)
        c=a["stratified_conditional_marginal_independence"]
        self.assertAlmostEqual(c["exact_conditional_lower_tail_inclusive"],
                               .045390436876504986,places=13)
        self.assertAlmostEqual(c["exact_conditional_upper_tail_inclusive"],
                               .9947827084049995,places=13)
        self.assertGreater(c["exact_two_inclusive_tails_doubled_min"],.05)
        self.assertTrue(c["independence_is_statistical_conditional_reference_not_measured"])

    def test_year_specific_one_off_not_a_second_species_replication(self):
        a=v.exact_route_analysis(self.coherent)["year_conditioned_sensitivity"]
        self.assertEqual(a["2007"]["n_mothers"],38)
        self.assertEqual(a["2008"]["n_mothers"],67)
        self.assertAlmostEqual(a["2007"]["exact_two_inclusive_tails_doubled_min"],
                               .2474747474747475,places=12)
        self.assertAlmostEqual(a["2008"]["exact_two_inclusive_tails_doubled_min"],
                               .4919160646714827,places=12)

    def test_original_2007_discrepancy_not_imputed_as_lateral_count(self):
        altered=copy.deepcopy(self.original)
        z=next(r for r in altered if r["year"]+"|"+r["plt"]=="2007|327")
        z["o.hds"]="20"
        with self.assertRaisesRegex(ValueError,
              "UNRESOLVED_SOURCE_MOTHER_COMPONENT_DISCREPANCY"):
            v.classify_original(altered)

    def test_marginals_and_strata_frozen(self):
        a=copy.deepcopy(self.coherent)
        z=next(r for r in a.values() if r["year"]=="2008" and r["arm"]=="di" and r["state"]=="both")
        z["state"]="apical_only"
        with self.assertRaisesRegex(ValueError,"ORIGINAL_FOUR_STATE_MOTHER_CELLS_CHANGED"):
            v.exact_route_analysis(a)

    def test_no_selection_or_genetic_modularity_claim(self):
        a=v.exact_route_analysis(self.coherent)
        x=a["scientific_claim_limits"]
        self.assertEqual(x["newly_observed_individual_plants_or_heads"],0)
        self.assertTrue(x["head_pair_null_does_not_establish_selected_bet_hedging"])
        self.assertTrue(x["no_independent_oviposition_or_effective_pollination_at_same_head"])
        self.assertTrue(x["two_observed_routes_are_not_lost_ancestral_genes"])


if __name__=="__main__":
    unittest.main()
