"""Original USDA maternal dual-route *coverage*, not an invented rescue process.

Each apical/lateral score refers to THE SAME mother Cirsium undulatum.
West & Louda (2021), AJB, 10.1002/ajb2.1744,
USDA 10.15482/USDA.ADC/1522649. Parent aza3 v2 is source locked.

The union advantage is algebraic counterfactual-free COVERAGE, not the
fitness consequence of organ removal, and not ancestral gene regain.
"""
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from audit_cirsium_undulatum_maternal_route_coverage_v3 import (
    coverage_from_counts, exact_discordance_test,
    audit_maternal_route_coverage, bh_correction
)


class ActualSourceMaternalCoverageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=audit_maternal_route_coverage()

    def test_real_original_USDA_nesting_and_bad_row_quarantine(self):
        s=self.report["source"]
        self.assertEqual(s["genuine_maternal_rows"],106)
        self.assertEqual(s["jointly_consistent_mothers"],105)
        self.assertEqual(s["original_individual_head_rows"],323)
        self.assertEqual(s["bad_joint_maternal_route_IDs"],["2007|327"])
        self.assertEqual(s["USDA_original_git_blob_mothers"],
                         "458800272c3d0a03c7653702027a8f8578d9e88c")

    def test_one_mother_four_way_seed_outcome_not_independent_heads(self):
        for year, arms in self.report["year_arm_observed_coverage"].items():
            for name, r in arms.items():
                counts=r["four_anatomical_routes"]
                self.assertEqual(r["n_mothers"],sum(counts.values()))
                self.assertAlmostEqual(
                    r["P_any_maternal_good_seed"],
                    r["P_apical_good_seed"]+r["P_lateral_good_seed"]
                    -r["P_both_apical_and_lateral"],places=14)
                self.assertAlmostEqual(
                    r["observed_union_surplus_over_better_single_route"],
                    min(counts["apical_only"],counts["lateral_only"])
                    /r["n_mothers"],places=14)
                self.assertGreaterEqual(r["observed_union_surplus_over_better_single_route"],0)
                self.assertFalse(r["organ_removal_causal_selection_estimated"])

    def test_2008_damage_insecticide_both_routes_fill_distinct_mother_failures(self):
        x=self.report["year_arm_observed_coverage"]["2008"]["di"]
        self.assertEqual(x["four_anatomical_routes"],
                         {"no_good_seed":0,"apical_only":4,
                          "lateral_only":8,"both_sites":7})
        self.assertEqual(x["n_mothers"],19)
        self.assertEqual(x["P_any_maternal_good_seed"],1)
        self.assertAlmostEqual(x["P_apical_good_seed"],11/19)
        self.assertAlmostEqual(x["P_lateral_good_seed"],15/19)
        self.assertAlmostEqual(x["observed_union_surplus_over_better_single_route"],4/19)
        self.assertAlmostEqual(x["observed_union_surplus_over_apical_route"],8/19)
        self.assertAlmostEqual(x["observed_union_surplus_over_lateral_route"],4/19)
        self.assertAlmostEqual(x["fraction_positive_mothers_using_exactly_one_route"],12/19)
        self.assertAlmostEqual(x["P_lateral_positive_given_apical_negative"],1)
        self.assertAlmostEqual(x["P_apical_positive_given_lateral_negative"],1)

    def test_2008_ambient_damage_joint_zeros_are_real_not_inaccessible_loci(self):
        x=self.report["year_arm_observed_coverage"]["2008"]["dni"]
        self.assertEqual(x["four_anatomical_routes"],
                         {"no_good_seed":5,"apical_only":0,
                          "lateral_only":7,"both_sites":0})
        self.assertAlmostEqual(x["P_joint_failure_both_routes"],5/12)
        self.assertAlmostEqual(x["P_lateral_positive_given_apical_negative"],7/12)
        self.assertAlmostEqual(x["observed_union_surplus_over_better_single_route"],0)
        self.assertFalse(self.report["biological_claim_limits"]["common_antagonist_proven_to_cause_both_failures"])

    def test_2008_intact_insecticide_source_success_and_rare_unique_lateral(self):
        x=self.report["year_arm_observed_coverage"]["2008"]["ndi"]
        self.assertEqual(x["four_anatomical_routes"],
                         {"no_good_seed":0,"apical_only":7,
                          "lateral_only":1,"both_sites":13})
        self.assertAlmostEqual(x["observed_union_surplus_over_better_single_route"],1/21)
        self.assertEqual(x["P_any_maternal_good_seed"],1)
        y=self.report["year_arm_observed_coverage"]["2008"]["ndni"]
        self.assertEqual(y["four_anatomical_routes"],
                         {"no_good_seed":12,"apical_only":1,
                          "lateral_only":1,"both_sites":1})
        self.assertAlmostEqual(y["P_joint_failure_both_routes"],12/15)

    def test_original_2007_groups_and_total_arithmetic(self):
        x=self.report["year_arm_observed_coverage"]["2007"]
        self.assertEqual([x[arm]["n_mothers"] for arm in
                          ("ndni","dni","ndi","di")],[9,8,10,11])
        self.assertEqual([x[arm]["four_anatomical_routes"]["lateral_only"]
                          for arm in ("ndni","dni","ndi","di")],[2,6,2,6])
        self.assertAlmostEqual(x["ndni"]["observed_union_surplus_over_better_single_route"],2/9)
        self.assertAlmostEqual(x["di"]["observed_union_surplus_over_better_single_route"],1/11)

    def test_exact_same_mother_discordance_mcnemar_not_false_independence(self):
        x=self.report["paired_route_discordance_exact_tests"]
        self.assertEqual(len(x),8)
        self.assertAlmostEqual(x["2008_di"]["p_two_sided"],0.3876953125)
        self.assertEqual(x["2008_di"]["number_of_discordant_mothers"],12)
        self.assertEqual(x["2008_di"]["apical_only_vs_lateral_only"],[4,8])
        self.assertGreater(x["2008_di"]["q_BH_eight"],.05)
        self.assertEqual(exact_discordance_test(0,0),1)
        self.assertAlmostEqual(exact_discordance_test(0,7),.015625)

    def test_known_ambient_damage_lateral_fisher_not_robust_after_multiplicity(self):
        p=self.report["original_alternative_route_presence_tests"]
        self.assertAlmostEqual(p["2008_damage_ambient"]["p_exact"],
                               0.03691838291380626,places=11)
        self.assertAlmostEqual(p["2008_damage_ambient"]["q_BH_eight"],
                               0.14767353165522504,places=11)
        self.assertLess(p["2008_insecticide_intact"]["q_BH_eight"],.05)

    def test_independence_counterfactual_is_not_a_measured_biological_process(self):
        x=self.report["year_arm_observed_coverage"]["2008"]["di"]
        self.assertAlmostEqual(
            x["conditional_independence_benchmark_P_union"],
            11/19+15/19-(11*15)/(19*19),places=14)
        self.assertGreater(x["P_any_maternal_good_seed"],
                           x["conditional_independence_benchmark_P_union"])
        self.assertFalse(x["conditional_independence_is_a_realistic_biological_null"])
        self.assertFalse(x["observed_correlated_routes_are_molecularly_independent"])

    def test_analytic_pure_cases_and_failure_preservation(self):
        z=coverage_from_counts({"no_good_seed":4,"apical_only":0,
                                "lateral_only":0,"both_sites":0})
        self.assertEqual(z["P_any_maternal_good_seed"],0)
        self.assertEqual(z["observed_union_surplus_over_better_single_route"],0)
        self.assertIsNone(z["P_lateral_positive_given_apical_negative"] if False else
                          coverage_from_counts({"no_good_seed":0,"apical_only":3,
                            "lateral_only":0,"both_sites":2})[
                              "P_lateral_positive_given_apical_negative"])
        with self.assertRaises(ValueError):
            coverage_from_counts({"no_good_seed":-1,"apical_only":1,
                                  "lateral_only":2,"both_sites":3})
        with self.assertRaises(ValueError):
            coverage_from_counts({"no_good_seed":0,"apical_only":0,
                                  "lateral_only":0,"both_sites":0})
        with self.assertRaises(ValueError):
            exact_discordance_test(-1,2)

    def test_bh_exact_multiplicity_is_deterministic(self):
        q=bh_correction({"x":0.03,"y":0.01,"z":.4})
        self.assertAlmostEqual(q["y"],.03)
        self.assertAlmostEqual(q["x"],.045)
        self.assertAlmostEqual(q["z"],.4)

    def test_no_free_genetic_memory_or_population_fitness_claims(self):
        r=self.report
        for name,value in r["biological_claim_limits"].items():
            self.assertFalse(value,name)
        self.assertEqual(r["conclusion"],
                         "OBSERVED_CONDITIONAL_MATERNAL_ROUTE_COVERAGE_NOT_CAUSAL_ARCHITECTURE_FITNESS_OR_NATIVE_REGAIN")
        self.assertTrue(r["biology"]["sibling_heads_are_present_before_this_season_not_evolutionarily_resurrected"])
        self.assertTrue(r["biology"]["same_mother_but_not_same_genetic_program_assayed"])


if __name__=="__main__":
    unittest.main()
