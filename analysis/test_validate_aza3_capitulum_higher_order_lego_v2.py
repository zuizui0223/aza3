"""Tests on aza3 outcome-blind module sampling design: no own plant data."""
from fractions import Fraction
from math import isclose,sqrt
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from validate_aza3_capitulum_higher_order_lego_v2 import (
    feasible_triple_bounds,no_complete_detection_probability,design_effect,
    exact_zero_success_two_sided_upper,run,
)


class Aza3LegoProspectiveTests(unittest.TestCase):
    def test_canonical_role_and_all_evidence_unmeasured(self):
        result=run()
        self.assertEqual(result["study_population"],
                         "Cirsium sieboldii; planned Phase-A n=40, not sampled for this analysis")
        self.assertEqual(result["registered_individuals_qualified_for_joint_resource_inference"],0)
        self.assertFalse(result["all_three_functional_resources_measured"])
        self.assertFalse(result["proven_genetic_LEGO_reuse"])
        self.assertFalse(result["new_biological_data_admitted"])
        self.assertFalse(result["validated_Ryukyu_flower_regain"])
        self.assertFalse(result["fitness_adaptation_demonstrated"])

    def test_sharp_pairwise_frechet_interval_in_design(self):
        self.assertEqual(feasible_triple_bounds(
            (Fraction(1,2),)*3,
            {(0,1):Fraction(1,4),(0,2):Fraction(1,4),(1,2):Fraction(1,4)}
        ),(Fraction(0),Fraction(1,4)))
        self.assertEqual(run()["formal_anto_moment_bound"],["0","1/4"])

    def test_40_independent_samples_do_not_guarantee_a_functional_recombination(self):
        x=run()["hypothetical_n40_draw_scenarios"]
        by_p={r["hypothetical_p_complete"]:r for r in x}
        self.assertEqual(by_p[0]["expected_complete_inventory_count_if_n_independent"],0)
        self.assertEqual(by_p[.125]["expected_complete_inventory_count_if_n_independent"],5)
        self.assertEqual(by_p[.25]["expected_complete_inventory_count_if_n_independent"],10)
        self.assertAlmostEqual(by_p[.125]["probability_zero_complete_in_40_independent_draws"],
                               (.875)**40)

    def test_null_observations_are_not_proof_zero_probability(self):
        upper=exact_zero_success_two_sided_upper(40)
        self.assertGreater(upper,.08)
        self.assertLess(upper,.09)

    def test_kinship_sensitivity_not_measured_icc(self):
        x=design_effect(40,10,.2)
        self.assertAlmostEqual(x["design_effect"],2.8)
        self.assertAlmostEqual(x["illustrative_effective_n"],40/2.8)
        r=run()["cluster_assumption_sensitivity"]
        self.assertEqual(len(r),4)
        self.assertTrue(all("heuristic" in z["warning"] for z in r))
        self.assertEqual(r[0]["illustrative_effective_n"],40)

    def test_no_observation_level_power_inference(self):
        r=run()
        self.assertFalse(r["real_empirical_power_study_completed"])
        self.assertIn("not sampled",r["study_population"])
        self.assertTrue(all(x["not_fitted"] for x in r["hypothetical_n40_draw_scenarios"]))

    def test_bad_probabilities_and_clustering_rejected(self):
        with self.assertRaises(ValueError):
            feasible_triple_bounds(
                (Fraction(1,10),)*3,
                {(0,1):Fraction(1,5),(0,2):Fraction(1,5),(1,2):Fraction(1,5)}
            )
        with self.assertRaises(ValueError):
            no_complete_detection_probability(-.1,40)
        with self.assertRaises(ValueError):
            no_complete_detection_probability(.2,0)
        with self.assertRaises(ValueError):
            design_effect(40,6,.2)
        with self.assertRaises(ValueError):
            design_effect(40,10,1.1)
        with self.assertRaises(ValueError):
            exact_zero_success_two_sided_upper(0)


if __name__=="__main__":
    unittest.main()
