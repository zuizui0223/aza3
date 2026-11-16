"""Unit checks for prospective inference limits; all fixture values are synthetic."""
import unittest
from analysis.phd_form_function_equivalence_gate_v1 import (
    ENDPOINTS, REQUIRED_GATES, classify_evidence, precision_halfwidth, demo
)

class FormFunctionGateTests(unittest.TestCase):
    def setUp(self):
        self.gates = dict.fromkeys(REQUIRED_GATES, True)
        self.margins = dict.fromkeys(ENDPOINTS, .1)
    def classify(self, intervals):
        return classify_evidence(gates=self.gates,
            distinct_morphologies=True, margins=self.margins,
            intervals=intervals, **kwargs)
    def test_perfect_identical_responses_with_narrow_intervals_can_pass(self):
        result=self.classify({k:(-.04,.04) for k in ENDPOINTS})
        self.assertEqual(result["status"], "CONDITIONAL_MULTI_ENDPOINT_EQUIVALENCE_COMPATIBLE")
    def test_same_final_output_with_distinct_pathways_is_not_equivalence(self):
        result=self.classify({
            "conspecific_pollen":(.2,.3),
            "direct_seed_damage":(.20,.30),
            "filled_achenes":(-.04,.04)})
        self.assertEqual(result["status"],"SAME_NET_FITNESS_COMPONENT_DIFFERENT_PATHWAYS")
    def test_nonsignificant_wide_effect_is_inconclusive(self):
        result=self.classify({k:(-.3,.3) for k in ENDPOINTS})
        self.assertEqual(result["status"],"INDETERMINATE_NOT_EQUIVALENT_BY_DEFAULT")
    def test_no_joint_records_fail_closed(self):
        bad=dict.fromkeys(REQUIRED_GATES,False)
        result=classify_evidence(gates=bad,distinct_morphologies=True,
          margins=self.margins,intervals={k:(-.01,.01) for k in ENDPOINTS})
        self.assertEqual(result["status"],"NOT_IDENTIFIABLE")
        self.assertEqual(len(result["missing_gates"]),len(REQUIRED_GATES))
    def test_not_morphologically_distinct_no_many_to_one(self):
        result=classify_evidence(gates=self.gates,distinct_morphologies=False,
          margins=self.margins,intervals={k:(-.01,.01) for k in ENDPOINTS})
        self.assertEqual(result["status"],"NO_DISTINCT_MORPHOLOGY_CONTRAST")
    def test_missing_endpoint_and_invalid_confidence(self):
        result=self.classify({"conspecific_pollen":(-.01,.01)})
        self.assertEqual(result["status"],"NOT_IDENTIFIABLE")
        with self.assertRaises(ValueError):
            precision_halfwidth(1,.4)
        with self.assertRaises(ValueError):
            precision_halfwidth(24,0)
        with self.assertRaises(ValueError):
            self.classify({k:(.3,-.3) for k in ENDPOINTS})
    def test_ci_precision_decreases_with_independent_mothers(self):
        self.assertGreater(precision_halfwidth(24,.4),precision_halfwidth(60,.4))
        self.assertGreater(precision_halfwidth(60,.4),precision_halfwidth(120,.4))
        self.assertTrue(demo()["not_empirical"])
        self.assertFalse(demo()["empirical_fitness_equivalence_authorized"])

if __name__=="__main__":
    unittest.main()
