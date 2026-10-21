"""Tests for logically sharp head-risk bounds; synthetic scenarios are not field estimates."""
from pathlib import Path
import copy
import json
import math
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT / "analysis"))
import audit_cirsium_pitcheri_head_risk_sharp_bounds_v1 as q

SOURCE=ROOT / "data/evidence/cirsium_pitcheri_2020_multi_enemy_rank_source_v1.json"


class TestCirsiumHeadInsuranceBounds(unittest.TestCase):
    def setUp(self):
        self.source=json.loads(SOURCE.read_text(encoding="utf-8"))

    def test_exact_frechet_pair(self):
        out=q.hypothetical_one_per_rank([.27,.50])
        self.assertEqual(out["joint_all_infested_sharp_bounds"],[0,.27])
        self.assertEqual(out["at_least_one_free_sharp_bounds"],[.73,1.0])
        self.assertAlmostEqual(out["independent_reference_at_least_one_free"],.865)
        self.assertAlmostEqual(out["independent_reference_extra_over_terminal_alone"],.135)

    def test_three_head_pair_all_attacked_overlap(self):
        out=q.hypothetical_one_per_rank([.27,.50,.80])
        self.assertEqual(out["joint_all_infested_sharp_bounds"],[0,.27])
        self.assertEqual(out["at_least_one_free_sharp_bounds"],[.73,1.0])
        self.assertAlmostEqual(out["independent_reference_at_least_one_free"],.892)

    def test_terminal_secondary_joint_can_span_perfect_redundancy_to_zero_overlap(self):
        # 27% of heads infested at rank T and 50% at rank S:
        # any coupling satisfying those marginals is admissible in the hypothetical.
        self.assertAlmostEqual(q.hypothetical_one_per_rank([.27,.50])[
            "sharp_extra_chance_over_terminal_alone"][0],0)
        self.assertAlmostEqual(q.hypothetical_one_per_rank([.27,.50])[
            "sharp_extra_chance_over_terminal_alone"][1],.27)

    def test_terminal_absence_2017_has_no_joint_secondary_seed_data(self):
        outcome=q.analyze(self.source)
        for site,p in (("SC",.27),("WFD",.33)):
            x=outcome["observed_terminal_absence_and_untested_recovery"][site]
            self.assertEqual(x["missing_terminal_fraction"],p)
            self.assertEqual(x["joint_missing_terminal_and_successful_secondary_sharp_bound_without_secondary_marginal"],[0,p])
            self.assertFalse(x["has_original_plant_level_joint_seed_counts"])
        self.assertTrue(outcome["inference_limit"]["no_plant_specific_compensation_claim"])

    def test_scenario_does_not_create_fecundity_finding(self):
        s=q.equal_effort_mechanism_scenarios(.27,.50,.80,.60)
        self.assertAlmostEqual(s["remaining_mature_seed_proxy_if_weevil_loss_rate_rank_invariant"]["terminal"],.838)
        self.assertAlmostEqual(s["remaining_mature_seed_proxy_if_weevil_loss_rate_rank_invariant"]["secondary"],.70)
        self.assertAlmostEqual(s["remaining_mature_seed_proxy_if_weevil_loss_rate_rank_invariant"]["tertiary"],.52)
        self.assertGreater(s["remaining_mature_seed_proxy_if_weevil_loss_rate_rank_invariant"]["terminal"],
                           s["remaining_mature_seed_proxy_if_weevil_loss_rate_rank_invariant"]["secondary"])
        self.assertFalse(s["assumptions_confirmed_in_original_study"])

    def test_catastrophic_loss_break_even_threshold_depends_on_side_fecundity(self):
        s=q.equal_effort_mechanism_scenarios(.27,.50,.80,.60)
        threshold=s["whole_terminal_loss_fraction_break_even_for_secondary_equal_cost_resource_ratios"]
        self.assertAlmostEqual(threshold["1.0"],1-.70/.838)
        self.assertAlmostEqual(threshold["0.5"],1-.35/.838)
        self.assertGreater(threshold["0.5"],threshold["1.0"])
        self.assertFalse(s["terminal_absence_2017_can_be_plugged_into_2016_weevil_scenario_as_estimate"])

    def test_forbid_unverified_deer_cause(self):
        src=copy.deepcopy(self.source)
        src["terminal_head_absence_by_site_year"]["specific_herbivore_causation_identified"]=True
        with self.assertRaisesRegex(ValueError,"DEER_CAUSE_UNVERIFIED"):
            q.analyze(src)

    def test_wrong_weevil_marginal_does_not_silently_generate_new_result(self):
        src=copy.deepcopy(self.source)
        src["rank_weevil_2016_WFD"]["infestation_evidence_fractions_approx"]["terminal"]=.10
        with self.assertRaisesRegex(ValueError,"ORIGINAL_SOURCE_WEevil_RATE_CHANGED"):
            q.analyze(src)

    def test_no_raw_data_availed_still_true(self):
        src=copy.deepcopy(self.source)
        src["material_limit"]["original_plant_and_head_rows_not_downloaded"]=False
        with self.assertRaisesRegex(ValueError,"ORIGINAL_PLANT_HEAD_DATA_NOT_HELD"):
            q.analyze(src)

    def test_probability_domain_and_one_sample(self):
        with self.assertRaisesRegex(ValueError,"INVALID_MARGINAL"):
            q.hypothetical_one_per_rank([.27,1.01])
        with self.assertRaisesRegex(ValueError,"NEED_AT_LEAST_TWO_HEAD_RANKS"):
            q.hypothetical_one_per_rank([.27])

    def test_head_not_infested_is_not_seeded_reproduction(self):
        r=q.analyze(self.source)
        self.assertTrue(r["inference_limit"]["not_infested_does_not_imply_filled_achene_output"])
        self.assertTrue(r["inference_limit"]["rank_is_not_orientation_or_phyllary_geometry"])
        self.assertFalse(r["inference_limit"]["adaptive_rank_insurance_established"])
        self.assertEqual(r["n_new_independent_biological_field_observations"],0)


if __name__=="__main__":
    unittest.main()
