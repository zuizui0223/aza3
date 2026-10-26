"""Real published *Cirsium* seed access data and stage-throughput identities.

Synthetic numbers in tests exercise mathematically valid/invalid conditions;
synthetic values are never inserted as field observations or fitness.
"""
import copy
import json
import math
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_q2_stage_gate_transport_v1 as gate

SOURCE_1988=ROOT/"data/evidence/cirsium_vulgare_1988_original_multistage_source_tables_v1.json"
SOURCE_2022=ROOT/"data/evidence/cirsium_arvense_four_arm_joint_response_independent_js_v1.json"
FROZEN=ROOT/"data/evidence/cirsium_q2_source_stage_transport_threshold_independent_JS_v1.json"


class ConditionalHeadPayoffTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=json.loads(SOURCE_1988.read_text(encoding="utf-8"))
        cls.unlinked=json.loads(SOURCE_2022.read_text(encoding="utf-8"))
        cls.frozen=json.loads(FROZEN.read_text(encoding="utf-8"))
        cls.computed=gate.analyze(cls.source,cls.unlinked)

    def test_actual_original_1981_1982_cage_study_not_1984_same_heads(self):
        p=self.computed["source"]
        self.assertEqual(p["source_doi"],"10.2307/2260602")
        self.assertEqual(p["exclosure_original_years"],[1981,1982])
        self.assertEqual(p["head_stem_years"],[1982,1984])
        self.assertEqual(p["original_number_of_exclosure_habitat_date_cells"],10)
        self.assertTrue(self.computed["definitions"]["no_same_seed_ID_link_across_1988_cohorts"])

    def test_real_cage_vs_open_seed_survival(self):
        p=self.computed["source_aggregates"]
        self.assertAlmostEqual(p["mean_seed_retention_percentage"]["cage"],93.8)
        self.assertAlmostEqual(p["mean_seed_retention_percentage"]["net"],34.6)
        self.assertAlmostEqual(p["mean_seed_retention_percentage"]["open"],36.5)
        self.assertAlmostEqual(p["mean_matched_cage_minus_open_survival_percentage_points"],57.3)
        self.assertEqual(p["within_habitat_date_cage_retention_exceeds_open"],10)
        self.assertEqual(p["within_habitat_date_cage_retention_exceeds_net"],10)
        self.assertEqual(p["n_source_cells_where_bird_net_retention_is_below_open"],6)

    def test_original_habitat_date_retention_ratio_source_consistency(self):
        p=self.computed["source_original_exclosure_cell_replay"]
        self.assertEqual(len(p),10)
        self.assertEqual(len({r["date"] for r in p}),5)
        self.assertEqual(len({r["habitat"] for r in p}),2)
        for r in p:
            q=r["measured_seed_retention_per_exposure_treatment"]
            for k,v in r["source_seed_predation_percentage"].items():
                self.assertAlmostEqual(q[k],1-v/100)
            self.assertAlmostEqual(
                r["head_seed_yield_ratio_threshold_for_hypothetical_cage_route_to_exceed_open"],
                q["open"]/q["cage"])
            self.assertAlmostEqual(r["source_recorded_cage_minus_open_retention"],
                                   q["cage"]-q["open"])

    def test_ratio_of_means_not_mean_of_ratios(self):
        p=self.computed["source_aggregates"]
        self.assertAlmostEqual(p["ratio_of_group_mean_retention_cage_divided_by_open"],
                               2.56986301369863)
        self.assertAlmostEqual(p["mean_of_cellwise_ratios_cage_divided_by_open"],
                               3.364887266101202)
        self.assertNotAlmostEqual(
            p["ratio_of_group_mean_retention_cage_divided_by_open"],
            p["mean_of_cellwise_ratios_cage_divided_by_open"])
        self.assertAlmostEqual(p["source_break_even_head_yield_ratio_min"],
                               0.1789473684210527)
        self.assertAlmostEqual(p["source_break_even_head_yield_ratio_max"],
                               0.8901098901098902)

    def test_same_late_survival_cannot_reverse_true_head_output_rank(self):
        for q in (0,.01,.365,.938,1.0):
            result=gate.constant_survival_difference(12,10,q)
            self.assertGreaterEqual(result,0)
            self.assertAlmostEqual(result,2*q)
            reversed_result=gate.constant_survival_difference(10,12,q)
            self.assertLessEqual(reversed_result,0)
            self.assertAlmostEqual(reversed_result,-2*q)
        self.assertEqual(gate.constant_survival_difference(0,1,0),0)

    def test_common_early_survival_cannot_reverse_rank_if_head_costs_in_output(self):
        for q in (.1,.2,.4,1):
            for early in (.01,.25,1):
                benefit=gate.constant_survival_difference(100,95,q)*early
                self.assertGreater(benefit,0)
        self.assertTrue(self.computed["counterfactual_algebra_only"]["stem_stage"][
            "original_1982_any_stem_attack_fraction_is_not_stem_death_probability"])

    def test_differential_late_survival_can_reverse_hypothetical_head_rank(self):
        self.assertLess(
            gate.conditional_survival_difference(.8,1,.20,.80),0)
        self.assertGreater(
            gate.conditional_survival_difference(.8,1,.90,.30),0)
        self.assertEqual(gate.conditional_survival_difference(1,1,.5,.5),0)

    def test_exact_break_even_threshold_and_strict_inequality(self):
        q1=.938
        q0=.365
        threshold=gate.rank_threshold_for_different_survival(q0,q1)
        self.assertAlmostEqual(threshold,.3891257995735608)
        self.assertAlmostEqual(
            gate.conditional_survival_difference(threshold,1,q1,q0),0)
        self.assertGreater(
            gate.conditional_survival_difference(threshold+.01,1,q1,q0),0)
        self.assertLess(
            gate.conditional_survival_difference(threshold-.01,1,q1,q0),0)

    def test_original_cage_is_not_an_architectural_allele(self):
        d=self.computed["definitions"]
        self.assertTrue(d["independent_cross_source_numeric_product_forbidden"])
        self.assertTrue(d["use_of_cage_as_if_genetic_module_could_change_postdispersal_route_is_scenario_not_empirical_trait_effect"])
        self.assertIn("NOT viable seedlings",d["units_of_yield"])

    def test_illustrative_seed_output_ratios_never_promoted_to_observed_trait_effect(self):
        ratios=self.computed["counterfactual_algebra_only"]["declared_illustrative_Hratios"]
        self.assertEqual([x["hypothetical_predispersal_head_seed_output_ratio_protected_over_open"]
                          for x in ratios],[.5,.75,.8,.9,1.0])
        self.assertEqual([x["n_habitat_dates_where_protected_head_route_would_yield_more_surviving_exposed_seeds"]
                          for x in ratios],[8,8,8,10,10])
        for x in ratios:
            self.assertFalse(x["head_output_ratio_is_observed"])
            self.assertFalse(x["trait_changes_seed_caging_or_deposition_observed"])
            self.assertEqual(
                x["n_habitat_dates_where_open_route_would_yield_more"]+
                x["n_habitat_dates_where_protected_head_route_would_yield_more_surviving_exposed_seeds"]+
                x["n_exact_ties"],10)

    def test_real_arvense_source_2022_unrelated_to_1988_dutch_source(self):
        b=self.computed["source"]["other_head_source"]
        self.assertEqual(b["paper_doi"],"10.21203/rs.3.rs-6634172/v1")
        self.assertEqual(b["original_head_rows"],615)
        self.assertEqual(b["matched_four_arm_blocks"],36)
        self.assertAlmostEqual(b["author_defined_bag_minus_control_good_seed_fraction"],
                               -.5230663124174955,places=6)
        self.assertTrue(b["not_a_count_H_or_genotype_matched"])
        self.assertTrue(b["no_original_individual_joint_source_link_to_1988"])

    def test_fake_joined_source_doi_causes_fail_closed(self):
        x=copy.deepcopy(self.unlinked)
        x["source"]["doi"]=self.source["source"]["doi"]
        with self.assertRaisesRegex(ValueError,"TRUE_SOURCE_COHORT_INTEGRITY_CHANGED"):
            gate.analyze(self.source,x)

    def test_1982_any_stem_attack_is_not_1982_stem_mortality(self):
        self.assertEqual(
            self.source["study_years"]["1982"]["stem_census"]
            ["percent_stems_attacked_any_herbivore_reported"],98)
        self.assertTrue(
            self.computed["counterfactual_algebra_only"]["stem_stage"]
            ["original_1982_any_stem_attack_fraction_is_not_stem_death_probability"])

    def test_mutated_original_exclosure_predation_fails_closed(self):
        source=copy.deepcopy(self.source)
        source["postdispersal_exclosure"]["rows"][3]["net"]=0
        with self.assertRaisesRegex(ValueError,
            "PUBLISHED_1981_1982_TABLE2_ORIGINAL_VALUES_DRIFT"):
            gate.analyze(source,self.unlinked)

    def test_invalid_stage_survival_rejected_before_thresholds(self):
        for q in (-1,1.1,float("nan")):
            with self.assertRaisesRegex(ValueError,"INVALID_LATE_STAGE_SURVIVAL"):
                gate.constant_survival_difference(10,8,q)
        with self.assertRaisesRegex(ValueError,"CANNOT_IDENTIFY_SURVIVAL_RATIO"):
            gate.rank_threshold_for_different_survival(.4,0)
        with self.assertRaisesRegex(ValueError,"INVALID_JOINT_STAGE_PRODUCT"):
            gate.conditional_survival_difference(-2,1,.9,.1)

    def test_no_fabricated_genetic_selection_or_recruitment_from_cage(self):
        p=self.computed["counterfactual_algebra_only"]
        self.assertTrue(p["costs_and_limits"]["no_claim_of_architecture_global_advantage_slk_Phi_or_SCH_conflict"])
        self.assertTrue(p["costs_and_limits"]["offspring_establishment_or_late_architecture_costs_need_separate_term"])
        self.assertTrue(p["costs_and_limits"]["1988_cage_effect_cannot_be_inferred_for_head_morphology"])
        status=self.computed["measurement_gate"]
        self.assertEqual(status["status"],
                         "SOURCE_RETENTION_QUANTIFIED_BUT_HEAD_MODULE_TRANSMISSION_NOT_IDENTIFIABLE")
        self.assertEqual(status["new_field_observations"],0)
        self.assertFalse(status["same_head_seed_batch_split_to_randomized_postdispersal_access"])

    def test_independently_calculated_JS_all_source_cells_and_scenarios(self):
        def match(a,b,path="root"):
            if isinstance(b,dict):
                self.assertEqual(set(a),set(b),path)
                for k in b:match(a[k],b[k],path+"."+str(k))
            elif isinstance(b,list):
                self.assertEqual(len(a),len(b),path)
                for i,(x,y) in enumerate(zip(a,b)):match(x,y,path+f"[{i}]")
            elif type(b) is bool:
                self.assertIs(type(a),bool,path)
                self.assertEqual(a,b,path)
            elif type(b) in (int,float):
                self.assertIsInstance(a,(int,float),path)
                self.assertTrue(math.isfinite(a),path)
                self.assertTrue(math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),
                                (path,a,b))
            else:
                self.assertEqual(a,b,path)
        match(self.computed,self.frozen)

if __name__=="__main__":
    unittest.main()
