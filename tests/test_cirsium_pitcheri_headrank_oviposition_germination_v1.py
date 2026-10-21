"""Scientific semantic gates for PUBLISHED table cells, not new seed experiments."""
import copy
import json
import pathlib
import sys
import unittest

ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_pitcheri_headrank_oviposition_germination_v1 as v

SOURCE=ROOT/"data/evidence/cirsium_pitcheri_published_rank_germination_cells_v1.json"


class TestGijsmanPublishedStages(unittest.TestCase):
    def setUp(self):
        self.s=json.loads(SOURCE.read_text(encoding="utf-8"))

    def test_exact_24_source_model_cells(self):
        c=v.validate(self.s["published_2021_model_rows"])
        self.assertEqual(c["n_source_model_cells"],24)
        self.assertEqual(c["n_per_collection_year"],{2016:12,2017:12})

    def test_2016_final_round_not_2017_individual_total(self):
        r=v.run(self.s)
        self.assertEqual(r["coverage"]["n_source_model_cells"],24)
        self.assertEqual(r["published_2016_and_2017_rank_germination"]["n_2016_terminal_better_out_of_six"],5)
        self.assertEqual(r["published_2016_and_2017_rank_germination"]["n_2017_terminal_better_out_of_six"],3)

    def test_published_terminal_advantage_is_context_sensitive_and_can_change_sign(self):
        a=v.run(self.s)["published_2016_and_2017_rank_germination"]
        by={(z["collection_year"],z["site"],z["day_temp_c"]):z for z in a["paired_model_rows"]}
        self.assertAlmostEqual(by[(2016,"SC",20)]["terminal_minus_secondary"],.631)
        self.assertAlmostEqual(by[(2016,"WFD",20)]["terminal_minus_secondary"],.314)
        self.assertAlmostEqual(by[(2017,"SC",20)]["terminal_minus_secondary"],.116)
        self.assertAlmostEqual(by[(2017,"WFD",20)]["terminal_minus_secondary"],-.033)
        self.assertTrue(all(by[(2017,"WFD",t)]["terminal_minus_secondary"]<0 for t in (20,25,30)))
        self.assertTrue(all(by[(2016,"WFD",t)]["terminal_minus_secondary"]>0 for t in (20,25,30)))
        self.assertTrue(a["not_an_inference_of_significant_between_year_site_by_rank_interaction"])

    def test_weevil_reported_percentages_approximate_only(self):
        z=v.run(self.s)["reported_weevil_seed_stage"]
        self.assertAlmostEqual(z["secondary_to_terminal_approx_infested_fraction_ratio"],.50/.27)
        self.assertAlmostEqual(z["tertiary_to_terminal_approx_infested_fraction_ratio"],.80/.27)
        self.assertEqual(z["mature_seeds_infested_vs_uninfested_observed_association_ratio"],.40)
        self.assertTrue(z["relative_risk_is_from_rounded_published_percentages_not_raw_binomial_estimator"])
        self.assertTrue(z["infestation_ratio_cannot_be_multipled_by_viable_seed_ratio_to_infer_per_rank_fitness"])

    def test_no_double_use_same_head_2020_2021(self):
        z=v.run(self.s)
        self.assertTrue(z["claim_ceiling"]["two_papers_do_not_yield_original_same_head_id_join"])
        self.assertTrue(z["claim_ceiling"]["weevil_infested_heads_are_not_randomized_to_weevils"])
        self.assertTrue(z["claim_ceiling"]["no_inferred_net_recruitment_reproductive_fitness_by_multiplication"])

    def test_duplicate_factorial_cell_rejected(self):
        s=copy.deepcopy(self.s)
        s["published_2021_model_rows"][1]=s["published_2021_model_rows"][0]
        with self.assertRaisesRegex(ValueError,"DUPLICATE_FITTED_SOURCE_CELL"):
            v.run(s)

    def test_false_cumulative_decline_rejected(self):
        s=copy.deepcopy(self.s)
        s["published_2021_model_rows"][0]["round2"]["fraction"]=.90
        with self.assertRaisesRegex(ValueError,"CUMULATIVE_GERMINATION_DECREASES"):
            v.run(s)

    def test_fitted_proportion_bounds_rejected(self):
        s=copy.deepcopy(self.s)
        s["published_2021_model_rows"][15]["total_fraction"]=1.2
        with self.assertRaisesRegex(ValueError,"INVALID_2017_GERMINATION_MODEL_RESULT"):
            v.run(s)

    def test_original_source_anchor_cell_locked(self):
        s=copy.deepcopy(self.s)
        s["published_2021_model_rows"][0]["round3"]["fraction"]=.829
        with self.assertRaisesRegex(ValueError,"PUBLISHED_PROVENANCE_ANCHOR_CELL_CHANGED"):
            v.run(s)

    def test_two_cohorts_noninterchangeable(self):
        s=copy.deepcopy(self.s)
        s["published_2021_model_rows"][0]["measurement"]="fitted_total_germination_after_three_rounds"
        with self.assertRaisesRegex(ValueError,"2016_SOURCE_ROUND3_DEFINITION_DRIFT"):
            v.run(s)

    def test_no_authoritative_unmeasured_head_angle(self):
        z=v.run(self.s)
        self.assertTrue(z["claim_ceiling"]["rank_related_to_size_and_anthesis_order_not_independent_orientation_or_spine_effect"])
        self.assertTrue(z["claim_ceiling"]["no_genetic_architecture_or_reexpression_measured"])

    def test_no_numerical_fabrication_from_unknown_cell_sizes(self):
        self.assertTrue(self.s["study_population_sample_info"]["samples_per_table_cell_missing"])
        self.assertTrue(self.s["claim_ceiling"]["fitted_mean_SE_is_not_individual_trial_data"])
        self.assertFalse(self.s["claim_ceiling"]["original_seed_level_data_downloaded_and_reanalysed"])


if __name__=="__main__":
    unittest.main()
