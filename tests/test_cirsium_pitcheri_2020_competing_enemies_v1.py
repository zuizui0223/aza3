"""2020 published Cirsium pitcheri head hierarchy: observational falsification guards."""
from pathlib import Path
import sys
import json
import unittest
import copy
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"analysis"))
import audit_cirsium_pitcheri_2020_competing_enemies_v1 as m

SOURCE=Path(__file__).resolve().parents[1]/"data/evidence/cirsium_pitcheri_2020_multi_enemy_rank_source_v1.json"

class TestCompetingHeadRisk(unittest.TestCase):
    def setUp(self):
        self.source=json.loads(SOURCE.read_text(encoding="utf-8"))

    def test_two_different_threat_processes_have_distinct_denominators(self):
        d=m.analyze(self.source)
        self.assertAlmostEqual(d["within_WFD_2016_weevil_head_rank"]["secondary_vs_terminal_infestation_ratio"],.50/.27)
        self.assertAlmostEqual(d["within_WFD_2016_weevil_head_rank"]["tertiary_vs_terminal_infestation_ratio"],.80/.27)
        self.assertTrue(d["inference_guards"]["cannot_compare_weevil_ranks_with_deer_caused_loss_as_one_same_cohort"])

    def test_terminal_missing_plant_fraction_increased_both_sites(self):
        d=m.analyze(self.source)["terminal_head_loss_cross_year_observations"]
        self.assertAlmostEqual(d["SC"]["difference_in_missing_terminal_plant_fraction"],.27)
        self.assertAlmostEqual(d["WFD"]["difference_in_missing_terminal_plant_fraction"],.241)
        self.assertFalse(d["WFD"]["original_cause_is_verified_deer_attack"])

    def test_seeds_quantity_and_quality_both_in_original_association(self):
        d=m.analyze(self.source)
        self.assertAlmostEqual(d["WFD_2016_head_seed_quantity"]["infested_to_uninfested_mature_seed_count_ratio"],.4)
        self.assertAlmostEqual(d["WFD_2016_head_seed_weight"]["terminal_infested_minus_uninfested_mg"],-1.83)
        self.assertEqual(d["WFD_2016_head_seed_weight"]["terminal_CI95_infested_minus_uninfested_mg"],[-2.73,-.93])
        self.assertTrue(d["WFD_2016_head_seed_weight"]["secondary_CI95_contains_zero"])
        self.assertFalse(d["WFD_2016_head_seed_weight"]["difference_of_terminal_and_secondary_effects_significance_known"])

    def test_site_comparison_not_equivalent_to_head_loss(self):
        d=m.analyze(self.source)
        self.assertAlmostEqual(d["site_fraction_of_attempted_seeds_mature_2016"]["WFD_minus_SC"],.33)
        self.assertFalse(d["site_fraction_of_attempted_seeds_mature_2016"]["WFD_site_difference_is_positive_effect_of_weevils"])

    def test_rank_2017_not_always_fecundity_recovery(self):
        d=m.analyze(self.source)["field_year_head_output"]
        self.assertEqual(d["terminal_mature_seed_count_relative_to_2016_approx"],.3)
        self.assertEqual(d["tertiary_mature_seeds_reported_at_both_sites_in_2017"],0)
        self.assertFalse(d["absence_of_detectable_secondary_year_difference_implies_causal_buffering"])

    def test_does_not_call_author_hypothesis_direct_deer_identification(self):
        s=copy.deepcopy(self.source)
        s["terminal_head_absence_by_site_year"]["specific_herbivore_causation_identified"]=True
        with self.assertRaisesRegex(ValueError,"DEER_CAUSE_MUST_BE_UNCERTAIN"):
            m.analyze(s)

    def test_rejects_invented_2020_raw_head_data(self):
        s=copy.deepcopy(self.source)
        s["source"]["primary_published_article_transcription_not_raw_head_rows"]=False
        with self.assertRaisesRegex(ValueError,"PUBLISHED_SOURCE_GRAIN_DRIFT"):
            m.analyze(s)

    def test_rejects_signed_secondary_seed_mass_without_estimate(self):
        s=copy.deepcopy(self.source)
        s["rank_specific_seed_mass_2016_WFD"]["secondary"]["infested_minus_uninfested_seed_mass_mg"]=-.5
        with self.assertRaisesRegex(ValueError,"SECONDARY_SEED_WEIGHT_ESTIMATE_NOT_REPORTED"):
            m.analyze(s)

    def test_refuses_cross_study_parent_seed_linkage(self):
        s=copy.deepcopy(self.source)
        s["material_limit"]["2020_and_2021_seeds_not_joined_by_original_head_id"]=False
        with self.assertRaisesRegex(ValueError,"DIFFERENT_STUDIES_HAVE_NO_SHARED_HEAD_IDS"):
            m.analyze(s)

    def test_another_population_model_not_reused_as_2007_biological_replicate(self):
        d=m.analyze(self.source)
        self.assertEqual(d["n_new_individual_head_field_observations"],0)
        self.assertTrue(d["inference_guards"]["head_rank_is_not_phyllary_spine_or_gravity_angle"])
        self.assertTrue(d["inference_guards"]["risk_diversification_is_prospective_not_established_adaptation"])

if __name__=="__main__":unittest.main()
