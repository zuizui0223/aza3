"""Original Cirsium vulgare 1988 stage-specific predation and denominator guards.

Tests replicate source Table 1 + 2 at their actual reporting level, never
claiming unobserved individual-level linked field outcomes.
"""
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_vulgare_1988_multistage_herbivore_gate_v1 as m

ORIG=ROOT/"data/evidence/cirsium_vulgare_1988_original_multistage_source_tables_v1.json"
FROZEN=ROOT/"data/evidence/cirsium_vulgare_1988_multistage_risk_independent_JS_v1.json"

class OriginalVulgareMultistagePredation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original=json.loads(ORIG.read_text(encoding="utf-8"))
        cls.frozen=json.loads(FROZEN.read_text(encoding="utf-8"))
        cls.python=m.analyze(cls.original)

    def test_one_real_peer_reviewed_study_not_synthetic_multi_paper_meta(self):
        self.assertEqual(self.original["source"]["doi"],"10.2307/2260602")
        self.assertEqual(self.original["source"]["paper_year"],1988)
        self.assertIn("Meijendel",self.original["source"]["study_locality"])
        self.assertFalse(self.original["source"]["original_publisher_pdf_bytes_verified"])
        self.assertTrue(self.original["source"]
                        ["reported_data_are_published_table_and_text_values_not_retrieved_individual_records"])

    def test_source_year_sample_sizes_match_original_sand_dunes(self):
        y=self.original["study_years"]
        self.assertEqual(y["1982"]["stem_census"]["n_random_flowering_individuals"],117)
        self.assertEqual(y["1984"]["stem_census"]["n_random_flowering_individuals"],67)
        self.assertEqual(self.python["year_biological_stage_contrasts"]["1982"]
                         ["n_head_samples"],268)
        self.assertEqual(self.python["year_biological_stage_contrasts"]["1984"]
                         ["n_head_samples"],97)

    def test_source_1982_undamaged_50_extra_plants_nonrandom_bias_preserved(self):
        x=self.original["study_years"]["1982"]["stem_census"]
        self.assertEqual(x["n_random_flowering_plants_undamaged_stems"],2)
        self.assertEqual(x["n_additional_nonstratified_undamaged_stem_plants_nonnrandom"],50)
        self.assertEqual(x["n_seed_sample_undamaged_stem_plants"],52)
        self.assertTrue(x["undamaged_plant_sample_not_random_population_sample"])
        self.assertTrue(self.python["warnings"]["nonnrandom_1982_undamaged_stem_comparator"])

    def test_source_1984_undamaged_comparator_is_not_1982_extra_sample(self):
        x=self.original["study_years"]["1984"]["stem_census"]
        self.assertEqual(x["n_random_flowering_plants_undamaged_stems"],51)
        self.assertEqual(x["n_additional_nonstratified_undamaged_stem_plants_nonnrandom"],0)
        self.assertEqual(x["n_seed_sample_undamaged_stem_plants"],51)
        self.assertFalse(x["undamaged_plant_sample_not_random_population_sample"])

    def test_stem_attacked_and_head_attacked_year_directions_are_opposite(self):
        d=self.python["year_change_1984_minus_1982"]
        self.assertEqual(d["stem_damaged_prevalence_pp"],-74)
        self.assertEqual(d["head_visible_damaged_prevalence_pp"],26)
        self.assertEqual(d["reported_stem_seed_reduction_pp"],-65)
        self.assertEqual(d["reported_predation_within_heads_seed_loss_pp"],14)
        self.assertTrue(d["year_change_is_different_field_cohorts_not_random_experimental_assignment"])

    def test_seed_means_are_at_original_plant_grain_not_head_grain(self):
        a=self.original["study_years"]
        self.assertEqual(a["1982"]["stem_census"]
                         ["mean_good_seeds_per_random_flowering_plant"],50)
        self.assertEqual(a["1982"]["stem_census"]
                         ["mean_good_seeds_per_undamaged_stem_plant"],246)
        self.assertEqual(a["1984"]["stem_census"]
                         ["mean_good_seeds_per_random_flowering_plant"],1470)
        self.assertEqual(a["1984"]["stem_census"]
                         ["mean_good_seeds_per_undamaged_stem_plant"],1713)
        self.assertAlmostEqual(self.python["year_biological_stage_contrasts"]
            ["1982"]["stem_seed_reduction_pct_replayed"],79.67479674796748)
        self.assertAlmostEqual(self.python["year_biological_stage_contrasts"]
            ["1984"]["stem_seed_reduction_pct_replayed"],14.185639229422065)

    def test_true_head_attacked_vs_unattacked_observed_good_seeds(self):
        z=self.python["year_biological_stage_contrasts"]
        self.assertEqual((z["1982"]["head_damage_count"],
                          z["1982"]["head_good_seeds_attacked_minus_unattacked_mean"]),
                         (67,-6))
        self.assertEqual((z["1984"]["head_damage_count"],
                          z["1984"]["head_good_seeds_attacked_minus_unattacked_mean"]),
                         (49,-69))
        self.assertAlmostEqual(z["1982"]["head_damage_fraction_from_counts"],.25)
        self.assertAlmostEqual(z["1984"]["head_damage_fraction_from_counts"],
                               49/97)
        self.assertEqual(z["1982"]["within_head_predation_mean_seed_loss_pct_reported"],3)
        self.assertEqual(z["1984"]["within_head_predation_mean_seed_loss_pct_reported"],17)

    def test_reproductive_rank_not_promoted_to_genetic_head_module_fitness(self):
        warnings=self.python["warnings"]
        self.assertTrue(warnings["not_legitimate_to_add_stem_seed_reduction_plus_flowerhead_predation_pct"])
        self.assertTrue(warnings["observed_head_groups_are_not_randomized_damage_treatment"])
        self.assertTrue(warnings["no_detailed_source_flowerhead_individual_ids"])
        self.assertTrue(warnings["no_genetically_inherited_capitulum_LEGO_selection_or_trait_reexpression_inferred"])
        self.assertEqual(warnings["new_own_observations"],0)

    def test_postdispersal_exact_original_5_dates_2_habitats_3_meshes(self):
        t=self.original["postdispersal_exclosure"]
        self.assertEqual(len(t["rows"]),10)
        self.assertEqual(t["treatments"]["cage"]["mesh_cm"],.6)
        self.assertEqual(t["treatments"]["net"]["mesh_cm"],2)
        self.assertIsNone(t["treatments"]["open"]["mesh_cm"])
        self.assertEqual(len({v["date"] for v in t["rows"]}),5)
        self.assertEqual(set(v["habitat"] for v in t["rows"]),
                         {"scrub","herbaceous_vegetation"})
        self.assertTrue(t["table2_values_are_group_level_percent_seed_loss_not_individual_container_rows"])

    def test_source_exclosure_qualitative_predator_isolation_is_rodent_consistent(self):
        g=self.python["postdispersal_exclosure"]
        self.assertEqual(g["n_source_habitat_dates"],10)
        self.assertEqual(g["original_measured_table2_cells"],30)
        self.assertAlmostEqual(g["equal_habitat_date_mean_predation_pct"]["cage"],6.2)
        self.assertAlmostEqual(g["equal_habitat_date_mean_predation_pct"]["net"],65.4)
        self.assertAlmostEqual(g["equal_habitat_date_mean_predation_pct"]["open"],63.5)
        self.assertAlmostEqual(g["equal_habitat_date_mean_open_minus_cage_pp"],57.3)
        self.assertAlmostEqual(g["equal_habitat_date_mean_net_minus_cage_pp"],59.2)
        self.assertAlmostEqual(g["equal_habitat_date_mean_net_minus_open_pp"],1.9)
        self.assertEqual(g["n_habitat_dates_open_greater_cage"],10)
        self.assertEqual(g["n_habitat_dates_net_greater_cage"],10)

    def test_habitat_specific_exclosure_differences_are_not_maternal_effects(self):
        h=self.python["postdispersal_exclosure"]["per_habitat"]
        self.assertAlmostEqual(h["scrub"]["open_minus_cage"],48)
        self.assertAlmostEqual(h["herbaceous_vegetation"]["open_minus_cage"],66.6)
        self.assertEqual(h["scrub"]["n_dates"],5)
        self.assertEqual(h["herbaceous_vegetation"]["n_dates"],5)
        self.assertTrue(self.python["postdispersal_exclosure"]
            ["descriptive_10_cells_not_ten_independent_origin_or_plant_fitness_measures"])

    def test_postdispersal_separate_source_studies_not_one_seed_fate_series(self):
        s=self.original["postdispersal_exclosure"]
        self.assertEqual(s["separate_seeds_sown_95_local_populations_1982"]
                         ["n_populations"],95)
        self.assertEqual(s["separate_seeds_sown_95_local_populations_1982"]
                         ["mean_seed_predation_percent_per_week"],39)
        self.assertEqual(s["separate_1981_3_week_sown_grid_experiment"]
                         ["n_grids"],9)
        self.assertEqual(s["separate_1981_3_week_sown_grid_experiment"]
                         ["author_reported_fraction_preyed_after_three_weeks_percent"],60)
        self.assertTrue(self.python["warnings"]["different_study_years_and_postdispersal_year"])

    def test_original_exclosure_table_falsification_detects_altered_seed_loss(self):
        o=copy.deepcopy(self.original)
        o["postdispersal_exclosure"]["rows"][0]["cage"]=80
        with self.assertRaisesRegex(ValueError,
              "PUBLISHED_1981_1982_TABLE2_ORIGINAL_VALUES_DRIFT"):
            m.analyze(o)

    def test_fabricated_random_1982_undamaged_stem_comparator_rejected(self):
        o=copy.deepcopy(self.original)
        o["study_years"]["1982"]["stem_census"][
            "n_additional_nonstratified_undamaged_stem_plants_nonnrandom"]=0
        with self.assertRaisesRegex(ValueError,"SOURCE_ORIGINAL_SAMPLE_GRAIN_DRIFT_1982"):
            m.analyze(o)

    def test_fabricated_1984_5_percent_head_damage_rejected(self):
        o=copy.deepcopy(self.original)
        o["study_years"]["1984"]["before_dispersal_flowerhead_census"][
            "pct_flowerheads_with_traces_of_seed_predation_reported"]=5
        with self.assertRaisesRegex(ValueError,
            "ORIGINAL_PAPER_TWO_STAGES_UNDER_TWO_YEARS_CHANGED"):
            m.analyze(o)

    def test_falsely_authenticated_raw_flowerhead_dataset_rejected(self):
        o=copy.deepcopy(self.original)
        o["source"]["original_publisher_pdf_bytes_verified"]=True
        with self.assertRaisesRegex(ValueError,"FICTITIOUS_VERIFIED_INDIVIDUAL_SOURCE"):
            m.analyze(o)

    def test_no_individual_container_replication_from_aggregate_table(self):
        o=copy.deepcopy(self.original)
        o["postdispersal_exclosure"][
            "separate_paired_container_level_raw_observations_retrieved"]=True
        with self.assertRaisesRegex(ValueError,
            "FICTITIOUS_INDEPENDENT_SEED_CONTAINER_REPLICATES"):
            m.analyze(o)

    def test_genuine_derived_receipt_matches_python_in_all_numerical_levels(self):
        def same(a,b,path="root"):
            if isinstance(b,dict):
                self.assertEqual(set(a),set(b),path)
                for k in b:same(a[k],b[k],path+"."+str(k))
            elif isinstance(b,list):
                self.assertEqual(len(a),len(b),path)
                for i,(x,y) in enumerate(zip(a,b)):
                    same(x,y,path+f"[{i}]")
            elif isinstance(b,bool):
                self.assertIs(type(a),bool,path)
                self.assertEqual(a,b,path)
            elif type(b) in (int,float):
                self.assertAlmostEqual(a,b,places=10,msg=path)
            else:
                self.assertEqual(a,b,path)
        same(self.python,self.frozen)

if __name__=="__main__":
    unittest.main()
