"""Test 134 published C. canescens means and 106 raw C. undulatum mother rows."""
import copy,json,math,pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"analysis"))
import audit_cirsium_cross_lifehistory_compensation_v1 as q

BASE=pathlib.Path(__file__).resolve().parents[1]
CAN=BASE/"data/evidence/cirsium_canescens_2018_original_table1_source_v1.json"
UND=BASE/"data/external/usda_wavyleaf_2007_2008/mother_plants.csv"


class TestOriginalAndPublishedGrain(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.can=json.loads(CAN.read_text(encoding="utf-8"))
        cls.und=q.read_undulatum_source(UND)

    def test_published_canescens_plants_are_NOT_available_raw(self):
        r=q.validate_canescens_table(self.can)
        self.assertEqual([r[k]["n_plants"] for k in q.C_GROUPS],[37,31,34,32])
        self.assertEqual(sum(r[k]["n_plants"] for k in q.C_GROUPS),134)
        self.assertFalse(self.can["provenance"]["original_raw_plant_rows_available_here"])

    def test_actual_undulatum_data_have_not_been_swapped_for_published_means(self):
        self.assertEqual(len(self.und),106)
        self.assertEqual(len(set(r["year"]+"|"+r["plt"] for r in self.und)),106)
        groups=q.means_from_undulatum(self.und)
        self.assertAlmostEqual(groups["year_by_arm"]["2007"]["reduced_intact"]["mean_good_seed_per_plant"],88.4)
        self.assertAlmostEqual(groups["year_by_arm"]["2008"]["ambient_intact"]["mean_good_seed_per_plant"],40/15)
        self.assertEqual(groups["original_sample_weighted_pooled_by_arm"]["reduced_intact"]["n_mothers"],31)

    def test_canescens_2018_original_author_apical_x_herbivory_interaction(self):
        c=q.analyze(self.can,self.und)
        r=c["compensation_contrasts"]["canescens_author_two_year_pooled"]
        self.assertAlmostEqual(r["damage_minus_intact_good_seeds_under_reduced_herbivory"],121.9)
        self.assertAlmostEqual(r["damage_minus_intact_good_seeds_under_ambient_herbivory"],-74.3)
        self.assertAlmostEqual(r["damage_differences_difference_reduced_minus_ambient"],196.2)
        self.assertAlmostEqual(r["damage_to_intact_seed_mean_ratio_by_herbivory"]["reduced"],477.8/355.9)
        self.assertAlmostEqual(r["damage_to_intact_seed_mean_ratio_by_herbivory"]["ambient"],112/186.3)
        self.assertEqual(c["study_registry"]["canescens"]["original_author_apical_damage_by_herbivory_model_p"],.04)

    def test_undulatum_annual_damage_response_is_not_constant(self):
        z=q.analyze(self.can,self.und)
        y=z["compensation_contrasts"]["undulatum_real_year_specific"]
        self.assertAlmostEqual(y["2007"]["damage_to_intact_seed_mean_ratio_by_herbivory"]["reduced"],60/88.4)
        self.assertAlmostEqual(y["2007"]["damage_to_intact_seed_mean_ratio_by_herbivory"]["ambient"],(373/9)/(539/9))
        self.assertAlmostEqual(y["2008"]["damage_to_intact_seed_mean_ratio_by_herbivory"]["reduced"],(2886/19)/(2906/21))
        self.assertAlmostEqual(y["2008"]["damage_to_intact_seed_mean_ratio_by_herbivory"]["ambient"],(107/12)/(40/15))
        self.assertLess(y["2007"]["damage_minus_intact_good_seeds_under_ambient_herbivory"],0)
        self.assertGreater(y["2008"]["damage_minus_intact_good_seeds_under_ambient_herbivory"],0)

    def test_pooled_undulatum_comes_from_original_maternal_counts(self):
        out=q.analyze(self.can,self.und)["compensation_contrasts"]
        pooled=out["undulatum_original_mother_sample_weighted_two_year_pooled"]
        r=pooled["mean_by_arm"]
        self.assertAlmostEqual(r["reduced_intact"],3790/31)
        self.assertAlmostEqual(r["reduced_damaged"],3546/30)
        self.assertAlmostEqual(r["ambient_intact"],579/24)
        self.assertAlmostEqual(r["ambient_damaged"],480/21)
        self.assertAlmostEqual(pooled["ratio_of_damage_ratios_reduced_to_ambient"],
            ((3546/30)/(3790/31))/((480/21)/(579/24)))

    def test_year_equal_is_distinct_from_n_weighted_pooled(self):
        a=q.analyze(self.can,self.und)["compensation_contrasts"]
        pooled=a["undulatum_original_mother_sample_weighted_two_year_pooled"]
        equal=a["undulatum_original_mother_equal_year_weighted"]
        self.assertNotAlmostEqual(pooled["damage_to_intact_seed_mean_ratio_by_herbivory"]["ambient"],
                                  equal["damage_to_intact_seed_mean_ratio_by_herbivory"]["ambient"])
        self.assertEqual(equal["mean_by_arm"]["reduced_intact"],
                         (88.4+2906/21)/2)

    def test_excluding_conflicting_original_2007_ids_changes_ambient_damage_sign(self):
        z=q.analyze(self.can,self.und)["compensation_contrasts"]
        before=z["undulatum_real_year_specific"]["2007"]["damage_minus_intact_good_seeds_under_ambient_herbivory"]
        after=z["undulatum_exclude_3_conflicting_mothers"]["by_year"]["2007"]["damage_minus_intact_good_seeds_under_ambient_herbivory"]
        self.assertAlmostEqual(before,(373-539)/9)
        self.assertAlmostEqual(after,353/8-326/8)
        self.assertLess(before,0)
        self.assertGreater(after,0)

    def test_canescens_published_value_tampering_rejected(self):
        bad=copy.deepcopy(self.can)
        bad["table1_plant_rows"][1]["undamaged_seeds_per_plant"]["mean"]=999
        with self.assertRaisesRegex(ValueError,"CANESCENS_TABLE1_PUBLISHED_VALUE_CHANGED"):
            q.analyze(bad,self.und)

    def test_canescens_raw_fiction_rejected(self):
        bad=copy.deepcopy(self.can)
        bad["provenance"]["original_raw_plant_rows_available_here"]=True
        with self.assertRaisesRegex(ValueError,"CANESCENS_PUBLISHED_NOT_RAW"):
            q.analyze(bad,self.und)

    def test_source_mother_treatment_fiction_rejected(self):
        bad=copy.deepcopy(self.und)
        bad[0]["trt"]="incorrect"
        # These original rows are loaded by frozen git blob; synthetic validator:
        with self.assertRaisesRegex(ValueError,"EMPTY_UNDULATUM_YEAR_TREATMENT_CELL|UNDULATUM"):
            q.means_from_undulatum([r for r in bad if r["trt"]!="incorrect"]
                                   if len([r for r in bad if r["trt"]=="incorrect"])>0
                                   else [])

    def test_cross_species_lifehistory_causation_disallowed(self):
        z=q.analyze(self.can,self.und)
        limits=z["all_scientific_inference_limits"]
        self.assertFalse(limits["between_species_life_history_effect_identified"])
        self.assertFalse(limits["the_two_species_represent_two_independently_replicated_life_history_classes"])
        self.assertFalse(limits["a_cross_species_significance_test_was_performed"])
        self.assertFalse(limits["canescens_134_original_mother_rows_recovered"])
        self.assertTrue(z["no_canescens_raw_bootstrap_or_constructed_p_from_SEs"])


if __name__=="__main__":
    unittest.main()
