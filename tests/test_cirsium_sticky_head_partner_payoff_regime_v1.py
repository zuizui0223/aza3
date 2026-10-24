"""Real named original studies, incompatible cohorts, and route-specific fitness caveats."""
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_sticky_head_partner_payoff_regime_v1 as v

SOURCE=ROOT/"data/evidence/cirsium_discolor_sticky_head_partner_access_payoff_source_study_registry_v1.json"
COMPARATOR=ROOT/"data/evidence/pedicularis_rex_two_access_operators_author_source_v1.json"

class OriginalCirsiumStickyPayoffTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw=json.loads(SOURCE.read_text(encoding="utf-8"))
        cls.comparator=json.loads(COMPARATOR.read_text(encoding="utf-8"))
        cls.gate=v.audit(cls.raw,cls.comparator)

    def test_actual_original_paper_sources_distinguished(self):
        self.assertEqual(self.gate["provenance"]["primary_publications"],3)
        self.assertEqual(self.gate["provenance"]["source_publication_dois"],
                         {"WAT1983":"10.2307/2425228",
                          "T2003":"10.22543/0090-0222.2085",
                          "T2007":"10.22543/0090-0222.2188"})
        self.assertTrue(self.gate["provenance"]
                        ["original_reports_not_raw_observation_matrix"])

    def test_sticky_original_cohorts_not_pooled_effect_meta(self):
        z=self.gate["original_seed_predation_support"]
        self.assertEqual(self.gate["provenance"]["original_source_organism_cohorts"],4)
        self.assertEqual(z["significant_source_reports"],1)
        self.assertEqual(z["reported_nondetections_not_null_effect_size"],3)
        self.assertIsNone(z["pooled_experimental_effect_size"])
        self.assertIsNone(z["confidence_interval_for_sticky_defence_benefit"])
        self.assertFalse(z["sign_reversal_demonstrated"])

    def test_actual_1980_cirsium_discolor_trial_has_positive_evidence(self):
        s=self.gate["original_seed_predation_support"]
        self.assertEqual(s["Cirsium discolor 1980"],
                         "DETECTED_GREATER_PREDATION_WHEN_STICKINESS_BLOCKED")
        self.assertEqual(s["Cirsium discolor 1982-1987"],
                         "NO_SUPPORT_FOR_SAME_EFFECT_IN_ORIGINAL_ABSTRACT")
        self.assertTrue(self.gate["provenance"]
                        ["1980_discolor_and_1982_1987_discolor_are_not_same_field_cohort"])

    def test_other_cirsium_not_assumed_to_share_positive_effect(self):
        s=self.gate["original_seed_predation_support"]
        self.assertEqual(s["Cirsium flodmanii 1980"],"NO_DETECTED_EFFECT")
        self.assertEqual(s["Cirsium flodmanii 1982"],"NO_DETECTED_EFFECT")
        self.assertTrue(self.gate["provenance"]
                        ["C_fodmani_1980_and_1982_are_not_two_independent_evolutionary_origins"])

    def test_sticky_surface_observational_roles_avoid_single_enemy_class(self):
        a=self.gate["independent_observed_visitor_access_2007"]
        self.assertEqual(len(a),4)
        self.assertEqual(a["reported_pollinators"],"DO_NOT_CONTACT_STICKY_TRAPS")
        self.assertEqual(a["specialist_seed_consumers"],
                         "SOME_SEED_PREDATORS_BYPASS_STICKY_TRAPS")
        self.assertEqual(a["ants_and_aphids"],"REPORTED_DETERRED_BY_STICKY_TRAPS")
        self.assertIn("SALTICID",a["small_prey_gleaning_arthropod_predators"])

    def test_absent_contact_not_equivalent_to_absent_pollination(self):
        a=self.gate["biological_inference"]
        self.assertTrue(a["pollinator_noncontact_not_zero_pollen_service"])
        self.assertFalse(a["trait_by_true_insect_contact_measured_in_same_head_to_fitness"])
        self.assertFalse(a["universal_functional_trait_payoff_transfers_to_other_Cirsium_lineages"])

    def test_pedicularis_is_original_external_functional_comparator(self):
        p=self.gate["positive_non_Cirsium_comparator"]
        self.assertEqual(p["organism"],"Pedicularis rex")
        self.assertEqual(p["source_populations_seed_outcome"],6)
        self.assertEqual(p["source_populations_visitor_frequency"],1)
        self.assertEqual(p["site_specific_predation_significant_after_conservative_BH"],4)
        self.assertTrue(p["not_homologous_to_sticky_Cirsium_head_intervention"])

    def test_three_publications_not_independent_ancestral_evolutionary_origins(self):
        x=self.raw["provenance_caveats"]
        self.assertTrue(x["papers_1983_2003_2007_are_not_three_independent_origins"])
        self.assertTrue(x["observation_2007_may_share_study_species_2003_but_no_head_links_verified"])
        self.assertTrue(self.gate["biological_inference"]["no_ancestral_regain_or_genetic_LEGO"])

    def test_published_author_1983_original_raw_p_missing_not_zero(self):
        for study in self.raw["primary_sources"]:
            for cohort in study.get("organismal_samples",[]):
                self.assertIsNone(cohort["source_p_value"])
                self.assertIsNone(cohort["source_mean_seed_loss"])
                self.assertFalse(cohort["source_individual_head_records_available"])

    def test_change_1980_sticky_sign_must_fail(self):
        x=copy.deepcopy(self.raw)
        x["primary_sources"][0]["organismal_samples"][0][
            "authors_reported_response"]="NO_EFFECT_DETECTED"
        with self.assertRaisesRegex(ValueError,"PUBLISHED_STICKY_SEED_OUTCOME_CODE_CHANGED"):
            v.audit(x,self.comparator)

    def test_change_2003_sticky_null_to_equivalence_must_fail(self):
        x=copy.deepcopy(self.raw)
        x["primary_sources"][1]["organismal_samples"][0][
            "authors_reported_response"]="NO_EFFECT_EQUIVALENCE_ESTABLISHED"
        with self.assertRaisesRegex(ValueError,"PUBLISHED_STICKY_SEED_OUTCOME_CODE_CHANGED"):
            v.audit(x,self.comparator)

    def test_change_sticky_effect_to_fabricated_numeric_throws(self):
        x=copy.deepcopy(self.raw)
        x["primary_sources"][0]["organismal_samples"][0]["source_p_value"]=.001
        with self.assertRaisesRegex(ValueError,
                "UNVERIFIED_STICKY_NUMERIC_EFFECT_OR_HEAD_RECORD_FABRICATED"):
            v.audit(x,self.comparator)

    def test_mutualist_role_egg_link_fabrication_fails(self):
        x=copy.deepcopy(self.raw)
        x["primary_sources"][2]["observed_interfaces"][0][
            "original_observation"]="NOTHING_POLLINATED"
        with self.assertRaisesRegex(ValueError,
                "PUBLISHED_ARTHROPOD_INTERFACE_ROLES_CHANGED"):
            v.audit(x,self.comparator)

    def test_false_cross_study_same_environment_causal_meta_fails(self):
        x=copy.deepcopy(self.raw)
        x["biology_contrast"]["source_study_environment_and_treatment_comparable"]=True
        with self.assertRaisesRegex(ValueError,
            "CAUSAL_TEMPORAL_REGIME_SHIFT_OR_META_EFFECT_FALSELY_IDENTIFIED"):
            v.audit(x,self.comparator)

    def test_wrong_comparator_bh_denominator_fails(self):
        x=copy.deepcopy(self.comparator)
        x["sources"]["drainage_2015"]["conservative_BH_significant_sites"]=5
        with self.assertRaisesRegex(ValueError,
            "PEDICULARIS_SOURCE_NOT_DISTINCT_OR_INCORRECT_GRAIN"):
            v.audit(self.raw,x)

    def test_no_false_original_head_fitness_and_new_field_fate(self):
        k=self.gate["biological_inference"]
        self.assertEqual(k["new_original_organismal_field_rows"],0)
        self.assertTrue(k["no_true_spines_capitulum_gravity_angle_measured"])
        self.assertFalse(k["sticky_exudate_protective_net_seed_fitness_in_all_environments_identified"])

if __name__=="__main__":
    unittest.main()
