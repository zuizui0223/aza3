"""Source-gated independent Pedicularis partner-specific water/exsertion interfaces."""
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_pedicularis_rex_two_access_operators_v1 as m

EVIDENCE=ROOT/"data/evidence/pedicularis_rex_two_access_operators_author_source_v1.json"

class PublishedPedicularisInterfaces(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw=json.loads(EVIDENCE.read_text(encoding="utf-8"))
        cls.r=m.audit(cls.raw)

    def test_two_exact_peer_reviewed_source_dois(self):
        self.assertEqual(self.r["source_papers"],{
            "wet_bract_intervention":"10.1093/aobpla/plv019",
            "floral_exsertion_observational_selection":"10.1093/aob/mcw097"})

    def test_2015_one_site_visitor_not_six_sites(self):
        paper=self.raw["sources"]["drainage_2015"]
        self.assertEqual(paper["seed_populations"],6)
        self.assertEqual(paper["pollinator_video_populations"],1)
        self.assertEqual(paper["pollinator_video_population"],"Shama")
        self.assertEqual(paper["visitation_2010_2011_randomized_treatment_plots"]["per_year"],10)
        self.assertEqual(paper["visitation_2010_2011_randomized_treatment_plots"]["water_drained"],5)
        self.assertEqual(paper["visitation_2010_2011_randomized_treatment_plots"]["intact"],5)

    def test_2015_published_treatment_endpoint_null_vs_positive_separation(self):
        r=self.r["2015"]
        self.assertAlmostEqual(r["legitimate_pollinator_Shama_treatment_p"],.958)
        self.assertAlmostEqual(r["nectar_robber_Shama_treatment_p"],.951)
        self.assertAlmostEqual(r["initial_seed_set_treatment_p"],.906)
        self.assertEqual(r["final_seed_set_treatment_p_upper"],.0001)
        self.assertEqual(r["seed_predation_treatment_p_upper"],.0001)

    def test_conservative_original_six_location_BH_review_not_exact_small_p(self):
        s=self.r["2015"]
        self.assertEqual(s["nominal_seed_predation_significant_sites"],5)
        self.assertEqual(s["conservative_six_site_BH_significant_sites"],4)
        q=s["per_site_p_and_q_upper"]
        self.assertAlmostEqual(q["Deqin"]["BH_q_conservative_upper"],.0552)
        self.assertFalse(q["Deqin"]["passes_conservative_BH_at_0p05"])
        self.assertAlmostEqual(q["Zhongdian"]["p_exact"],.789)
        self.assertEqual(sum(x["passes_conservative_BH_at_0p05"] for x in q.values()),4)
        self.assertTrue(all(x["BH_q_conservative_upper"]<=.0002+1e-12
                            for k,x in q.items() if k in
                            ("Baishuitai","Sanba","Daxueshan")))
        self.assertAlmostEqual(q["Shama"]["BH_q_conservative_upper"],.0075)

    def test_conservative_source_p_lt_not_fabricated_as_exact_zero(self):
        q=self.r["2015"]["per_site_p_and_q_upper"]
        self.assertIsNone(q["Sanba"]["p_exact"])
        self.assertAlmostEqual(q["Sanba"]["p_reported_less_than"],.0001)
        self.assertTrue(self.raw["sources"]["drainage_2015"][
            "p_values_less_0p0001_treated_as_conservative_upper_limit"])

    def test_original_2016_stigmatic_pollen_grain_is_not_visits(self):
        p=self.r["2016"]
        self.assertEqual(p["trait_pollen_populations"],14)
        self.assertEqual(p["individual_pollen_plants"],299)
        self.assertEqual(p["individual_pollen_flowers"],598)
        self.assertTrue(p["corolla_exsertion_positively_associated_with"][
            "stigmatic_pollen_receipt"])

    def test_original_2016_flower_traits_do_not_join_all_population_seed_data(self):
        p=self.r["2016"]
        self.assertEqual(p["linked_predation_trait_populations"],7)
        self.assertEqual(p["linked_population_codes"],[1,3,5,8,9,10,11])
        self.assertEqual(self.raw["sources"]["exsertion_2016"][
            "lost_morphology_labels_populations"],[2,4,6,7,12])

    def test_truly_different_water_and_exsertion_operators(self):
        p=self.raw["synthesis"]
        self.assertTrue(p["same_species_different_trials"])
        self.assertTrue(p["not_same_individual_factorial_water_x_exsertion"])
        self.assertTrue(p["study_2015_water_cut_not_sham_water_only"])
        self.assertTrue(p["study_2015_visitor_non_significance_not_equivalence"])
        self.assertFalse(p["not_same_individual_factorial_water_x_exsertion"] is False)

    def test_2016_predation_pollen_positive_but_final_context_interactive(self):
        z=self.r["2016"]
        self.assertTrue(z["corolla_exsertion_positively_associated_with"][
            "predispersal_seed_predation"])
        self.assertEqual(z["model_best_AICc"],{
            "pollen_arrival":1700.1,
            "predispersal_seed_predation":-156.11,
            "initial_seed_set":-493.44,
            "final_viable_seed":-1000.03,
        })

    def test_2016_source_xls_not_claimed_downloaded(self):
        s=self.raw["sources"]["exsertion_2016"]
        self.assertTrue(s["individual_supplement_xls_url"].endswith(".xls"))
        self.assertFalse(s["individual_source_bytes_verified"])
        self.assertTrue(self.raw["synthesis"][
            "study_2016_original_raw_individual_reanalysis_not_performed"])

    def test_wrong_source_doi_is_rejected(self):
        d=copy.deepcopy(self.raw)
        d["sources"]["drainage_2015"]["doi"]="10.1234/invalid"
        with self.assertRaisesRegex(ValueError,"TWO_PRIMARY_PAPERS_MISSING"):
            m.audit(d)

    def test_wrong_per_site_significance_false_claim_rejected(self):
        d=copy.deepcopy(self.raw)
        p=d["sources"]["drainage_2015"]["site_specific_seed_predation"]
        next(x for x in p if x["name"]=="Deqin")["passes_BH_conservative"]=True
        with self.assertRaisesRegex(ValueError,"PUBLISHED_BH_SITE_CLASSIFICATION_CHANGED"):
            m.audit(d)

    def test_wrong_morphology_seed_population_join_rejected(self):
        d=copy.deepcopy(self.raw)
        d["sources"]["exsertion_2016"]["n_linked_trait_predation_populations"]=14
        with self.assertRaisesRegex(ValueError,"SOURCE_2016_DIFFERENT_POPULATION_DENOMINATORS_CHANGED"):
            m.audit(d)

    def test_randomized_factorial_without_actual_field_work_is_rejected(self):
        d=copy.deepcopy(self.raw)
        d["synthesis"]["not_same_individual_factorial_water_x_exsertion"]=False
        with self.assertRaisesRegex(ValueError,"FALSE_SAME_HEAD_FACTORIAL_CAUSAL_LOSS"):
            m.audit(d)
        self.assertEqual(self.r["biology"]["new_Cirsium_field_events"],0)

    def test_population_q_formula_exact_upper_threshold(self):
        values=[.0001,.0001,.789,.046,.0001,.005]
        q=m.conservative_bh(values)
        self.assertEqual(len(q),6)
        self.assertAlmostEqual(q[3],.0552)
        with self.assertRaisesRegex(ValueError,"INVALID_OR_INCOMPLETE_ORIGINAL_SIX_LOCATION_P_VALUES"):
            m.conservative_bh([.01])

if __name__=="__main__":
    unittest.main()
