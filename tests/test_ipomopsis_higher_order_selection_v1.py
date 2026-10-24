"""Original Ipomopsis 2010–18 higher-order trait × life-stage audit tests.

No synthetic data are entered into real-data repository outputs.
"""
import copy
import json
import math
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_ipomopsis_higher_order_selection_v1 as v

RECEIPT=ROOT/"data/evidence/ipomopsis_source_higher_order_selection_independent_js_v1.json"


class IpomopsisAuthenticNonlinearSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=json.loads(RECEIPT.read_text(encoding="utf-8"))
        cls.m={x["model_id"]:x for x in cls.r["models"]}

    def test_author_study_not_an_unverified_new_thistle_trial(self):
        self.assertEqual(self.r["source"]["paper"],"10.1086/716740")
        self.assertEqual(self.r["source"]["author_file_md5"],
                         "12ed3dbe6d269c5f808f9c6b230d468c")
        self.assertEqual(self.r["source"]["source_original_records"],437)
        self.assertTrue(self.r["source"]["source_original_plant_ids_unavailable"])
        self.assertTrue(self.r["limits"]["different_organism_than_Cirsium"])

    def test_real_life_stage_does_not_use_a_fictitious_seed_product(self):
        self.assertEqual(self.r["design"]["primary_stage_responses"],
            ["seed_initiation_log1p","fruit_fly_escape_fraction"])
        self.assertTrue(self.r["limits"]["no_verified_unique_plant_ids"])
        self.assertTrue(self.r["limits"][
            "whole_outcome_observational_and_exposure_not_randomized"])

    def test_exact_models_are_one_quadratic_plus_two_floral_pairs(self):
        self.assertEqual([x["id"] for x in v.MODELS],
            ["width_quadratic","width_by_color","width_by_sepal"])
        self.assertEqual(self.r["design"]["models"],
            [x["id"] for x in v.MODELS])
        self.assertEqual(self.r["multiplicity"]["n_tests"],6)
        self.assertEqual(v.PERMUTATIONS,5000)
        self.assertEqual(v.BOOTSTRAPS,1500)

    def test_source_original_no_plant_id_and_valid_label_covariates(self):
        self.assertEqual(v.HEADERS,
            ("Site","Year","Type","width","color","nectar","sepal",
             "totpinene","indoleN","seeds","seedsinit","flyescape"))
        self.assertTrue(self.r["design"]
            ["original_site_labels_are_not_distinct_physical_sites"])
        self.assertEqual(self.r["design"]["year_levels"],list(range(2010,2019)))

    def test_source_head_and_sepal_missingness_reduce_trait_cohort(self):
        self.assertEqual(
            [(m["model_id"],m["n_original_joint_records"],
              m["n_eligible_original_source_rows"],
              m["n_original_site_year_type_groups"]) for m in self.r["models"]],
            [("width_quadratic",367,354,34),
             ("width_by_color",301,288,29),
             ("width_by_sepal",96,88,14)])

    def test_actual_inner_corolla_curvature_small_and_unstable(self):
        w=self.m["width_quadratic"]
        init=w["endpoints"]["seed_initiation_log1p"]
        escape=w["endpoints"]["fruit_fly_escape_fraction"]
        self.assertAlmostEqual(init["partial_R2"],.011900458290506002)
        self.assertAlmostEqual(escape["partial_R2"],.0042590421485211435)
        self.assertAlmostEqual(init["extra_standardized_beta"],
                               .07775729151264825)
        self.assertAlmostEqual(init["freedman_lane_within_source_group_p"],
                               .0643871225754849)
        self.assertGreater(init["bootstrap_beta_ci95"][0]*init["bootstrap_beta_ci95"][1],-1e6)
        self.assertLess(init["bootstrap_beta_ci95"][0],0)
        self.assertGreater(init["bootstrap_beta_ci95"][1],0)

    def test_true_pairwise_width_color_does_not_robustly_infer_LEGO(self):
        w=self.m["width_by_color"]
        self.assertAlmostEqual(w["endpoints"]["seed_initiation_log1p"]
            ["partial_R2"],.005054805728157206)
        self.assertGreater(w["endpoints"]["seed_initiation_log1p"]
            ["freedman_lane_within_source_group_p"],.2)
        self.assertGreater(w["endpoints"]["fruit_fly_escape_fraction"]
            ["freedman_lane_within_source_group_p"],.6)

    def test_true_width_sepal_source_subset_is_sparse(self):
        s=self.m["width_by_sepal"]
        self.assertEqual(s["n_eligible_original_source_rows"],88)
        self.assertEqual(s["n_original_site_year_type_groups"],14)
        self.assertAlmostEqual(s["endpoints"]["seed_initiation_log1p"]
                               ["partial_R2"],.026229981135367997)
        self.assertGreater(s["endpoints"]["seed_initiation_log1p"]
                           ["freedman_lane_within_source_group_p"],.05)
        self.assertTrue(self.r["limits"][
            "original_authors_multivariate_selection_results_not_refuted"])

    def test_six_exploratory_conditional_comparisons_are_multiplicity_guarded(self):
        q=self.r["multiplicity"]["six_source_exploratory_tests_BH_q"]
        self.assertEqual(len(q),6)
        self.assertTrue(all(.05 < x <= 1 for x in q.values()))
        self.assertAlmostEqual(
            q["width_quadratic|seed_initiation_log1p"],.3392921415716857)

    def test_actual_year_holdout_improvement_mostly_negative(self):
        stats={
          (model["model_id"],stage):detail["heldout_year_relative_SSE_gain"]
          for model in self.r["models"]
          for stage,detail in model["leave_one_year_out_summary"].items()
        }
        self.assertEqual(len(stats),6)
        self.assertEqual(sum(x>0 for x in stats.values()),1)
        self.assertAlmostEqual(
            stats[("width_quadratic","seed_initiation_log1p")],
            -.05769615258402491)
        self.assertAlmostEqual(
            stats[("width_by_color","fruit_fly_escape_fraction")],
            -.00663522918300958)
        self.assertAlmostEqual(
            stats[("width_by_sepal","seed_initiation_log1p")],
            .005116776664043072)

    def test_actual_group_rows_are_eligible_within_original_source_strata(self):
        for m in self.r["models"]:
            g=m["model_source_group_rows"]
            self.assertEqual(sum(x["n"] for x in g),
                             m["n_eligible_original_source_rows"])
            self.assertEqual(len(g),m["n_original_site_year_type_groups"])
            ids=[i for z in g for i in z["source_file_rows"]]
            self.assertEqual(len(ids),len(set(ids)))
            self.assertTrue(all(2<=i<=438 for i in ids))
            self.assertTrue(all(x["n"]>=2 for x in g))

    def test_math_solver_and_null_identity_negative_control(self):
        x=[[1,1],[1,2],[1,3],[1,4]]
        y=[1,2,3,4]
        a=v.fit(x,y,1)
        b=v.fit(x,y,2)
        self.assertIsNotNone(a)
        self.assertIsNotNone(b)
        self.assertLess(abs(b["sse"]),1e-12)
        self.assertGreater(a["sse"],0)
        self.assertIsNone(v.solve([[0,0],[0,0]],[0,0]))

    def test_zero_rng_seed_rejected(self):
        with self.assertRaisesRegex(ValueError,"INVALID_OR_ZERO_SEED"):
            next(v.xorshift32(0))
        g=v.xorshift32(20261025)
        h=v.xorshift32(20261025)
        self.assertEqual([next(g) for _ in range(12)],
                         [next(h) for _ in range(12)])

    def test_invalid_fake_author_file_rejected(self):
        with self.assertRaisesRegex(ValueError,"ORIGINAL_AUTHOR_ZENODO_CSV_MD5_CHANGED"):
            v.read_author_csv(
                b"Site,Year,Type,width,color,nectar,sepal,totpinene,indoleN,seeds,seedsinit,flyescape\n")

    def test_source_bootstrap_nonnegative_partial_r_squared_is_not_significance(self):
        for model in self.r["models"]:
            for val in model["endpoints"].values():
                self.assertTrue(
                    val["n_valid_bootstrap"]<=v.BOOTSTRAPS)
                self.assertTrue(
                    0<=val["partial_R2"]<=1)
                self.assertIn("not a test",val["note"])
        self.assertTrue(self.r["limits"]
            ["bootstrap_partial_R2_interval_not_evidence_against_zero"])

    def test_source_context_prevents_false_evolutionary_reexpression(self):
        self.assertTrue(self.r["limits"]
            ["quadratic_extra_term_is_not_an_evolutionary_reversion"])
        self.assertTrue(self.r["limits"]
            ["interaction_geometry_not_proven_Pedicularis_or_Cirsium_mechanism"])
        self.assertEqual(self.r["limits"]["no_new_field_observations"],0)

if __name__=="__main__":
    unittest.main()
