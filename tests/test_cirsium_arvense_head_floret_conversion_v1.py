"""Original-data safeguards for Cirsium head florets × mature seed accounting."""
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_arvense_head_floret_conversion_v1 as v

RECEIPT=ROOT/"data/evidence/cirsium_arvense_head_floret_conversion_independent_js_v1.json"

def mk(F,V,N,D,M,heads=1):
    assert F==V+N+D+M
    return [{
      "Flowers":str(F),"Viable_Seeds":str(V),
      "NViable_Seeds":str(N),"Damaged_Seeds":str(D),
      "Missing_Seeds":str(M),
    } for i in range(heads)]

def plant(F,V):
    n=F-V
    return v.plant_head_counts(mk(F,V,n,0,0))

def example():
    out=[]
    for site in v.SITES:
        out.append({"site":site,"block":"B1","arms":{
            "Control":plant(100,80),
            "Pollinator Exclusion":plant(100,10),
            "Predator Exclusion":plant(90,70),
            "Predator Model":plant(95,75),
        }})
    return out

class CirsiumHeadStageTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt=json.loads(RECEIPT.read_text(encoding="utf-8"))

    def test_exact_two_factor_product(self):
        z=v.simulate_control()
        self.assertEqual(z,"PASS_SOURCE_FLORET_PRODUCT_AND_TWO_CATEGORY_IDENTITIES")
        a=plant(100,70)
        b=plant(80,40)
        q=v.compare_one(a,b,"campus","B1")
        self.assertEqual(q["good_seeds_per_head_delta"],30)
        self.assertAlmostEqual(q["floret_quantity_component"],12)
        self.assertAlmostEqual(q["good_fraction_component"],18)
        self.assertAlmostEqual(q["floret_quantity_component"]+
                               q["good_fraction_component"],
                               q["good_seeds_per_head_delta"])

    def test_zero_and_negative_source_missing_are_not_automatically_imputed(self):
        a=v.plant_head_counts(mk(10,7,5,0,-2))
        self.assertEqual(a["signed_missing"],-2)
        self.assertAlmostEqual(a["signed_missing_fraction"],-.2)
        self.assertAlmostEqual(a["florets_per_head"]*a["good_fraction"],
                               a["good_seeds_per_head"])

    def test_bad_sum_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"SOURCE_FLORET_CATEGORY_IDENTITY"):
            v.plant_head_counts([{"Flowers":"10","Viable_Seeds":"7",
              "NViable_Seeds":"1","Damaged_Seeds":"1","Missing_Seeds":"0"}])

    def test_one_inferred_plant_group_cannot_be_empty(self):
        with self.assertRaisesRegex(ValueError,"UNOBSERVED_HEAD_GROUP"):
            v.plant_head_counts([])

    def test_source_grain_not_confused_with_615_independent_experiments(self):
        a=self.receipt["source_counts"]
        self.assertEqual((a["raw_heads"],a["eligible_heads"],
                          a["plant_proxies"],a["site_blocks"],
                          a["complete_blocks"]),(615,612,164,45,36))
        self.assertEqual(a["complete_sites"],
                         {"campus":12,"construction":9,"ditch":8,"rapseed":7})
        self.assertEqual(self.receipt["source"]["seeds_csv_git_blob"],
                         v.SOURCE_GIT_BLOB)

    def test_bag_good_seed_rate_dominates_floral_count_component(self):
        m=self.receipt["comparisons"]["bag_minus_control"]["measurements"]
        self.assertAlmostEqual(m["good_seeds_per_head_delta"]["mean"],
                               -41.58818342151675)
        self.assertAlmostEqual(m["good_fraction_component"]["mean"],
                               -41.25479554408028)
        self.assertAlmostEqual(m["floret_quantity_component"]["mean"],
                               -.33338787743648635)
        self.assertAlmostEqual(self.receipt[
            "bag_fraction_component_share_of_mean_good_seed_difference"],
            0.9919835912509706)
        self.assertEqual(m["good_fraction_delta"]["sign"],
                         {"negative":36,"positive":0,"zero":0})

    def test_four_sites_differ_in_floret_contribution_sign(self):
        m=self.receipt["comparisons"]["bag_minus_control"]["measurements"]
        q=m["floret_quantity_component"]["site_means"]
        self.assertLess(q["campus"],0)
        self.assertGreater(q["construction"],0)
        self.assertLess(q["ditch"],0)
        self.assertGreater(q["rapseed"],0)
        self.assertTrue(all(v<0 for v in
            m["good_fraction_component"]["site_means"].values()))

    def test_original_category_identity_is_exact_not_causal_partition(self):
        for comp in self.receipt["comparisons"].values():
            self.assertEqual(comp["n_blocks"],36)
            self.assertAlmostEqual(comp["total_identity_mean_residual"],0,places=9)
            self.assertAlmostEqual(
                comp["good_fraction_categories_identity_mean_residual"],0,
                places=9)
            for b in comp["paired_block_rows"]:
                self.assertAlmostEqual(
                    b["good_seeds_per_head_delta"],
                    b["floret_quantity_component"]+
                    b["good_fraction_component"],places=9)
                self.assertAlmostEqual(
                    b["good_fraction_delta"]+
                    b["nonviable_fraction_delta"]+
                    b["direct_damage_fraction_delta"]+
                    b["signed_missing_fraction_delta"],0,places=9)

    def test_bag_category_components_not_direct_seed_feeding(self):
        m=self.receipt["comparisons"]["bag_minus_control"]["measurements"]
        self.assertAlmostEqual(m["nonviable_fraction_delta"]["mean"],
                               0.38871814419478146)
        self.assertAlmostEqual(m["direct_damage_fraction_delta"]["mean"],
                               -0.035998803120399506)
        self.assertAlmostEqual(m["signed_missing_fraction_delta"]["mean"],
                               0.17034697134311352)
        self.assertTrue(self.receipt["limits"][
            "signed_missing_cannot_be_called_direct_consumption"])

    def test_synthetic_fixed_sites_with_no_effect_bootstrap(self):
        blocks=example()
        out=v.summarize_blocks(blocks)
        out=v.bootstrap(out,blocks,reps=100)
        q=out["bag_minus_control"]["measurements"][
            "good_seeds_per_head_delta"]["bootstrap_ci95"]
        self.assertEqual(q,[-70,-70])

    def test_missing_study_site_is_not_imputed_with_zero(self):
        with self.assertRaisesRegex(ValueError,"FOUR_ORIGINAL_SITES_REQUIRED"):
            v.summarize_blocks(example()[:-1])

    def test_bootstrap_budget_and_rng_guards(self):
        with self.assertRaisesRegex(ValueError,"INVALID_BOOTSTRAP_DRAWS"):
            v.bootstrap(v.summarize_blocks(example()),example(),reps=0)
        with self.assertRaisesRegex(ValueError,"INVALID_RNG_SEED"):
            next(v.rng32(0))

    def test_both_plant_metrics_are_observed_post_exposure(self):
        lim=self.receipt["limits"]
        self.assertTrue(lim["not_true_pre_treatment_floret_count"])
        self.assertTrue(lim["flower_count_is_author_Flowers_not_calibrated_head_geometry"])
        self.assertTrue(lim["ratio_of_components_not_separate_causal_path_effects"])
        self.assertTrue(lim["bag_blocks_both_pollinators_and_some_antagonists"])
        self.assertEqual(lim["zero_new_fieldwork_or_anatomical_measurements"],0)

    def test_original_factorial_four_arm_comparison_count(self):
        self.assertEqual(len(self.receipt["comparisons"]),5)
        self.assertEqual(len(v.METRICS),8)
        self.assertEqual(v.N_BOOT,10000)
        self.assertEqual(v.SEED,20261023)

if __name__=="__main__":
    unittest.main()
