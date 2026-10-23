"""Standalone regression tests: source-grain, joint response, and causal limits."""
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import audit_cirsium_arvense_joint_access_response_v1 as v


def head(viable=.8, insects=2, damaged=.1, perhead=50):
    return {
        "viable_fraction": viable,
        "emerged_insects_per_head": insects,
        "observed_damaged_fraction": damaged,
        "viable_per_head": perhead,
    }


def sample_complete():
    arr = []
    for j, site in enumerate(v.SITES):
        arr.append({
            "site": site, "block": "B" + str(j+1),
            "arms": {
                "Control": head(.8, 2, .1, 50),
                "Pollinator Exclusion": head(.3, .5, .03, 18),
                "Predator Exclusion": head(.7, 1.5, .1, 45),
                "Predator Model": head(.65, 2, .1, 44),
            },
        })
    return arr


class JointAccessSourceTest(unittest.TestCase):
    def test_source_known_contrasts_and_site_order(self):
        self.assertEqual(len(v.COMPARISONS), 5)
        self.assertEqual(len(v.METRICS), 4)
        self.assertEqual(tuple(v.SITES),
                         ("campus", "construction", "ditch", "rapseed"))
        self.assertEqual(v.SEED, 20261023)
        self.assertEqual(v.DRAWS, 10000)

    def test_two_outcome_sign_does_not_equate_insect_loss_to_seed_gain(self):
        b = sample_complete()
        c = v.COMPARISONS[0]
        self.assertEqual(
            v.joint_counts(b, c),
            {"viable_down__emerged_down": 4},
        )
        c = v.COMPARISONS[1]
        self.assertEqual(
            v.joint_counts(b, c),
            {"viable_down__emerged_down": 4},
        )

    def test_site_and_original_block_comparisons(self):
        d, by_site = v.summarize(sample_complete())
        self.assertEqual(len(by_site), 4)
        self.assertEqual(d["bag_minus_control"]["n_original_matched_blocks"], 4)
        z = d["bag_minus_control"]["measures"]["viable_fraction"]
        self.assertAlmostEqual(z["mean"], -.5)
        self.assertEqual(z["sign"],
                         {"positive": 0, "negative": 4, "zero": 0})
        p = d["pred_model_minus_pred_exclusion"]["measures"]["viable_fraction"]
        self.assertAlmostEqual(p["mean"], -.05)

    def test_reproducible_site_stratified_rng(self):
        a, by_site = v.summarize(sample_complete())
        b, by_site2 = v.summarize(sample_complete())
        x = v.bootstrap(a, by_site, draws=200)
        y = v.bootstrap(b, by_site2, draws=200)
        self.assertEqual(x, y)
        self.assertEqual(
            x["bag_minus_control"]["measures"]["viable_fraction"]
            ["site_stratified_bootstrap_ci95"], [-.5, -.5])

    def test_missing_site_must_fail_closed(self):
        vals = [x for x in sample_complete() if x["site"] != "ditch"]
        with self.assertRaisesRegex(ValueError, "ALL_FOUR_SITES"):
            v.summarize(vals)

    def test_at_least_two_sites_required(self):
        with self.assertRaisesRegex(ValueError, "MULTIPLE_ORIGINAL_SITES"):
            v.summarize(sample_complete()[:1])

    def test_missing_fourth_arm_does_not_become_zero(self):
        a = sample_complete()
        del a[0]["arms"]["Predator Model"]
        with self.assertRaises(KeyError):
            v.summarize(a)

    def test_invalid_bootstrap_replicates_fail_closed(self):
        data, s = v.summarize(sample_complete())
        with self.assertRaisesRegex(ValueError, "INVALID_BOOTSTRAP_DRAWS"):
            v.bootstrap(data, s, draws=0)

    def test_zero_rng_seed_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "INVALID_XORSHIFT_SEED"):
            next(v.xorshift32(0))

    def test_synthetic_structural_and_sign_controls(self):
        self.assertEqual(
            v.synthetic_controls(),
            "PASS_MUTUALLY_EXCLUSIVE_JOINT_SIGNS_AND_PREDATOR_NEGATIVE_CONTROL",
        )
        self.assertEqual(v.sign_counts([-1, 0, 2]),
                         {"positive": 1, "negative": 1, "zero": 1})

    def test_independent_original_source_receipt_has_36_not_615_independent_blocks(self):
        p = ROOT / "data/evidence/cirsium_arvense_four_arm_joint_response_independent_js_v1.json"
        z = json.loads(p.read_text(encoding="utf-8"))
        self.assertEqual(z["source"]["seeds_git_blob_sha"], v.SOURCE_GIT_BLOB)
        self.assertEqual(z["source_counts"]["heads_all"], 615)
        self.assertEqual(z["source_counts"]["plant_proxy_groups"], 164)
        self.assertEqual(z["source_counts"]["complete_four_arm_blocks"], 36)
        self.assertEqual(z["source_counts"]["site_complete_counts"],
                         {"campus":12,"construction":9,"ditch":8,"rapseed":7})

    def test_original_observed_joint_quadrants_and_both_partner_limits(self):
        p = ROOT / "data/evidence/cirsium_arvense_four_arm_joint_response_independent_js_v1.json"
        z = json.loads(p.read_text(encoding="utf-8"))
        bag = z["comparisons"]["bag_minus_control"]
        self.assertEqual(bag["joint_response_counts"], {
            "viable_down__emerged_down": 25,
            "viable_down__emerged_equal": 8,
            "viable_down__emerged_up": 3,
        })
        self.assertEqual(
            bag["measures"]["viable_fraction"]["sign"],
            {"positive": 0, "negative": 36, "zero": 0})
        self.assertTrue(z["limits"][
            "pollinator_bag_compound_access_and_microclimate_intervention"])
        self.assertTrue(z["limits"][
            "no_true_cirsium_phyllary_spine_or_head_orientation_measurement"])
        self.assertTrue(z["limits"][
            "predator_exclusion_stem_adhesive_plus_manual_removal"])
        self.assertEqual(z["limits"]["new_biological_field_observations"], 0)


if __name__ == "__main__":
    unittest.main()
