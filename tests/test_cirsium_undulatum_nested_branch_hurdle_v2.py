"""Same-original-branch Cirsium seed hurdle: provenance and inference guards."""
from pathlib import Path
import copy
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import audit_cirsium_undulatum_nested_branch_hurdle_v2 as v

SRC = ROOT / "data/external/usda_wavyleaf_2007_2008"


class NestedBranchHurdleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mothers = v.read_locked(
            SRC / "mother_plants.csv", "mother_plants.csv",
            ("plt", "year", "dam", "ins", "trt")
        )
        cls.heads = v.read_locked(
            SRC / "flower_heads.csv", "flower_heads.csv",
            ("plt", "year", "dam", "ins", "hd", "av.dam", "hd.gsd")
        )
        cls.dictionary = v.read_locked(SRC / "source_dictionary.txt",
                                       "source_dictionary.txt")
        cls.results = v.analyze(cls.mothers, cls.heads, cls.dictionary, 5000)

    def test_original_source_population_and_20_real_branch_pairs(self):
        r = self.results
        self.assertEqual(r["source_audit"]["n_original_mothers"], 106)
        self.assertEqual(r["source_audit"]["n_original_head_rows"], 323)
        self.assertEqual(r["source_audit"]["n_clean_head_rows"], 297)
        self.assertEqual(r["n_exact_branch_pairs"], 20)
        self.assertEqual(r["n_original_mothers_with_exact_branch_pairs"], 14)
        self.assertEqual(r["unmatched_intermediate_branch_ids"],
                         ["2007|355|3"])

    def test_each_pair_has_actual_terminal_and_intermediate_on_same_branch(self):
        self.assertEqual(len({p["id"] for p in self.results["paired_branches"]}),
                         20)
        for p in self.results["paired_branches"]:
            k = int(p["id"].split("|")[2])
            self.assertEqual(int(p["terminal_head_code"]), 100 * k)
            self.assertEqual(int(p["intermediate_head_code"]), 100 * k + 10)
            self.assertEqual(p["good_seed_delta"],
                             p["intermediate_good_seed"]
                             - p["terminal_good_seed"])
            self.assertTrue(0 <= p["terminal_damage"] <= 1)
            self.assertTrue(0 <= p["intermediate_damage"] <= 1)

    def test_original_year_route_hurdle_and_branch_signs(self):
        r = self.results["by_year"]
        self.assertEqual(r["2007"]["n_paired_branches"], 11)
        self.assertEqual(r["2008"]["n_paired_branches"], 9)
        self.assertEqual(r["2007"]["branch_seed_routes"], {
            "neither": 0, "terminal_only": 2,
            "intermediate_only": 1, "both": 8})
        self.assertEqual(r["2008"]["branch_seed_routes"], {
            "neither": 1, "terminal_only": 5,
            "intermediate_only": 0, "both": 3})
        self.assertEqual(r["2007"]["paired_branch_seed_directions"],
                         {"negative": 6, "positive": 4, "equal": 1})
        self.assertEqual(r["2008"]["paired_branch_seed_directions"],
                         {"negative": 7, "positive": 1, "equal": 1})
        self.assertEqual(r["2008"]["paired_branch_damage_directions"],
                         {"negative": 2, "positive": 3, "equal": 4})
        self.assertEqual(r["2008"]["n_terminal_good_seed_positive"], 8)
        self.assertEqual(r["2008"]["n_intermediate_good_seed_positive"], 3)

    def test_branch_mean_not_naively_used_as_independent_mother_test(self):
        y = self.results["by_year"]["2008"]
        self.assertEqual(y["n_paired_mothers"], 7)
        self.assertAlmostEqual(y["paired_branch_seed_mean_delta"], -310/9)
        self.assertAlmostEqual(y["paired_mother_mean_seed_delta"], -283/7)
        self.assertEqual(y["paired_mother_seed_directions"],
                         {"negative": 7, "positive": 0, "equal": 0})
        self.assertEqual(y["posthoc_mother_exact_two_sided_sign_p"],
                         0.015625)
        self.assertEqual(self.results["by_year"]["2007"]
                         ["paired_mother_seed_directions"],
                         {"negative": 4, "positive": 2, "equal": 1})

    def test_positive_exception_in_real_branch_not_erased_by_mother_average(self):
        p = {x["id"]: x for x in self.results["paired_branches"]}
        self.assertEqual(p["2008|WR56|3"]["terminal_good_seed"], 68)
        self.assertEqual(p["2008|WR56|3"]["intermediate_good_seed"], 85)
        self.assertEqual(p["2008|WR56|3"]["good_seed_delta"], 17)
        self.assertEqual(p["2008|WR56|3"]["damage_delta"], 0.0)
        self.assertEqual(p["2008|SD05|2"]["route"], "neither")
        self.assertAlmostEqual(
            self.results["by_year"]["2008"]
            ["mean_seed_delta_conditioned_on_both_positive"], -70/3
        )

    def test_outlier_leave_one_original_mother_not_new_fake_individual(self):
        deltas = self.results["by_year"]["2008"]["paired_mother_seed_deltas"]
        self.assertEqual(len(deltas), 7)
        self.assertEqual(deltas["2008|WG25"], -130)
        remainder = [v for k, v in deltas.items() if k != "2008|WG25"]
        self.assertAlmostEqual(sum(remainder)/len(remainder), -25.5)

    def test_source_does_not_identify_head_age_or_genetic_reexpression(self):
        lim = self.results["limits"]
        self.assertTrue(lim["head_code_not_flowering_age_or_gravity_orientation"])
        self.assertTrue(lim["cannot_separate_pollen_resource_phenology_insect_causes"])
        self.assertTrue(lim["no_genetic_lego_or_evolutionary_regain_test"])
        self.assertEqual(lim["new_field_observations"], 0)

    def test_negative_control_dictionary_change(self):
        with self.assertRaisesRegex(ValueError, "SOURCE_HEAD_CODE_ONTOLOGY"):
            v.build_pairs(self.mothers, self.heads,
                          self.dictionary.replace(
                              "between the terminals and stem",
                              "not independently verified"))

    def test_negative_control_off_dictionary_head_digit(self):
        h = copy.deepcopy(self.heads)
        target = next(x for x in h if x["year"] == "2008"
                      and x["plt"] == "WR56" and x["hd"] == "310")
        target["hd"] = "320"
        with self.assertRaisesRegex(ValueError, "HEAD_CODE_NOT_SINGLE_BRANCH"):
            v.build_pairs(self.mothers, h, self.dictionary)

    def test_negative_control_falsified_good_seed_count(self):
        h = copy.deepcopy(self.heads)
        target = next(x for x in h if x["year"] == "2008"
                      and x["plt"] == "WR56" and x["hd"] == "310")
        target["hd.gsd"] = "-9"
        with self.assertRaisesRegex(ValueError, "GOOD_SEED_COUNT"):
            v.build_pairs(self.mothers, h, self.dictionary)

    def test_negative_control_branch_replication_count(self):
        with self.assertRaisesRegex(ValueError, "ORIGINAL_20_SAME_BRANCH"):
            h = copy.deepcopy(self.heads)
            target = next(x for x in h if x["year"] == "2008"
                          and x["plt"] == "WR56" and x["hd"] == "310")
            # Keep raw row count and intermediate count unchanged, but
            # reassign to a branch with no matching terminal head.
            target["hd"] = "410"
            v.build_pairs(self.mothers, h, self.dictionary)

    def test_negative_control_bootstrap_budget(self):
        with self.assertRaisesRegex(ValueError, "INVALID_RESAMPLING_BUDGET"):
            v.yearly_results(self.results["paired_branches"], 0)


if __name__ == "__main__":
    unittest.main()
