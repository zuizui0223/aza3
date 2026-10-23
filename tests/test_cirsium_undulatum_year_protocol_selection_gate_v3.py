"""USDA 106-mother selection denominator and 2012 external protocol guard."""
import copy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import audit_cirsium_undulatum_year_protocol_selection_gate_v3 as v

RAW = ROOT / "data/external/usda_wavyleaf_2007_2008"


class YearProtocolSelectionGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = v.read_locked(
            RAW / "mother_plants.csv", "mother_plants.csv",
            ("year", "plt", "trt", "plt.gsd", "r.diam"))
        cls.h = v.read_locked(
            RAW / "flower_heads.csv", "flower_heads.csv",
            ("year", "plt", "dam", "ins", "hd", "hd.gsd", "av.dam"))
        cls.dictionary = v.read_locked(
            RAW / "source_dictionary.txt", "source_dictionary.txt")
        cls.r = v.analyze(cls.m, cls.h, cls.dictionary)

    def test_original_106_mothers_only_fourteen_selected(self):
        d = self.r["original_denominator_selection"]
        self.assertEqual((d["2007"]["all_original_mothers"],
                          d["2008"]["all_original_mothers"]), (39, 67))
        self.assertEqual((d["2007"]["mothers_with_both_branch_code_positions"],
                          d["2008"]["mothers_with_both_branch_code_positions"]),
                         (7, 7))
        self.assertAlmostEqual(d["2007"]["fraction_with_both"], 7/39)
        self.assertAlmostEqual(d["2008"]["fraction_with_both"], 7/67)
        self.assertEqual(self.r["original_source_audit"]["n_clean_head_rows"],
                         297)

    def test_seed_production_selection_is_only_descriptive(self):
        d = self.r["original_denominator_selection"]
        for y, selected, other in (
                ("2007", 102.57142857142857, 54.3125),
                ("2008", 147.71428571428572, 81.75)):
            self.assertAlmostEqual(
                d[y]["selected"]["author_good_seed_total"]["mean"], selected)
            self.assertAlmostEqual(
                d[y]["not_selected"]["author_good_seed_total"]["mean"], other)
            self.assertGreater(selected, other)
        self.assertTrue(self.r["limits"][
            "conditioning_on_observed_both_head_types_is_post_outcome_selection"])

    def test_rosette_nonmissing_and_no_missing_equals_zero(self):
        d = self.r["original_denominator_selection"]
        self.assertEqual(d["2007"]["selected"]
                         ["source_rosette_diameter_nonmissing"]["n_observed"], 7)
        self.assertEqual(d["2007"]["not_selected"]
                         ["source_rosette_diameter_nonmissing"]["n_observed"], 29)
        self.assertEqual(d["2008"]["selected"]
                         ["source_rosette_diameter_nonmissing"]["n_observed"], 7)
        self.assertEqual(d["2008"]["not_selected"]
                         ["source_rosette_diameter_nonmissing"]["n_observed"], 49)
        self.assertAlmostEqual(d["2007"]["selected"]
                               ["source_rosette_diameter_nonmissing"]["mean"],
                               37.142857142857146)
        self.assertAlmostEqual(d["2008"]["selected"]
                               ["source_rosette_diameter_nonmissing"]["mean"],
                               30.142857142857142)
        self.assertTrue(self.r["limits"][
            "rosette_timing_in_current_export_not_ascertained"])

    def test_four_2007_randomization_arms_never_replaced_by_selected_heads(self):
        z = self.r["original_denominator_selection"]["2007"]
        a = z["by_year_treatment_arm"]
        self.assertEqual(
            [a[x]["n_original_mothers"] for x in v.ARMS], [10, 9, 11, 9])
        self.assertEqual(
            [a[x]["n_with_matching_branch_codes"] for x in v.ARMS],
            [1, 1, 4, 1])
        self.assertEqual(sum(a[x]["n_original_mothers"] for x in v.ARMS), 39)

    def test_four_2008_randomization_arms_never_replaced_by_selected_heads(self):
        z = self.r["original_denominator_selection"]["2008"]
        a = z["by_year_treatment_arm"]
        self.assertEqual(
            [a[x]["n_original_mothers"] for x in v.ARMS], [21, 15, 19, 12])
        self.assertEqual(
            [a[x]["n_with_matching_branch_codes"] for x in v.ARMS],
            [1, 1, 3, 2])
        self.assertEqual(sum(a[x]["n_original_mothers"] for x in v.ARMS), 67)

    def test_2008_one_pair_per_mother_downshifts_apparent_binary_result(self):
        z = self.r["nested_branch_binary_route_sensitivity"]["2008"]
        a = z["lowest_code_one_per_mother"]
        self.assertEqual(z["all_actual_recorded_branch_code_pairs"], 9)
        self.assertEqual(a["n_distinct_mothers"], 7)
        self.assertEqual(a["route_counts"],
                         {"neither": 1, "terminal_only": 4,
                          "intermediate_only": 0, "both": 2})
        self.assertEqual(
            a["two_sided_exact_discordant_symmetry_reference"], 0.125)
        b = z["highest_code_one_per_mother"]
        self.assertEqual(b["route_counts"],
                         {"neither": 0, "terminal_only": 5,
                          "intermediate_only": 0, "both": 2})
        self.assertEqual(
            b["two_sided_exact_discordant_symmetry_reference"], 0.0625)
        self.assertEqual(z["all_pair_route_counts_descriptive_not_independent"],
                         {"neither": 1, "terminal_only": 5,
                          "intermediate_only": 0, "both": 3})

    def test_2007_one_pair_per_mother_can_change_discordant_type(self):
        z = self.r["nested_branch_binary_route_sensitivity"]["2007"]
        self.assertEqual(z["all_actual_recorded_branch_code_pairs"], 11)
        a = z["lowest_code_one_per_mother"]
        self.assertEqual(a["route_counts"],
                         {"neither": 0, "terminal_only": 2,
                          "intermediate_only": 1, "both": 4})
        self.assertEqual(
            a["two_sided_exact_discordant_symmetry_reference"], 1.0)
        b = z["highest_code_one_per_mother"]
        self.assertEqual(b["route_counts"],
                         {"neither": 0, "terminal_only": 2,
                          "intermediate_only": 0, "both": 5})
        self.assertEqual(
            b["two_sided_exact_discordant_symmetry_reference"], 0.5)

    def test_2012_documented_year_protocol_is_not_ecology_only(self):
        p = self.r["manual_external_source_method_provenance"]
        self.assertEqual(p["2007"]["insecticide_active_ingredient"],
                         "acephate")
        self.assertEqual(p["2008"]["insecticide_active_ingredient"],
                         "bifenthrin")
        self.assertEqual(p["2007"]["application_interval_days"], "7-10")
        self.assertEqual(p["2008"]["application_interval_days"], "14")
        self.assertEqual(p["2007"]["flowerhead_census_interval_weeks"], 2)
        self.assertEqual(p["2008"]["flowerhead_census_interval_weeks"], 4)
        self.assertTrue(self.r["limits"][
            "year_environment_protocol_is_perfectly_conflated"])
        self.assertTrue(self.r["limits"][
            "2012_method_metadata_is_externally_manually_verified_not_CSV_encoded"])

    def test_original_dissertation_location_does_not_fake_same_head_image(self):
        p = self.r["original_sources"]
        self.assertTrue(
            p["figure_2_2_subordinate_heads_below_terminal_on_same_branch"])
        self.assertTrue(self.r["limits"][
            "dissertation_figure_code_correspondence_not_per_head_photo_join"])
        self.assertIn("bioscidiss/41", p["dissertation_2012_chapter2"])

    def test_exact_reference_small_discordant_cell_counts(self):
        self.assertEqual(v.exact_discordant_symmetry(4, 0), 0.125)
        self.assertEqual(v.exact_discordant_symmetry(5, 0), 0.0625)
        self.assertEqual(v.exact_discordant_symmetry(2, 1), 1.0)
        self.assertEqual(v.exact_discordant_symmetry(0, 0), 1.0)
        self.assertEqual(v.exact_discordant_symmetry(6, 1), 0.125)
        with self.assertRaisesRegex(ValueError, "INVALID_CODE_CHOICE_RULE"):
            v.select_one_branch_per_mother([], "2008", "largest_seed")

    def test_negative_mother_good_seed_rejected(self):
        m = copy.deepcopy(self.m)
        m[0]["plt.gsd"] = "-5"
        with self.assertRaisesRegex(ValueError, "MOTHER_GOOD_SEED_COUNT"):
            v.analyze(m, self.h, self.dictionary)

    def test_negative_mother_rosette_rejected(self):
        m = copy.deepcopy(self.m)
        m[0]["r.diam"] = "-10"
        with self.assertRaisesRegex(ValueError, "INVALID_ROSETTE_DIAMETER"):
            v.analyze(m, self.h, self.dictionary)

    def test_negative_duplicate_mother_identity_rejected(self):
        m = copy.deepcopy(self.m)
        m[0]["plt"] = m[1]["plt"]
        with self.assertRaisesRegex(ValueError, "DUPLICATE_ORIGINAL_MOTHER"):
            v.analyze(m, self.h, self.dictionary)

    def test_new_biological_data_or_causal_selection_not_claimed(self):
        limit = self.r["limits"]
        self.assertEqual(limit["new_field_observations"], 0)
        self.assertTrue(limit["no_supported_evolutionary_adaptation_genetic_lego_or_trait_regain"])
        self.assertTrue(limit["separate_mechanisms_pollen_injury_resources_not_identified"])


if __name__ == "__main__":
    unittest.main()
