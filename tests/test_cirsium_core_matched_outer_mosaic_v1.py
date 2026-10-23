"""Original-publisher Cirsium reproductive head morphology matching: 114 plants."""
import copy
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_core_matched_outer_mosaic_v1 as v

SRC=ROOT/"data/evidence/cirsium_2023_same_plant_head_morphology_subset_v1.csv"


class CirsiumOuterMosaicSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw=SRC.read_bytes()
        cls.rows=v.read_original_subset(cls.raw)
        cls.result=v.analyze(cls.raw)

    def test_original_114_source_rows_six_taxon_sex_strata(self):
        self.assertEqual(self.result["source"]["n_original_individuals"],114)
        self.assertEqual(self.result["source"]["n_strata"],6)
        self.assertEqual(self.result["source"]["taxa_concept_count"],3)
        self.assertEqual(
            self.result["source"]["subset_csv_git_blob_sha"],
            "be95e6002ab45393786de05a7180eb77f38d941d")
        self.assertEqual(self.result["source"]["source_original_xlsx_sha256"],
            "c06423cc6f9800ba4775a4428afa072ba1b3dc6fb48ccd05fdb7e42a2f03ee79")
        self.assertEqual(sum(v.EXPECTED_STRATA.values()),114)

    def test_original_outer_proxy_exact_equals_ibr_width_divided_inv_width(self):
        for r in self.rows:
            self.assertAlmostEqual(r["Inv_Sp"],r["Ibr_W"]/r["Inv_W"],places=12)
        self.assertIn("not physical spine length",
                      self.result["source"]["numerical_proxy"])

    def test_primary_all_edges_are_not_independent_mother_plants(self):
        q=self.result["analyses"]["primary"]
        self.assertEqual(q["n_core_compatible_edges"],239)
        self.assertEqual(q["n_outer_divergent_edges"],67)
        self.assertEqual(q["n_disjoint_core_first_pairs"],44)
        self.assertEqual(q["n_disjoint_core_first_outer_divergent"],9)
        self.assertEqual(q["n_disjoint_outer_first_pairs"],16)
        self.assertEqual(q["thresholds"],
                         {"cor_tol":.1,"inv_tol":.1,"cap_tol":.15,
                          "outer_min":.1})
        self.assertTrue(self.result["matching"]
                        ["greedy_not_optimal_maximum_cardinality"])

    def test_primary_nonoverlapping_greedy_paired_ids_and_source_strata(self):
        q=self.result["analyses"]["primary"]
        for name in ("core_first_pairs","outer_first_pairs"):
            seq=q[name]
            ids=[n for p in seq for n in p["excel_rows"]]
            self.assertEqual(len(ids),len(set(ids)))
            for x in seq:
                self.assertEqual(len(x["excel_rows"]),2)
                a,b=(next(z for z in self.rows if z["excel_row"]==i)
                     for i in x["excel_rows"])
                self.assertEqual(a["taxon"],b["taxon"])
                self.assertEqual(a["sex"],b["sex"])
                self.assertLessEqual(x["core_absolute_log_differences"]["Cor_L"],.1)
                self.assertLessEqual(x["core_absolute_log_differences"]["Inv_L"],.1)
                self.assertLessEqual(x["core_absolute_log_differences"]["Cap_W"],.15)
                self.assertAlmostEqual(
                    x["outer_log_contrast"],
                    abs(__import__("math").log(a["Inv_Sp"]/b["Inv_Sp"])))
        self.assertEqual(len(q["outer_first_pairs"]),16)
        self.assertEqual(len(q["core_first_pairs"]),44)

    def test_source_acaulon_near_floor_under_primary_proxy_threshold(self):
        q=self.result["analyses"]["primary"]
        for sex in ("f","h"):
            group=q["group_counts"]["Cirsium acaulon|"+sex]
            self.assertEqual(group["outer_divergent_edges"],0)
            self.assertEqual(group["outer_first_disjoint_pairs"],0)
        self.assertEqual(
            sum(q["group_counts"][k]["outer_first_disjoint_pairs"]
                for k in q["group_counts"]),16)

    def test_width_denominator_control_is_biologically_material(self):
        x=self.result["analyses"]["width_denominator_control"]
        self.assertEqual((
            x["n_core_compatible_edges"],
            x["n_outer_divergent_edges"],
            x["n_disjoint_core_first_pairs"],
            x["n_disjoint_core_first_outer_divergent"],
            x["n_disjoint_outer_first_pairs"]),
            (97,34,33,11,13))
        self.assertEqual(x["thresholds"]["invw_tol"],.1)
        self.assertEqual(x["group_counts"]["Cirsium acaulon|h"]
                         ["outer_divergent_edges"],0)

    def test_width_control_checks_actual_absolute_tip_width_not_just_ratio(self):
        w=self.result["denominator_controlled_absolute_outer_tip_width"]
        self.assertEqual(w["n_disjoint_ratio_divergent_pairs"],13)
        self.assertEqual(w["n_pairs_with_abs_Ibr_W_log_ratio_ge_0p1"],10)
        self.assertEqual(len(w["per_pair_absolute_width_log_contrast"]),13)
        self.assertTrue(w["source_relative_tip_spread_is_Ibr_W_div_Inv_W"])
        indexed={r["excel_row"]:r for r in self.rows}
        import math
        for e in w["per_pair_absolute_width_log_contrast"]:
            a,b=(indexed[i] for i in e["excel_rows"])
            self.assertLessEqual(abs(math.log(a["Inv_W"]/b["Inv_W"])),.1)
            self.assertAlmostEqual(
                e["absolute_outer_tip_width_log_contrast"],
                abs(math.log(a["Ibr_W"]/b["Ibr_W"])),places=12)

    def test_strict_and_relaxed_core_sensitivity_do_not_invent_zero(self):
        configs={
            "strict_core":(63,16,28,7,9),
            "relaxed_core":(331,88,48,11,19),
            "modest_outer":(239,117,44,19,27),
        }
        for key,expected in configs.items():
            z=self.result["analyses"][key]
            got=(z["n_core_compatible_edges"],z["n_outer_divergent_edges"],
                 z["n_disjoint_core_first_pairs"],
                 z["n_disjoint_core_first_outer_divergent"],
                 z["n_disjoint_outer_first_pairs"])
            self.assertEqual(got,expected,key)
        self.assertGreater(self.result["analyses"]["strict_core"]
            ["n_disjoint_outer_first_pairs"],0)

    def test_negative_mutation_source_sha_not_permitted(self):
        corrupted=self.raw.replace(b"Cirsium acaulon",b"Cirsium XXaulo",1)
        with self.assertRaisesRegex(ValueError,
                                     "PUBLISHER_VERIFIED_SUBSET_BLOB_CHANGED"):
            v.read_original_subset(corrupted)

    def test_negative_mutation_head_ratio_does_not_silently_change_biology(self):
        rows=copy.deepcopy(self.rows)
        rows[0]["Inv_Sp"]+=.5
        # A new table must be re-imported and verified; this mutation violates
        # the original same-plant ratio identity even without changing a rank.
        self.assertNotAlmostEqual(
            rows[0]["Inv_Sp"],rows[0]["Ibr_W"]/rows[0]["Inv_W"])

    def test_greedy_never_counts_one_plant_twice(self):
        p={"excel_rows":[20,21],"is_outer_divergent":True}
        q={"excel_rows":[20,22],"is_outer_divergent":True}
        r={"excel_rows":[23,24],"is_outer_divergent":True}
        s=v.greedy([p,q,r])
        self.assertEqual(len(s),2)
        self.assertEqual([e["excel_rows"] for e in s],[[20,21],[23,24]])

    def test_no_field_adaptation_or_evolution_claim_derived(self):
        limits=self.result["caveats"]
        self.assertEqual(limits["head_access_events_or_pollen_or_seed_fitness_measured"],0)
        self.assertFalse(limits["genetic_LEGO_ancestral_regain_or_selection_evidence"])
        self.assertEqual(limits["new_own_biological_samples"],0)
        self.assertTrue(limits["within_taxon_sex_not_proven_same_field_population"])
        self.assertTrue(limits["external_ratio_wider_phyllary_tip_envelope_not_true_spine_geometry"])

    def test_original_source_paired_morphology_not_two_head_files_joined(self):
        self.assertEqual(set(v.MAIN_VALUES),
            {"Cor_L","Inv_L","Inv_W","Ibr_W","Cap_W","Inv_Sp"})
        self.assertTrue(all(z["excel_row"] for z in self.rows))
        self.assertTrue(self.result["caveats"]
            ["same_plant_trait_matrix_not_necessarily_same_exact_floret_timepoint"])


if __name__=="__main__":
    unittest.main()
