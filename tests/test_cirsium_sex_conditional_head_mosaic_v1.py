"""Falsification, real-source evidence and reproducible sex-specific head morphology."""
import math
from pathlib import Path
import sys
import unittest
import copy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_sex_conditional_head_mosaic_v1 as v

SRC=ROOT/"data/evidence/cirsium_2023_same_plant_head_morphology_subset_v1.csv"


class OriginalCirsiumSexHeadModules(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw=SRC.read_bytes()
        cls.groups,cls.source_sha=v.read_source(cls.raw)

    def test_true_original_individual_sample_grain(self):
        self.assertEqual(cls_len:=sum(len(g["all"]) for g in self.groups),114)
        self.assertEqual(len(self.groups),3)
        self.assertEqual([(g["n_f"],g["n_h"]) for g in self.groups],
                         [(10,30),(10,24),(10,30)])
        self.assertEqual(self.source_sha,v.EXPECTED_BLOB)
        self.assertEqual(v.AUTHOR_ORIGINAL_XLSX_SHA256,
                         "c06423cc6f9800ba4775a4428afa072ba1b3dc6fb48ccd05fdb7e42a2f03ee79")

    def test_source_external_proxy_is_ratio_not_anatomical_spine_length(self):
        self.assertAlmostEqual(
            v.source_features({
                "Cor_L":"30","Inv_L":"29","Inv_W":"20",
                "Ibr_W":"30","Cap_W":"40","Inv_Sp":"1.5"
            })["outer_ratio_log"],math.log(1.5))
        with self.assertRaisesRegex(ValueError,"ORIGINAL_OUTER_TIP_RATIO_INVALID"):
            v.source_features({
                "Cor_L":"30","Inv_L":"29","Inv_W":"20",
                "Ibr_W":"30","Cap_W":"40","Inv_Sp":"2"
            })

    def test_original_group_original_ids_hard_locked(self):
        bad=self.raw.replace(b"Cirsium erisithales",b"Cirsium mutisithal",1)
        with self.assertRaisesRegex(ValueError,"SOURCE_ORIGINAL_114_GIT_BLOB_DRIFT"):
            v.read_source(bad)

    def test_negative_nonpositive_measurement_is_not_log_transformed(self):
        x={"Cor_L":"0","Inv_L":"29","Inv_W":"20",
           "Ibr_W":"30","Cap_W":"40","Inv_Sp":"1.5"}
        with self.assertRaisesRegex(ValueError,"SOURCE_INVALID_HEAD_LENGTH_Cor_L"):
            v.source_features(x)

    def test_measured_sex_corolla_shortening_three_original_taxa(self):
        exp={
            "Cirsium acaulon":-.12468196150717636,
            "Cirsium bertolonii - Monte Sagro":-.06235879411812961,
            "Cirsium erisithales":-.12422772117668801,
        }
        for group in self.groups:
            got=v.contrast(group["female"],group["herm"])
            self.assertLess(got["corolla_log"],0)
            self.assertAlmostEqual(got["corolla_log"],exp[group["taxon"]],
                                   places=12)

    def test_outer_ratio_is_not_same_sex_contrast_in_each_species(self):
        exp={
            "Cirsium acaulon":.0045185882283212166,
            "Cirsium bertolonii - Monte Sagro":-.06291449552781503,
            "Cirsium erisithales":-.020073736600489583,
        }
        for group in self.groups:
            got=v.contrast(group["female"],group["herm"])
            self.assertAlmostEqual(got["outer_ratio_log"],exp[group["taxon"]],
                                   places=12)
        self.assertGreater(exp["Cirsium acaulon"],0)
        self.assertLess(exp["Cirsium bertolonii - Monte Sagro"],0)

    def test_rng_fixed_across_repeat_calls_and_no_zero_seed(self):
        a=v.xorshift32(20261024)
        b=v.xorshift32(20261024)
        self.assertEqual([next(a) for i in range(12)],
                         [next(b) for i in range(12)])
        with self.assertRaisesRegex(ValueError,"INVALID_OR_ZERO_DETERMINISTIC_SEED"):
            next(v.xorshift32(0))

    def test_pooled_source_contrast_not_genus_phylogenetic_sample(self):
        observed={g["taxon"]:v.contrast(g["female"],g["herm"])
                  for g in self.groups}
        pooled={k:v.mean([observed[t][k] for t in v.TAXA])
                for k in v.FEATURES}
        self.assertAlmostEqual(pooled["corolla_log"],-.10375615893399799)
        self.assertAlmostEqual(pooled["outer_ratio_log"],-.026156547966661132)
        self.assertAlmostEqual(100*math.expm1(pooled["corolla_log"]),
                               -9.85549200466496)

    def test_relative_corolla_length_has_all_three_negative_signs(self):
        for group in self.groups:
            traits=v.contrast(group["female"],group["herm"])
            self.assertLess(traits["relative_floret_log"],0)

    def test_sex_label_permutation_is_not_reproductive_role_causal_trial(self):
        self.assertEqual(v.PRIMARY,["corolla_log","outer_ratio_log"])
        self.assertEqual(v.PERMUTATIONS,20000)
        self.assertEqual(v.BOOTSTRAPS,5000)
        self.assertEqual(v.FEATURES[0],"corolla_log")

    def test_taxon_specific_study_does_not_give_measured_insect_fitness(self):
        src=(ROOT/"data/intake/capitulum_guild_access_event_ledger_v1.csv")
        self.assertEqual(len(src.read_text(encoding="utf-8").splitlines()),1)
        for g in self.groups:
            self.assertTrue(all(x["values"]["corolla_log"]>0 for x in g["all"]))

    def test_positive_study_conclusion_falsifiable_by_permutation(self):
        d={g["taxon"]:v.contrast(g["female"],g["herm"]) for g in self.groups}
        self.assertEqual(sum(d[t]["corolla_log"]<0 for t in v.TAXA),3)
        self.assertLess(abs(d["Cirsium acaulon"]["outer_ratio_log"]),
                        abs(d["Cirsium bertolonii - Monte Sagro"]["outer_ratio_log"]))

    def test_source_is_not_new_experiment_or_genetically_assigned(self):
        self.assertEqual(len(v.TAXA),3)
        self.assertEqual(sum(v.F_COUNTS),30)
        self.assertEqual(sum(v.H_COUNTS),84)
        self.assertTrue(all(isinstance(g["n_f"],int) for g in self.groups))

if __name__=="__main__":
    unittest.main()
