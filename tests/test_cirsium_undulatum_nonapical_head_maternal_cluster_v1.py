"""Tests distinguish original lateral HEAD records from actual mother replicates."""
import copy
import math
import pathlib
import sys
import unittest

ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_undulatum_nonapical_head_maternal_cluster_v1 as v

DATA=ROOT/"data/external/usda_wavyleaf_2007_2008"


class TestActualLateralHeadAccessAndSeedStage(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m=v.read_locked(DATA/"mother_plants.csv","mother_plants.csv",
                            ("plt","year","dam","ins","trt","o.hds","plt.gsd"))
        cls.h=v.read_locked(DATA/"flower_heads.csv","flower_heads.csv",
                            ("plt","year","dam","ins","hd","av.dam","hd.gsd"))

    def test_full_source_orphan_and_conflicting_head_label_rows(self):
        ix,clean,q=v.source_audit(self.m,self.h)
        self.assertEqual(len(self.m),106)
        self.assertEqual(len(self.h),323)
        self.assertEqual(len(clean),297)
        self.assertEqual(q["n_head_with_mother"],309)
        self.assertEqual(q["n_head_without_mother"],14)
        self.assertEqual(q["n_head_treatment_label_mismatch"],12)
        self.assertEqual(q["head_conflicting_source_mother_ids"],
                         ["2007|323","2007|375"])
        self.assertEqual(q["n_matched_heads_with_nonmissing_damage_score"],290)
        self.assertEqual(len(ix),106)

    def test_lateral_observations_n_are_not_mother_n(self):
        _,clean,_=v.source_audit(self.m,self.h)
        groups=v.year_arm_records(v.real_mother_nonapical(self.m,clean))
        self.assertEqual(
            [len(groups["2008"][a]) for a in v.ARMS],[11,12,17,15])
        self.assertEqual(sum(len(groups[y][a]) for y in v.YEARS for a in v.ARMS),
                         87)
        records=v.head_arm_summary(clean,self.m)
        self.assertEqual(sum(records[y][a]["n_lateral_head_records"]
                             for y in v.YEARS for a in v.ARMS),194)

    def test_2008_nonapical_damage_and_seed_positive_are_original_measured(self):
        z=v.analyze(self.m,self.h,bootstrap_reps=200)
        h=z["head_record_group_by_year_treatment"]["2008"]
        self.assertEqual([h[a]["n_lateral_head_records"] for a in v.ARMS],
                         [22,23,28,38])
        self.assertEqual([h[a]["n_highest_recorded_damage_score"] for a in v.ARMS],
                         [18,17,6,4])
        self.assertEqual([h[a]["n_any_good_seed_heads"] for a in v.ARMS],
                         [3,7,22,33])
        self.assertTrue(h["ndni"]["high_score_is_max_observed_0p88_not_exact_fraction_destroyed"])

    def test_parent_weighted_lateral_damage_not_head_weighted(self):
        _,clean,_=v.source_audit(self.m,self.h)
        groups=v.year_arm_records(v.real_mother_nonapical(self.m,clean))
        a=v.group_mother_stats(groups["2008"]["ndni"])
        b=v.group_mother_stats(groups["2008"]["ndi"])
        self.assertAlmostEqual(a["mean_scored_damage_per_mother"],
                               .8672159090909092)
        self.assertAlmostEqual(b["mean_scored_damage_per_mother"],
                               .5004656862745098)
        self.assertAlmostEqual(b["mean_scored_damage_per_mother"]-
                               a["mean_scored_damage_per_mother"],
                               -.36675022281639935)
        self.assertEqual(a["n_mothers"],11)
        self.assertEqual(b["n_mothers"],17)

    def test_2008_parent_level_lateral_success_exact_test(self):
        z=v.analyze(self.m,self.h,bootstrap_reps=100)
        p=z["exact_maternal_any_lateral_seed_comparisons"]["2008"]
        self.assertEqual(p["intact"]["n_any_lateral_seeds_reduced_ambient"],
                         [14,2])
        self.assertEqual(p["intact"]["n_mothers_reduced_ambient"],[21,15])
        self.assertAlmostEqual(p["intact"]["fisher_p"],
                               0.0021770639606882775)
        self.assertAlmostEqual(p["intact"]["BH_q_four_year_damage_strata"],
                               0.00870825584275311)
        self.assertAlmostEqual(p["damaged"]["fisher_p"],
                               .25344980115401355)
        self.assertGreater(p["damaged"]["BH_q_four_year_damage_strata"],.05)

    def test_2008_damage_and_seed_parent_cluster_bootstrap_has_finite_intervals(self):
        _,clean,_=v.source_audit(self.m,self.h)
        groups=v.year_arm_records(v.real_mother_nonapical(self.m,clean))
        boot=v.cluster_bootstrap(groups,reps=300,seed=v.BOOT_SEED)
        self.assertAlmostEqual(
            boot["2008"]["damaged"]["mean_damage"]["estimated_reduced_minus_ambient"],
            -.3185)
        self.assertGreater(
            boot["2008"]["damaged"]["fraction_positive"]["estimated_reduced_minus_ambient"],
            .57)
        for y in v.YEARS:
            for d in ("intact","damaged"):
                for k in v.METRICS:
                    lower,upper=boot[y][d][k]["maternal_percentile_CI95"]
                    self.assertTrue(math.isfinite(lower) and math.isfinite(upper))
                    self.assertLessEqual(lower,upper)

    def test_within_mother_apical_lateral_pair_is_paired_not_treated_as_independent(self):
        z=v.analyze(self.m,self.h,bootstrap_reps=100)
        q=z["within_mother_apical_vs_lateral_damage"]["2008"]["di"]
        self.assertEqual(q["n_mothers_with_both_position_damage_scores"],15)
        self.assertEqual([q["n_lateral_less_damaged"],
                          q["n_lateral_more_damaged"],q["n_equal"]],[12,2,1])
        self.assertLess(q["nonapical_minus_apical_damage_mother_mean"],0)
        self.assertEqual(z["source_audit"]["n_apical_lateral_damage_pairs"],82)

    def test_missing_av_damage_remains_missing(self):
        self.assertIsNone(v.value({"av.dam":"NA"},"av.dam"))
        self.assertIsNone(v.value({"av.dam":""},"av.dam"))
        self.assertNotEqual(v.value({"av.dam":"0"},"av.dam"),None)

    def test_source_conflict_is_not_reassigned_to_whichever_treatment_seems_better(self):
        forged=copy.deepcopy(self.h)
        for r in forged:
            if r["plt"]=="323" and r["year"]=="2007":
                r["dam"]="nd" if r["dam"]=="d" else "d"
        with self.assertRaisesRegex(ValueError,"ORIGINAL_MOTHER_HEAD_SOURCE_JOIN_CHANGED"):
            v.source_audit(self.m,forged)

    def test_plant_without_lateral_not_imputed_an_undamaged_lateral(self):
        _,clean,_=v.source_audit(self.m,self.h)
        observed=v.real_mother_nonapical(self.m,clean)
        self.assertEqual(len(observed),87)
        self.assertEqual(len(self.m)-len(observed),19)
        self.assertTrue(all(r["n_lateral"]>0 for r in observed))

    def test_no_mechanical_structure_or_evolutionary_reexpression_claim(self):
        z=v.analyze(self.m,self.h,bootstrap_reps=100)
        self.assertTrue(z["inference_limits"]["head_rank_not_orientation_or_phyllary_spine_geometry"])
        self.assertTrue(z["inference_limits"]["head_damage_score_is_measured_after_access_treatment_and_not_causal_independent_covariate"])
        self.assertTrue(z["inference_limits"]["individual_source_head_rows_are_not_independent_mother_replicates"])
        self.assertEqual(z["inference_limits"]["new_biological_field_samples_added"],0)


if __name__=="__main__":
    unittest.main()
