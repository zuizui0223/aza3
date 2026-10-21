"""Original-plant analysis and scientific fail-closed tests, not fabricated field observations."""
from __future__ import annotations
import copy
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_cirsium_undulatum_mother_head_compensation_v1 as m
SOURCE=ROOT/"data/external/usda_wavyleaf_2007_2008"


class TestSourceLinkedMaternalCompensation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plants=m.load(SOURCE/"mother_plants.csv","mother_plants.csv",m.MOTHER_HEADERS)
        cls.heads=m.load(SOURCE/"flower_heads.csv","flower_heads.csv",m.HEAD_HEADERS)

    def test_real_original_usda_source_and_groups(self):
        m.validate_mothers(self.plants)
        groups=m.year_groups(self.plants)
        self.assertEqual(m.SOURCE_BLOBS["mother_plants.csv"],
                         "458800272c3d0a03c7653702027a8f8578d9e88c")
        self.assertEqual(len(self.plants),106)
        self.assertEqual(len(self.heads),323)
        self.assertEqual({k:len(v) for k,v in groups["2007"].items()},
                         {"ndni":9,"dni":9,"ndi":10,"di":11})
        self.assertEqual({k:len(v) for k,v in groups["2008"].items()},
                         {"ndni":15,"dni":12,"ndi":21,"di":19})

    def test_join_has_two_conflicting_id_groups_not_pseudoreplicated_heads(self):
        a=m.quality(self.plants,self.heads)
        self.assertEqual(a["n_heads_with_matching_mother"],309)
        self.assertEqual(a["n_head_records_without_mother"],14)
        self.assertEqual(a["n_head_rows_treatment_label_conflicts"],12)
        self.assertEqual(a["treatment_label_conflict_mother_ids"],
                         ["2007|323","2007|375"])
        self.assertEqual(a["component_mismatch_mother_ids"],["2007|327"])
        self.assertEqual(a["mother_ids_with_no_head_records"],
                         ["2007|329","2007|367"])
        self.assertEqual(a["n_matched_parent_head_sum_mismatches"],0)

    def test_real_mean_plant_seed_outputs_by_original_treatment(self):
        gs=m.year_groups(self.plants)
        self.assertAlmostEqual(m.group_stats(gs["2007"]["ndi"])["total"],88.4)
        self.assertAlmostEqual(m.group_stats(gs["2007"]["di"])["total"],60.)
        self.assertAlmostEqual(m.group_stats(gs["2008"]["ndni"])["total"],40/15)
        self.assertAlmostEqual(m.group_stats(gs["2008"]["di"])["total"],2886/19)
        self.assertEqual(m.group_stats(gs["2008"]["ndni"])["zero_seed_mothers"],12)

    def test_lateral_share_can_rise_despite_absolute_lateral_seed_decline(self):
        g=m.year_groups(self.plants)["2007"]
        nd=m.group_stats(g["ndi"])
        d=m.group_stats(g["di"])
        self.assertGreater(d["lateral_seed_share_of_pooled_good_seed_total"],
                           nd["lateral_seed_share_of_pooled_good_seed_total"])
        self.assertLess(d["lateral"],nd["lateral"])
        self.assertGreater(d["flowered"],nd["flowered"])

    def test_real_seed_presence_hurdle_2008(self):
        g=m.year_groups(self.plants)["2008"]
        p={key:m.group_stats(v)["has_seeds"] for key,v in g.items()}
        self.assertEqual(p["di"],p["ndi"])
        self.assertEqual(p["di"],1)
        self.assertAlmostEqual(p["dni"],7/12)
        self.assertAlmostEqual(p["ndni"],3/15)
        self.assertGreater(p["dni"]-p["ndni"],0)

    def test_reconciliation_flags_are_not_mutated_to_force_alignment(self):
        h=copy.deepcopy(self.heads)
        h[0]["plt"]="invented"
        with self.assertRaisesRegex(ValueError,"USDA_ORIGINAL_FLOWERHEAD_N_CHANGED|USDA_EXTERNAL_HEAD_JOIN_OR_ACCOUNTING_DRIFT"):
            m.quality(self.plants,h)

    def test_bad_treatment_code_rejected(self):
        p=copy.deepcopy(self.plants)
        p[0]["trt"]="bad"
        with self.assertRaisesRegex(ValueError,"MOTHER_TREATMENT_CODING_CHANGED"):
            m.validate_mothers(p)

    def test_duplicate_mother_id_rejected(self):
        p=copy.deepcopy(self.plants)
        p[1]["plt"]=p[0]["plt"]
        p[1]["year"]=p[0]["year"]
        with self.assertRaisesRegex(ValueError,"DUPLICATE_MOTHER_YEAR_PLANT_ID"):
            m.validate_mothers(p)

    def test_actual_bootstrap_seed_and_finite_signs(self):
        rng=m.xs32(20261021)
        self.assertEqual([next(rng) for _ in range(10)],[
            553487182,2736026372,2213040766,3176989453,3131645939,
            3014744152,2415432869,820197491,2856518437,2396979446
        ])
        g=m.year_groups(self.plants)
        sample=m.bootstrap_contrasts(g,reps=200,seed=20261021)
        result=sample["2008"]["total"]["insecticide_under_apical_damage"]
        self.assertAlmostEqual(result["estimate"],2886/19-107/12)
        self.assertTrue(all(isinstance(x,float) for x in result["bootstrap_percentile_CI95"]))

    def test_sensitivity_removing_source_conflicts_is_published_explicitly(self):
        z=m.analyze(self.plants,self.heads,bootstrap_reps=100)
        y=z["sensitivity_drop_two_treatment_mismatch_and_one_component_mismatch"]["2007"]
        self.assertEqual([y[a]["n"] for a in m.ARMS],[8,8,9,11])
        self.assertAlmostEqual(y["damage_effects"]["total"]["ambient"],44.125-40.75)
        self.assertTrue(z["claims"]["source_treatment_conflict_ids_are_not_silently_recoded"])
        self.assertFalse(z["claims"]["same_head_pollinator_oviposition_and_seed_fitness_measured"])

    def test_source_plant_is_statistical_unit_not_323_independent_replicates(self):
        z=m.analyze(self.plants,self.heads,bootstrap_reps=100)
        self.assertEqual(z["cohort"]["n_mothers"],106)
        self.assertEqual(z["cohort"]["n_heads"],323)
        self.assertTrue(z["cohort"]["head_rows_are_not_independent_maternal_replicates"])
        self.assertEqual(z["claims"]["external_independent_plant_samples_newly_collected"],0)


if __name__=="__main__":
    unittest.main()
