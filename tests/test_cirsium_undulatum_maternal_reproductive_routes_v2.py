"""Independent tests for REAL USDA Cirsium undulatum maternal reproductive routes.

One maternal plant, NOT one flowerhead, is the unit. Apical damage differs
from apical functional abolition; 2007 maternal 327 has mismatching totals.
Source: USDA DOI 10.15482/USDA.ADC/1522649, West & Louda 2021
DOI 10.1002/ajb2.1744 (the experiment and broad compensation effect are prior art).
"""
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from audit_cirsium_undulatum_maternal_reproductive_routes_v2 import (
    classify_route, summarize_year_arm, maternal_routes,
    fixed_fisher_eight_contrasts
)
from audit_cirsium_undulatum_maternal_hurdle_size_v1 import read_locked

BASE=ROOT/"data/external/usda_wavyleaf_2007_2008"


def synthetic(apical=0,lateral=7,total=7,year="2008",
              treatment="dni",plant="T1"):
    return {"plt":plant,"year":year,"trt":treatment,
            "dam":"d" if treatment.startswith("d") else "nd",
            "ins":"i" if treatment.endswith("i") and not treatment.endswith("ni") else "ni",
            "plt.gsd":str(total),"ap.hd":str(apical),"o.hds":str(lateral)}


class MotherRouteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows=read_locked(BASE/"mother_plants.csv","mother_plants.csv",
            ("plt","year","trt","dam","ins","plt.gsd","ap.hd","o.hds"))
        cls.r=maternal_routes(cls.rows,strict=True)

    def test_source_is_real_106_mothers_and_exactly_one_mismatched_apical_seed_sum(self):
        x=self.r["source"]
        self.assertEqual(x["original_mother_rows"],106)
        self.assertEqual(x["coherent_mother_rows"],105)
        self.assertEqual(x["inconsistent_mother_identifiers"],["2007|327"])
        self.assertEqual(x["original_heads_in_related_USDA_source"],323)
        self.assertTrue(x["missing_fraction_rejected_not_imputed"])

    def test_all_four_anatomical_outcome_states_exist_in_2008_intact_ambient(self):
        y=self.r["maternal_route_states"]["2008"]
        self.assertEqual(y["ndni"]["n"],15)
        self.assertEqual(y["ndni"]["states"],{
            "no_good_seed":12,"apical_only":1,"lateral_only":1,"both_sites":1
        })
        self.assertEqual(y["dni"]["states"],{
            "no_good_seed":5,"apical_only":0,"lateral_only":7,"both_sites":0
        })
        self.assertEqual(y["ndi"]["states"],{
            "no_good_seed":0,"apical_only":7,"lateral_only":1,"both_sites":13
        })
        self.assertEqual(y["di"]["states"],{
            "no_good_seed":0,"apical_only":4,"lateral_only":8,"both_sites":7
        })

    def test_damage_does_not_mean_apical_production_always_disappears(self):
        y=self.r["maternal_route_states"]["2008"]
        self.assertEqual(y["di"]["apical_positive_mothers"],11)
        self.assertEqual(y["dni"]["apical_positive_mothers"],0)
        self.assertEqual(y["di"]["lateral_positive_mothers"],15)
        self.assertEqual(y["ndi"]["lateral_positive_mothers"],14)
        self.assertFalse(self.r["scientific_limits"]["apical_damage_means_apical_function_completely_removed"])

    def test_worst_enemy_conditions_have_lateral_only_or_zero(self):
        y=self.r["maternal_route_states"]["2008"]
        self.assertEqual(y["dni"]["states"]["lateral_only"],7)
        self.assertEqual(y["dni"]["states"]["no_good_seed"],5)
        self.assertEqual(y["ndni"]["states"]["no_good_seed"],12)
        self.assertAlmostEqual(y["dni"]["lateral_success_fraction"],7/12)
        self.assertAlmostEqual(y["ndni"]["lateral_success_fraction"],2/15)
        self.assertAlmostEqual(y["dni"]["lateral_mean_good_seed_count"],107/12)
        self.assertAlmostEqual(y["ndni"]["lateral_mean_good_seed_count"],21/15)

    def test_year2007_bad_record_does_not_force_false_lateral_count(self):
        y=self.r["maternal_route_states"]["2007"]
        self.assertEqual(y["dni"]["n"],8)
        self.assertEqual(y["dni"]["excluded_incoherent_mothers"],1)
        self.assertEqual(y["dni"]["states"]["lateral_only"],6)
        self.assertEqual(y["dni"]["states"]["no_good_seed"],2)
        self.assertEqual(y["ndni"]["states"],{
            "no_good_seed":0,"apical_only":2,"lateral_only":2,"both_sites":5
        })
        self.assertEqual(y["ndi"]["states"],{
            "no_good_seed":0,"apical_only":3,"lateral_only":2,"both_sites":5
        })
        self.assertEqual(y["di"]["states"],{
            "no_good_seed":0,"apical_only":1,"lateral_only":6,"both_sites":4
        })

    def test_eight_fisher_tests_and_BH_adjustment_do_not_inflate_one_pvalue(self):
        t=self.r["lateral_presence_fisher_tests"]
        self.assertEqual(len(t),8)
        self.assertAlmostEqual(
            t["2008_damage_ambient"]["p_exact"],0.03691838291380626,places=12)
        self.assertAlmostEqual(
            t["2008_damage_ambient"]["q_BH_eight"],0.14767353165522504,places=12)
        self.assertAlmostEqual(
            t["2008_insecticide_intact"]["p_exact"],0.0021770639606882775,places=12)
        self.assertAlmostEqual(
            t["2008_insecticide_intact"]["q_BH_eight"],0.01741651168550622,places=12)
        self.assertAlmostEqual(t["2007_damage_ambient"]["p_exact"],1,places=12)
        self.assertFalse(self.r["scientific_limits"]["new_treatment_effect_is_randomized_between_plant_genotypes_proven"])

    def test_parental_source_hurdle_consistent_with_reproductive_states(self):
        for year in ("2007","2008"):
            for name,state in self.r["maternal_route_states"][year].items():
                self.assertEqual(sum(state["states"].values()),state["n"])
                self.assertEqual(state["apical_positive_mothers"],
                                 state["states"]["apical_only"]+state["states"]["both_sites"])
                self.assertEqual(state["lateral_positive_mothers"],
                                 state["states"]["lateral_only"]+state["states"]["both_sites"])
                self.assertEqual(state["seed_positive_mothers"],
                                 state["n"]-state["states"]["no_good_seed"])

    def test_synthetic_apical_only_vs_lateral_only_not_pseudoregain(self):
        self.assertEqual(classify_route(synthetic(apical=5,lateral=0,total=5)),"apical_only")
        self.assertEqual(classify_route(synthetic(apical=0,lateral=5,total=5)),"lateral_only")
        self.assertEqual(classify_route(synthetic(apical=5,lateral=6,total=11)),"both_sites")
        self.assertEqual(classify_route(synthetic(apical=0,lateral=0,total=0)),"no_good_seed")
        self.assertIsNone(classify_route(synthetic(apical=0,lateral=1,total=20)))
        with self.assertRaises(ValueError):
            classify_route(synthetic(apical=-1,lateral=2,total=1))

    def test_no_genetic_memory_or_same_head_fitness_claims(self):
        for flag in self.r["scientific_limits"].values():
            self.assertFalse(flag)
        self.assertEqual(self.r["inference"]["observed_entity"],"maternal_good_seed_anatomical_route")
        self.assertEqual(self.r["inference"]["biological_conclusion"],
                         "OBSERVED_REAL_SPATIAL_REPRODUCTIVE_FALLBACK_NOT_ADAPTIVE_SELECTION_OR_GENETIC_REGAIN")


if __name__=="__main__":
    unittest.main()
