"""RED/GREEN source-join audits: visitor bouts vs inferred seed parent unit.

Do not test or assert same-head effects or causal pollinator-mediated fitness.
"""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"analysis"))
from audit_slk_anto_arvense_bout_seed_join_v1 import (
    normalize_arm, join_source_labels, pearson, controlled_correlation,
    summarize_source_biological_join,
)


def head(site="campus",block="B1",arm="Control",viable=8,
         damaged=1,nonviable=1):
    return {"site":site,"block":block,"treatment":arm,"vial":"1",
            "Viable_Seeds":viable,"NViable_Seeds":nonviable,
            "Damaged_Seeds":damaged,"Missing_Seeds":0,
            "Flowers":viable+damaged+nonviable,"Insects":0}


def bout(key,site="campus",block="b1",arm="Control",day="2/7/2022",
         start="9:00",end="9:15",openheads="2"):
    return {"Start":day+" "+start,"End":day+" "+end,
            "Site_Code":site,"Block_Code":block,"Treatment_Code":arm,
            "Number_openFlowerhead":openheads,"Key":key}


def visitor(key,order="Hymenoptera",n="1"):
    return {"Key":key,"Order":order,"Visited_Flowers":n}


class SourceJoinTests(unittest.TestCase):
    def test_treatment_alias_and_known_arm_scope(self):
        self.assertEqual(normalize_arm("Exclude predators"),"predator exclusion")
        self.assertEqual(normalize_arm("Predator Exclusion"),"predator exclusion")
        self.assertEqual(normalize_arm("Predator model"),"predator model")
        self.assertEqual(normalize_arm("Control"),"control")
        self.assertNotEqual(normalize_arm("Pollinator Exclusion"),"predator exclusion")

    def test_genuine_label_group_is_not_one_head(self):
        heads=[head(viable=8),head(viable=4,damaged=2,nonviable=4)]
        bouts=[bout("v1")]
        visits=[visitor("v1","Hymenoptera","3"),
                visitor("v1","Diptera","1")]
        j=join_source_labels(heads,bouts,visits,strict=False)
        self.assertEqual(j["source_join"]["matched_label_groups_before_effort_filter"],1)
        self.assertEqual(j["source_join"]["matched_label_groups_after_effort_filter"],1)
        self.assertEqual(j["source_join"]["valid_nonzero_head_bouts"],1)
        self.assertEqual(j["source_join"]["excluded_ambiguous_key_bouts"],0)
        self.assertEqual(j["group_records"][0]["n_heads"],2)
        self.assertAlmostEqual(j["group_records"][0]["viable_fraction"],.6)
        self.assertAlmostEqual(j["group_records"][0]["hym_visit_share"],.75)
        self.assertAlmostEqual(j["group_records"][0]["rate_visit_per_nominal_head_hour"],8)
        self.assertEqual(j["group_records"][0]["n_bouts"],1)
        self.assertTrue(j["source_join"]["group_label_join_not_verified_same_head"])

    def test_ambiguous_bout_key_does_not_turn_into_two_observations(self):
        h=[head()]
        bouts=[bout("duplicate"),bout("duplicate",start="9:30",end="9:45"),
               bout("good",start="10:00",end="10:15")]
        visits=[visitor("duplicate","Hymenoptera","8"),
                visitor("good","Diptera","1")]
        j=join_source_labels(h,bouts,visits,strict=False)
        self.assertEqual(j["source_join"]["excluded_ambiguous_key_bouts"],2)
        self.assertEqual(j["source_join"]["valid_nonzero_head_bouts"],1)
        self.assertEqual(j["group_records"][0]["hym_count"],0)
        self.assertEqual(j["group_records"][0]["dip_count"],1)

    def test_zero_open_head_or_zero_duration_out_of_rate_denominator(self):
        b=[bout("a",openheads="0"),bout("b",start="10:00",end="10:00"),
           bout("c",start="10:30",end="10:45",openheads="2")]
        j=join_source_labels([head()],b,[visitor("a"),visitor("b"),
                                                 visitor("c")],strict=False)
        self.assertEqual(j["source_join"]["excluded_zero_open_head_bouts"],1)
        self.assertEqual(j["source_join"]["excluded_zero_duration_bouts"],1)
        self.assertEqual(j["source_join"]["valid_nonzero_head_bouts"],1)

    def test_site_and_treatment_demeaning_matches_hand_synthetic(self):
        data=[{"site":"a","treatment":"control","x":0,"y":0},
              {"site":"a","treatment":"control","x":2,"y":4},
              {"site":"b","treatment":"model","x":10,"y":5},
              {"site":"b","treatment":"model","x":12,"y":9}]
        p=controlled_correlation(data,"x","y")
        self.assertEqual(p["n"],4)
        self.assertEqual(p["n_strata"],2)
        self.assertAlmostEqual(p["r"],1.0)
        self.assertIsNone(pearson([1,1,1],[0,1,2]))

    def test_no_species_inference_from_pollinator_guild_label(self):
        r=summarize_source_biological_join(
           [head(),head(block="B2")],
           [bout("p1"),bout("p2",block="b2")],
           [visitor("p1"),visitor("p2","Diptera")],strict=False)
        gates=r["scientific_boundaries"]
        self.assertFalse(gates["experimentally_verified_pollination_effectiveness"])
        self.assertFalse(gates["same_biological_head_contact_and_seed_fate"])
        self.assertFalse(gates["anatomical_SLK_division_of_labour_profit_estimated"])
        self.assertFalse(gates["functional_ANTO_reexpression_competence_measured"])
        self.assertFalse(gates["intergenerational_environment_markov_rate_inferred"])
        self.assertEqual(r["interpretation"]["analysis_status"],
                         "EXPLORATORY_POST_HOC_SOURCE_LABEL_ASSOCIATION")

    def test_observed_zero_exposure_is_not_missing_head_or_genetic_loss(self):
        j=join_source_labels([head()],
                    [bout("v",openheads="2")],
                    [visitor("v",order="",n="")],strict=False)
        self.assertEqual(j["group_records"][0]["hym_count"],0)
        self.assertIsNone(j["group_records"][0]["hym_visit_share"])
        self.assertEqual(j["group_records"][0]["n_heads"],1)

    def test_invalid_seed_component_and_unmatched_code_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"SEED_MASS_BALANCE"):
            join_source_labels([dict(head(), Flowers=500)],
                               [bout("v")],[visitor("v")],strict=False)
        with self.assertRaisesRegex(ValueError,"VISITOR_NEGATIVE_COUNT"):
            join_source_labels([head()],[bout("v")],
                               [visitor("v",n="-3")],strict=False)

    def test_public_json_never_contains_original_bout_key(self):
        r=summarize_source_biological_join([head()],
                            [bout("source-sensitive-uuid")],
                            [visitor("source-sensitive-uuid")],strict=False)
        self.assertNotIn("source-sensitive-uuid",str(r))
        self.assertNotIn("b1",str(r))
        self.assertEqual(r["results"]["actual_plant_genomic_data_records"],0)


if __name__=="__main__":
    unittest.main()
