"""Synthetic unit tests for source-verified Cirsium arvense visitor clock audit."""
import csv
import io
import sys
import unittest
from datetime import date
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"analysis"))
from audit_arvense_original_visitor_clock_v1 import (
    git_blob_sha1, parse_source_csv, parse_bout_time, audit
)


BOUT_HEADER = ("Start","End","Site_Code","Block_Code","Treatment_Code",
               "Number_openFlowerhead","Remarks","Instance_ID","Instance_Name","Key")
VISITOR_HEADER = ("Common_Name","Visited_Flowers","Key","Feeding_Guild",
                  "Pollination","Order","Family","Genus","Species","Species_Code")


def as_bytes(columns,rows):
    output=io.StringIO()
    writer=csv.DictWriter(output,fieldnames=columns,lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode("utf-8")


def bout(key,day,site="campus",block="b1",treatment="Control",note=""):
    return dict(zip(BOUT_HEADER,[
        f"{day} 09:00",f"{day} 09:15",site,block,treatment,"2",note,
        key,"7/2/2022",key,
    ]))


def insect(key,order="Hymenoptera"):
    return dict(zip(VISITOR_HEADER,[
        "visitor","1",key,"Herbivore","Yes",order,"Family","","","code"
    ]))


class TestArvenseClock(unittest.TestCase):
    def test_quoted_comments_roundtrip_and_source_hash(self):
        source=as_bytes(BOUT_HEADER,[bout("k1","2/7/2022",note='With bees, flies and "visitors"')])
        self.assertEqual(len(parse_source_csv(source, BOUT_HEADER)),1)
        self.assertEqual(parse_source_csv(source,BOUT_HEADER)[0]["Remarks"],
                         'With bees, flies and "visitors"')
        self.assertEqual(len(git_blob_sha1(source)),40)
        self.assertEqual(parse_bout_time("13/07/22 09:11").date(),date(2022,7,13))
        self.assertEqual(parse_bout_time("2/7/2022 9:11").date(),date(2022,7,2))
        with self.assertRaises(ValueError):
            parse_bout_time("2022-07-02")

    def test_source_drift_fails_closed(self):
        raw=as_bytes(BOUT_HEADER,[bout("k1","2/7/2022")])
        with self.assertRaisesRegex(ValueError,"SOURCE_SHA"):
            parse_source_csv(raw,BOUT_HEADER,expected_git_sha1="0"*40)
        with self.assertRaisesRegex(ValueError,"SOURCE_HEADER"):
            parse_source_csv(raw,VISITOR_HEADER)

    def test_ambiguous_bout_join_excluded_not_arbitrarily_matched(self):
        bouts=[bout("a","2/7/2022"),bout("b","16/07/22"),
               bout("dup","16/07/22"),bout("dup","16/07/22")]
        visitors=[insect("a","Hymenoptera"),insect("a","Diptera"),
                  insect("b","Diptera"),insect("dup","Hymenoptera"),
                  insect("not-a-key","Hymenoptera"),insect("b","")]
        result=audit(bouts,visitors,require_original_counts=False)
        q=result["matched_taxonomic_observations"]
        self.assertEqual(result["source_rows"]["bouts"],4)
        self.assertEqual(result["source_rows"]["visitor_lines"],6)
        self.assertEqual(q["named_order_visitor_lines"],3)
        self.assertEqual(q["named_order_bouts"],2)
        self.assertEqual(q["mixed_Hymenoptera_Diptera_bouts"],1)
        self.assertEqual(q["ambiguous_bout_key_visitor_lines"],1)
        self.assertEqual(q["unmatched_bout_key_visitor_lines"],1)
        self.assertEqual(q["blank_order_visitor_lines"],1)
        self.assertEqual(result["source_rows"]["distinct_calendar_dates"],2)

    def test_bout_equal_weight_differs_from_insect_record_weight(self):
        bouts=[bout("a","2/7/2022"),bout("b","2/7/2022"),
               bout("c","16/07/22"),bout("d","16/07/22")]
        # Early: giant H bout and one small D; later two moderate mixed.
        visitors=[insect("a","Hymenoptera") for _ in range(9)]
        visitors += [insect("b","Diptera")]
        visitors += [insect("c","Hymenoptera"),insect("c","Diptera"),
                     insect("d","Hymenoptera"),insect("d","Diptera")]
        result=audit(bouts,visitors,require_original_counts=False)
        site=result["site_time_comparisons"][0]
        self.assertGreater(site["record_weighted_final_minus_first"],-1)
        self.assertLess(site["record_weighted_final_minus_first"],0)
        self.assertEqual(site["equal_bout_final_minus_first"],0)
        self.assertEqual(result["source_rows"]["distinct_calendar_dates"],2)

    def test_unmeasured_intergenerational_capacity_not_imputed(self):
        result=audit([bout("x","2/7/2022")],[insect("x")],
                     require_original_counts=False)
        self.assertFalse(result["scientific_boundaries"]["Markov_environment_transition_probability_identified"])
        self.assertFalse(result["scientific_boundaries"]["same_plant_heritable_functional_modules_linked"])
        self.assertFalse(result["scientific_boundaries"]["same_head_reproductive_fitness_from_visits_identified"])
        self.assertFalse(result["scientific_boundaries"]["genealogy_regain_or_fixation_inferred"])

    def test_actual_public_count_requirement_rejects_partial_sources(self):
        with self.assertRaisesRegex(ValueError,"SOURCE_ROW_COUNT"):
            audit([bout("x","2/7/2022")],[insect("x")],require_original_counts=True)

    def test_no_raw_keys_or_locations_in_output(self):
        r=audit([bout("sensitive-individual-key","2/7/2022")],
                [insect("sensitive-individual-key")],require_original_counts=False)
        text=repr(r)
        self.assertNotIn("sensitive-individual-key",text)
        self.assertNotIn("Instance_ID",text)
        self.assertNotIn("coordinates",text)


if __name__=="__main__":
    unittest.main()
