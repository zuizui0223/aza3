"""Negative-control sensitivity for five published source protein candidate IDs."""
import json
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from run_anto_cirsium_reference_jackknife_v4 import main

class ReferenceJackknifeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=main()

    def test_all_five_are_same_initial_strong_positive_under_full_panel(self):
        self.assertEqual(self.r["number_of_screened_public_transcripts"],5)
        self.assertEqual(len({row["source_assembly_group"] for row in self.r["records"]}),4)
        self.assertTrue(all(row["comparisons"]["all_seven"]["strong_and_same_family"]
                            for row in self.r["records"]))

    def test_only_two_of_five_stay_stable_under_all_reference_omissions(self):
        self.assertTrue(self.r["frozen_published_source_reproducibility_verified"])
        self.assertEqual(self.r["stable_family_and_strong_in_all_scenarios"],2)
        self.assertEqual(self.r["strong_when_Saussurea_removed"],3)
        ids={r["transcript_id"] for r in self.r["records"]
             if r["stable_across_all_reference_omissions"]}
        self.assertEqual(ids,{
            "TR35526|c0_g1_i1","TR39766|c0_g1_i2"})
        byid={r["transcript_id"]:r for r in self.r["records"]}
        for id in ("TR38875|c0_g4_i5","TR38875|c0_g4_i7"):
            self.assertEqual(byid[id]["comparisons"]["leave_out_Saussurea_medusa"]["best_group"],
                             "ANR_negative_control")
        self.assertEqual(byid["TR58302|c0_g1_i1"]["comparisons"]["leave_out_Arabidopsis_positive_only"]["best_group"],
                         "FLS_negative_control")

    def test_no_reference_scenario_removes_all_same_family_proteins(self):
        for item in self.r["records"]:
            self.assertEqual(set(item["comparisons"]),{
                "all_seven","leave_out_Saussurea_medusa",
                "leave_out_Gerbera_hybrida",
                "leave_out_Arabidopsis_positive_only",
            })
            self.assertTrue(all(z["reference_count"]>=5 for z in item["comparisons"].values()))

    def test_sensitivity_outcome_does_not_certify_orthology_or_colour_regain(self):
        self.assertEqual(self.r["claim_limits"]["orthology_certified"],False)
        self.assertEqual(self.r["claim_limits"]["floral_pigment_reexpression_demonstrated"],False)
        self.assertTrue(all(not row["gene_function_verified"] for row in self.r["records"]))

    def test_jackknife_result_is_self_consistent(self):
        expected=sum(all(c["strong_and_same_family"] for c in row["comparisons"].values())
                     for row in self.r["records"])
        self.assertEqual(self.r["stable_family_and_strong_in_all_scenarios"],expected)
        expected_saus=sum(row["comparisons"]["leave_out_Saussurea_medusa"]["strong_and_same_family"]
                          for row in self.r["records"])
        self.assertEqual(self.r["strong_when_Saussurea_removed"],expected_saus)

    def test_every_reference_group_winner_is_annotated_not_functional(self):
        for row in self.r["records"]:
            for z in row["comparisons"].values():
                self.assertIn(z["best_group"],{
                    "DFR_reference","ANS_reference",
                    "ANR_negative_control","FLS_negative_control"
                })
                self.assertIsNotNone(z["best_reference_id"])

if __name__=="__main__":
    unittest.main()
