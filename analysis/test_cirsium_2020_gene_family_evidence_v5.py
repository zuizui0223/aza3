"""Negative controls for source-audited gene-family phylogenetic record."""
from copy import deepcopy
import json
import sys
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from validate_cirsium_2020_gene_family_evidence_v5 import (
    check_source_record,compare_tree_analysis,main
)

class RealSourceFamilyEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=check_source_record()

    def test_five_distinct_source_transcripts_seven_tree_evaluations(self):
        r=self.data["reported_results"]
        self.assertEqual(len(r),7)
        self.assertEqual(len({x["source_transcript"] for x in r}),5)
        self.assertEqual({x["case"] for x in r},
                         {"DFR_ANR_primary","DFR_ANR_alt_isoform","ANS_FLS"})

    def test_no_biological_orthology_or_ryukyu_reexpression_claim(self):
        ceiling=self.data["source_claim_ceiling"]
        self.assertFalse(ceiling["orthology_verified"])
        self.assertFalse(ceiling["genome_structural_toolkit_retention_in_Ryukyu"])
        self.assertFalse(ceiling["confirmed_ancestral_white_to_coloured_reexpression"])

    def test_every_relative_identity_interval_includes_zero(self):
        for r in self.data["reported_results"]:
            lo,hi=r["bootstrap_95pct_interval"]
            self.assertLessEqual(lo,0)
            self.assertGreaterEqual(hi,0)
            self.assertTrue(r["interval_contains_zero"])
            self.assertFalse(r["enzyme_function_certified"])

    def test_actual_source_motif_is_not_enzymatic_disruption_proof(self):
        r=next(r for r in self.data["reported_results"]
               if r["case"]=="DFR_ANR_primary"
               and r["source_transcript"]=="TR39766|c0_g1_i2")
        x=r["biochemical_motif_evidence"]
        self.assertEqual((x["covered"],x["matching"],x["substituted"]),(5,4,1))
        self.assertFalse(x["function_proven"])
        other=next(r for r in self.data["reported_results"]
                   if r["case"]=="ANS_FLS"
                   and r["source_transcript"]=="TR58302|c0_g1_i1")
        self.assertEqual(other["biochemical_motif_evidence"]["covered"],2)

    def test_cannot_silently_make_better_identity_classification(self):
        mutated=deepcopy(self.data)
        mutated["reported_results"][0]["bootstrap_95pct_interval"]=[.02,.08]
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"mutated.json"
            path.write_text(json.dumps(mutated),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"zero"):
                check_source_record(path)

    def test_cannot_silently_upgrade_motif_to_function(self):
        mutated=deepcopy(self.data)
        mutated["reported_results"][0]["enzyme_function_certified"]=True
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"mutated.json"
            path.write_text(json.dumps(mutated),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"function"):
                check_source_record(path)

    def test_report_status_explicit(self):
        report=main()
        self.assertEqual(report["source_evaluations"],7)
        self.assertEqual(report["complete_NJ_family_classifications_supported"],0)
        self.assertEqual(report["noisy_relative_aa_identity_intervals_contain_zero"],7)
        self.assertFalse(report["gene_function_verified"])
        self.assertFalse(report["floral_phenotype_regain_demonstrated"])

if __name__=="__main__":
    unittest.main()
