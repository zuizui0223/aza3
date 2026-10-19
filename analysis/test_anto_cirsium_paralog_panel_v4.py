"""Tests for DFR/ANS versus ANR/FLS control panel with real public sequences."""
from pathlib import Path
from tempfile import TemporaryDirectory
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
sys.path.insert(0,str(ROOT/"analysis/anto_empirical"))
from cirsium_paralog_specificity_v4 import (
    GROUPS,load_panel,match_protein_segments,analyze_transcript,run_panel
)
from cirsium_source_homology_v3 import read_fasta
from run_anto_cirsium_paralog_panel_v4 import main

DATA=ROOT/"data/anto_source_handoff/molecular"

class ParalogReferencePanelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panel=load_panel(
            DATA/"extended_DFR_ANS_ANR_FLS_reference_proteins_v4.fasta",
            DATA/"extended_reference_protein_provenance_v4.json")

    def test_seven_independently_annotated_positive_and_negative_references(self):
        self.assertEqual(len(self.panel),7)
        self.assertEqual({r["group"] for r in self.panel},set(GROUPS))
        for r in self.panel:
            self.assertTrue(r["protein"])
            self.assertGreater(len(r["protein"]),300)

    def test_every_reference_self_classifies_to_registered_family(self):
        # This validates *reference specificity*, not correct Cirsium function.
        for ref in self.panel:
            with self.subTest(accession=ref["accession"]):
                p=match_protein_segments([{
                    "aa_segment":ref["protein"],
                    "strand":"+","reading_frame":0,
                    "segment_length_aa":len(ref["protein"]),
                }],self.panel)
                winner=max(p,key=lambda group:p[group]["normalized_alignment_score"])
                self.assertEqual(winner,ref["group"])

    def test_negative_controls_not_trivially_equivalent_to_positives(self):
        for ref in self.panel:
            if ref["group"] not in ("ANR_negative_control","FLS_negative_control"):
                continue
            p=match_protein_segments([{
                "aa_segment":ref["protein"],"strand":"+",
                "reading_frame":0,"segment_length_aa":len(ref["protein"])
            }],self.panel)
            self.assertGreater(p[ref["group"]]["normalized_alignment_score"],
                max(p[g]["normalized_alignment_score"] for g in
                    ("DFR_reference","ANS_reference")))

    def test_26_pub_source_sequences_remain_distinct_from_2026_ryukyu(self):
        r=main()
        self.assertEqual(r["n_publisher_cirsium_transcripts"],26)
        self.assertEqual(r["n_reference_proteins"],7)
        self.assertEqual(r["2026_Ryukyu_focal_genotypes_included"],False)
        self.assertFalse(r["flower_colour_regain_demonstrated"])
        self.assertFalse(r["functional_orthology_demonstrated"])
        self.assertEqual(len(r["records"]),26)
        self.assertEqual(sum(r["family_top_hit_counts"].values()),26)
        self.assertTrue(all(not row["orthology_verified"] for row in r["records"]))

    def test_uninformative_short_sequence_does_not_prove_gene_absence(self):
        result=analyze_transcript("ATGTAA",self.panel)
        self.assertEqual(result["screen"],"NO_TRANSLATED_SEGMENT_LONG_ENOUGH")
        self.assertFalse(result["orthology_verified"])
        self.assertFalse(result["enzyme_function_verified"])

    def test_mutated_source_candidate_fails_reproducibility(self):
        src=DATA/"cirsium_japonicum_2020_26_source_sequences_v6.fasta"
        s=src.read_text(encoding="utf-8")
        ix=s.index("\n")+1
        with TemporaryDirectory() as temp:
            changed=Path(temp)/"swapped.fa"
            changed.write_text(s[:ix]+("A" if s[ix]!="A" else "C")+s[ix+1:])
            with self.assertRaisesRegex(ValueError,"Corrupted source transcript"):
                run_panel(changed,
                    DATA/"cirsium_japonicum_2020_26_sequence_provenance_v6.json",
                    DATA/"extended_DFR_ANS_ANR_FLS_reference_proteins_v4.fasta",
                    DATA/"extended_reference_protein_provenance_v4.json")

    def test_positive_family_support_does_not_imply_flower_function(self):
        r=main()
        for x in r["records"]:
            self.assertFalse(x["enzyme_function_verified"])
            self.assertFalse(x["floral_pathway_function_verified"])
        self.assertFalse(r["functional_orthology_demonstrated"])

if __name__=="__main__":
    unittest.main()
