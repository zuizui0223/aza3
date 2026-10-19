"""Actual source-verified, complete 104890-transcript competition result tests.

No raw FASTQ, full RNA library or plant phenotypic assays are required;
tests are strict about source registration and conservative claim ceilings.
"""
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from copy import deepcopy
import json
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from validate_ryukyu_full_104890_unigene_mapping_v7 import audit,DEFAULT


class FullReferenceCompetitionSourceTests(unittest.TestCase):
    def test_all_six_cirsium_source_libraries_and_full_reference(self):
        r=audit()
        self.assertEqual(r["source_Cirsium_2020_full_transcripts_validated"],104890)
        self.assertEqual(r["focal_Ryukyu_2026_R1_biological_libraries"],6)
        self.assertEqual(r["young_leaf_R1_prefix_reads"],747234)

    def test_shared_two_sources_survive_full_2020_competition(self):
        r=audit()
        self.assertEqual(r["full_104890_reference_competition_selected_26_candidate_matches"],215)
        self.assertEqual(r["full_104890_reference_competition_TR39766_fragments"],66)
        self.assertEqual(r["full_104890_reference_competition_TR41106_fragments"],106)
        self.assertEqual(r["full_104890_reference_competition_shuffled_null_fragments"],0)
        self.assertTrue(r["matching_TR39766_and_TR41106_each_in_all_six_R1"])

    def test_matches_are_not_molecular_reexpression_or_enzyme_function(self):
        r=audit()
        self.assertFalse(r["focal_gene_orthology_verified"])
        self.assertFalse(r["functional_flower_anthocyanin_proven"])
        self.assertFalse(r["reexpression_re_evolution_supported"])

    def test_swapping_taxon_label_fails(self):
        d=json.loads(DEFAULT.read_text(encoding="utf-8"))
        d["source_specimens"][0]["taxon_as_source_title"]="Cirsium wrong"
        with TemporaryDirectory() as folder:
            file=Path(folder)/"spurious.json"
            file.write_text(json.dumps(d))
            with self.assertRaisesRegex(ValueError,"identity"):
                audit(report=file)

    def test_adding_candidate_read_without_matching_total_fails(self):
        d=json.loads(DEFAULT.read_text(encoding="utf-8"))
        d["source_specimens"][0]["strict_source_candidates_from_original_26"]+=1
        with TemporaryDirectory() as folder:
            file=Path(folder)/"spurious.json"
            file.write_text(json.dumps(d))
            with self.assertRaisesRegex(ValueError,"original source"):
                audit(report=file)

    def test_reference_is_not_a_complete_ryukyu_genome(self):
        d=json.loads(DEFAULT.read_text(encoding="utf-8"))
        d["reference_source_validation"]["not_complete_Ryukyu_genome"]=False
        with TemporaryDirectory() as folder:
            file=Path(folder)/"spurious.json"
            file.write_text(json.dumps(d))
            with self.assertRaisesRegex(ValueError,"reference"):
                audit(report=file)

    def test_gene_catalysis_evidence_upgrade_fails(self):
        d=json.loads(DEFAULT.read_text(encoding="utf-8"))
        d["claim_ceiling"]["catalytic_DFR_ANS_activity_assayed"]=True
        with TemporaryDirectory() as folder:
            file=Path(folder)/"spurious.json"
            file.write_text(json.dumps(d))
            with self.assertRaisesRegex(ValueError,"Unsupported phenotype"):
                audit(report=file)

    def test_reuse_of_same_six_biological_libraries_is_not_independent_replication(self):
        d=json.loads(DEFAULT.read_text(encoding="utf-8"))
        d["importantly_not_independent_biological_validation_of_26_transcript_run"]=False
        with TemporaryDirectory() as folder:
            file=Path(folder)/"spurious.json"
            file.write_text(json.dumps(d))
            with self.assertRaisesRegex(ValueError,"not independent"):
                audit(report=file)


if __name__=="__main__":
    unittest.main()
