"""Regression tests for actual 2026 Ryukyu six-source 8MiB full-reference signal.

Biological evidence remains limited to public young-leaf RNA fragments;
positive DFR/ANS orthology and ancestral phenotype regain are NOT supported.
"""
from copy import deepcopy
import json
import sys
import tempfile
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from validate_ryukyu_104890_full_reference_v7 import validate,PRIMARY


def mutate_and_check(mutator):
    doc=json.loads(PRIMARY.read_text(encoding="utf-8"))
    mutator(doc)
    with tempfile.TemporaryDirectory() as folder:
        p=Path(folder)/"tampered.json"
        p.write_text(json.dumps(doc),encoding="utf-8")
        validate(primary=p)


class RyukyuSixSourceAllTranscriptsTests(unittest.TestCase):
    def test_actual_six_source_result_is_audited(self):
        x=validate()
        self.assertEqual(x["full_source_2020_unigene_background_verified"],104890)
        self.assertEqual(x["six_public_young_leaf_R1_specimens"],6)
        self.assertEqual(x["full_source_competition_validated_complete_reads"],747234)
        self.assertEqual(x["strict_26_source_candidate_reads_under_full_background"],215)
        self.assertEqual(x["candidate_counts_under_full_reference"],{
            "TR39766|c0_g1_i2":66,"TR41106|c0_g1_i5":106})

    def test_8_and_16_source_prefixes_are_not_independent_replicates(self):
        x=validate()
        self.assertTrue(x["source_8MiB_16MiB_limited_reference_presence_agree"])
        self.assertTrue(x["source_prefixes_are_nested_not_independent_replication"])
        self.assertTrue(x["filter_thresholds_differ_so_counts_not_expression_fold_changes"])

    def test_no_orthology_floret_regain_claim(self):
        x=validate()
        self.assertFalse(x["confirmed_functional_DFR_ANS_orthology"])
        self.assertFalse(x["confirmed_floret_anthocyanin_expression"])
        self.assertFalse(x["confirmed_ancestral_phenotype_reexpression"])
        self.assertEqual(x["read_prefix_only_R1_bytes_per_library"],8388608)

    def test_changed_BioSample_join_fails(self):
        with self.assertRaisesRegex(ValueError,"BioSample|Run accession"):
            mutate_and_check(lambda d:d["samples"][0].update(
                {"biosample_accession":"SAMN000000000"}))

    def test_changed_cand_per_run_fails(self):
        with self.assertRaisesRegex(ValueError,"read count|totals"):
            mutate_and_check(lambda d:d["samples"][0].update(
                {"source_2020_26_candidate_strict_unique_fragments":0}))

    def test_functional_claim_injected_into_source_fails(self):
        with self.assertRaisesRegex(ValueError,"function"):
            mutate_and_check(lambda d:d["assay"].update(
                {"source_2020_candidates_functionally_assayed_for_DFR_ANS":True}))

    def test_invalid_source_background_fails(self):
        with self.assertRaisesRegex(ValueError,"version/provenance"):
            mutate_and_check(lambda d:d.update({"reference_transcripts":26}))

    def test_gene_occupancy_in_all_six_must_have_real_positive_count(self):
        with self.assertRaisesRegex(ValueError,"read count"):
            mutate_and_check(lambda d:d["samples"][0].update(
                {"TR39766_cirsium_japonicum_like_strict_fragments":0}))

    def test_no_transcript_prefix_lost_or_added(self):
        with self.assertRaisesRegex(ValueError,"exactly six"):
            mutate_and_check(lambda d:d["samples"].pop())


if __name__=="__main__":
    unittest.main()
