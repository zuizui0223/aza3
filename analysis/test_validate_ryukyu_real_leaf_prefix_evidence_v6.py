"""Regression checks on real 2026 Cirsium young-leaf public read fragments."""
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from validate_ryukyu_real_leaf_prefix_evidence_v6 import audit, PRIMARY,SECONDARY

class RealRyukyuLeafReadEvidenceTests(unittest.TestCase):
    def test_real_six_source_libraries_and_source_like_fragments(self):
        result=audit()
        self.assertEqual(result["n_independent_source_libraries"],6)
        self.assertEqual(result["primary_16MiB_total_R1_reads"],1499192)
        self.assertEqual(result["primary_16MiB_stringent_2020_source_like_fragments"],636)
        self.assertEqual(result["primary_16MiB_shuffled_decoy_fragments"],0)
        self.assertEqual(result["primary_16MiB_distinct_source_transcript_IDs_with_matches"],18)

    def test_competing_smaller_prefix_and_distinct_thresholds_not_independent(self):
        r=audit()
        self.assertEqual(r["supporting_8MiB_total_R1_reads"],747995)
        self.assertEqual(r["supporting_8MiB_stringent_source_fragments"],378)
        self.assertEqual(r["supporting_8MiB_shuffled_decoy_fragments"],0)
        self.assertTrue(r["same_biological_libraries_not_independent_experiments"])

    def test_dfr_like_and_ans_blast_source_detected_in_all_six_leaf_libraries(self):
        r=audit()
        self.assertTrue(r["TR39766_DFR_related_source_in_six_leaf_libraries"])
        self.assertTrue(r["TR41106_ANS_Blast_annotated_source_in_six_leaf_libraries"])
        self.assertEqual(r["TR39766_16MiB_fragments"],161)
        self.assertEqual(r["TR41106_16MiB_fragments"],270)

    def test_cannot_conclude_flower_colour_reexpression_or_function(self):
        r=audit()
        self.assertFalse(r["functional_orthology_assayed"])
        self.assertFalse(r["floral_expression_assayed"])
        self.assertFalse(r["genetic_evolutionary_reacquisition_demonstrated"])
        self.assertFalse(r["between_species_differential_expression_permitted"])

    def test_swapped_sample_ids_fail_closed(self):
        d=json.loads(PRIMARY.read_text(encoding="utf-8"))
        d["per_run"][0]["biosample_accession"]="SAMN00000000"
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/"bogus.json"
            p.write_text(json.dumps(d),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"BioSample"):
                audit(primary=p)

    def test_unsupported_per_candidate_added_counts_fail_closed(self):
        d=json.loads(PRIMARY.read_text(encoding="utf-8"))
        d["per_run"][0]["per_published_2020_Cirsium_source_transcript_strict_counts"]["TR39766|c0_g1_i2"]+=1
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/"bogus.json"
            p.write_text(json.dumps(d),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"source counts"):
                audit(primary=p)

    def test_arbitrary_family_not_on_original_2020_26_reference_rejected(self):
        d=json.loads(PRIMARY.read_text(encoding="utf-8"))
        d["per_run"][0]["per_published_2020_Cirsium_source_transcript_strict_counts"]["INVENTED_DFR"]=2
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/"bogus.json"
            p.write_text(json.dumps(d),encoding="utf-8")
            with self.assertRaises(ValueError):
                audit(primary=p)

    def test_8MiB_data_cannot_be_compressed_independent_biological_replication(self):
        d=json.loads(SECONDARY.read_text(encoding="utf-8"))
        d["per_run"][0]["strict_decoy_fragments"]=1
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/"bogus.json"
            p.write_text(json.dumps(d),encoding="utf-8")
            with self.assertRaises(ValueError):
                audit(secondary=p)


if __name__=="__main__":
    unittest.main()
