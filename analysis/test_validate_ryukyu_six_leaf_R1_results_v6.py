"""Negative controls for actual first six Ryukyu 2026 leaf-R1 candidate matches."""
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from validate_ryukyu_six_leaf_R1_results_v6 import validate,EVIDENCE


class RyukyuSixEmpiricalEvidenceAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=json.loads(EVIDENCE.read_text(encoding="utf-8"))

    def _changed_source(self,patch):
        with TemporaryDirectory() as d:
            data=deepcopy(self.source)
            patch(data)
            path=Path(d)/"source.json"
            path.write_text(json.dumps(data),encoding="utf-8")
            with self.assertRaises(ValueError):
                validate(evidence=path)

    def test_true_six_source_run_results(self):
        v=validate()
        self.assertEqual(v["independent_source_accessions"],6)
        self.assertEqual(v["real_leaf_R1_reads_processed"],1499192)
        self.assertEqual(v["candidate_source_transcript_R1_fragment_alignments"],636)
        self.assertEqual(v["TR39766_mapping_reads_across_all_six"],161)
        self.assertEqual(v["shuffled_decoy_matches"],0)
        self.assertEqual(v["original_2020_candidate_ids_detected_in_selected_reference_panel"],18)

    def test_three_by_three_biological_source_not_reversal_replication(self):
        counts={x["taxon_as_study_title"] for x in self.source["rows"]}
        self.assertEqual(counts,{"Cirsium brevicaule","Cirsium irumtiense"})
        self.assertEqual(len({x["specimen_code"] for x in self.source["rows"]}),6)
        self.assertTrue(all(x["R1_fragments_to_2020_TR39766"]>0 for x in self.source["rows"]))
        self.assertFalse(self.source["claim_ceiling"]["confirmed_orthology"])
        self.assertFalse(self.source["claim_ceiling"]["historical_ancestral_reexpression_event_inferred"])

    def test_false_gene_function_claim_rejected(self):
        self._changed_source(
            lambda doc:doc["claim_ceiling"].__setitem__(
                "pigment_structural_gene_function_verified",True))

    def test_swapped_biosample_or_taxon_rejected(self):
        self._changed_source(
            lambda doc:doc["rows"][0].__setitem__("sample_accession","SAMN00000000"))
        self._changed_source(
            lambda doc:doc["rows"][0].__setitem__("taxon_as_study_title","Cirsium bogus"))

    def test_fake_source_hit_count_rejected(self):
        self._changed_source(
            lambda doc:doc["rows"][0].__setitem__(
                "R1_fragments_to_2020_TR39766",doc["rows"][0]["R1_fragments_to_2020_TR39766"]+1))

    def test_nonrandom_r1_cannot_be_promoted_to_whole_transcriptome(self):
        self.assertFalse(self.source["read_processing"]["read_order_random_sample"])
        self.assertTrue(self.source["read_processing"]["R1_only"])
        self._changed_source(
            lambda doc:doc["read_processing"].__setitem__("read_order_random_sample",True))

if __name__=="__main__":
    unittest.main()
