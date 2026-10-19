"""No source-backed family method concordance may be called orthology."""
from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from validate_cirsium_2020_method_concordance_v5 import audit_frozen,compare_analysis,main

class FrozenMethodConcordanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=audit_frozen()

    def test_seven_query_configurations_five_source_transcripts(self):
        self.assertEqual(self.source["n_query_configurations"],7)
        self.assertEqual(self.source["distinct_source_transcripts"],5)
        self.assertEqual(self.source["n_tree_configurations"],3)

    def test_nj_and_approxml_are_equally_unresolved_not_supported_function(self):
        self.assertTrue(all(r["NJ_unrooted_reference_family_placement"]=="BETWEEN_FAMILY_ANCHORS"
                            and r["FastTree_WAG_unrooted_reference_family_placement"]=="BETWEEN_FAMILY_ANCHORS"
                            for r in self.source["records"]))
        self.assertFalse(any(r["functional_orthology_verified"] for r in self.source["records"]))
        self.assertEqual(self.source["n_candidates_confidently_identified_as_functional_dfr_or_ans"],0)

    def test_no_ryukyu_coding_data_included(self):
        d=main()
        self.assertFalse(d["Ryukyu_coloured_flower_regain_tested"])
        self.assertFalse(d["biochemical_activity_assayed"])
        self.assertEqual(d["confident_DFR_ANS_family_labels"],0)

    def test_relabeling_ambiguous_tree_fails_closed(self):
        d=deepcopy(self.source)
        d["records"][0]["FastTree_WAG_unrooted_reference_family_placement"]="POSITIVE_FAMILY_SIDE"
        with TemporaryDirectory() as tmp:
            file=Path(tmp)/"mutated.json"
            file.write_text(json.dumps(d),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"Family placement"):
                audit_frozen(file)

    def test_unsupported_function_claim_is_not_admitted(self):
        d=deepcopy(self.source)
        d["n_candidates_confidently_identified_as_functional_dfr_or_ans"]=2
        with TemporaryDirectory() as tmp:
            file=Path(tmp)/"mutated.json"
            file.write_text(json.dumps(d),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"Unsupported"):
                audit_frozen(file)

    def test_original_four_vs_three_isoform_reuse_is_not_independent_event(self):
        records=self.source["records"]
        dfr={r["source_transcript"] for r in records if r["case"].startswith("DFR")}
        self.assertEqual(len(dfr),4)
        self.assertEqual(len({x.rsplit("_i",1)[0] for x in dfr}),3)
        self.assertIn("TR38875|c0_g4_i5",dfr)
        self.assertIn("TR38875|c0_g4_i7",dfr)

if __name__=="__main__":
    unittest.main()
