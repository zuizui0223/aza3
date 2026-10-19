"""Regression tests: source-backed Cirsium empirical algorithms now live in aza3."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
sys.path.insert(0,str(ROOT/"analysis/anto_empirical"))
from run_anto_cirsium_empirical_sources_v3 import main
from cirsium_source_homology_v3 import ref_proteins,read_fasta,run_study,triage_sequence
from cirsium_run_intake_v3 import load_rows
from cirsium_organ_evidence_v3 import validate_registry

DAT=ROOT/"data/anto_source_handoff"


class MigratedEmpiricalSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=main(output=Path("/tmp/aza3_cirsium_empirical_test_result_v3.json"))

    def test_canonical_source_owner(self):
        self.assertEqual(self.result["source_owner"],"zuizui0223/aza3")
        self.assertEqual(self.result["upstream_theory_owner"],"zuizui0223/anto")

    def test_real_published_not_independent_Ryukyu_data(self):
        self.assertEqual(self.result["cirsium_public_2026_sra"]["manifest_n_runs"],25)
        self.assertEqual(self.result["cirsium_public_2026_sra"]["ryukyu_n_runs"],6)
        self.assertEqual(self.result["public_asteraceae_reference_preflight"]["registered_screen_count"],18)
        self.assertEqual(self.result["public_asteraceae_reference_preflight"]["n_candidate_sequences_supplied"],0)
        self.assertEqual(self.result["cirsium_public_2020_protein_homology_triage"]["n_transcripts"],26)
        self.assertFalse(self.result["cirsium_public_2020_protein_homology_triage"]["Ryukyu_2026_sequences_included"])

    def test_distinct_transcripts_not_assumed_gene_copy_count(self):
        p=self.result["cirsium_public_2020_protein_homology_triage"]
        self.assertLessEqual(p["n_distinct_assembly_groups"],26)
        self.assertGreaterEqual(p["n_distinct_assembly_groups"],1)
        self.assertTrue(all(r["orthology_confirmed"] is False for r in p["records"]))
        self.assertTrue(all(r["enzyme_function_confirmed"] is False
                            for r in p["records"] if "enzyme_function_confirmed" in r))

    def test_protein_reference_controls(self):
        refs=ref_proteins(DAT/"molecular/asteraceae_reference_cds_v5.fasta")
        self.assertEqual(len(refs["DFR"]),342)
        self.assertEqual(len(refs["ANS"]),356)
        nt=read_fasta(DAT/"molecular/asteraceae_reference_cds_v5.fasta")
        smdfr=next(v for k,v in nt.items() if k.startswith("EF600682.1|"))
        classification=triage_sequence(smdfr,refs)
        self.assertEqual(classification["best_reference_family"],"DFR")
        self.assertTrue(classification["orthology_confirmed"] is False)

    def test_registered_leaf_runs_cannot_become_flower_assays(self):
        rows=load_rows(DAT/"public_sra/prjna1311153_run_manifest_20261017.csv")
        self.assertTrue(all(x["individual_flower_colour"]=="NOT_MEASURED_IN_THIS_DATASET" for x in rows))
        self.assertTrue(all(x["tissue_from_project_description"]=="young leaf" for x in rows))

    def test_organ_expression_remains_Cj2020_not_Ryukyu(self):
        doc=validate_registry(DAT/"molecular/cirsium_japonicum_organ_evidence_v6.json")
        self.assertEqual(len(doc["records"]),11)
        self.assertFalse(doc["six_ryukyu_specimens_linked"])
        self.assertFalse(doc["trait_regain_event_observed"])

    def test_accepted_claim_bounds(self):
        for key,val in self.result["interpretation_gates"].items():
            self.assertIs(val,False,f"Inappropriate biological inference promoted: {key}")

    def test_output_materializes_and_is_reproducible(self):
        with TemporaryDirectory() as d:
            output=Path(d)/"source_report.json"
            r=main(output)
            a=json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(a["analysis_version"],r["analysis_version"])
            self.assertEqual(a["cirsium_public_2020_protein_homology_triage"]["n_transcripts"],26)


if __name__=="__main__":
    unittest.main()
