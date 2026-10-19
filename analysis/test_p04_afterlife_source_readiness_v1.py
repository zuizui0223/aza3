"""Test source-only evidence qualification; never turn missing assays into negative ecology."""
from __future__ import annotations
from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"/"anto_empirical"))
import p04_afterlife_source_readiness_v1 as p


class AfterlifeP04SourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest=p.csv_rows("manifest")
        cls.qc=p.json_data("qc")
        cls.morph=p.csv_rows("takaoense")
        cls.family=p.json_data("family")
        cls.panel=p.csv_rows("panel")

    def build(self, manifest=None, qc=None, morph=None, family=None, panel=None):
        return p.group_readiness(
            deepcopy(self.manifest if manifest is None else manifest),
            deepcopy(self.qc if qc is None else qc),
            deepcopy(self.morph if morph is None else morph),
            deepcopy(self.family if family is None else family),
            deepcopy(self.panel if panel is None else panel))

    def test_actual_source_ledger_roundtrip(self):
        result=p.run()
        frozen=json.loads((ROOT/"data/anto_source_handoff/p04_afterlife_public_source_readiness_v1.json").read_text(encoding="utf-8"))
        self.assertEqual(result,frozen)
        s=result["source_summary"]
        self.assertEqual((s["public_run_count"],s["young_leaf_runs"],s["floral_RNA_runs_in_this_source"]),(25,25,0))
        self.assertEqual((s["Ryukyu_brevicaule_leaf_runs"],s["Ryukyu_irumtiense_leaf_runs"]),(3,3))
        self.assertEqual((s["takaoense_W_leaf_libraries"],s["takaoense_BP_leaf_libraries"]),(3,3))
        self.assertEqual((s["unresolved_2020_NJ_and_FastTree_evaluations"],s["distinct_2020_candidate_transcripts"]),(7,5))
        self.assertEqual(s["phaseA_target_individuals"],298)
        self.assertEqual(len(result["per_taxon"]),12)

    def test_no_source_run_confused_with_functionally_verified_floral_resource(self):
        r=self.build()
        for x in r["per_taxon"]:
            self.assertEqual(x["public_flower_stage_RNA_seq_runs"],0)
            self.assertEqual(x["same_plant_verified_functional_toolkit"],p.UNMEASURED)
            self.assertEqual(x["heritable_repeated_reexpression"],p.UNMEASURED)
            self.assertEqual(x["reexpression_origin_history"],p.UNKNOWN)
            self.assertEqual(x["actual_phase_A_individuals_collected_by_this_work"],0)

    def test_duplicate_public_run_fails_closed(self):
        rows=deepcopy(self.manifest)
        rows[1]["run_accession"]=rows[0]["run_accession"]
        with self.assertRaisesRegex(ValueError,"Repeated source accession"):self.build(manifest=rows)

    def test_tissue_drift_fails_closed(self):
        rows=deepcopy(self.manifest)
        rows[0]["tissue_from_project_description"]="flower"
        with self.assertRaisesRegex(ValueError,"Tissue not uniformly young leaf"):self.build(manifest=rows)

    def test_nonadmitted_flower_colour_is_not_filled_from_species_label(self):
        rows=deepcopy(self.manifest)
        rows[0]["individual_flower_colour"]="WHITE"
        with self.assertRaisesRegex(ValueError,"phenotype status changed"):self.build(manifest=rows)

    def test_morph_crosswalk_one_to_one_provenance(self):
        rows=deepcopy(self.morph)
        rows[0]["specimen_code_from_experiment_title"]="IMAGINARY"
        with self.assertRaisesRegex(ValueError,"wrong biological specimen"):self.build(morph=rows)

    def test_no_upgrading_of_gene_family_from_matching_methods(self):
        family=deepcopy(self.family)
        family["records"][0]["functional_orthology_verified"]=True
        with self.assertRaisesRegex(ValueError,"Different family placement"):self.build(family=family)

    def test_no_silent_sampling_panel_change(self):
        panel=deepcopy(self.panel)
        for row in panel:
            if row["phase"]=="A":
                row["recommended_n"]="50"
                break
        with self.assertRaisesRegex(ValueError,"Phase-A denominator"):self.build(panel=panel)

    def test_source_file_pin_matches_and_wrong_pin_fails(self):
        b=p.load_source("manifest")
        self.assertTrue(len(b)>8000)
        self.assertEqual(p.git_blob_sha(b),p.SOURCE["manifest"][1])
        with self.assertRaisesRegex(ValueError,"source blob identity changed"):
            p.load_source("manifest",root=ROOT/"not_existing_placeholder") if False else self._simulate_bad_source(b)

    def _simulate_bad_source(self,b):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as td:
            root=Path(td)
            relative=p.SOURCE["manifest"][0]
            dest=root/relative
            dest.parent.mkdir(parents=True)
            dest.write_bytes(b+b"X")
            return p.load_source("manifest",root=root)


if __name__=="__main__":unittest.main()
