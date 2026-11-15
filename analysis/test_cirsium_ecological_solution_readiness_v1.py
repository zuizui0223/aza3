#!/usr/bin/env python3
"""Synthetic tests: no actual field or RAD data in these fixtures."""
from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT=Path(__file__).with_name("audit_cirsium_ecological_solution_readiness_v1.py")
spec=importlib.util.spec_from_file_location("head_readiness",SCRIPT)
gate=importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(gate)


def populate(root:Path, **tables):
    for kind,relpath in gate.FILES.items():
        path=root/relpath
        path.parent.mkdir(parents=True,exist_ok=True)
        fields=sorted(gate.SCHEMA[kind] | {"observation_bout_id"})
        with path.open("w",encoding="utf-8",newline="") as handle:
            writer=csv.DictWriter(handle,fieldnames=fields,extrasaction="ignore")
            writer.writeheader()
            writer.writerows(tables.get(kind,[]))


def complete():
    k={"population_id":"POP","individual_id":"PLANT","capitulum_id":"HEAD"}
    return {
        "dna":[{**k,"voucher_id":"V1","leaf_dna_tissue_id":"DNA1","admission_status":"admitted"}],
        "shape":[{**k,"field_access_permission_id":"K","head_gravity_angle_deg":"80",
                  "minimum_phyllary_entry_gap_mm":"1.4",
                  "gravity_reference_visible":"yes","scale_reference_visible":"true"}],
        "pollen":[{**k,"floret_id":"F1","observation_bout_id":"B1",
                   "conspecific_pollen_grains":"12","pollen_assay_method":"counted",
                   "head_link_verified":"yes","admission_status":"admitted",
                   "authorization_key":"K","evidence_uri":"SOURCE1"}],
        "seed":[{**k,"achenes_assessed":"100","achenes_filled":"60",
                 "achenes_damaged":"20","head_link_verified":"yes",
                 "admission_status":"admitted","authorization_key":"K","evidence_uri":"SOURCE2"}],
        "effort":[{**k,"evaluable_duration_s":"120","valid_rate_denominator":"true",
                   "head_zone_fully_covered":"true"}],
        "events":[{**k,"record_id":"E1","oviposition_confirmed":"yes","evidence_uri":"SOURCE3"}],
    }


class AuditTests(unittest.TestCase):
    def test_blank_public_intakes_hold(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            populate(root)
            out=gate.audit(root)
            self.assertEqual(out["documented_dna_shape_pollen_seed_effort_joint_head_candidates"],0)
            self.assertTrue(out["status"].startswith("HOLD_"))
            self.assertFalse(out["functional_equivalence_or_pareto_inference_authorized"])

    def test_linked_synthetic_rows_still_not_ecological_results(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            populate(root,**complete())
            out=gate.audit(root)
            self.assertEqual(out["documented_dna_shape_pollen_seed_effort_joint_head_candidates"],1)
            self.assertEqual(out["documented_positive_oviposition_heads"],1)
            self.assertFalse(out["genotype_sequence_qc_verified_by_this_audit"])
            self.assertFalse(out["functional_equivalence_or_pareto_inference_authorized"])
            self.assertNotIn("PLANT",json.dumps(out))

    def test_different_biological_head_not_accepted_as_joint(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            data=complete()
            data["seed"][0]["capitulum_id"]="OTHER_HEAD"
            populate(root,**data)
            self.assertEqual(gate.audit(root)["documented_dna_shape_pollen_seed_effort_joint_head_candidates"],0)

    def test_visits_without_direct_pollen_or_final_seed_not_enough(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            data=complete()
            data["pollen"]=[]
            data["seed"]=[]
            populate(root,**data)
            out=gate.audit(root)
            self.assertEqual(out["documented_dna_shape_pollen_seed_effort_joint_head_candidates"],0)

    def test_seed_counts_cannot_exceed_examined(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            data=complete()
            data["seed"][0]["achenes_filled"]="101"
            populate(root,**data)
            with self.assertRaises(ValueError):
                gate.audit(root)

    def test_duplicate_pollen_floret_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            data=complete()
            data["pollen"].append(dict(data["pollen"][0]))
            populate(root,**data)
            with self.assertRaises(ValueError):
                gate.audit(root)


if __name__=="__main__":
    unittest.main()
