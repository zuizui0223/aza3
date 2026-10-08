#!/usr/bin/env python3
"""Fail-closed, real-anatomy head contact ontology for Cirsium access.

Only a verified insect touching actual phyllary/spine tissue AND verified
attempt-level physical blockage can be called a physical barrier. Candidate
seed feeders and pollinators are independently classified before their
approach outcomes. Repeated attempt rows do not create independent heads.

No field events have been committed; tests are explicitly SYNTHETIC.
"""
from __future__ import annotations
import argparse
import copy
import csv
import json
import math
from pathlib import Path

COLUMNS=[
 "record_id","individual_id","population_id","capitulum_id",
 "observation_bout_id","approach_episode_id","attempt_index",
 "video_clip_id","head_stage","contact_surface","contact_resolved",
 "spine_touch_verified","physical_blockage_verified",
 "entry_not_achieved_at_attempt","actual_spine_length_mm",
 "minimum_phyllary_gap_mm","anatomy_scale_evidence_uri",
 "video_contact_evidence_uri","pre_entry_taxon_evidence_uri",
 "scorer_blinded_to_oviposition_result"
]
KEYS=("individual_id","population_id","capitulum_id",
      "observation_bout_id","approach_episode_id")
ATTEMPT_KEY=KEYS+("attempt_index",)
SURFACES={
 "spine_tip","spine_shaft","phyllary_lamina","phyllary_gap",
 "floret_disc","stem","no_surface_contact","unresolved"
}
ARMATURE_SURFACES={"spine_tip","spine_shaft","phyllary_lamina","phyllary_gap"}
VERIFIED_GUILDS={
 "pre_entry_taxon_and_local_life_history_verified",
 "independently_assigned_before_access"
}
FLAGS={"0","1","NA"}

def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        return list(reader.fieldnames or []),list(reader)

def validate(rows, attempts, events):
    if not rows:
        return {"status":"NO_REAL_TRUE_ARMATURE_CONTACT_ANNOTATIONS",
                "n_contact_rows":0, "n_verified_spine_touches":None,
                "n_verified_physical_blockages":None,
                "n_independent_heads":0,
                "empirical_inference":"NOT_IDENTIFIABLE"}

    at_index={}
    for a in attempts:
        key=tuple(a.get(k,"") for k in ATTEMPT_KEY)
        if not all(key) or key in at_index:
            raise ValueError("DUPLICATE_OR_MISSING_PARENT_CONTACT_ATTEMPT")
        at_index[key]=a
    ev_index={}
    for ev in events:
        key=tuple(ev.get(k,"") for k in KEYS)
        if not all(key) or key in ev_index:
            raise ValueError("DUPLICATE_OR_MISSING_PARENT_APPROACH_EPISODE")
        ev_index[key]=ev
    counts={"touches":0,"blockages":0,"unknown":0}
    seen=set()
    heads=set()
    for i,row in enumerate(rows,1):
        if set(row)!=set(COLUMNS):
            raise ValueError("ANNOTATION_HEADER_OR_COLUMN_DRIFT")
        key=tuple(row[k].strip() for k in ATTEMPT_KEY)
        if not all(key) or key in seen:
            raise ValueError("ANNOTATION_EMPTY_OR_DUPLICATE_ATTEMPT")
        seen.add(key)
        a=at_index.get(key)
        e=ev_index.get(key[:len(KEYS)])
        if a is None or e is None:
            raise ValueError("ANNOTATION_MUST_MATCH_PARENT_ATTEMPT_AND_EPISODE")
        if row["video_clip_id"]!=a["video_clip_id"] or row["head_stage"]!=a["head_stage"]:
            raise ValueError("ANNOTATION_VIDEO_OR_STAGE_MISMATCH")
        if row["head_stage"]!=e["phenological_stage"]:
            raise ValueError("ANNOTATION_EPISODE_STAGE_MISMATCH")
        if (e["pre_entry_guild"]=="unresolved" or
            e["pre_entry_guild_evidence"] not in VERIFIED_GUILDS or
            not row["pre_entry_taxon_evidence_uri"].strip()):
            raise ValueError("INSECT_ROLE_CANNOT_BE_ASSIGNED_AFTER_EGG_OR_POLLEN_SUCCESS")
        if a["continuous_coverage"]!="1":
            raise ValueError("UNCOVERED_CONTACT_SEQUENCE_CANNOT_PROVE_BARRIER")
        if row["contact_surface"] not in SURFACES:
            raise ValueError("UNRECOGNIZED_TRUE_BOTANICAL_CONTACT_SURFACE")
        for col in ("contact_resolved","spine_touch_verified",
                    "physical_blockage_verified","entry_not_achieved_at_attempt"):
            if row[col] not in FLAGS:
                raise ValueError("INVALID_TRI_STATE_CONTACT_FLAG")
        if row["scorer_blinded_to_oviposition_result"] not in {"0","1"}:
            raise ValueError("REVIEW_BLINDING_FLAG_MUST_BE_EXPLICIT")
        if not row["video_contact_evidence_uri"].strip():
            raise ValueError("VIDEO_CITATION_REQUIRED_FOR_EVERY_REAL_ANNOTATION")
        surface=row["contact_surface"]
        spine=row["spine_touch_verified"]
        if surface=="unresolved":
            if row["contact_resolved"]!="0" or spine!="NA":
                raise ValueError("UNRESOLVED_SURFACE_CANNOT_BE_NO_SPINE_TOUCH")
            counts["unknown"]+=1
        else:
            if row["contact_resolved"]!="1":
                raise ValueError("CLASSIFIED_SURFACE_NEEDS_REVIEWED_CONTACT")
            if surface in {"spine_tip","spine_shaft"} and spine!="1":
                raise ValueError("TRUE_SPINE_SURFACE_REQUIRES_VERIFIED_TOUCH")
            if surface not in {"spine_tip","spine_shaft"} and spine!="0":
                raise ValueError("NONSPINE_OR_NONE_CANNOT_BE_SPINE_TOUCH")
        if spine=="1":
            counts["touches"]+=1
        if row["entry_not_achieved_at_attempt"]=="1" and a["access_observed"]!="0":
            raise ValueError("FAILED_ACCESS_CONTRADICTS_PARENT_ATTEMPT")
        if (row["physical_blockage_verified"]=="1" or
            row["entry_not_achieved_at_attempt"]=="1"):
            if row["contact_resolved"]!="1":
                raise ValueError("UNASSESSED_CONTACT_CANNOT_SUPPORT_CONFIRMED_FAILURE")
        if row["physical_blockage_verified"]=="1":
            if (surface not in ARMATURE_SURFACES or
                row["entry_not_achieved_at_attempt"]!="1" or
                a["obstruction_observed"]!="1" or
                a["access_observed"]!="0"):
                raise ValueError("PHYSICAL_BLOCKAGE_REQUIRES_REAL_BRACT_TOUCH_AND_PARENT_OBSTRUCTION")
            if not row["anatomy_scale_evidence_uri"].strip():
                raise ValueError("UNSCALED_ARMATURE_CANNOT_IDENTIFY_MECHANICAL_GATE")
            counts["blockages"]+=1
        for col in ("actual_spine_length_mm","minimum_phyllary_gap_mm"):
            value=row[col]
            if value not in ("NA",""):
                try:x=float(value)
                except ValueError:raise ValueError("ANATOMICAL_SCALE_NOT_NUMERIC") from None
                if not math.isfinite(x) or x<0:
                    raise ValueError("ANATOMICAL_DIMENSION_MUST_BE_NONNEGATIVE")
        if row["physical_blockage_verified"]=="1":
            required_dimension=("actual_spine_length_mm" if surface in
                                {"spine_tip","spine_shaft"} else
                                "minimum_phyllary_gap_mm")
            if row[required_dimension] in ("","NA"):
                raise ValueError("CONFIRMED_MECHANICAL_OBSTRUCTION_LACKS_RELEVANT_MEASURED_ANATOMY")
            if surface in {"spine_tip","spine_shaft"} and float(row[required_dimension])<=0:
                raise ValueError("POSITIVE_TRUE_SPINE_TOUCH_REQUIRES_NONZERO_SPINE_LENGTH")
        heads.add(tuple(row[k] for k in
                        ("population_id","individual_id","capitulum_id")))
    return {
        "status":"DESCRIPTIVE_TRUE_CONTACT_GATES_ONLY",
        "n_contact_rows":len(rows),
        "n_verified_spine_touches":counts["touches"],
        "n_verified_physical_blockages":counts["blockages"],
        "n_unresolved_contacts":counts["unknown"],
        "n_independent_heads":len(heads),
        "empirical_inference":"CONTACT_BEHAVIOUR_ONLY_NOT_CAUSAL_SELECTION_OR_FITNESS",
        "cannot_infer":["treatment effects", "seed set", "genotypic fitness",
                        "selection on trait covariance", "historical coevolution"],
    }

def synthetic_tests():
    # Independent pre-entry candidate guild, parent attempt recording
    # confirmed involucre obstruction with complete coverage.
    a={
        "individual_id":"I1","population_id":"POP1","capitulum_id":"H1",
        "observation_bout_id":"B1","approach_episode_id":"E1",
        "attempt_index":"1","video_clip_id":"synthetic://clip",
        "head_stage":"bud","obstruction_observed":"1",
        "access_observed":"0","continuous_coverage":"1",
    }
    e={
        **{k:a[k] for k in KEYS},
        "pre_entry_guild":"seed_feeder_candidate",
        "pre_entry_guild_evidence":"independently_assigned_before_access",
        "phenological_stage":"bud",
    }
    x=dict(zip(COLUMNS,[
        "R1",a["individual_id"],a["population_id"],a["capitulum_id"],
        a["observation_bout_id"],a["approach_episode_id"],"1",
        a["video_clip_id"],"bud","spine_tip","1","1","1","1",
        "3.4","NA","synthetic://caliper","synthetic://clip",
        "synthetic://taxon","1"
    ]))
    def check(row=x, attempt=a, episode=e):
        return validate([copy.deepcopy(row)],[copy.deepcopy(attempt)],
                        [copy.deepcopy(episode)])
    assert check()["n_verified_physical_blockages"]==1
    assert check()["n_verified_spine_touches"]==1
    assert validate([],[],[])["n_verified_spine_touches"] is None
    bad=[
        ({**x,"spine_touch_verified":"0"},a,e,"TRUE_SPINE_SURFACE"),
        ({**x,"anatomy_scale_evidence_uri":""},a,e,"UNSCALED_ARMATURE"),
        ({**x,"actual_spine_length_mm":"NA"},a,e,"LACKS_RELEVANT_MEASURED"),
        ({**x,"entry_not_achieved_at_attempt":"0"},a,e,"PHYSICAL_BLOCKAGE"),
        ({**x,"video_contact_evidence_uri":""},a,e,"VIDEO_CITATION"),
        (x,{**a,"continuous_coverage":"0"},e,"UNCOVERED_CONTACT"),
        (x,a,{**e,"pre_entry_guild_evidence":"unresolved"},"INSECT_ROLE"),
        ({**x,"contact_surface":"unresolved"},a,e,"UNRESOLVED_SURFACE"),
        ({**x,"approach_episode_id":"E_DOES_NOT_EXIST"},a,e,"ANNOTATION_MUST_MATCH"),
        ({**x,"actual_spine_length_mm":"-2"},a,e,"POSITIVE_TRUE_SPINE_TOUCH"),
    ]
    for rr,aa,ee,needle in bad:
        try:check(rr,aa,ee)
        except ValueError as er:
            assert needle in str(er),(needle,str(er))
        else:
            raise AssertionError(f"INVALID_REAL_CONTACT_INFERRED:{needle}")
    return {
        "status":"PASS_SYNTHETIC_REAL_ARMATURE_AND_PARENT_LINK_GATES",
        "positive_fixtures":1,"rejected_bad_scientific_claims":len(bad),
        "real_organismal_field_records":0,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--contact",type=Path,default=Path(
        "data/intake/capitulum_true_armature_contact_v1.csv"))
    p.add_argument("--attempts",type=Path,default=Path(
        "data/intake/capitulum_contact_attempt_sequence_v1.csv"))
    p.add_argument("--events",type=Path,default=Path(
        "data/intake/capitulum_guild_access_event_ledger_v1.csv"))
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    ch,contacts=read_csv(a.contact)
    ah,attempts=read_csv(a.attempts)
    eh,events=read_csv(a.events)
    if ch!=COLUMNS or not set(ATTEMPT_KEY).issubset(ah) or not set(KEYS).issubset(eh):
        raise ValueError("CONTACT_ATTEMPT_OR_EVENT_CSV_CONTRACT_DRIFT")
    result={
        "version":"capitulum_true_armature_contact_gate_v1",
        "date":"2026-10-08",
        "synthetic_tests":synthetic_tests(),
        "actual_intake":validate(contacts,attempts,events),
        "science":"True structural barrier = independently identified before access, documented anatomical touch, verified blockage, unsuccessful attempt. No phylogenetic or plant-fitness effect from this alone.",
    }
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",
                         encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
