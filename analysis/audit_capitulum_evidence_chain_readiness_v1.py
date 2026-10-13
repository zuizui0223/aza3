#!/usr/bin/env python3
"""Cross-ledger evidence-chain readiness: measured Cirsium heads, not invented effects.

The intended unit for reproductive effects is the preregistered biological
head nested in a plant. An observed event/attempt/parasitoid emergence is NOT
a substitute for filled seeds. This script reports availability and join
consistency, not correlations, trait effects or selection.

Zero-row EAzami and aza3 CSVs yield explicit NO_REAL_FIELD_COHORT.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import math
from pathlib import Path

from validate_capitulum_video_effort_denominator_v1 import (
    K, read_csv, validate as validate_video_effort
)
from validate_capitulum_true_armature_contact_v1 import (
    validate as validate_true_armature, COLUMNS as CONTACT_COLUMNS,
    ATTEMPT_KEY, KEYS
)

HEAD_KEYS=("individual_id","population_id","capitulum_id")
BOUT_KEYS=K
ORIGINAL_HEAD_REQUIRED=set(HEAD_KEYS)|{
    "taxon","phenological_stage","natural_orientation_deg",
    "spine_length_mm","involucre_spread_deg",
    "total_achenes","filled_achenes"
}
ORIGINAL_BOUT_REQUIRED=set(BOUT_KEYS)|{
    "observation_date","start_time_local","end_time_local",
    "phenological_stage"
}

def index_rows(rows, keys, name):
    index={}
    for i,row in enumerate(rows,1):
        if not set(keys).issubset(row):
            raise ValueError(f"{name}_MISSING_COLUMNS")
        key=tuple((row.get(k) or "").strip() for k in keys)
        if not all(key):
            raise ValueError(f"{name}_EMPTY_ID_AT_{i}")
        if key in index:
            raise ValueError(f"{name}_DUPLICATE_ID_AT_{i}")
        index[key]=row
    return index

def numeric_or_none(value, field):
    if value is None or value.strip() in {"","NA","not_assessed"}:
        return None
    try:
        x=float(value)
    except ValueError:
        raise ValueError(f"{field}_INVALID_NUMERIC") from None
    if not math.isfinite(x) or x<0:
        raise ValueError(f"{field}_NONFINITE_OR_NEGATIVE")
    return x

def analyze(heads, original_bouts, efforts, events, attempts, contacts,
            video_spec):
    if any(not ORIGINAL_HEAD_REQUIRED.issubset(row) for row in heads):
        raise ValueError("NATIVE_EAZAMI_HEAD_SCHEMA_DRIFT")
    if any(not ORIGINAL_BOUT_REQUIRED.issubset(row) for row in original_bouts):
        raise ValueError("NATIVE_EAZAMI_BOUT_SCHEMA_DRIFT")
    h=index_rows(heads,HEAD_KEYS,"NATIVE_HEAD")
    b=index_rows(original_bouts,BOUT_KEYS,"NATIVE_BOUT")
    e=index_rows(efforts,BOUT_KEYS,"VIDEO_EFFORT")
    events_index=index_rows(events,KEYS,"APPROACH_EVENT")
    attempts_index=index_rows(attempts,ATTEMPT_KEY,"CONTACT_ATTEMPT")
    contacts_index=index_rows(contacts,ATTEMPT_KEY,"TRUE_CONTACT")

    for key,row in b.items():
        if key[:3] not in h:
            raise ValueError("ORIGINAL_BOUT_WITHOUT_REGISTERED_HEAD")
    for key,row in e.items():
        if key not in b or key[:3] not in h:
            raise ValueError("CAMERA_BOUT_WITHOUT_PREREGISTERED_SOURCE_HEAD")
        if row["head_stage"]!=b[key]["phenological_stage"]:
            raise ValueError("CAMERA_STAGE_CONTRADICTS_SOURCE_BOUT")
    for key,ev in events_index.items():
        if key[:4] not in e:
            raise ValueError("EVENT_WITHOUT_VIDEO_EFFORT")
    for key,a in attempts_index.items():
        if key[:5] not in events_index:
            raise ValueError("ATTEMPT_WITHOUT_INDEPENDENT_APPROACH_EPISODE")
    for key,c in contacts_index.items():
        if key not in attempts_index:
            raise ValueError("ANATOMICAL_CONTACT_WITHOUT_PARENT_ATTEMPT")

    effort_check=validate_video_effort(efforts,events,video_spec)
    true_check=validate_true_armature(contacts,attempts,events)
    valid_bouts={key for key,x in e.items() if x["valid_rate_denominator"]=="1"}
    registered_valid_footage_heads={key[:3] for key in valid_bouts}
    if not h and any((b,e,events_index,attempts_index,contacts_index)):
        raise ValueError("DOWNSTREAM_CIRSIUM_RECORDS_WITHOUT_ANY_NATIVE_HEAD")

    true_morphology=set()
    filled_counts=set()
    damaged_head_info=collections.Counter()
    for key,row in h.items():
        orient=numeric_or_none(row["natural_orientation_deg"],"HEAD_ORIENTATION")
        spine=numeric_or_none(row["spine_length_mm"],"TRUE_SPINE_LENGTH")
        phyllary=numeric_or_none(row["involucre_spread_deg"],"PHYL_SPREAD")
        filled=numeric_or_none(row["filled_achenes"],"FILLED_ACHENES")
        total=numeric_or_none(row["total_achenes"],"TOTAL_ACHENES")
        if filled is not None and (
            not float(filled).is_integer() or
            (total is not None and filled>total)):
            raise ValueError("FILLED_ACHENES_EXCEED_TOTAL_OR_NOT_INTEGER")
        if all(v is not None for v in (orient,spine,phyllary)):
            true_morphology.add(key)
        if filled is not None:
            filled_counts.add(key)
        if row.get("seed_predator_damage_fraction","") not in ("","NA"):
            damaged_head_info["head_with_raw_damage_field"]+=1

    triple=registered_valid_footage_heads & true_morphology & filled_counts
    statuses={
        "registered_native_head_cohort":(
            "PRESENT" if h else "HOLD_NO_PREREGISTERED_HEADS"),
        "same_head_valid_continuous_video":(
            "PRESENT" if valid_bouts else "HOLD_NO_VALID_FILMED_HEAD_BOUTS"),
        "independently_identified_approaches":(
            "PRESENT" if events_index else "HOLD_NO_OBSERVED_INSECT_APPROACHES"),
        "anatomically_scored_eligible_attempts":true_check["status"],
        "measured_orientation_spine_phyllary":(
            "PRESENT" if true_morphology else "HOLD_NO_TRUE_HEAD_TRAITS"),
        "filled_achene_endpoint":(
            "PRESENT" if filled_counts else "HOLD_NO_FILLED_ACHENE_COUNTS"),
        "head_with_joint_video_morphology_seed_data":(
            "PRESENT" if triple else "HOLD_NO_JOINT_HEAD_LEVEL_ENDPOINTS"),
    }
    return {
        "status":("NO_REAL_FIELD_COHORT" if not h and not e and not events
                  and not attempts and not contacts
                  else "PARTIAL_DATA_READY_FOR_DESCRIPTIVE_AUDIT"),
        "counts":{
            "original_EAzami_heads":len(h),
            "original_EAzami_bouts":len(b),
            "video_effort_bouts":len(e),
            "verified_continuous_video_bouts":len(valid_bouts),
            "independently_recorded_approach_episodes":len(events_index),
            "source_contact_attempts":len(attempts_index),
            "true_anatomical_contact_annotations":len(contacts_index),
            "heads_with_true_orientation_and_armature":len(true_morphology),
            "heads_with_filled_achene_counts":len(filled_counts),
            "heads_with_joint_verified_video_morphology_filled_achenes":len(triple),
        },
        "readiness":statuses,
        "video_denominator_status":effort_check["status"],
        "real_contact_denominator_status":true_check["status"],
        "limitations":[
            "The absence of filled-achene counts is not zero seed production.",
            "The presence of an anatomically scored attempt is not proof of adaptation.",
            "Filled achenes are not independently measured germination.",
            "A measurable association is not an experimental treatment effect.",
            "Insect-order automatic classification does not identify a guild.",
            "Recorded approach episodes are nested in biological heads and plants.",
            "No propensity to approach is inferred without complete video effort.",
            "No evolutionary origin from trait covariance or descriptive access rates.",
        ],
        "only_current_claim":"SOURCE_AVAILABILITY_AND_JOIN_INTEGRITY_NOT_NEW_BIOLOGICAL_FINDINGS",
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--native-heads",type=Path,required=True)
    ap.add_argument("--native-bouts",type=Path,required=True)
    ap.add_argument("--efforts",type=Path,default=Path(
        "data/intake/capitulum_video_effort_denominator_v1.csv"))
    ap.add_argument("--events",type=Path,default=Path(
        "data/intake/capitulum_guild_access_event_ledger_v1.csv"))
    ap.add_argument("--attempts",type=Path,default=Path(
        "data/intake/capitulum_contact_attempt_sequence_v1.csv"))
    ap.add_argument("--contacts",type=Path,default=Path(
        "data/intake/capitulum_true_armature_contact_v1.csv"))
    ap.add_argument("--effort-contract",type=Path,default=Path(
        "data/contracts/capitulum_video_effort_denominator_v1.json"))
    ap.add_argument("--out",type=Path)
    args=ap.parse_args()
    heads_header,heads=read_csv(args.native_heads)
    bouts_header,bouts=read_csv(args.native_bouts)
    eh,efforts=read_csv(args.efforts)
    gh,events=read_csv(args.events)
    ah,attempts=read_csv(args.attempts)
    ch,contacts=read_csv(args.contacts)
    spec=json.loads(args.effort_contract.read_text(encoding="utf-8"))
    if (not ORIGINAL_HEAD_REQUIRED.issubset(heads_header) or
        not ORIGINAL_BOUT_REQUIRED.issubset(bouts_header) or
        eh!=spec["columns"] or len(gh)!=32 or len(ah)!=17 or
        ch!=CONTACT_COLUMNS):
        raise ValueError("SOURCE_CSV_HEADER_OR_CONTRACT_DRIFT")
    report={
        "version":"cirsium_head_evidence_chain_readiness_v1",
        "date":"2026-10-13",
        "analysis":analyze(heads,bouts,efforts,events,attempts,contacts,spec),
    }
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",
                            encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
