#!/usr/bin/env python3
"""Enforce visible head+bout denominators before analysing insect access rates.

This source-only validator is deliberately zero-safe:
  - no video effort means NO REAL DATA (not zero arrival)
  - no event records with independently verified full video CAN be zero
  - event-triggered recordings cannot estimate unconditional arrival rates.
Synthetic examples are never added to real observation CSVs.
"""
from __future__ import annotations
import argparse
import copy
import csv
import json
from pathlib import Path

K=("individual_id","population_id","capitulum_id","observation_bout_id")
EPS=1e-7

def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        return list(reader.fieldnames or []), list(reader)

def fail(message):
    raise ValueError(message)

def validate(efforts,events,spec):
    index={}
    for i,b in enumerate(efforts,1):
        if any(x not in b for x in spec["columns"]):
            fail("EFFORT_SCHEMA_MISSING_FIELD")
        # Alternative blooming Cirsium resources can redistribute flower-head
        # weevil oviposition; do not turn an unmeasured resource into zero.
        census=b["other_cirsium_resource_census_status"]
        if census not in spec["allowed"]["other_cirsium_resource_census_status"]:
            fail(f"EFFORT_{i}_UNKNOWN_OTHER_CIRSIUM_CENSUS_STATUS")
        req=("other_cirsium_buds_count","other_cirsium_open_heads_count",
             "other_cirsium_census_radius_m","other_cirsium_census_reference")
        if census=="not_assessed":
            if any(b[k].strip() for k in req):
                fail(f"EFFORT_{i}_UNMEASURED_OTHER_CIRSIUM_MUST_NOT_BE_ZERO_OR_GUESSED")
        else:
            try:
                bcount=int(b[req[0]])
                ocount=int(b[req[1]])
                radius=float(b[req[2]])
            except (ValueError,TypeError):
                fail(f"EFFORT_{i}_INVALID_OTHER_CIRSIUM_COUNTS")
            if (not b[req[0]].isdigit() or not b[req[1]].isdigit() or
                bcount<0 or ocount<0 or not 0<radius<float("inf") or
                not b[req[3]].strip()):
                fail(f"EFFORT_{i}_INVALID_OTHER_CIRSIUM_CENSUS")
        key=tuple(b[k].strip() for k in K)
        if not all(key):fail(f"EFFORT_{i}_MISSING_JOIN_ID")
        if key in index:fail(f"EFFORT_{i}_DUPLICATE_BOUT")
        index[key]=b
        for field in ("video_recording_id","evidence_uri"):
            if not b[field].strip() and b["valid_rate_denominator"]=="1":
                fail(f"EFFORT_{i}_VALID_BUT_MISSING_{field}")
        for field,values in spec["allowed"].items():
            if field in ("binary","optional_binary"):continue
            if b[field] not in values:fail(f"EFFORT_{i}_UNKNOWN_{field}")
        for flag in ("approach_zone_fully_covered","head_zone_fully_covered",
                     "approach_detection_recall_verified","zero_approaches_confirmed",
                     "valid_rate_denominator","reviewer_blinded_to_head_traits"):
            if b[flag] not in ("0","1"):fail(f"EFFORT_{i}_INVALID_{flag}")
        try:
            start=float(b["video_window_start_s"])
            end=float(b["video_window_end_s"])
            sec=float(b["evaluable_duration_s"])
        except (ValueError,TypeError):
            fail(f"EFFORT_{i}_INVALID_VIDEO_TIMES")
        if not(0<=start<end and 0<=sec<=end-start+EPS):
            fail(f"EFFORT_{i}_INVALID_VIDEO_WINDOW")
        b["_start"]=start;b["_end"]=end;b["_sec"]=sec
        valid=(sec>0 and abs(sec-(end-start))<EPS
               and b["approach_zone_fully_covered"]=="1"
               and b["head_zone_fully_covered"]=="1"
               and b["video_review_mode"] in
                  ("full_continuous_manual","continuous_algorithm_audited")
               and (b["video_review_mode"]=="full_continuous_manual"
                    or b["approach_detection_recall_verified"]=="1"))
        if b["valid_rate_denominator"]=="1" and not valid:
            fail(f"EFFORT_{i}_FALSE_VALID_RATE_DENOMINATOR")
        if (b["valid_rate_denominator"]=="1" and
            b["video_review_mode"]=="continuous_algorithm_audited" and
            not b["recall_validation_reference"].strip()):
            fail(f"EFFORT_{i}_AUDITED_RECALL_REFERENCE_MISSING")
        if b["zero_approaches_confirmed"]=="1" and b["valid_rate_denominator"]!="1":
            fail(f"EFFORT_{i}_FALSE_ZERO_FROM_INCOMPLETE_FOOTAGE")
    visits_by_key={key:0 for key in index}
    for i,row in enumerate(events,1):
        key=tuple(str(row.get(k,"")).strip() for k in K)
        if key not in index:
            fail(f"EVENT_{i}_UNMATCHED_EFFORT_BOUT")
        b=index[key]
        try:sec=float(row["time_from_bout_start_s"])
        except (KeyError,ValueError,TypeError):fail(f"EVENT_{i}_INVALID_TIME")
        if not(b["_start"]-EPS<=sec<=b["_end"]+EPS):
            fail(f"EVENT_{i}_OUTSIDE_REVIEWED_EFFORT_WINDOW")
        if (b["head_stage"]!=row.get("phenological_stage","") or
            (b["reproductive_sex_state"]!="not_determined" and
             row.get("reproductive_sex_state") not in
             (b["reproductive_sex_state"],"not_determined"))):
            fail(f"EVENT_{i}_HEAD_STAGE_SEX_CONTRADICTION")
        if row.get("approached")!="1":fail(f"EVENT_{i}_NOT_A_REAL_APPROACH")
        visits_by_key[key]+=1
    for key,b in index.items():
        if b["valid_rate_denominator"]=="1" and visits_by_key[key]==0 and b["zero_approaches_confirmed"]!="1":
            fail("VALID_COMPLETE_EFFORT_WITH_NO_EVENTS_REQUIRES_VERIFIED_ZERO")
        if b["zero_approaches_confirmed"]=="1" and visits_by_key[key]>0:
            fail("FALSE_ZERO_CONFIRMED_WITH_OBSERVED_APPROACH")
    confirmed_minutes=sum(b["_sec"] for b in index.values()
                          if b["valid_rate_denominator"]=="1")/60
    visits=sum(visits_by_key[k] for k,b in index.items()
               if b["valid_rate_denominator"]=="1")
    return {
      "n_measured_other_cirsium_bouts":sum(
          b["other_cirsium_resource_census_status"]=="measured"
          for b in index.values()),
      "n_unassessed_other_cirsium_bouts":sum(
          b["other_cirsium_resource_census_status"]=="not_assessed"
          for b in index.values()),
      "status":"VALIDATED_EVENT_AND_EFFORT_LINK_ONLY" if index
               else "NOT_IDENTIFIABLE_NO_REAL_VIDEO_BOUTS",
      "n_effort_bouts":len(index),"n_events_total":len(events),
      "n_valid_approach_denominator_bouts":sum(b["valid_rate_denominator"]=="1" for b in index.values()),
      "n_verified_zero_approach_bouts":sum(b["zero_approaches_confirmed"]=="1" for b in index.values()),
      "eligible_total_minutes":confirmed_minutes,
      "eligible_approach_episode_count":visits,
      "observed_raw_arrival_rate_per_minute":visits/confirmed_minutes if confirmed_minutes>0 else None,
      "scientific_limit":"Raw observed episode/min only, NOT independent biological replicates, effective pollination, adaptive fitness or treatment effect. Detection bias and stage, plant, guild effects must be modelled separately.",
    }

def make_effort():
    return dict(
        individual_id="TEST_PLANT",population_id="TEST_POP",
        capitulum_id="TEST_HEAD",observation_bout_id="TEST_BOUT_1",
        video_recording_id="synthetic/continuous.mp4",
        head_stage="full_anthesis",reproductive_sex_state="hermaphroditic",
        observed_head_orientation_deg="90",head_rank="terminal",
        video_window_start_s="0",video_window_end_s="120",
        evaluable_duration_s="120",approach_zone_fully_covered="1",
        head_zone_fully_covered="1",video_review_mode="full_continuous_manual",
        approach_detection_recall_verified="0",recall_validation_reference="",
        other_cirsium_resource_census_status="not_assessed",
        other_cirsium_buds_count="",other_cirsium_open_heads_count="",
        other_cirsium_census_radius_m="",other_cirsium_census_reference="",
        zero_approaches_confirmed="1",valid_rate_denominator="1",
        incomplete_reason="",reviewer_blinded_to_head_traits="1",
        evidence_uri="synthetic/continuous.mp4"
    )

def make_event():
    return dict(individual_id="TEST_PLANT",population_id="TEST_POP",
       capitulum_id="TEST_HEAD",observation_bout_id="TEST_BOUT_1",
       time_from_bout_start_s="35",phenological_stage="full_anthesis",
       reproductive_sex_state="hermaphroditic",approached="1")

def tests(spec):
    zero=make_effort()
    v=validate([copy.deepcopy(zero)],[],spec)
    assert v["status"]=="VALIDATED_EVENT_AND_EFFORT_LINK_ONLY"
    assert v["n_verified_zero_approach_bouts"]==1
    assert v["observed_raw_arrival_rate_per_minute"]==0.0
    good=make_effort()
    good["zero_approaches_confirmed"]="0"
    event=make_event()
    x=validate([copy.deepcopy(good)],[event],spec)
    assert x["observed_raw_arrival_rate_per_minute"]==0.5
    congener_measured=make_effort()
    congener_measured.update(
        other_cirsium_resource_census_status="measured",
        other_cirsium_buds_count="0",
        other_cirsium_open_heads_count="0",
        other_cirsium_census_radius_m="5",
        other_cirsium_census_reference="synthetic/census/001",
    )
    q=validate([copy.deepcopy(congener_measured)],[],spec)
    assert q["n_measured_other_cirsium_bouts"]==1
    cases=[
     ("UNMEASURED_CONGENER_IS_NOT_ZERO",
       [{**zero,"other_cirsium_open_heads_count":"0"}],[],
       "UNMEASURED_OTHER_CIRSIUM_MUST_NOT_BE_ZERO_OR_GUESSED"),
     ("MEASURED_CONGENER_MUST_HAVE_RADIUS",
       [{**congener_measured,"other_cirsium_census_radius_m":""}],[],
       "INVALID_OTHER_CIRSIUM_COUNTS"),
     ("UNMATCHED_EFFORT",[],[event],"UNMATCHED_EFFORT_BOUT"),
     ("NO_FALSE_ZERO",[zero],[event],"FALSE_ZERO_CONFIRMED_WITH_OBSERVED_APPROACH"),
     ("MISSING_FULL_HEAD",[{**good,"head_zone_fully_covered":"0"}],[event],"FALSE_VALID_RATE_DENOMINATOR"),
     ("NO_EVENT_TRIGGERED_DENOM",[{
         **good,"video_review_mode":"event_triggered_only"}],[event],"FALSE_VALID_RATE_DENOMINATOR"),
     ("NO_PARTIAL_INTERVAL_RATE",[{
         **good,"evaluable_duration_s":"100"}],[event],"FALSE_VALID_RATE_DENOMINATOR"),
     ("NO_UNREVIEWED_ZERO",[{
         **zero,"valid_rate_denominator":"0"}],[],"FALSE_ZERO_FROM_INCOMPLETE_FOOTAGE"),
     ("NO_MISSING_ZERO",[good],[],"REQUIRES_VERIFIED_ZERO"),
     ("NO_TIME_OUTSIDE_EFFORT",[good],[{
         **event,"time_from_bout_start_s":"125"}],"OUTSIDE_REVIEWED_EFFORT_WINDOW"),
     ("NO_STAGE_MISMATCH",[good],[{
         **event,"phenological_stage":"bud"}],"HEAD_STAGE_SEX_CONTRADICTION"),
     ("NO_UNCALIBRATED_ALGO",[{
         **good,"video_review_mode":"continuous_algorithm_audited",
         "approach_detection_recall_verified":"0"}],[event],"FALSE_VALID_RATE_DENOMINATOR"),
     ("NO_MISSING_ALGO_CALIBRATION",[{
         **good,"video_review_mode":"continuous_algorithm_audited",
         "approach_detection_recall_verified":"1"}],[event],"AUDITED_RECALL_REFERENCE_MISSING")
    ]
    for label,es,vs,error in cases:
        try:validate([copy.deepcopy(i) for i in es],vs,spec)
        except ValueError as e:
            assert error in str(e),(label,str(e))
        else:raise AssertionError(f"TEST_FAILURE_{label}")
    return {"status":"PASS_SYNTHETIC_DENOMINATOR_FIREWALL",
            "positive_control_bouts":3,"negative_cases":len(cases),
            "no_real_organismal_field_data_submitted":True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--contract",type=Path,required=True)
    ap.add_argument("--effort",type=Path,required=True)
    ap.add_argument("--events",type=Path,required=True)
    ap.add_argument("--self-test",action="store_true")
    a=ap.parse_args()
    spec=json.loads(a.contract.read_text(encoding="utf-8"))
    h,rows=read_csv(a.effort)
    assert h==spec["columns"],"EFFORT_HEADER_CONTRACT_DRIFT"
    eh,events=read_csv(a.events)
    assert set(K).issubset(eh)
    result=validate(rows,events,spec)
    if a.self_test:result["synthetic_tests"]=tests(spec)
    print(json.dumps(result,indent=2))

if __name__=="__main__":main()
