#!/usr/bin/env python3
"""Exact-concurrent-bout head-guild route audit, with zero-safe intake.

Use stage/rank/sex matching PLUS exactly concurrent calendar-time intervals.
Only validated continuous full-space effort and independently identified
pre-entry guilds are eligible. Route scores are descriptive, never fitness.
Current actual sidecar has zero rows; synthetic examples remain tests only.
"""
from __future__ import annotations
import argparse
import collections
import copy
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from validate_capitulum_video_effort_denominator_v1 import K, make_effort, read_csv, validate

SIDE_FIELDS=[
    "individual_id","population_id","capitulum_id","observation_bout_id",
    "video_start_local_iso","video_end_local_iso","recording_evidence_uri",
    "contemporaneous_weather_reference","floral_scent_status",
    "floral_scent_sample_reference",
]
GUILDS=("legitimate_pollinator_candidate","seed_feeder_candidate")
ROUTES=("above","side","below","stem_crawl")
PRE_ENTRY_EVIDENCE={
    "pre_entry_taxon_and_local_life_history_verified",
    "independently_assigned_before_access",
}
MIN_APPROACHES=10
MIN_PLANTS=3
EPS=1e-5

def overlap(a,b):
    na,nb=sum(a.values()),sum(b.values())
    return (sum(min(a[r]/na,b[r]/nb) for r in ROUTES)
            if na>0 and nb>0 else None)

def parse_time(value):
    try:
        d=datetime.fromisoformat(value)
    except (ValueError,TypeError):
        raise ValueError("CLOCK_NOT_ISO8601") from None
    if d.tzinfo is None or d.utcoffset() is None:
        raise ValueError("CLOCK_WITHOUT_EXPLICIT_TIMEZONE")
    return d

def make_side_index(sidecar,efforts):
    known={tuple(e[k] for k in K):e for e in efforts}
    index={}
    for row in sidecar:
        if set(row)!=set(SIDE_FIELDS):
            raise ValueError("CLOCK_SIDECAR_SCHEMA_DRIFT")
        key=tuple(row[k].strip() for k in K)
        if not all(key) or key not in known or key in index:
            raise ValueError("CLOCK_ORPHAN_EMPTY_ID_OR_DUPLICATE")
        if not row["recording_evidence_uri"].strip():
            raise ValueError("CLOCK_SOURCE_VIDEO_REFERENCE_MISSING")
        if row["floral_scent_status"] not in {"not_assessed","measured"}:
            raise ValueError("CLOCK_SCENT_STATUS_INVALID")
        if row["floral_scent_status"]=="measured" and not row["floral_scent_sample_reference"].strip():
            raise ValueError("CLOCK_MEASURED_SCENT_WITHOUT_SAMPLE")
        if row["floral_scent_status"]=="not_assessed" and row["floral_scent_sample_reference"].strip():
            raise ValueError("CLOCK_UNASSESSED_SCENT_HAS_FALSE_SAMPLE")
        if not row["contemporaneous_weather_reference"].strip():
            raise ValueError("CLOCK_WEATHER_CONTEXT_MUST_BE_SOURCE_OR_NOT_ASSESSED")
        start=parse_time(row["video_start_local_iso"])
        end=parse_time(row["video_end_local_iso"])
        if start.utcoffset()!=end.utcoffset():
            raise ValueError("CLOCK_OFFSET_SWITCH_WITHIN_BOUT")
        seconds=(end-start).total_seconds()
        if not 0<seconds<=900+EPS:
            raise ValueError("CLOCK_EXACT_BOUT_MUST_BE_MAX_15_MINUTES")
        e=known[key]
        if e["valid_rate_denominator"]=="1" and abs(float(e["evaluable_duration_s"])-seconds)>EPS:
            raise ValueError("CLOCK_EFFORT_DURATION_DISAGREES_WITH_ABSOLUTE_TIME")
        index[key]={
            "start_utc":start.astimezone(timezone.utc).isoformat(),
            "end_utc":end.astimezone(timezone.utc).isoformat(),
            "scored_scent":row["floral_scent_status"]=="measured",
        }
    return index

def analyze(efforts,events,sidecar,contract):
    # Shared source contracts reject false zeros, incomplete video and head-stage mismatch.
    receipt=validate(efforts,events,contract)
    if not efforts and not sidecar:
        return {
            "status":"NO_REAL_HEAD_VIDEO_OR_CLOCK_EFFORT",
            "n_bouts":0,"n_events":len(events),
            "n_clock_matched_strata":0,
            "pooled_overlap":None,"exact_concurrent_overlap":None,
            "claim_ceiling":"NOT_IDENTIFIABLE",
        }
    index=make_side_index(sidecar,efforts)
    valid={
        tuple(e[k] for k in K):e for e in efforts if e["valid_rate_denominator"]=="1"
    }
    missing=set(valid)-set(index)
    if missing:
        return {
            "status":"HOLD_VALID_HEAD_BOUTS_LACK_ABSOLUTE_CLOCK_SIDECAR",
            "n_bouts":len(efforts),"n_events":len(events),
            "n_valid_bouts_missing_clock":len(missing),
            "pooled_overlap":None,"exact_concurrent_overlap":None,
            "claim_ceiling":"NOT_IDENTIFIABLE",
        }
    strata=collections.defaultdict(lambda:collections.defaultdict(
        lambda:{"n":0,"plants":set(),"heads":set(),
                "routes":collections.Counter(),"unknown":0,
                "taxa":collections.Counter(),"scent_sampled":False,
                "access_known":0,"access_yes":0,"access_unknown":0}
    ))
    pooled={g:collections.Counter() for g in GUILDS}
    unresolved_guild=0
    for ev in events:
        key=tuple(ev[k] for k in K)
        bout=valid.get(key)
        if bout is None:
            continue
        guild=ev.get("pre_entry_guild")
        if guild=="unresolved":
            unresolved_guild+=1
            continue
        if guild not in GUILDS:
            continue
        if (ev.get("pre_entry_guild_evidence") not in PRE_ENTRY_EVIDENCE or
            not ev.get("evidence_uri","").strip() or
            not ev.get("visitor_taxon_label","").strip()):
            raise ValueError("POST_ACCESS_GUILD_OR_MISSING_PRECONTACT_TAXON")
        route=ev.get("entry_route_world","undetermined")
        if route not in set(ROUTES)|{"undetermined"}:
            raise ValueError("INVALID_WORLD_ROUTE")
        source=index[key]
        group=(
            ev["population_id"],bout["head_stage"],bout["head_rank"],
            bout["reproductive_sex_state"],source["start_utc"],source["end_utc"]
        )
        datum=strata[group][guild]
        datum["n"]+=1
        datum["plants"].add((ev["population_id"],ev["individual_id"]))
        datum["heads"].add((ev["population_id"],ev["individual_id"],ev["capitulum_id"]))
        datum["taxa"][ev["visitor_taxon_label"]]+=1
        datum["scent_sampled"] |= source["scored_scent"]
        access=ev.get("reproductive_zone_reached","NA")
        if access not in {"0","1","NA"}:
            raise ValueError("REPRODUCTIVE_ZONE_ACCESS_IS_NOT_ASSESSED_OR_BINARY")
        if access=="NA":
            datum["access_unknown"]+=1
        else:
            datum["access_known"]+=1
            datum["access_yes"]+=int(access=="1")
        if route=="undetermined":
            datum["unknown"]+=1
        else:
            datum["routes"][route]+=1
            pooled[guild][route]+=1

    eligible=[]
    held=[]
    for group, data in sorted(strata.items()):
        failures={}
        for g in GUILDS:
            item=data[g]
            reasons=[]
            if item["n"]<MIN_APPROACHES:reasons.append("TOO_FEW_PRE_ENTRY_APPROACHES")
            if len(item["plants"])<MIN_PLANTS:reasons.append("TOO_FEW_INDEPENDENT_PLANTS")
            if item["unknown"]:reasons.append("UNRESOLVED_WORLD_ROUTES")
            if reasons:failures[g]=reasons
        if failures:
            held.append({
                "stratum":list(group),
                "failure_by_guild":failures,
                "approaches_by_guild":{g:data[g]["n"] for g in GUILDS},
            })
            continue
        eligible.append({
            "stratum":list(group),
            "overlap":overlap(data[GUILDS[0]]["routes"],data[GUILDS[1]]["routes"]),
            "approaches":{g:data[g]["n"] for g in GUILDS},
            "independent_plants":{g:len(data[g]["plants"]) for g in GUILDS},
            "taxon_labels":{g:dict(data[g]["taxa"]) for g in GUILDS},
            "access_given_independently_identified_approach":{
                g:{
                    "known":data[g]["access_known"],
                    "unknown":data[g]["access_unknown"],
                    "zone_reached":data[g]["access_yes"],
                    "observed_proportion_if_all_assessed":(
                        data[g]["access_yes"]/data[g]["n"]
                        if data[g]["access_unknown"]==0 else None
                    ),
                } for g in GUILDS
            },
            "scent_measured_in_any_included_bout":
                any(data[g]["scent_sampled"] for g in GUILDS),
        })
    naive=overlap(pooled[GUILDS[0]],pooled[GUILDS[1]])
    concurrent=(sum(x["overlap"] for x in eligible)/len(eligible)
                if eligible else None)
    return {
        "status": ("EXACT_CLOCK_MATCHED_DESCRIPTIVE_ONLY" if eligible and not held
                   else "PARTIAL_EXACT_CLOCK_MATCHED_DESCRIPTIVE_ONLY" if eligible
                   else "NO_IDENTIFIABLE_CONCURRENT_GUILD_CONTRAST"),
        "n_bouts":len(efforts),
        "n_valid_continuous_bouts":len(valid),
        "n_events":len(events),
        "n_pre_entry_unresolved_guild":unresolved_guild,
        "n_clock_matched_strata":len(eligible),
        "n_held_strata":len(held),
        "pooled_overlap":naive,
        "exact_concurrent_overlap":concurrent,
        "eligible_strata":eligible,
        "held_strata":held,
        "scent_not_controlled_by_clock_matching":True,
        "warning":"Even exact concurrent clock windows do not identify head odour,"
                  " true phyllary-spine anatomy, effective pollination, causal"
                  " entry filtering, evolutionary adaptation, or viable seeds.",
        "effort_validation":receipt["status"],
        "claim_ceiling":"DESCRIPTIVE_HEAD_STAGE_AND_CLOCK_ROUTE_OVERLAP_ONLY",
    }

def synthetic_tests(contract):
    efforts,events,sidecar=[],[],[]
    for suffix,counts,route,start in [
        ("morning",(90,10),"above","2026-06-02T08:00:00+09:00"),
        ("afternoon",(10,90),"below","2026-06-02T15:00:00+09:00"),
    ]:
        start_dt=parse_time(start)
        from datetime import timedelta
        end=(start_dt+timedelta(minutes=15)).isoformat()
        for p in range(3):
            b=make_effort()
            b.update(individual_id=f"{suffix}_p{p}",
                     capitulum_id=f"{suffix}_h{p}",
                     observation_bout_id=f"{suffix}_b{p}",
                     video_recording_id=f"synthetic://{suffix}/{p}",
                     head_stage="full_anthesis",head_rank="terminal",
                     video_window_end_s="900",evaluable_duration_s="900",
                     zero_approaches_confirmed="0",
                     evidence_uri=f"synthetic://{suffix}/{p}")
            efforts.append(b)
            sidecar.append(dict(zip(SIDE_FIELDS,[
                b["individual_id"],b["population_id"],b["capitulum_id"],
                b["observation_bout_id"],start,end,
                f"synthetic://{suffix}/{p}","not_assessed",
                "not_assessed",""
            ])))
        for guild,n in zip(GUILDS,counts):
            for i in range(n):
                p=i%3
                events.append({
                    "population_id":"TEST_POP",
                    "individual_id":f"{suffix}_p{p}",
                    "capitulum_id":f"{suffix}_h{p}",
                    "observation_bout_id":f"{suffix}_b{p}",
                    "approach_episode_id":f"{suffix}_{guild}_{i}",
                    "time_from_bout_start_s":str(i*8),
                    "phenological_stage":"full_anthesis",
                    "reproductive_sex_state":"hermaphroditic",
                    "approached":"1",
                    "pre_entry_guild":guild,
                    "pre_entry_guild_evidence":"independently_assigned_before_access",
                    "visitor_taxon_label":f"synthetic:{guild}",
                    "evidence_uri":"synthetic://clip",
                    "entry_route_world":route,
                    "reproductive_zone_reached":"1",
                })
    def check(x,y,z):
        return analyze(copy.deepcopy(x),copy.deepcopy(y),copy.deepcopy(z),contract)
    out=check(efforts,events,sidecar)
    assert out["status"]=="EXACT_CLOCK_MATCHED_DESCRIPTIVE_ONLY"
    assert abs(out["pooled_overlap"]-.2)<1e-12
    assert abs(out["exact_concurrent_overlap"]-1)<1e-12
    assert out["n_clock_matched_strata"]==2
    assert analyze([],[],[],contract)["exact_concurrent_overlap"] is None
    missing=check(efforts,events,sidecar[:-1])
    assert missing["status"]=="HOLD_VALID_HEAD_BOUTS_LACK_ABSOLUTE_CLOCK_SIDECAR"
    for bad,expected in [
        ({**sidecar[0],"video_start_local_iso":"2026-06-02T08:00:00"},"CLOCK_WITHOUT_EXPLICIT_TIMEZONE"),
        ({**sidecar[0],"floral_scent_status":"measured"},"CLOCK_MEASURED_SCENT_WITHOUT_SAMPLE"),
        ({**sidecar[0],"video_end_local_iso":"2026-06-02T09:15:00+09:00"},"CLOCK_EXACT_BOUT_MUST_BE_MAX_15_MINUTES"),
    ]:
        altered=copy.deepcopy(sidecar)
        altered[0]=bad
        try:check(efforts,events,altered)
        except ValueError as e:
            assert expected in str(e),(str(e),expected)
        else:raise AssertionError("Invalid clock data accepted")
    invalid=copy.deepcopy(events)
    invalid[0]["pre_entry_guild_evidence"]="unresolved"
    try:check(efforts,invalid,sidecar)
    except ValueError as e:
        assert "POST_ACCESS_GUILD" in str(e)
    else:raise AssertionError("Outcome-based guild identification accepted")
    # No valid contrast if the same interval has only one pre-entry guild.
    only=[e for e in events if e["pre_entry_guild"]==GUILDS[0]]
    one=check(efforts,only,sidecar)
    assert one["exact_concurrent_overlap"] is None
    return {
        "status":"PASS_SYNTHETIC_EXACT_CONCURRENT_CLOCK_MATCHING",
        "naive_same_stage_overlap":out["pooled_overlap"],
        "within_exact_bout_window_overlap":out["exact_concurrent_overlap"],
        "matched_windows":out["n_clock_matched_strata"],
        "negative_guards":5,
        "real_observation_rows":0,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--contract",type=Path,default=Path("data/contracts/capitulum_video_effort_denominator_v1.json"))
    ap.add_argument("--effort",type=Path,default=Path("data/intake/capitulum_video_effort_denominator_v1.csv"))
    ap.add_argument("--events",type=Path,default=Path("data/intake/capitulum_guild_access_event_ledger_v1.csv"))
    ap.add_argument("--sidecar",type=Path,default=Path("data/intake/capitulum_clock_matched_effort_sidecar_v1.csv"))
    ap.add_argument("--out",type=Path)
    a=ap.parse_args()
    spec=json.loads(a.contract.read_text(encoding="utf-8"))
    eh,eff=read_csv(a.effort)
    vh,events=read_csv(a.events)
    sh,side=read_csv(a.sidecar)
    if eh!=spec["columns"] or len(vh)!=32 or sh!=SIDE_FIELDS:
        raise ValueError("CLOCK_INPUT_HEADERS_DRIFT_FROM_V1_CONTRACT")
    result={
        "version":"capitulum_clock_matched_guild_access_v1",
        "date":"2026-10-08",
        "synthetic_tests":synthetic_tests(spec),
        "real_data":analyze(eff,events,side,spec),
    }
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                         encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
