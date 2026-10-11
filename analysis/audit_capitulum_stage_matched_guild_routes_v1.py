#!/usr/bin/env python3
"""Phenophase-matched Cirsium guild routes; guard against pooled Simpson artifacts.

Rates and overlap are only descriptive for independently pre-identified
pollinator and seed-feeder candidates with valid CONTINUOUS video effort.
Same-site, same-head-stage, same-head-rank matching does not ensure matching
calendar dates/weather, and no causal/evolutionary inference is licensed.
"""
from __future__ import annotations

import argparse
import collections
import copy
import json
from pathlib import Path
from validate_capitulum_video_effort_denominator_v1 import (
    K, make_effort, read_csv, validate,
)

GUILDS=("legitimate_pollinator_candidate", "seed_feeder_candidate")
ROUTES=("above", "side", "below", "stem_crawl")
MIN_APPROACHES=10
MIN_INDEPENDENT_PLANTS=3


def overlap(a, b):
    ta,tb=sum(a.values()),sum(b.values())
    if ta==0 or tb==0:
        return None
    return sum(min(a[x]/ta,b[x]/tb) for x in ROUTES)


def analyze(efforts, events, effort_contract):
    # Delegate the exact verified-denominator/zero/matching/time/stage gates
    # to the existing shared validator, before any route calculations.
    effort_validation=validate(efforts,events,effort_contract)
    if not efforts:
        return {
            "status":"NO_REAL_CONTINUOUS_BOUTS",
            "n_real_bouts":0,"n_events":len(events),
            "naive_pooled_overlap":None,
            "stage_matched_equal_stratum_overlap":None,
            "ecological_claim":"NOT_IDENTIFIABLE",
        }

    valid={
        tuple(b[k] for k in K):b
        for b in efforts
        if b["valid_rate_denominator"]=="1"
    }
    strata=collections.defaultdict(
        lambda:collections.defaultdict(
            lambda:{"routes":collections.Counter(),
                    "plants":set(),"heads":set(),"n":0,"unknown_routes":0}
        )
    )
    raw={g:collections.Counter() for g in GUILDS}
    raw_events={g:0 for g in GUILDS}
    unresolved=0
    for e in events:
        key=tuple(str(e.get(k,"")).strip() for k in K)
        bout=valid.get(key)
        if bout is None:
            continue
        g=e.get("pre_entry_guild","unresolved")
        if g=="unresolved":
            unresolved+=1
            continue
        if g not in GUILDS:
            continue
        if (e.get("pre_entry_guild_evidence") not in
            {"pre_entry_taxon_and_local_life_history_verified",
             "independently_assigned_before_access"} or
            not e.get("evidence_uri","").strip()):
            raise ValueError("POST_ACCESS_ROLE_IS_NOT_PRE_ENTRY_GUILD")
        route=e.get("entry_route_world","undetermined")
        if route not in set(ROUTES)|{"undetermined"}:
            raise ValueError("INVALID_APPROACH_WORLD_ROUTE")
        if e.get("phenological_stage")!=bout["head_stage"]:
            raise ValueError("EVENT_STAGE_MUST_MATCH_EFFORT")
        # These are precise matched *observation strata*, not time-matched
        # calendar windows; full date/time is not in current data contract.
        stratum=(e["population_id"],bout["head_stage"],bout["head_rank"])
        v=strata[stratum][g]
        v["n"]+=1
        raw_events[g]+=1
        v["plants"].add((e["population_id"],e["individual_id"]))
        v["heads"].add((e["population_id"],e["individual_id"],e["capitulum_id"]))
        if route=="undetermined":
            v["unknown_routes"]+=1
        else:
            v["routes"][route]+=1
            raw[g][route]+=1

    matched=[]
    holds=[]
    for stratum, gs in sorted(strata.items()):
        failures={}
        for g in GUILDS:
            x=gs[g]
            why=[]
            if x["n"]<MIN_APPROACHES: why.append("FEWER_THAN_10_APPROACHES")
            if len(x["plants"])<MIN_INDEPENDENT_PLANTS:
                why.append("FEWER_THAN_3_INDEPENDENT_PLANTS")
            if x["unknown_routes"]>0:
                why.append("UNKNOWN_WORLD_ROUTE_NOT_A_ZERO")
            if why:
                failures[g]=why
        if failures:
            holds.append({
                "stratum":list(stratum),
                "reason_by_guild":failures,
                "pre_entry_approaches":{g:gs[g]["n"] for g in GUILDS},
            })
            continue
        score=overlap(gs[GUILDS[0]]["routes"],gs[GUILDS[1]]["routes"])
        matched.append({
            "stratum":list(stratum),
            "n_approaches":{g:gs[g]["n"] for g in GUILDS},
            "n_independent_plants":{g:len(gs[g]["plants"]) for g in GUILDS},
            "n_heads":{g:len(gs[g]["heads"]) for g in GUILDS},
            "route_counts":{g:{r:gs[g]["routes"][r] for r in ROUTES}
                            for g in GUILDS},
            "within_stratum_overlap":score,
        })

    naive=overlap(raw[GUILDS[0]],raw[GUILDS[1]])
    mean=sum(v["within_stratum_overlap"] for v in matched)/len(matched) if matched else None
    coverage_complete=not holds and unresolved==0
    status=(
        "STAGE_MATCHED_DESCRIPTIVE_ONLY"
        if matched and coverage_complete
        else "PARTIAL_STAGE_COVERAGE_DESCRIPTIVE_ONLY"
        if matched
        else "NO_IDENTIFIABLE_SAME_STAGE_GUILD_CONTRAST"
    )
    return {
        "status":status,
        "n_real_bouts":len(efforts),
        "n_valid_continuous_bouts":len(valid),
        "n_events":len(events),
        "pre_entry_classification_unresolved_events":unresolved,
        "n_matched_population_stage_rank_strata":len(matched),
        "n_unmatched_or_low_coverage_strata":len(holds),
        "naive_pooled_overlap":naive,
        "naive_pooled_is_not_an_adaptive_interface_test":True,
        "stage_matched_equal_stratum_overlap":mean,
        "matched_strata":matched,
        "held_strata":holds,
        "denominator_status":effort_validation["status"],
        "interpretation":"Route overlap measures coarse approach direction only; "
            "there is no anatomical spine-contact, causal defence, pollen "
            "delivery, viable achene or phylogenetic adaptation effect.",
        "unresolved_calendar_confounding":
            "Same population/stage/rank need not share date, time, weather, "
            "host-race or insect pressure; prospective exact-bout-time matching required.",
    }


def synthetic_tests(contract):
    # Pure Simpson-style compositional counterexample:
    # EACH stage's pollinator and seed-feeder route distribution is IDENTICAL.
    # Guilds occur at very different relative frequencies across stages.
    efforts=[]
    events=[]
    for stage, poll_n, feeder_n, route in (
        ("bud",90,10,"above"),
        ("full_anthesis",10,90,"below"),
    ):
        for plant in range(3):
            b=make_effort()
            suffix=f"{stage}_{plant}"
            b.update(individual_id=f"P_{suffix}",
                     capitulum_id=f"H_{suffix}",
                     observation_bout_id=f"B_{suffix}",
                     video_recording_id=f"synthetic://video/{suffix}",
                     evidence_uri=f"synthetic://video/{suffix}",
                     head_stage=stage, head_rank="terminal",
                     zero_approaches_confirmed="0",
                     video_window_end_s="3600",
                     evaluable_duration_s="3600")
            efforts.append(b)
        for guild,n in ((GUILDS[0],poll_n),(GUILDS[1],feeder_n)):
            for i in range(n):
                plant=i%3
                suffix=f"{stage}_{plant}"
                events.append({
                    "individual_id":f"P_{suffix}",
                    "population_id":"TEST_POP",
                    "capitulum_id":f"H_{suffix}",
                    "observation_bout_id":f"B_{suffix}",
                    "approach_episode_id":f"E_{stage}_{guild}_{i}",
                    "time_from_bout_start_s":str((i+1)*10),
                    "phenological_stage":stage,
                    "reproductive_sex_state":"hermaphroditic",
                    "approached":"1",
                    "pre_entry_guild":guild,
                    "pre_entry_guild_evidence":"independently_assigned_before_access",
                    "evidence_uri":f"synthetic://video/{suffix}/{i}",
                    "entry_route_world":route,
                })
    r=analyze(copy.deepcopy(efforts),copy.deepcopy(events),contract)
    assert r["status"]=="STAGE_MATCHED_DESCRIPTIVE_ONLY"
    assert r["n_matched_population_stage_rank_strata"]==2
    assert abs(r["naive_pooled_overlap"]-.2)<1e-12
    assert abs(r["stage_matched_equal_stratum_overlap"]-1.0)<1e-12
    assert analyze([],[],contract)["stage_matched_equal_stratum_overlap"] is None

    altered=copy.deepcopy(events)
    # Route unknown invalidates the affected stratum, never silently
    # counted as missing guild and never produces artificial separation.
    altered[0]["entry_route_world"]="undetermined"
    x=analyze(copy.deepcopy(efforts),altered,contract)
    assert x["n_matched_population_stage_rank_strata"]==1
    assert x["n_unmatched_or_low_coverage_strata"]==1
    assert x["status"]=="PARTIAL_STAGE_COVERAGE_DESCRIPTIVE_ONLY"
    only_one_guild=[e for e in events if e["pre_entry_guild"]==GUILDS[0]]
    h=analyze(copy.deepcopy(efforts),only_one_guild,contract)
    assert h["status"]=="NO_IDENTIFIABLE_SAME_STAGE_GUILD_CONTRAST"
    # If pre-entry label was assigned after outcome, fail closed.
    invalid=copy.deepcopy(events)
    invalid[0]["pre_entry_guild_evidence"]="unresolved"
    try:
        analyze(copy.deepcopy(efforts),invalid,contract)
    except ValueError as e:
        assert "POST_ACCESS_ROLE" in str(e)
    else:
        raise AssertionError("Post-outcome guild selection was accepted")
    return {
        "status":"PASS_SYNTHETIC_STAGE_SIMPS0N_FIREWALL",
        "example_pooled_overlap":r["naive_pooled_overlap"],
        "example_within_stage_overlap":r["stage_matched_equal_stratum_overlap"],
        "negative_tests":["unknown_route","only_one_guild",
                          "post_access_guild_label","zero_real_effort"],
        "real_insect_data":0,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--contract",type=Path,
                    default=Path("data/contracts/capitulum_video_effort_denominator_v1.json"))
    ap.add_argument("--effort",type=Path,
                    default=Path("data/intake/capitulum_video_effort_denominator_v1.csv"))
    ap.add_argument("--events",type=Path,
                    default=Path("data/intake/capitulum_guild_access_event_ledger_v1.csv"))
    ap.add_argument("--out",type=Path)
    a=ap.parse_args()
    contract=json.loads(a.contract.read_text(encoding="utf-8"))
    eh,efforts=read_csv(a.effort)
    vh,events=read_csv(a.events)
    if eh!=contract["columns"]:
        raise ValueError("EFFORT_HEADER_CONTRACT_MISMATCH")
    if len(vh)!=32:
        raise ValueError("EVENT_HEADER_CONTRACT_MISMATCH")
    result={
        "version":"capitulum_stage_matched_route_audit_v1",
        "date":"2026-10-08",
        "synthetic_tests":synthetic_tests(contract),
        "real_data":analyze(efforts,events,contract),
    }
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                         encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
