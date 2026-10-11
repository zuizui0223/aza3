#!/usr/bin/env python3
"""Summarize head-entry routes only when video efforts are fully enumerated.

Conceptual endpoint: route overlap between *independently pre-identified*
legitimate pollinator candidates and seed-feeder candidates. Same-guild
visits after outcome-selected identification do not count. Overlap may
indicate shared exposure to the same geometry, not actual selected defence.
"""
from __future__ import annotations
import argparse
import collections
import json
from pathlib import Path
from validate_capitulum_video_effort_denominator_v1 import read_csv,validate

GUILDS=("legitimate_pollinator_candidate","seed_feeder_candidate")
ROUTES=("above","side","below","stem_crawl")
MIN_EVENTS=10
MIN_HEADS=3

def overlap(count_a,count_b):
    a=sum(count_a.values());b=sum(count_b.values())
    if not a or not b: return None
    return sum(min(count_a.get(k,0)/a,count_b.get(k,0)/b) for k in ROUTES)

def tabulate(efforts,events,contract):
    joint=validate(efforts,events,contract)
    if not efforts:
        return {"status":"NO_REAL_VIDEO_EVENTS_OR_EFFORT",
                "conditional_route_overlap":None,
                "n_valid_efforts":0,"n_real_events":len(events),
                "claim_boundary":"NOT_IDENTIFIABLE__NO_FIELD_DATA"}
    valid={tuple(b[k] for k in ("individual_id","population_id","capitulum_id","observation_bout_id"))
           for b in efforts if b["valid_rate_denominator"]=="1"}
    by_guild={g:collections.Counter() for g in GUILDS}
    success=collections.defaultdict(lambda:{"all":0,"reached":0,"heads":set(),"unknown_routes":0})
    for e in events:
        key=tuple(e[k] for k in ("individual_id","population_id","capitulum_id","observation_bout_id"))
        if key not in valid:continue
        g=e["pre_entry_guild"]
        if g not in GUILDS:continue
        r=e["entry_route_world"]
        h=(e["individual_id"],e["population_id"],e["capitulum_id"])
        stat=success[g]
        stat["all"]+=1
        stat["heads"].add(h)
        if e["reproductive_zone_reached"]=="1":stat["reached"]+=1
        if r not in ROUTES:stat["unknown_routes"]+=1
        else:by_guild[g][r]+=1
    constraints={}
    for g in GUILDS:
        stat=success[g]
        n=stat["all"]
        constraints[g]={
            "eligible_approaches":n,
            "eligible_heads":len(stat["heads"]),
            "reproductive_zone_reached":stat["reached"],
            "observed_zone_reach_given_approach":stat["reached"]/n if n else None,
            "known_world_route_n":sum(by_guild[g].values()),
            "unknown_world_route_n":stat["unknown_routes"],
            "world_route_counts":{r:by_guild[g][r] for r in ROUTES},
            "valid_conditional_route_estimation":n>=MIN_EVENTS and len(stat["heads"])>=MIN_HEADS and
                stat["unknown_routes"]==0,
        }
    good=all(v["valid_conditional_route_estimation"] for v in constraints.values())
    val=overlap(by_guild[GUILDS[0]],by_guild[GUILDS[1]]) if good else None
    status="DESCRIPTIVE_GUILD_ROUTE_OVERLAP_ONLY" if good else "INSUFFICIENT_OR_UNKNOWN_PRE_ENTRY_GUILD"
    return {
        "status":status,
        "n_effort_bouts":len(efforts),
        "n_valid_rate_bouts":len(valid),
        "total_events_in_all_bouts":len(events),
        "guilds":constraints,
        "conditional_route_overlap":val,
        "definition":"sum_{r in [above,side,below,stem_crawl]} min[P(route r|pre-entry pollinator candidate),P(route r|pre-entry seed-feeder candidate)]",
        "interpretation_if_zero":"Different world-centred approach-route usage, not proof no floral-defence tradeoff; anatomy/contact is unmeasured.",
        "interpretation_if_one":"Same gross world-centred approach routes, not evidence real spine contact or an evolutionary tradeoff.",
        "not_modelled":["stage and sex","species/genotype","biotic density","actual head-facing approach transforms","time-to-egg laying","pollen receipt","parasitoid access","filled seeds"],
        "inferential_limit":"No phylogenetic or causal treatment comparison; specimens/approach events are nested within plant and head. No adaptative origin or selection estimate."
    }

def synthetic_tests():
    a=collections.Counter({"above":5,"side":5})
    b=collections.Counter({"below":8,"stem_crawl":2})
    assert overlap(a,b)==0.0
    assert overlap(a,a)==1.0
    c=collections.Counter({"above":6,"below":4})
    assert abs(overlap(a,c)-.5)<1e-12
    return {"status":"PASS_ROUTE_OVERLAP_ARITHMETIC_ONLY",
            "n_synthetic_tests":3,"not_field_data":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--contract",type=Path,required=True)
    p.add_argument("--effort",type=Path,required=True)
    p.add_argument("--events",type=Path,required=True)
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    spec=json.loads(a.contract.read_text())
    eh,eff=read_csv(a.effort)
    vh,events=read_csv(a.events)
    assert eh==spec["columns"]
    result=tabulate(eff,events,spec)
    if a.self_test:result["synthetic_tests"]=synthetic_tests()
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
