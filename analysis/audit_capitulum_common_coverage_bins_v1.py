#!/usr/bin/env python3
"""Five-minute common-coverage bins for Cirsium capitulum access.

Repair to exact-start/end analysis: compare two complete recordings during
the SAME predeclared UTC five-minute bin even if video bout boundaries differ.
Every included head must have full continuous, already audited footage over
the ENTIRE bin; partial, missing and event-triggered windows cannot certify
zero activity or exposure. Candidate guild is fixed BEFORE access outcome.
Observational route overlap is NOT adaptive structural access or plant fitness.
"""
from __future__ import annotations
import argparse
import collections
import copy
import json
from datetime import datetime, timezone, timedelta
from pathlib import Path

from audit_capitulum_concurrent_clock_guilds_v1 import (
    SIDE_FIELDS, GUILDS, ROUTES, PRE_ENTRY_EVIDENCE,
    make_side_index, crosswalk_original_bouts, overlap
)
from validate_capitulum_video_effort_denominator_v1 import (
    K, make_effort, read_csv, validate,
)

BIN_SECONDS=300
MIN_PER_GUILD=10
MIN_PLANTS_PER_GUILD=3
MIN_DUAL_GUILD_HEADS=3
EPS=1e-6

def bin_start(dt):
    if dt.tzinfo is None:
        raise ValueError("UNZONED_EVENT_TIME")
    epoch=int(dt.timestamp())
    return datetime.fromtimestamp((epoch//BIN_SECONDS)*BIN_SECONDS,
                                  tz=timezone.utc)

def interval_bins(start,end):
    """Only bins fully covered by a continuous validated video interval."""
    cursor=bin_start(start)
    if cursor<start:
        cursor+=timedelta(seconds=BIN_SECONDS)
    ans=[]
    while cursor+timedelta(seconds=BIN_SECONDS)<=end:
        ans.append(cursor)
        cursor+=timedelta(seconds=BIN_SECONDS)
    return ans

def audit(efforts, events, sidecar, contract):
    receipt=validate(efforts,events,contract)
    if not efforts and not sidecar:
        return {
            "status":"NO_REAL_VALIDATED_CLOCK_VIDEO",
            "n_bouts":0,"n_events":len(events),
            "n_shared_covered_bins":0,"n_qualified_comparison_bins":0,
            "naive_pooled_route_overlap":None,
            "matched_binned_route_overlap":None,
            "claim":"NOT_IDENTIFIABLE",
        }
    details=make_side_index(sidecar,efforts)
    valid={tuple(x[k] for k in K):x for x in efforts
           if x["valid_rate_denominator"]=="1"}
    if set(valid)-set(details):
        return {
            "status":"HOLD_MISSING_CLOCK_FOR_VALID_VIDEO_BOUT",
            "n_missing_bouts":len(set(valid)-set(details)),
            "n_shared_covered_bins":0,
            "n_qualified_comparison_bins":0,
            "naive_pooled_route_overlap":None,
            "matched_binned_route_overlap":None,
            "claim":"NOT_IDENTIFIABLE",
        }

    coverage=collections.defaultdict(set)
    # grouping is population × head stage × rank × sex × UTC 5min bin
    group_key={}
    for key,b in valid.items():
        src=details[key]
        st=datetime.fromisoformat(src["start_utc"])
        en=datetime.fromisoformat(src["end_utc"])
        g=(b["population_id"],b["head_stage"],b["head_rank"],
           b["reproductive_sex_state"])
        group_key[key]=g
        for bin_utc in interval_bins(st,en):
            coverage[(g,bin_utc)].add(key)

    # A second overlapping bout/camera on the SAME biological head cannot
    # be treated as another replicate. Without independently reconciled insect
    # identities it could also double-count a single arrival.
    for (_group, _time_bin), covered_keys in coverage.items():
        seen_heads=set()
        for key in covered_keys:
            head=(key[1],key[0],key[2])
            if head in seen_heads:
                raise ValueError("DUPLICATE_HEAD_COVERAGE_BOUTS_NEED_DEDUPLICATION")
            seen_heads.add(head)

    stats=collections.defaultdict(lambda:collections.defaultdict(lambda:{
        "n":0, "route_counts":collections.Counter(), "unknown":0,
        "plants":set(), "heads":set(), "access_yes":0, "access_known":0,
        "access_unknown":0, "taxa":collections.Counter(),
        "per_head_routes":collections.defaultdict(collections.Counter),
        "per_head_unknown":collections.Counter(),
    }))
    pooled={g:collections.Counter() for g in GUILDS}
    n_excluded_in_partial_bin=0
    unresolved_role=0

    for ev in events:
        key=tuple(ev.get(k,"") for k in K)
        bout=valid.get(key)
        if bout is None:
            continue
        src=details[key]
        st=datetime.fromisoformat(src["start_utc"])
        en=datetime.fromisoformat(src["end_utc"])
        try:
            relative=float(ev["time_from_bout_start_s"])-float(bout["video_window_start_s"])
        except (KeyError,ValueError,TypeError):
            raise ValueError("EVENT_TIME_CANNOT_MAP_TO_ABSOLUTE_CLOCK") from None
        if not 0<=relative<(en-st).total_seconds():
            raise ValueError("EVENT_OUTSIDE_VALID_SIDE_CAR_INTERVAL")
        when=st+timedelta(seconds=relative)
        bin_utc=bin_start(when)
        stratum=(group_key[key],bin_utc)
        if key not in coverage.get(stratum,set()):
            n_excluded_in_partial_bin+=1
            continue
        guild=ev.get("pre_entry_guild","unresolved")
        if guild=="unresolved":
            unresolved_role+=1
            continue
        if guild not in GUILDS:
            continue
        if (ev.get("pre_entry_guild_evidence") not in PRE_ENTRY_EVIDENCE or
            not ev.get("evidence_uri","").strip() or
            not ev.get("visitor_taxon_label","").strip()):
            raise ValueError("POST_ACCESS_GUILD_OR_MISSING_TAXON")
        route=ev.get("entry_route_world")
        if route not in set(ROUTES)|{"undetermined"}:
            raise ValueError("INVALID_WORLD_ROUTE")
        access=ev.get("reproductive_zone_reached","NA")
        if access not in {"0","1","NA"}:
            raise ValueError("INVALID_ACCESS_OR_UNKNOWN_STATE")
        datum=stats[stratum][guild]
        datum["n"]+=1
        plant=(bout["population_id"],bout["individual_id"])
        head=(bout["population_id"],bout["individual_id"],bout["capitulum_id"])
        datum["plants"].add(plant)
        datum["heads"].add(head)
        datum["taxa"][ev["visitor_taxon_label"]]+=1
        if route=="undetermined":
            datum["unknown"]+=1
            datum["per_head_unknown"][head]+=1
        else:
            datum["route_counts"][route]+=1
            datum["per_head_routes"][head][route]+=1
            pooled[guild][route]+=1
        if access=="NA":
            datum["access_unknown"]+=1
        else:
            datum["access_known"]+=1
            datum["access_yes"]+=int(access=="1")

    eligible=[]
    held=[]
    comparable_cover=0
    bins_with_no_guild=0
    bins_with_one_guild=0
    bins_with_both_guilds=0
    for (g, bin_utc), covered_heads in sorted(coverage.items()):
        if len({(x[1],x[0],x[2]) for x in covered_heads})<MIN_DUAL_GUILD_HEADS:
            continue
        comparable_cover+=1
        counts=stats[(g,bin_utc)]
        guilds_present=sum(counts[z]["n"]>0 for z in GUILDS)
        if guilds_present==0:
            bins_with_no_guild+=1
        elif guilds_present==1:
            bins_with_one_guild+=1
        else:
            bins_with_both_guilds+=1
        failures=[]
        for guild in GUILDS:
            x=counts[guild]
            if x["n"]<MIN_PER_GUILD:
                failures.append(guild+":INSUFFICIENT_APPROACHES")
            if len(x["plants"])<MIN_PLANTS_PER_GUILD:
                failures.append(guild+":INSUFFICIENT_PLANTS")
            if x["unknown"]>0:
                failures.append(guild+":UNKNOWN_WORLD_ROUTE")
        dual=counts[GUILDS[0]]["heads"] & counts[GUILDS[1]]["heads"]
        if len(dual)<MIN_DUAL_GUILD_HEADS:
            failures.append("FEWER_THAN_3_SAME_HEAD_DUAL_GUILD_EXPOSURES")
        bin_info={
            "population_stage_rank_sex":list(g),
            "bin_start_utc":bin_utc.isoformat(),
            "duration_s":BIN_SECONDS,
            "n_heads_full_video":len(covered_heads),
            "n_same_head_both_guilds":len(dual),
            "n_candidates":{guild:counts[guild]["n"] for guild in GUILDS},
        }
        if failures:
            held.append({**bin_info,"reasons":failures})
            continue
        value=overlap(counts[GUILDS[0]]["route_counts"],
                      counts[GUILDS[1]]["route_counts"])
        head_scores=[]
        for head in sorted(dual):
            a=counts[GUILDS[0]]["per_head_routes"][head]
            b=counts[GUILDS[1]]["per_head_routes"][head]
            n_a=sum(a.values())
            n_b=sum(b.values())
            if n_a==0 or n_b==0:
                raise ValueError("SAME_HEAD_ROUTES_MISSING_AFTER_GUILD_MATCH")
            head_scores.append({
                "head_id":list(head),
                "route_overlap":overlap(a,b),
                "approaches_by_guild":{GUILDS[0]:n_a,GUILDS[1]:n_b},
            })
        equal_head_mean=sum(x["route_overlap"] for x in head_scores)/len(head_scores)
        eligible.append({
            **bin_info, "route_overlap":value,
            "mean_equal_head_weight_route_overlap":equal_head_mean,
            "head_conditioned_route_overlap":head_scores,
            "n_independent_plants":{
                z:len(counts[z]["plants"]) for z in GUILDS
            },
            "access_given_approach":{
                z:{
                    "yes":counts[z]["access_yes"],
                    "known":counts[z]["access_known"],
                    "unknown":counts[z]["access_unknown"],
                    "proportion_when_fully_assessed":(
                        counts[z]["access_yes"]/counts[z]["n"]
                        if counts[z]["access_unknown"]==0 else None
                    )
                } for z in GUILDS
            },
            "taxa":{z:dict(counts[z]["taxa"]) for z in GUILDS},
        })
    naive=overlap(pooled[GUILDS[0]],pooled[GUILDS[1]])
    matched=(sum(x["route_overlap"] for x in eligible)/len(eligible)
             if eligible else None)
    head_matched=(
        sum(x["mean_equal_head_weight_route_overlap"] for x in eligible)/len(eligible)
        if eligible else None
    )
    return {
        "status":("COMPLETE_BINNED_DESCRIPTIVE_OVERLAP_ONLY"
                  if eligible and not held else
                  "PARTIAL_BINNED_DESCRIPTIVE_OVERLAP_ONLY"
                  if eligible else
                  "NO_COMPARABLE_SAME_TIME_AND_HEAD_GUILDS"),
        "n_bouts":len(efforts),
        "n_valid_full_coverage_bouts":len(valid),
        "n_events":len(events),
        "n_shared_covered_bins":comparable_cover,
        "n_covered_bins_with_neither_identified_guild":bins_with_no_guild,
        "n_covered_bins_with_only_one_identified_guild":bins_with_one_guild,
        "n_covered_bins_with_both_identified_guilds":bins_with_both_guilds,
        "n_qualified_comparison_bins":len(eligible),
        "n_held_bins":len(held),
        "n_approach_events_in_partial_time_bins_excluded":n_excluded_in_partial_bin,
        "n_unresolved_pre_entry_guild_in_covered_bins":unresolved_role,
        "naive_pooled_route_overlap":naive,
        "matched_binned_route_overlap":matched,
        "matched_same_head_equal_weight_route_overlap":head_matched,
        "qualified_bins":eligible,
        "held_bins":held,
        "denominator_validation":receipt["status"],
        "claim":"DESCRIPTIVE_ONLY_NO_ANATOMICAL_OR_CAUSAL_FITNESS_EFFECT",
        "important_limitations":[
            "Same full-coverage bin controls clock opportunity only, not flower scent.",
            "World approach direction does not establish contact with a genuine spine.",
            "Repeated same-head arrivals are not independent trait selection replicates.",
            "Repeated approach episodes do not identify unique insect individuals; repeated visits by one insect cannot be treated as independent insects.",
            "A duplicated concurrent video bout for the same biological head is rejected to prevent double-counted exposure.",
            "Pooling insects across differently preferred heads can create an apparent route contrast even when routes match within every head.",
            "Equal-head-weighted overlap conditions on both guilds visiting the same head; this does not estimate the unconditional head population effect.",
            "Unknown guilds or incomplete video are not zeros.",
            "Restricting to dual-guild encounter bins conditions on arrivals; "
            "these selected bins do NOT represent the whole head population.",
            "Between-site, host-race, species and weather confounding remains.",
            "Without real viable achenes no plant fitness or adaptation is measurable.",
        ],
    }


def synthetic_tests(contract):
    from datetime import timedelta
    a0=datetime.fromisoformat("2026-06-02T09:00:00+09:00")
    efforts,sidecars,events=[],[],[]
    # Different video boundaries; overlap only 09:05-09:15.
    # Inside each 5min overlap period both guilds have IDENTICAL route.
    for i in range(6):
        start=a0+timedelta(minutes=5 if i>=3 else 0)
        end=start+timedelta(minutes=15)
        b=make_effort()
        b.update(
            individual_id=f"P{i}",capitulum_id=f"H{i}",
            observation_bout_id=f"B{i}",
            video_recording_id=f"synthetic://cam/{i}",
            evidence_uri=f"synthetic://cam/{i}",
            video_window_end_s="900",evaluable_duration_s="900",
            zero_approaches_confirmed="0",
            head_stage="full_anthesis",head_rank="terminal",
        )
        efforts.append(b)
        sidecars.append(dict(zip(SIDE_FIELDS,[
            b["individual_id"],b["population_id"],
            b["capitulum_id"],b["observation_bout_id"],
            start.isoformat(),end.isoformat(),"synthetic://cam",
            "not_assessed","not_assessed","",
        ])))
    for pos,(when,poll,seed,route) in enumerate([
        (a0+timedelta(minutes=5),90,10,"above"),
        (a0+timedelta(minutes=10),10,90,"below"),
    ]):
        for guild,n in zip(GUILDS,(poll,seed)):
            for j in range(n):
                i=j%6
                start=a0+timedelta(minutes=5 if i>=3 else 0)
                instant=when+timedelta(seconds=1+(j%90)*2)
                rel=(instant-start).total_seconds()
                events.append({
                    "population_id":"TEST_POP","individual_id":f"P{i}",
                    "capitulum_id":f"H{i}","observation_bout_id":f"B{i}",
                    "approach_episode_id":f"E_{pos}_{guild}_{j}",
                    "time_from_bout_start_s":str(rel),
                    "phenological_stage":"full_anthesis",
                    "reproductive_sex_state":"hermaphroditic",
                    "approached":"1", "pre_entry_guild":guild,
                    "pre_entry_guild_evidence":"independently_assigned_before_access",
                    "visitor_taxon_label":f"synthetic:{guild}",
                    "evidence_uri":"synthetic://episode",
                    "entry_route_world":route,
                    "reproductive_zone_reached":"1",
                })
    def check(es=efforts,vs=events,ss=sidecars):
        return audit(copy.deepcopy(es),copy.deepcopy(vs),copy.deepcopy(ss),
                     contract)
    result=check()
    assert result["status"]=="PARTIAL_BINNED_DESCRIPTIVE_OVERLAP_ONLY"
    assert result["n_qualified_comparison_bins"]==2
    assert result["n_held_bins"]==2
    assert result["n_covered_bins_with_neither_identified_guild"]==2
    assert result["n_covered_bins_with_only_one_identified_guild"]==0
    assert result["n_covered_bins_with_both_identified_guilds"]==2
    assert abs(result["naive_pooled_route_overlap"]-.2)<EPS
    assert abs(result["matched_binned_route_overlap"]-1.0)<EPS
    assert abs(result["matched_same_head_equal_weight_route_overlap"]-1.0)<EPS
    assert result["n_approach_events_in_partial_time_bins_excluded"]==0
    assert audit([],[],[],contract)["status"]=="NO_REAL_VALIDATED_CLOCK_VIDEO"
    partial=check(es=efforts,ss=sidecars[:-1])
    assert partial["status"]=="HOLD_MISSING_CLOCK_FOR_VALID_VIDEO_BOUT"
    # Duplicated concurrent cameras on one biological head cannot double
    # the independent head count or the approach exposure.
    double=copy.deepcopy(efforts[0])
    double["observation_bout_id"]="B_EXTRA_CAMERA"
    double["video_recording_id"]="synthetic://extra_cam"
    double["evidence_uri"]="synthetic://extra_cam"
    double["zero_approaches_confirmed"]="1"
    side_duplicate=copy.deepcopy(sidecars[0])
    side_duplicate["observation_bout_id"]="B_EXTRA_CAMERA"
    side_duplicate["recording_evidence_uri"]="synthetic://extra_cam"
    try:
        check(es=efforts+[double],ss=sidecars+[side_duplicate])
    except ValueError as e:
        assert "DUPLICATE_HEAD_COVERAGE_BOUTS_NEED_DEDUPLICATION" in str(e)
    else:
        raise AssertionError("Duplicate camera-head accepted as independent")
    # Moving arrival outside the evaluated sidecar cannot enter denominator.
    invalid=copy.deepcopy(events)
    invalid[0]["time_from_bout_start_s"]="900"
    try:check(vs=invalid)
    except ValueError as e:
        assert "EVENT_OUTSIDE_VALID_SIDE_CAR_INTERVAL" in str(e)
    else:raise AssertionError("Event past end accepted")
    unknown=copy.deepcopy(events)
    unknown[0]["entry_route_world"]="undetermined"
    u=check(vs=unknown)
    assert u["n_qualified_comparison_bins"]==1
    # A bare event-triggered video can never be a continuous denominator.
    bad_eff=copy.deepcopy(efforts)
    bad_eff[0]["video_review_mode"]="event_triggered_only"
    try:check(es=bad_eff)
    except ValueError as e:
        assert "FALSE_VALID_RATE_DENOMINATOR" in str(e)
    else:raise AssertionError("Event-triggered video accepted")
    # Second independent falsification: SAME 5-minute clock interval,
    # identical insect-guild routes WITHIN each head, but both guilds
    # prefer different heads. Pooling across heads fakes route divergence.
    eh,sh,ev=[],[],[]
    start=a0+timedelta(minutes=5)
    end=start+timedelta(minutes=5)
    for plant in range(3):
        b=make_effort()
        b.update(individual_id=f"HC{plant}",capitulum_id=f"HH{plant}",
                 observation_bout_id=f"HB{plant}",
                 video_recording_id=f"synthetic://head-mix/{plant}",
                 evidence_uri=f"synthetic://head-mix/{plant}",
                 video_window_end_s="300",evaluable_duration_s="300",
                 zero_approaches_confirmed="0",
                 head_stage="full_anthesis",head_rank="terminal")
        eh.append(b)
        sh.append(dict(zip(SIDE_FIELDS,[
            b["individual_id"],b["population_id"],b["capitulum_id"],
            b["observation_bout_id"],start.isoformat(),end.isoformat(),
            "synthetic://head-mix","not_assessed","not_assessed",""
        ])))
        route="below" if plant==1 else "above"
        n_poll,n_seed=(10,90) if plant==1 else (90,10)
        for guild,n in zip(GUILDS,(n_poll,n_seed)):
            for j in range(n):
                ev.append({
                    "population_id":"TEST_POP",
                    "individual_id":f"HC{plant}",
                    "capitulum_id":f"HH{plant}",
                    "observation_bout_id":f"HB{plant}",
                    "approach_episode_id":f"HC{plant}_{guild}_{j}",
                    "time_from_bout_start_s":str(1+j*2.5),
                    "phenological_stage":"full_anthesis",
                    "reproductive_sex_state":"hermaphroditic",
                    "approached":"1",
                    "pre_entry_guild":guild,
                    "pre_entry_guild_evidence":"independently_assigned_before_access",
                    "visitor_taxon_label":f"synthetic:{guild}",
                    "evidence_uri":"synthetic://head-mix",
                    "entry_route_world":route,
                    "reproductive_zone_reached":"1",
                })
    headmix=audit(eh,ev,sh,contract)
    assert headmix["n_qualified_comparison_bins"]==1
    assert abs(headmix["matched_binned_route_overlap"]-(20/110+10/190))<EPS
    assert abs(headmix["matched_same_head_equal_weight_route_overlap"]-1)<EPS
    assert all(abs(x["route_overlap"]-1)<EPS
               for x in headmix["qualified_bins"][0]["head_conditioned_route_overlap"])
    return {
        "status":"PASS_SYNTHETIC_PARTIAL_OVERLAP_TIME_BIN_REPAIR",
        "naive_pooled_overlap":result["naive_pooled_route_overlap"],
        "same_5min_overlap":result["matched_binned_route_overlap"],
        "n_qualified_bins":result["n_qualified_comparison_bins"],
        "head_composition_counterexample_pooled":headmix["matched_binned_route_overlap"],
        "head_composition_counterexample_within_head":headmix["matched_same_head_equal_weight_route_overlap"],
        "first_bout":"09:00-09:15",
        "second_bout":"09:05-09:20",
        "shared_exposure":"09:05-09:15",
        "real_biological_observations":0,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--contract",type=Path,
                    default=Path("data/contracts/capitulum_video_effort_denominator_v1.json"))
    ap.add_argument("--effort",type=Path,
                    default=Path("data/intake/capitulum_video_effort_denominator_v1.csv"))
    ap.add_argument("--events",type=Path,
                    default=Path("data/intake/capitulum_guild_access_event_ledger_v1.csv"))
    ap.add_argument("--sidecar",type=Path,
                    default=Path("data/intake/capitulum_clock_matched_effort_sidecar_v1.csv"))
    ap.add_argument("--original-bouts",type=Path,
                    help="Native EAzami Aim2 bouts file, fixed version")
    ap.add_argument("--out",type=Path)
    a=ap.parse_args()
    contract=json.loads(a.contract.read_text(encoding="utf-8"))
    eh,eff=read_csv(a.effort)
    vh,events=read_csv(a.events)
    sh,sidecar=read_csv(a.sidecar)
    if eh!=contract["columns"] or len(vh)!=32 or sh!=SIDE_FIELDS:
        raise ValueError("DATA_LEDGER_SCHEMA_DRIFT")
    crosswalk=None
    if a.original_bouts is not None:
        oh,bouts=read_csv(a.original_bouts)
        if not {"observation_date","start_time_local",
                "end_time_local","phenological_stage"}.issubset(oh):
            raise ValueError("AUTHORITATIVE_ORIGINAL_BOUT_SCHEMA_DRIFT")
        crosswalk=crosswalk_original_bouts(bouts,sidecar,eff)
    result={
        "version":"capitulum_five_minute_common_coverage_v1",
        "date":"2026-10-08",
        "original_bout_crosswalk":crosswalk,
        "synthetic_tests":synthetic_tests(contract),
        "real_data":audit(eff,events,sidecar,contract),
    }
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",
                         encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
