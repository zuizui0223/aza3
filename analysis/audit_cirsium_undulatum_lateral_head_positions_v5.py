#!/usr/bin/env python3
"""Independent within-plant spatial positional opportunities in Cirsium undulatum.

Use original-authenticated West & Louda (2021) USDA CC0 mother/head records
(10.15482/USDA.ADC/1522649), via the exact 2008 mother–head crosswalk.
'Head position' codes 200, 210, 300, ... denote the original source branch
positions, not validated time order or developmental fate.
Source original authors own the experiment and compensate-potential result;
we neither infer its ecological cause nor re-use its sample as independent.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from audit_cirsium_undulatum_scored_head_failure_v4 import (
    source_locked_tables, head_mother_crosswalk, good_seed_int
)

ARMS=("ndni","dni","ndi","di")
ROUTE_CLASSES=(
    "first_zero_later_positive",
    "first_positive_later_positive",
    "first_positive_later_zero",
    "first_zero_later_zero",
)


def require(ok,why):
    if not ok:
        raise ValueError(why)


def source_prefix_group(plant):
    match=re.fullmatch(r"([A-Z]+)(?:[0-9]+(?:\.[0-9]+)?)",plant)
    require(match is not None,"UNRECOGNIZED_SOURCE_PLANT_ID_PREFIX")
    return match.group(1)


def classify_lateral_position_routes(heads):
    """Lowest numeric source head position, versus higher coded positions.

    This is intentionally *not* chronological head phenology.
    """
    if not heads:
        return {"class":"none_scored_lateral",
                "first_scored_lateral_position":None,
                "first_scored_lateral_good_seed":None,
                "any_higher_position_good_seed":None,
                "n_scored_lateral_heads":0}
    parsed=[]
    codes=set()
    for h in heads:
        raw=str(h.get("hd",""))
        require(raw.isdigit() and int(raw)>100,
                "INVALID_HEAD_POSITION")
        pos=int(raw)
        require(pos not in codes,"DUPLICATE_LATERAL_HEAD_POSITION")
        codes.add(pos)
        count=good_seed_int(h,"hd.gsd")
        parsed.append((pos,count))
    parsed.sort()
    first,first_count=parsed[0]
    if len(parsed)==1:
        status="one_scored_lateral"
        higher=None
    else:
        higher=any(seed>0 for _,seed in parsed[1:])
        first_positive=first_count>0
        status=("first_positive_" if first_positive else "first_zero_")+(
            "later_positive" if higher else "later_zero")
    return {
        "class":status,
        "first_scored_lateral_position":first,
        "first_scored_lateral_good_seed":first_count,
        "any_higher_position_good_seed":higher,
        "n_scored_lateral_heads":len(parsed),
    }


def analyse_lateral_positions(mothers=None, heads=None,strict=True):
    if mothers is None and heads is None:
        mothers,heads=source_locked_tables()
    require(mothers is not None and heads is not None,
            "NEED_BOTH_MOTHER_AND_HEAD_SOURCE")
    joined,source2007=head_mother_crosswalk(
        mothers,heads,strict2008=strict)
    arm_report={}
    by_prefix={}
    for arm in ARMS:
        arm_items=[x for x in joined if x["arm"]==arm]
        n2=0
        classified={key:0 for key in ROUTE_CLASSES}
        first_zero_higher_positive=[]
        apical_zero_for_that=0
        lat_pos=0
        for p in arm_items:
            m=p["mother"]
            code=source_prefix_group(m["plt"])
            require(code in ("ES","SD","WG","WR"),
                    "NONCANONICAL_SOURCE_PREFIX")
            prefix=by_prefix.setdefault(code,{
                "total_mothers":0,
                "original_arms":{a:{
                    "n":0,"both_seed_routes_zero":0,
                    "lateral_good_seed_positive":0,
                    "mother_has_scored_lateral":0}
                    for a in ARMS}})
            prefix["total_mothers"]+=1
            row=prefix["original_arms"][arm]
            row["n"]+=1
            ap=good_seed_int(m,"ap.hd")
            lateral=good_seed_int(m,"o.hds")
            row["both_seed_routes_zero"]+=int(
                ap==0 and lateral==0)
            row["lateral_good_seed_positive"]+=int(lateral>0)
            row["mother_has_scored_lateral"]+=int(bool(p["lateral"]))
            lat_pos+=int(lateral>0)
            state=classify_lateral_position_routes(p["lateral"])
            if state["n_scored_lateral_heads"]<2:
                continue
            n2+=1
            status=state["class"]
            require(status in classified,
                    "LATERAL_TWOHEAD_CLASSIFICATION_INVALID")
            classified[status]+=1
            if status=="first_zero_later_positive":
                first_zero_higher_positive.append(m["plt"])
                apical_zero_for_that+=int(ap==0)
        first_zero_higher_positive.sort()
        arm_report[arm]={
            "original_mothers":len(arm_items),
            "lateral_seed_positive_mothers":lat_pos,
            "mothers_with_2plus_scored_lateral_heads":n2,
            "lateral_rank_presence_2plus":{
                "mothers":n2,
                "first_zero_later_positive":
                    classified["first_zero_later_positive"],
                "first_positive_later_positive":
                    classified["first_positive_later_positive"],
                "first_positive_later_zero":
                    classified["first_positive_later_zero"],
                "first_zero_later_zero":
                    classified["first_zero_later_zero"],
            },
            "first_zero_later_positive_mothers":
                first_zero_higher_positive,
            "first_zero_later_positive_mothers_with_apical_good_seed_zero":
                apical_zero_for_that,
            "first_zero_later_positive_in_positive_mothers":
                len(first_zero_higher_positive),
        }
    order=("ES","SD","WG","WR")
    by_prefix={x:by_prefix[x] for x in order}
    report={
        "schema":"aza3_USDA_Cirsium_undulatum_lateral_position_opportunity_v5",
        "source":{
            "dataset_doi":"10.15482/USDA.ADC/1522649",
            "paper_doi":"10.1002/ajb2.1744",
            "original_mothers":len(mothers),
            "original_heads":len(heads),
            "jointly_validated_2008_mothers":len(joined),
            "original_2008_scored_lateral_heads":
                sum(len(x["lateral"]) for x in joined),
            "mothers_with_2plus_scored_lateral_heads":
                sum(x["mothers_with_2plus_scored_lateral_heads"]
                    for x in arm_report.values()),
            "all_2008_head_good_seed_sums_match_maternal_source":True,
            "2007_head_archive_full_identity_not_verified":
                not source2007["full_2007_mother_head_joint_state_authenticated"],
            "study_source_not_a_new_independent_sample":True,
        },
        "by_2008_treatment":arm_report,
        "original_id_prefix_strata":by_prefix,
        "biology":{
            "original_nonapical_heads_do_have_multiple_source_position_codes":True,
            "first_position_reproductive_failure_can_coexist_with_other_position_positive":True,
            "head_position_coding_is_spatial_not_observed_event_timing":True,
            "original_joint_insecticide_effect_not_a_new_discovery":True,
        },
        "claim_limits":{
            "head_position_codes_are_proven_chronological_order":False,
            "higher_numbered_head_is_proven_later_flowering":False,
            "source_id_prefix_verified_as_geographic_site":False,
            "first_lower_code_head_failed_due_to_specific_arthropod":False,
            "higher_position_head_caused_fitness_compensation":False,
            "endogenous_head_number_safe_as_causal_adjustment":False,
            "observed_route_mosaic_is_genetic_division_of_labour":False,
            "native_ancestral_organ_regain_from_F2_or_Cirsium_shown":False,
            "direct_SLK_rare_invasion_or_complete_WF_architecture_measured":False,
            "ANTO_regained_molecular_competence_identified":False,
        },
        "conclusion":
          "ACTUAL_NONAPICAL_HEAD_CODE_ALTERNATIVES_WITH_POSITIONAL_FALLBACK_NOT_EVOLUTIONARY_REGENERATION",
    }
    if strict:
        expect={
            "ndni":(6,1,0,0,5),
            "dni":(6,2,0,3,1),
            "ndi":(8,1,6,1,0),
            "di":(13,0,12,1,0),
        }
        require(len(joined)==67 and
                report["source"]["original_2008_scored_lateral_heads"]==111,
                "PINNED_2008_MOTHER_HEAD_COUNT_DRIFT")
        for arm,(n,v1,v2,v3,v4) in expect.items():
            t=arm_report[arm]["lateral_rank_presence_2plus"]
            require(tuple(t.values())==(n,v1,v2,v3,v4),
                    "ORIGINAL_2008_POSITIONAL_ROUTE_COUNTS_DRIFT_"+arm)
        require(arm_report["dni"]["first_zero_later_positive_mothers"]
                ==["ES68","SD05"],
                "ORIGINAL_DAMAGED_AMBIENT_MULTILEVEL_FALLBACK_DRIFT")
        require({
            prefix:by_prefix[prefix]["original_arms"]["dni"]
                    ["both_seed_routes_zero"] for prefix in order
        }=={"ES":0,"SD":2,"WG":2,"WR":1},
                "ORIGINAL_PREFIX_FAILURE_DISTRIBUTION_DRIFT")
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    outcome=analyse_lateral_positions(strict=True)
    output=json.dumps(outcome,indent=2,sort_keys=True,
                      ensure_ascii=False,allow_nan=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(output,encoding="utf-8")
    else:
        print(output)


if __name__=="__main__":
    main()
