#!/usr/bin/env python3
"""Source-authenticated actual Cirsium undulatum within-mother head injury and good seed ranks.

Observational same-year same-MOTHER contrasts of original reproductive heads:
West & Louda 2021 DOI 10.1002/ajb2.1744, USDA CC0 10.15482/USDA.ADC/1522649.

Biological question: does having a LOWER observed receptacle-damage score
necessarily yield MORE author-classified good seeds on the SAME mother's
other lateral, specifically another branch-terminal head? Source authentic
negative answers do not prove a positive or negative injury selection
gradient. Head-level good-seed number is also driven by source-unmeasured
size, pollen, timing, stage, resources and developmental position.

The main head-site-comparable diagnostic uses multiple terminal positions
(number x00 except main apical 100) per mother; the secondary analysis uses
all lateral positions (x00 and x10). Select source minimum/maximum injury
levels WITHIN each source mother; when several heads tie for either level,
average their good-seed counts, rather than cherry-pick seed outcomes.
Such an injury-extrema conditional comparison is descriptive and not an
instrument or a randomized injury manipulation.

Source av.dam is an averaged QUARTER-coded receptacle injury SCORE (original
maximum 0.88), NOT literal 88% tissue killed, insect identification or
confirmed seed predation. "Good seed" is author-scored visually developed
undamaged achene, NOT independent viability or germination assay.

No new field observations or Cirsium genetic trait reexpression studies.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path

SHA={
    "mother_plants.csv":"458800272c3d0a03c7653702027a8f8578d9e88c",
    "flower_heads.csv":"637a4dcd5194d0bf5b26e7f9cf0dde56a50d7616",
    "source_dictionary.txt":"5e20c9d4609abc982287bd12f739c76d7c6d1d2f",
}
YEARS=("2007","2008")
ARMS=("ndni","dni","ndi","di")
PAPER="10.1002/ajb2.1744"
DATASET="10.15482/USDA.ADC/1522649"

def ensure(c,reason):
    if not c:raise ValueError(reason)

def source_blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode("ascii")+b"\0"+data).hexdigest()

def parse_locked(path,basename):
    data=Path(path).read_bytes()
    ensure(source_blob(data)==SHA[basename],"SOURCE_BLOB_CHANGED_"+basename)
    if basename=="source_dictionary.txt":
        s=data.decode("utf-8")
        ensure('Apical head is ""100""' in s and
               "four quarters of the receptacle" in s,
               "AUTHOR_DICTIONARY_POSITION_OR_DAMAGE_CHANGED")
        return s
    reader=csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
    require=({"plt","year","trt","dam","ins"} if basename=="mother_plants.csv"
             else {"plt","year","dam","ins","hd","av.dam","hd.gsd"})
    ensure(require.issubset(reader.fieldnames or ()),
           "ORIGINAL_MOTHER_HEAD_SCHEMA_DRIFT")
    return list(reader)

def mother_id(r):
    return r["year"]+"|"+r["plt"]

def good(r):
    v=float(r["hd.gsd"])
    ensure(math.isfinite(v) and v>=0 and v.is_integer(),
           "INVALID_AUTHOR_GOOD_SEEDS")
    return int(v)

def damage(r):
    try:v=float(r["av.dam"])
    except (ValueError,TypeError):return None
    ensure(math.isfinite(v) and 0<=v<=.88,
           "INVALID_CODED_RECEPTACLE_DAMAGE")
    return v

def audit_sources(mothers, heads, dictionary):
    ensure('Apical head is ""100""' in dictionary and
           "four quarters of the receptacle" in dictionary,
           "AUTHOR_DICTIONARY_POSITION_OR_DAMAGE_CHANGED")
    ensure(len(mothers)==106 and len(heads)==323,
           "ORIGINAL_106_MOTHER_323_HEAD_SOURCE_CHANGED")
    mother={mother_id(m):m for m in mothers}
    ensure(len(mother)==106,"DUPLICATE_SOURCE_ORIGINAL_MOTHER")
    ensure(all(m["year"] in YEARS and m["trt"] in ARMS and
               m["trt"]==m["dam"]+m["ins"] for m in mothers),
           "ORIGINAL_MOTHER_TREATMENT_ASSIGNMENT_CHANGED")
    original_head_ids=[mother_id(h)+"|"+h["hd"] for h in heads]
    ensure(len(original_head_ids)==len(set(original_head_ids)),
           "DUPLICATE_ORIGINAL_HEAD_RECORD")
    match=[h for h in heads if mother_id(h) in mother]
    conflict=[h for h in match if
              (h["dam"]!=mother[mother_id(h)]["dam"] or
               h["ins"]!=mother[mother_id(h)]["ins"])]
    clean=[h for h in match if
           (h["dam"]==mother[mother_id(h)]["dam"] and
            h["ins"]==mother[mother_id(h)]["ins"])]
    ensure(len(match)==309 and len(conflict)==12 and len(clean)==297 and
           sorted({mother_id(r) for r in conflict})==
               ["2007|323","2007|375"],
           "WRONG_ORIGINAL_PARENT_HEAD_JOIN_OR_TREATMENT_CONFLICTS")
    lateral=[h for h in clean if h["hd"]!="100"]
    scored=[h for h in lateral if damage(h) is not None]
    terminal=[h for h in scored if int(h["hd"])%100==0]
    intermediate=[h for h in scored if int(h["hd"])%100!=0]
    ensure(len(lateral)==194 and len(scored)==192 and
           len(terminal)==171 and len(intermediate)==21,
           "ORIGINAL_HEAD_SITE_AND_DAMAGE_SCORE_COUNTS_DRIFT")
    ensure(all(int(h["hd"])>100 for h in lateral) and
           all(int(h["hd"])%100==10 for h in intermediate),
           "UNEXPECTED_SOURCE_HEAD_POSITION_DIGITS")
    # Some original source heads have missing damage; never change to 0.
    # The original 0.88 is top recorded coded score, not percent tissue.
    ensure(max(damage(x) for x in scored)==.88,"CODED_DAMAGE_CEILING_CHANGED")
    for h in scored:good(h)
    return mother,scored,terminal,intermediate,{
        "original_mother_rows":len(mothers),
        "original_head_rows":len(heads),
        "head_with_source_mother":len(match),
        "head_without_source_mother":len(heads)-len(match),
        "treatment_mismatch_head_rows":len(conflict),
        "consistent_heads":len(clean),
        "consistent_lateral_head_rows":len(lateral),
        "scored_lateral_head_rows":len(scored),
        "scored_terminal_head_rows":len(terminal),
        "scored_intermediate_head_rows":len(intermediate),
    }

def make_stage(rows,mothers):
    grouped=defaultdict(list)
    for h in rows:
        grouped[mother_id(h)].append(h)
    records=[]
    for key,allheads in grouped.items():
        if len(allheads)<2:continue
        inj=[damage(x) for x in allheads]
        low,high=min(inj),max(inj)
        if low==high:continue
        lesser=[h for h in allheads if damage(h)==low]
        greater=[h for h in allheads if damage(h)==high]
        def head_summary(heads):
            return sorted([{
                "head_code":h["hd"],
                "damage_score":damage(h),
                "good_seed_count":good(h),
            } for h in heads],key=lambda z:int(z["head_code"]))
        lo=head_summary(lesser);hi=head_summary(greater)
        lseed=sum(x["good_seed_count"] for x in lo)/len(lo)
        hseed=sum(x["good_seed_count"] for x in hi)/len(hi)
        records.append({
            "mother":key,"year":allheads[0]["year"],
            "treatment":mothers[key]["trt"],
            "eligible_scored_heads":len(allheads),
            "least_damaged_head_count":len(lo),
            "most_damaged_head_count":len(hi),
            "least_damaged_score":low,"most_damaged_score":high,
            "damage_score_gap":high-low,
            "least_damaged_mean_good_seed":lseed,
            "most_damaged_mean_good_seed":hseed,
            "lower_damage_minus_higher_damage_good_seed":lseed-hseed,
            "lower_damage_positive_rate":sum(x["good_seed_count"]>0 for x in lo)/len(lo),
            "higher_damage_positive_rate":sum(x["good_seed_count"]>0 for x in hi)/len(hi),
            "severe_success_vs_low_failure":(
                all(x["good_seed_count"]==0 for x in lo)
                and any(x["good_seed_count"]>0 for x in hi)),
            "low_heads":lo,"high_heads":hi,
        })
    return sorted(records,key=lambda z:z["mother"])

def sign_counts(data):
    vals=[x["lower_damage_minus_higher_damage_good_seed"] for x in data]
    return {
        "lower_damage_has_more_good_seeds":sum(x>0 for x in vals),
        "higher_damage_has_more_good_seeds":sum(x<0 for x in vals),
        "equal_good_seed_counts":sum(x==0 for x in vals),
    }

def exact_sign(pos,neg):
    ensure(all(type(x) is int and x>=0 for x in (pos,neg)),
           "INVALID_EXACT_SIGNS")
    n=pos+neg
    if not n:return 1.0
    return min(1.0,2*sum(math.comb(n,k)
                        for k in range(min(pos,neg)+1))/2**n)

def summary(rows):
    s=sign_counts(rows)
    n=len(rows)
    return {
        "n_original_mothers":n,
        "directions":s,
        "exact_two_sided_p_for_symmetry_exploratory":exact_sign(
            s["lower_damage_has_more_good_seeds"],
            s["higher_damage_has_more_good_seeds"]),
        "mean_lower_minus_higher_good_seed":(
            sum(x["lower_damage_minus_higher_damage_good_seed"]
                for x in rows)/n if n else None),
        "n_higher_damage_head_produces_any_good_seed":sum(
            any(z["good_seed_count"]>0 for z in x["high_heads"])
            for x in rows),
        "n_lower_damage_head_fails_while_higher_damage_succeeds":sum(
            x["severe_success_vs_low_failure"] for x in rows),
    }

def strata(records):
    return {
        year:{
            "total":summary([x for x in records if x["year"]==year]),
            "by_original_arm":{
                arm:summary([x for x in records if
                            x["year"]==year and x["treatment"]==arm])
                for arm in ARMS
            },
        } for year in YEARS
    }

def same_damage_groups(terminal,mothers):
    groups=defaultdict(list)
    for h in terminal:
        groups[mother_id(h)+"|"+h["av.dam"]].append(h)
    result=[]
    for key,hs in groups.items():
        if len(hs)<2:continue
        mkey=mother_id(hs[0])
        seeds=[good(x) for x in hs]
        result.append({
            "mother":mkey,
            "year":hs[0]["year"],
            "treatment":mothers[mkey]["trt"],
            "damage_score":damage(hs[0]),
            "n_heads":len(hs),
            "min_good_seed":min(seeds),"max_good_seed":max(seeds),
            "different_good_seed_counts":max(seeds)!=min(seeds),
            "mixed_zero_and_positive":(0 in seeds and any(x>0 for x in seeds)),
            "heads":sorted([{
                "head_code":h["hd"],"good_seed_count":good(h)
            } for h in hs],key=lambda z:int(z["head_code"])),
        })
    return sorted(result,key=lambda z:z["mother"])

def same_summary(records):
    return {
        "n_same_source_damage_score_terminal_groups":len(records),
        "n_distinct_source_mothers":len({r["mother"] for r in records}),
        "n_same_damage_but_different_good_seed_counts":sum(
            x["different_good_seed_counts"] for x in records),
        "n_same_damage_but_seed_positive_vs_zero":sum(
            x["mixed_zero_and_positive"] for x in records),
        "by_year":{year:{
            "eligible_groups":sum(x["year"]==year for x in records),
            "variable_good_seed_groups":sum(
                x["year"]==year and x["different_good_seed_counts"]
                for x in records),
        } for year in YEARS},
    }

def nested_pairs(scored):
    by_code=defaultdict(dict)
    for h in scored:
        n=int(h["hd"])
        block=n//100
        which="terminal" if n%100==0 else "intermediate"
        key=mother_id(h)+"|"+str(block)
        ensure(which not in by_code[key],"DUPLICATE_SOURCE_NESTED_BRANCH_LOCATION")
        by_code[key][which]=h
    matches=[x for x in by_code.values() if set(x)=={"terminal","intermediate"}]
    ensure(len(matches)==20,"ORIGINAL_NESTED_20_PAIRS_CHANGED")
    out={}
    def direction(z):
        return "positive" if z>1e-10 else "negative" if z< -1e-10 else "equal"
    for year in YEARS:
        # Source code family key still includes "YYYY|"
        rs=[v for k,v in by_code.items() if k.startswith(year+"|")
            and set(v)=={"terminal","intermediate"}]
        ct={}
        for pair in rs:
            term,inner=pair["terminal"],pair["intermediate"]
            d=damage(inner)-damage(term)
            s=good(inner)-good(term)
            key=direction(d)+"_damage__"+direction(s)+"_good_seed"
            ct[key]=ct.get(key,0)+1
        out[year]={
            "n":len(rs),
            "damage_and_seed_direction_cross_tab":ct,
        }
    return {
        "n_parent_source_code_matched_pairs":len(matches),
        "by_year":out,
    }

def analyze(mothers,heads,dictionary):
    ix,scored,terminal,intermediate,audit=audit_sources(
        mothers,heads,dictionary)
    a=make_stage(scored,ix)
    b=make_stage(terminal,ix)
    e=same_damage_groups(terminal,ix)
    result={
        "schema":"aza3_real_undulatum_same_mother_lateral_damage_seed_rank_original_JS_v1",
        "as_of":"2026-10-24",
        "source":{
            "author":"West & Louda 2021",
            "paper_doi":PAPER,
            "usda_doi":DATASET,
            "original_mother_git_blob_sha":SHA["mother_plants.csv"],
            "original_head_git_blob_sha":SHA["flower_heads.csv"],
            "original_dictionary_git_blob_sha":SHA["source_dictionary.txt"],
        },
        "grain":(
          "Within one true original year x mother, compare minimum vs "
          "maximum author receptacle damage score among eligible scored "
          "lateral heads; mean seeds across ties; compare terminal x00 "
          "separately to hold source head position class constant; two "
          "independent mothers may not be same randomized trait treatment."),
        "source_audit":audit,
        "damage_ontology":{
            "maximum_recorded_score":.88,
            "max_score_is_coded_not_exact_fraction_of_tissue_lost":True,
            "author_good_seeds_are_not_verified_germination":True,
        },
        "analyses":{
            "lateral_all_classes":{
                "summary":summary(a),"by_year_arm":strata(a),"per_mother":a,
            },
            "terminal_only":{
                "summary":summary(b),"by_year_arm":strata(b),"per_mother":b,
            },
        },
        "equal_scored_damage_distinct_seed_outputs":{
            "summary":same_summary(e),"groups":e,
        },
        "prior_source_nested_branch_germination_not_measured":
            nested_pairs(scored),
        "limits":{
            "only_mothers_with_min_two_scored_head_same_class_and_different_damage_eligible":True,
            "min_max_selection_on_observed_posttreatment_damage_not_randomized":True,
            "head_size_and_flowering_age_not_identified":True,
            "source_damaged_appearance_and_good_seed_are_related_but_different_scoring":True,
            "exclusions_conflicting_treatments_and_orphans_not_imputed_zero":True,
            "no_year_environment_only_conclusion_due_to_changed_protocol":True,
            "no_true_spine_phyllary_angle_or_verified_pollinator_service":True,
            "no_causal_insect_seed_fitness_or_genetic_LEGO_support":True,
            "n_new_field_observations":0,
        },
    }
    ensure(len(a)==48 and len(b)==46 and
           result["analyses"]["terminal_only"]["summary"]["directions"]=={
             "lower_damage_has_more_good_seeds":23,
             "higher_damage_has_more_good_seeds":18,
             "equal_good_seed_counts":5,
           } and len(e)==13 and
           result["equal_scored_damage_distinct_seed_outputs"]
                 ["summary"]["n_same_damage_but_seed_positive_vs_zero"]==8,
           "ORIGINAL_REAL_HEAD_DISCORDANCE_SOURCE_CHANGED")
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mothers",required=True,type=Path)
    p.add_argument("--heads",required=True,type=Path)
    p.add_argument("--dictionary",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    mothers=parse_locked(args.mothers,"mother_plants.csv")
    heads=parse_locked(args.heads,"flower_heads.csv")
    dictionary=parse_locked(args.dictionary,"source_dictionary.txt")
    result=analyze(mothers,heads,dictionary)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,
                                    allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "source":result["source_audit"],
        "lateral":result["analyses"]["lateral_all_classes"]["summary"],
        "terminal":result["analyses"]["terminal_only"]["summary"],
        "identical_damage":result["equal_scored_damage_distinct_seed_outputs"]
                          ["summary"],
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
