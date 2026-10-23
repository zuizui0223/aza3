#!/usr/bin/env python3
"""Observed branch-terminal vs stemward-intermediate flower heads within mothers.

West & Louda (2021), DOI 10.1002/ajb2.1744, USDA CC0
DOI 10.15482/USDA.ADC/1522649. Source has 106 mothers and 323 real heads.
The original USDA head-position dictionary describes hd=100 as APICAL,
other x00 codes as branch-terminal, and inner-branch x10 codes as
between terminal ends and stem. These categorical placement codes are not
angles, spine strength, head age, habitat exposure, randomized assignments
or evidence of genetic LEGO components.

Only 14 mothers have both lateral head placement categories: 7 per year.
All positional contrasts condition on those actual mothers, rather than
pretending 194 lateral heads or 106 mothers are independent matches.
Experimental damage and insecticide treatment and mother/year differences
remain confounded with placement and appearance. This is post hoc,
source-descriptive and exploratory, not proof of position/structural
access causation or replicated life-history evolution.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import json
import math
import re
from pathlib import Path
from audit_cirsium_undulatum_maternal_hurdle_size_v1 import (
    read_locked, source_groups, value, BLOBS
)

YEARS=("2007","2008")
CATS=("apical","branch_terminal","intermediate")
METRICS=("good_seed_per_head","seed_positive_fraction","damage_score")
RNG_SEED=20261023
REPS=5000
MASK=0xffffffff

def check(v,msg):
    if not v: raise ValueError(msg)

def identity(h):
    return h["year"]+"|"+h["plt"]

def position(head):
    code=str(head.get("hd",""))
    check(bool(re.fullmatch(r"[1-9]\d*0",code)),"INVALID_HEAD_POSITION_CODE")
    n=int(code)
    if n==100:return "apical"
    if n%100==0:return "branch_terminal"
    if n%10==0:return "intermediate"
    raise ValueError("INVALID_HEAD_POSITION_CODE")

def avg(xs):
    check(bool(xs),"EMPTY_ORIGINAL_HEAD_GROUP")
    return sum(xs)/len(xs)

def head_stats(heads):
    check(bool(heads),"EMPTY_ORIGINAL_HEAD_GROUP")
    scored=[value(h,"av.dam") for h in heads if value(h,"av.dam") is not None]
    return{
      "n_heads":len(heads),
      "n_scored_damage":len(scored),
      "n_positive_seed_heads":sum(value(h,"hd.gsd")>0 for h in heads),
      "mean_good_seed":avg([value(h,"hd.gsd") for h in heads]),
      "seed_positive_fraction":
         avg([float(value(h,"hd.gsd")>0) for h in heads]),
      "mean_damage":avg(scored) if scored else None,
    }

def raw_audit(mothers,heads):
    parents={identity(m):m for m in mothers}
    check(len(mothers)==106 and len(parents)==106 and len(heads)==323,
          "ORIGINAL_USDA_106_MOTHERS_323_HEADS_CHANGED")
    original_hid=[identity(h)+"|"+h["hd"] for h in heads]
    check(len(set(original_hid))==len(heads),"DUPLICATE_SOURCE_HEAD_ID")
    matched=[h for h in heads if identity(h) in parents]
    conflicting=[h for h in matched
                 if parents[identity(h)]["dam"]!=h["dam"]
                 or parents[identity(h)]["ins"]!=h["ins"]]
    clean=[h for h in matched if h not in conflicting]
    check(len(matched)==309 and len(conflicting)==12 and len(clean)==297
          and sorted({identity(h) for h in conflicting})==
          ["2007|323","2007|375"],"UNRESOLVED_SOURCE_TREATMENT_JOIN_CHANGED")
    counts=Counter(position(h) for h in clean)
    check(dict(counts)=={"apical":103,"branch_terminal":173,"intermediate":21},
          "ORIGINAL_BRANCHED_HEAD_COUNTS_DRIFT")
    return parents,clean,{
       "n_original_mothers":len(mothers),
       "n_original_head_rows":len(heads),
       "n_head_rows_matching_mother":len(matched),
       "n_original_rows_with_conflicting_treatment":len(conflicting),
       "n_clean_head_rows":len(clean),
       "source_position_counts":{c:counts[c] for c in CATS},
       "source_conflict_mother_ids":["2007|323","2007|375"],
       "head_counts_not_independent_replicate_mothers":True,
    }

def pair_mothers(mothers,clean):
    ix=defaultdict(lambda:{c:[] for c in CATS})
    for h in clean:ix[identity(h)][position(h)].append(h)
    out=[]
    for ident in sorted(ix):
        categories=ix[ident]
        if not (categories["branch_terminal"] and categories["intermediate"]):
            continue
        branch=head_stats(categories["branch_terminal"])
        within=head_stats(categories["intermediate"])
        check(branch["mean_damage"] is not None
              and within["mean_damage"] is not None,
              "UNSCORED_HEAD_POSITIONS_CANNOT_BE_ZERO_IMPUTED")
        m=next(m for m in mothers if identity(m)==ident)
        out.append({
          "id":ident,"year":m["year"],"arm":m["trt"],
          "n_branch_terminal":branch["n_heads"],
          "n_intermediate":within["n_heads"],
          "delta":{
            "good_seed_per_head":within["mean_good_seed"]-branch["mean_good_seed"],
            "seed_positive_fraction":
                within["seed_positive_fraction"]-branch["seed_positive_fraction"],
            "damage_score":within["mean_damage"]-branch["mean_damage"],
          },
        })
    check(len(out)==14 and
          Counter(x["year"] for x in out)=={"2007":7,"2008":7},
          "PAIRED_INTERMEDIATE_HEAD_MOTHER_COHORT_CHANGED")
    return out

def choose(n,k):
    if k<0 or k>n:return 0
    k=min(k,n-k)
    out=1
    for j in range(1,k+1):
        out=out*(n-k+j)//j
    return out

def exact_sign_p(plus,minus):
    n=plus+minus
    if n==0:return 1.0
    return min(1.,2*sum(choose(n,k) for k in
                        range(min(plus,minus)+1))/(2**n))

def bh(ps):
    # exact, conservative family of 3 exploratory outcomes x 2 years
    order=sorted(ps,key=lambda z:ps[z])
    k=len(ps)
    vals={}
    bound=1.
    for i in range(k,0,-1):
        name=order[i-1]
        bound=min(bound,k*ps[name]/i)
        vals[name]=bound
    return vals

def xorshift(seed):
    check(0<seed<2**32,"BOOT_SEED_OUT_OF_BOUNDS")
    state=seed
    while True:
        state^=(state<<13)&MASK
        state^=state>>17
        state^=(state<<5)&MASK
        state&=MASK
        yield state

def quantile(a,p):
    z=sorted(a)
    index=(len(z)-1)*p
    j=int(math.floor(index))
    frac=index-j
    return z[j]*(1-frac)+z[min(j+1,len(z)-1)]*frac

def analyze(mothers,heads,source_dictionary,reps=REPS):
    check("Flowerhead position on plant" in source_dictionary
          and 'Apical head is ""100""' in source_dictionary
          and "between the terminals and stem" in source_dictionary,
          "SOURCE_HEAD_POSITION_ONTOLOGY_UNVERIFIED")
    check(100<=reps<=15000,"REPLICATE_BUDGET_INVALID")
    source_groups(mothers)
    parents,clean,audit=raw_audit(mothers,heads)
    pairs=pair_mothers(mothers,clean)
    ledger={}
    by_year={}
    tests={}
    for y in YEARS:
        r=[x for x in pairs if x["year"]==y]
        out={"n_paired_mothers":len(r),
             "n_branch_terminal_heads_in_pairs":sum(x["n_branch_terminal"] for x in r),
             "n_intermediate_heads_in_pairs":sum(x["n_intermediate"] for x in r),
             "position_differences":{}}
        for m in METRICS:
            d=[z["delta"][m] for z in r]
            sorted_d=sorted(d)
            plus=sum(z>1e-9 for z in d)
            minus=sum(z< -1e-9 for z in d)
            ties=len(r)-plus-minus
            pv=exact_sign_p(plus,minus)
            tests[y+"|"+m]=pv
            out["position_differences"][m]={
               "mean_intermediate_minus_branch_terminal":avg(d),
               "median_intermediate_minus_branch_terminal":
                   (sorted_d[3]+sorted_d[4])/2,
               "n_increase":plus,
               "n_decrease":minus,
               "n_equal":ties,
               "posthoc_exact_two_sided_sign_p":pv,
               "paired_mother_delta_ordered":d,
            }
        by_year[y]=out
    adjusted=bh(tests)
    for y in YEARS:
        for metric in METRICS:
            by_year[y]["position_differences"][metric][
                "BH_q_six_exploratory_sign_tests"]=adjusted[y+"|"+metric]
    rng=xorshift(RNG_SEED)
    draws={y:{m:[] for m in METRICS} for y in YEARS}
    for _ in range(reps):
        for y in YEARS:
            r=[x for x in pairs if x["year"]==y]
            z={m:0. for m in METRICS}
            for _i in range(len(r)):
                x=r[next(rng)%len(r)]
                for m in METRICS:z[m]+=x["delta"][m]
            for m in METRICS:draws[y][m].append(z[m]/len(r))
    for y in YEARS:
        for metric in METRICS:
            v=by_year[y]["position_differences"][metric]
            v["paired_mother_bootstrap_ci95"]=[
                quantile(draws[y][metric],.025),
                quantile(draws[y][metric],.975),
            ]
    all_rows=[]
    for y in YEARS:
        all_rows.append({
           "year":y,
           "head_position":{c:head_stats([
                h for h in clean if h["year"]==y and position(h)==c
           ]) for c in CATS}
        })
    assert set(parents)=={identity(m) for m in mothers}
    return{
        "version":"aza3_cirsium_undulatum_within_branch_position_audit_v1",
        "source":{
           "dataset_doi":"10.15482/USDA.ADC/1522649",
           "paper_doi":"10.1002/ajb2.1744",
           "git_blob_sources":BLOBS,
           "head_dictionary_position_rule":"100 apex; other x00 branch-terminal; x10 between branch terminal and stem",
        },
        "source_audit":audit,
        "n_mothers_with_any_intermediate_head":len({
            identity(h) for h in clean if position(h)=="intermediate"
        }),
        "n_mothers_with_both_lateral_positions":len(pairs),
        "head_position_by_year":all_rows,
        "paired_mothers":pairs,
        "paired_within_mother_year_results":by_year,
        "bootstrap":{
            "n_draws":reps,"seed":RNG_SEED,
            "resampling_unit":"MOTHER_WITH_BOTH_LATERAL_POSITION_TYPES",
            "sampling_order":"replicate > year 2007 then 2008 > mothers lexicographically",
            "ci_scope":"descriptive bootstrap conditional on 14 paired mothers, 7 per year",
        },
        "limits":{
           "not_original_head_position_randomization":True,
           "no_explicit_videotaped_oviposition_or_pollination":True,
           "posthoc_six_tests_do_not_confirm_architectural_adaptation":True,
           "conditional_on_position_both_types_recorded":True,
           "external_tissue_damage_measurement_not_a_direct_florivore_causal_path":True,
           "very_small_seven_paired_mothers_in_2008":True,
           "head_code_not_anatomical_angle_or_spine_stiffness":True,
           "new_individual_head_observations_collected":0,
        },
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mothers",required=True,type=Path)
    p.add_argument("--heads",required=True,type=Path)
    p.add_argument("--dictionary",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--bootstrap",type=int,default=REPS)
    q=p.parse_args()
    mothers=read_locked(q.mothers,"mother_plants.csv",
        ("plt","year","dam","ins","trt"))
    heads=read_locked(q.heads,"flower_heads.csv",
        ("plt","year","dam","ins","hd","av.dam","hd.gsd"))
    dictionary=read_locked(q.dictionary,"source_dictionary.txt")
    z=analyze(mothers,heads,dictionary,q.bootstrap)
    q.out.parent.mkdir(parents=True,exist_ok=True)
    q.out.write_text(json.dumps(z,indent=2,ensure_ascii=False,
                                allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "source_audit":z["source_audit"],
        "n_mothers_paired":z["n_mothers_with_both_lateral_positions"],
        "paired_year_results":{
            y:{m:{k:v for k,v in payload.items() if k!="paired_mother_delta_ordered"}
               for m,payload in d["position_differences"].items()}
            for y,d in z["paired_within_mother_year_results"].items()}
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
