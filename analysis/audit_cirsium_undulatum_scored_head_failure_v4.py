#!/usr/bin/env python3
"""Observed head inventory versus current-season seed failure in Cirsium undulatum.

The focal positive biological evidence is **a scored nonapical head exists**
in the original *flower_heads.csv* even on maternal plants with no good seeds.
The recorded head is not automatically verified to have reached anthesis.

Original experiment West & Louda (2021), 10.1002/ajb2.1744; USDA CC0
10.15482/USDA.ADC/1522649. Existing aza3 source Git-SHA checks reused.
Only 2008 mothers/head records perfectly reconcile; 2007 head-level
reconciliation has missing mothers and treatment inconsistencies.
The damaging agent, pollen limitation, and genetic afterlife are unmeasured.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

from audit_cirsium_undulatum_maternal_hurdle_size_v1 import (
    read_locked, BLOBS, value
)
from audit_cirsium_undulatum_maternal_reproductive_routes_v2 import BASE

ARMS=("ndni","dni","ndi","di")
EXPECTED_SOURCE_MOTHERS=106
EXPECTED_SOURCE_HEADS=323


def ensure(ok, message):
    if not ok:
        raise ValueError(message)


def source_locked_tables():
    mothers=read_locked(
        BASE/"mother_plants.csv", "mother_plants.csv",
        ("plt","year","dam","ins","trt","plt.gsd","ap.hd",
         "o.hds","n.hds","no.flw"))
    heads=read_locked(
        BASE/"flower_heads.csv","flower_heads.csv",
        ("plt","year","dam","ins","hd","av.dam","hd.gsd"))
    return mothers,heads


def good_seed_int(row,field):
    x=value(row,field)
    ensure(x is not None and x>=0 and x.is_integer(),
           "INVALID_OR_MISSING_ORIGINAL_GOOD_SEED_"+field)
    return int(x)


def score_damage_as_number(raw):
    if raw is None or str(raw).strip().upper() in ("NA","","N/A","."):
        return None
    try:
        x=float(raw)
    except (ValueError,TypeError):
        raise ValueError("INVALID_DAMAGE_SCORE_VALUE") from None
    ensure(math.isfinite(x) and 0<=x<=1,"INVALID_DAMAGE_SCORE_VALUE")
    return x


def head_mother_crosswalk(mothers,heads,strict2008=True):
    """Link original rows by year + plant + position without repairing source.

    All 2008 head positions and seed sums reconcile. The 2007 archival
    limitations are OUTPUT, never silently propagated into an analysis.
    """
    ensure(len(mothers)>0 and len(heads)>0, "EMPTY_USDA_HEAD_OR_MOTHER_SOURCE")
    mothers_by_key={}
    for row in mothers:
        ident=(str(row["year"]),str(row["plt"]))
        ensure(ident not in mothers_by_key,"DUPLICATE_MOTHER_ID")
        mothers_by_key[ident]=row

    heads_by_key=defaultdict(list)
    positions=set()
    for row in heads:
        key=(str(row["year"]),str(row["plt"]),str(row["hd"]))
        ensure(key not in positions,"DUPLICATE_HEAD_POSITION")
        positions.add(key)
        ensure(key[2].isdigit() and int(key[2])>=100,
               "INVALID_SOURCE_HEAD_POSITION")
        _=good_seed_int(row,"hd.gsd")
        _=score_damage_as_number(row["av.dam"])
        heads_by_key[key[:2]].append(row)

    mismatch_lat=[]
    mismatch_ap=[]
    missing=[]
    wrong_treatment=[]
    verified=[]
    all_2008=[]
    for (year,plant),mother in mothers_by_key.items():
        source=heads_by_key.get((year,plant),[])
        ap=[h for h in source if h["hd"]=="100"]
        lat=[h for h in source if h["hd"]!="100"]
        sum_a=sum(good_seed_int(h,"hd.gsd") for h in ap)
        sum_l=sum(good_seed_int(h,"hd.gsd") for h in lat)
        expect_a=good_seed_int(mother,"ap.hd")
        expect_l=good_seed_int(mother,"o.hds")
        if not source and year=="2007":
            missing.append(f"{year}|{plant}")
        if expect_l != sum_l and year=="2007":
            mismatch_lat.append(f"{year}|{plant}")
        if expect_a != sum_a and year=="2007":
            mismatch_ap.append(f"{year}|{plant}")
        conflicting=[h for h in source if h["dam"]!=mother["dam"]
                      or h["ins"]!=mother["ins"]]
        if conflicting and year=="2007":
            wrong_treatment.append(f"{year}|{plant}")
        if year=="2008":
            ensure(len(ap)==1, "SOURCE_2008_APICAL_HEAD_NOT_UNIQUE")
            ensure(bool(source), "SOURCE_2008_MOTHER_WITHOUT_SCORED_HEAD")
            ensure(not conflicting, "SOURCE_2008_TREATMENT_CONFLICT")
            ensure(expect_a==sum_a and expect_l==sum_l,
                   "MOTHER_HEAD_GOOD_SEED_SUM_CONFLICT")
            ensure(good_seed_int(mother,"plt.gsd")==expect_a+expect_l,
                   "SOURCE_2008_MOTHER_TOTAL_GOOD_SEED_CONFLICT")
            ensure(mother["trt"] in ARMS,
                   "UNRECOGNIZED_ORIGINAL_USDA_TREATMENT")
            all_2008.append({
                "mother":mother,"apical":ap[0],"lateral":lat,
                "year":"2008","arm":mother["trt"]
            })
            verified.append(f"{year}|{plant}")

    if strict2008:
        ensure(len(mothers)==EXPECTED_SOURCE_MOTHERS and
               len(heads)==EXPECTED_SOURCE_HEADS,
               "ORIGINAL_USDA_SOURCE_COUNTS_CHANGED")
        ensure(len(all_2008)==67,
               "ORIGINAL_2008_MATERNAL_COHORT_CHANGED")
        ensure(sum(1 for h in heads if h["year"]=="2008")==178,
               "ORIGINAL_2008_HEAD_COUNT_CHANGED")
        ensure(len(verified)==67,
               "ORIGINAL_2008_MOTHER_HEAD_CROSSWALK_INCOMPLETE")

    report2007={
        "year":"2007",
        "mother_rows":sum(m["year"]=="2007" for m in mothers),
        "head_rows":sum(h["year"]=="2007" for h in heads),
        "mother_ids_without_scored_head_rows":sorted(missing),
        "mother_ids_with_lateral_seed_sum_disagreement":sorted(mismatch_lat),
        "mother_ids_with_apical_seed_sum_disagreement":sorted(mismatch_ap),
        "mother_ids_with_treatment_label_conflicts":sorted(wrong_treatment),
        "full_2007_mother_head_joint_state_authenticated":False,
    }
    return all_2008,report2007


def summarize_2008_group(mothers):
    """Aggregate one original mother per row, lateral heads separate scope."""
    r={
        "mothers":len(mothers),
        "mothers_with_lateral_head_record":0,
        "mothers_without_scored_lateral_head":0,
        "mothers_with_two_plus_lateral_head_records":0,
        "mothers_with_positive_lateral_good_seed":0,
        "mothers_with_scored_lateral_but_zero_lateral_seed":0,
        "apical_seed_positive_mothers":0,
        "mothers_both_seed_routes_zero":0,
        "joint_failure_mothers_with_scored_lateral_head":0,
        "joint_failure_mothers_with_two_plus_scored_lateral_heads":0,
        "joint_failure_mothers_with_3plus_source_buds":0,
        "joint_failure_mothers_with_lateral_scored_good_seed_zero":0,
        "scored_lateral_head_rows":0,
        "scored_lateral_good_seed_positive_head_rows":0,
        "scored_lateral_good_seed_zero_head_rows":0,
        "lateral_head_rows_damage_score_missing":0,
        "lateral_head_rows_damage_ge_080":0,
        "lateral_positive_seed_heads_damage_ge_080":0,
        "lateral_head_rows_at_max_observed_damage_088":0,
        "lateral_positive_seed_heads_at_max_damage_088":0,
        "lateral_zero_seed_heads_at_max_damage_088":0,
    }
    for t in mothers:
        m=t["mother"]
        lateral=t["lateral"]
        ap=good_seed_int(m,"ap.hd")>0
        lat=good_seed_int(m,"o.hds")>0
        zero=not (ap or lat)
        n=len(lateral)
        r["mothers_with_lateral_head_record"]+=int(n>0)
        r["mothers_without_scored_lateral_head"]+=int(n==0)
        r["mothers_with_two_plus_lateral_head_records"]+=int(n>=2)
        r["mothers_with_positive_lateral_good_seed"]+=int(lat)
        r["mothers_with_scored_lateral_but_zero_lateral_seed"]+=int(n>0 and not lat)
        r["apical_seed_positive_mothers"]+=int(ap)
        r["mothers_both_seed_routes_zero"]+=int(zero)
        r["joint_failure_mothers_with_scored_lateral_head"]+=int(zero and n>0)
        r["joint_failure_mothers_with_two_plus_scored_lateral_heads"]+=int(
            zero and n>=2)
        buds=value(m,"n.hds")
        r["joint_failure_mothers_with_3plus_source_buds"]+=int(
            zero and buds is not None and buds>=3)
        r["joint_failure_mothers_with_lateral_scored_good_seed_zero"]+=int(
            zero and n>0 and all(good_seed_int(h,"hd.gsd")==0 for h in lateral))
        for h in lateral:
            r["scored_lateral_head_rows"]+=1
            seed=good_seed_int(h,"hd.gsd")
            r["scored_lateral_good_seed_positive_head_rows"]+=int(seed>0)
            r["scored_lateral_good_seed_zero_head_rows"]+=int(seed==0)
            dam=score_damage_as_number(h["av.dam"])
            r["lateral_head_rows_damage_score_missing"]+=int(dam is None)
            if dam is None:continue
            r["lateral_head_rows_damage_ge_080"]+=int(dam>=.8)
            r["lateral_positive_seed_heads_damage_ge_080"]+=int(
                dam>=.8 and seed>0)
            r["lateral_head_rows_at_max_observed_damage_088"]+=int(
                math.isclose(dam,.88,rel_tol=0,abs_tol=1e-12))
            r["lateral_positive_seed_heads_at_max_damage_088"]+=int(
                math.isclose(dam,.88,rel_tol=0,abs_tol=1e-12) and seed>0)
            r["lateral_zero_seed_heads_at_max_damage_088"]+=int(
                math.isclose(dam,.88,rel_tol=0,abs_tol=1e-12) and seed==0)
    assert r["mothers_with_lateral_head_record"]+r["mothers_without_scored_lateral_head"]==r["mothers"]
    assert r["scored_lateral_good_seed_positive_head_rows"]+r["scored_lateral_good_seed_zero_head_rows"]==r["scored_lateral_head_rows"]
    return r


def summarize_within_mother_damage(scored_mothers):
    """Same mother having scored 0 and positive lateral-seed heads.

    Mean damage positive vs negative siblings, not a causal effect.
    Missing damage excludes only those heads and is counted separately.
    This report intentionally DOES NOT treat 111 heads as independent plants.
    """
    by_arm={}
    out={"n_2008_mothers_with_both_positive_and_zero_lateral_heads":0,
         "successful_lateral_heads_lower_damage_mothers":0,
         "successful_lateral_heads_higher_damage_mothers":0,
         "successful_lateral_heads_same_damage_mothers":0}
    for arm in ARMS:
        counts={"n":0,"lower":0,"higher":0,"same":0}
        for parent in (x for x in scored_mothers if x["arm"]==arm):
            available=[h for h in parent["lateral"]
                       if score_damage_as_number(h["av.dam"]) is not None]
            pos=[h for h in available if good_seed_int(h,"hd.gsd")>0]
            zero=[h for h in available if good_seed_int(h,"hd.gsd")==0]
            if not (pos and zero):
                continue
            def d(xs):
                return sum(score_damage_as_number(r["av.dam"]) for r in xs)/len(xs)
            delta=d(zero)-d(pos)
            counts["n"]+=1
            if delta>1e-10:counts["lower"]+=1
            elif delta < -1e-10:counts["higher"]+=1
            else:counts["same"]+=1
        by_arm[arm]=counts
        out["n_2008_mothers_with_both_positive_and_zero_lateral_heads"]+=counts["n"]
        out["successful_lateral_heads_lower_damage_mothers"]+=counts["lower"]
        out["successful_lateral_heads_higher_damage_mothers"]+=counts["higher"]
        out["successful_lateral_heads_same_damage_mothers"]+=counts["same"]
    out["by_treatment"]=by_arm
    out["one_signed_damage_pattern_in_all_mothers_inferred"]=False
    return out


def audit_scored_heads(mothers=None,heads=None,strict=True):
    if mothers is None or heads is None:
        ensure(mothers is None and heads is None,
               "MUST_SUPPLY_BOTH_ORIGINAL_TABLES")
        mothers,heads=source_locked_tables()
    rec,report2007=head_mother_crosswalk(mothers,heads,strict2008=strict)
    arms={arm:summarize_2008_group([r for r in rec if r["arm"]==arm])
          for arm in ARMS}
    failures=sum(t["mothers_both_seed_routes_zero"] for t in arms.values())
    maxdamage=[score_damage_as_number(h["av.dam"]) for h in heads
               if score_damage_as_number(h["av.dam"]) is not None]
    observed_max=max(maxdamage) if maxdamage else None
    if strict:
        ensure(observed_max==.88,
               "SOURCE_MAX_OBSERVED_DAMAGE_SCORE_CHANGED")
        ensure(failures==17 and
               arms["dni"]["joint_failure_mothers_with_scored_lateral_head"]==5,
               "SOURCE_2008_JOINT_FAILURE_WITH_LATERAL_HEAD_DRIFT")
        ensure(arms["ndni"]["joint_failure_mothers_with_scored_lateral_head"]==9,
               "SOURCE_2008_UNDAMAGED_FAILURE_WITH_LATERAL_HEAD_DRIFT")
    semantics={
        "source_mothers_with_zero_no_flw_but_positive_total_seed":
          sum(value(m,"no.flw")==0 and good_seed_int(m,"plt.gsd")>0
              for m in mothers),
        "source_2008_mothers_zero_no_flw_positive_seed":
          sum(m["year"]=="2008" and value(m,"no.flw")==0
              and good_seed_int(m,"plt.gsd")>0 for m in mothers),
        "no_flw_is_complete_verified_plant_head_anthesis_inventory":False,
        "score_hd_gsd_zero_guarantees_no_flower_opened":False,
    }
    return {
        "schema":"aza3_USDA_Cirsium_undulatum_scored_head_failure_v4",
        "source":{
            "original_USDA_doi":"10.15482/USDA.ADC/1522649",
            "original_paper_doi":"10.1002/ajb2.1744",
            "mother_source_sha1":BLOBS["mother_plants.csv"],
            "head_source_sha1":BLOBS["flower_heads.csv"],
            "original_mother_rows":len(mothers),
            "original_head_rows":len(heads),
            "original_2008_mothers":len(rec),
            "original_2008_head_rows":sum(h["year"]=="2008" for h in heads),
            "original_2008_lateral_head_rows":
                sum(len(t["lateral"]) for t in rec),
            "original_2008_apical_head_rows":len(rec),
            "mother_head_reconciliation_all_2008":True,
            "max_recorded_head_damage_score":observed_max,
            "original_source_and_prior_compensation_are_not_new_studies":True,
        },
        "year2008_by_original_treatment":arms,
        "within_mother_lateral_head_damage_comparison":
            summarize_within_mother_damage(rec),
        "other_year_archive_quality":report2007,
        "source_semantics":semantics,
        "scientific_claim_limits":{
            "scored_lateral_head_means_confirmed_anthesis":False,
            "damage_score_deterministically_classifies_reproductive_outcome":False,
            "within_mother_damage_association_is_causal":False,
            "one_insect_antagonist_proven_joint_failure_agent":False,
            "apical_damage_means_apex_was_removed_completely":False,
            "structural_head_availability_implies_pollen_service":False,
            "new_floral_SLK_division_of_labour_payoff_estimated":False,
            "genetic_ANTO_ancestral_functional_re_expression_estimated":False,
            "head_records_independent_genotypes":False,
            "insecticide_selectively_affects_seed_predators_only":False,
            "two_original_years_share_the_same_maternal_genotypes":False,
            "good_seed_equals_established_seedling":False,
            "observed_damage_correlation_establishes_pre_exposure_cause":False,
        },
        "conclusion":
            "REAL_LATERAL_HEADS_RECORDED_DURING_JOINT_SEED_FAILURE_CAUSE_UNIDENTIFIED",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    report=audit_scored_heads(strict=True)
    t=json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n"
    if args.output is None:print(t)
    else:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(t,encoding="utf-8")


if __name__=="__main__":
    main()
