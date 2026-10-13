#!/usr/bin/env python3
"""Observed head-level seed-predator incidence versus severity decomposition.

Uses the ORIGINAL EAzami Aim2 head registry without adding or imputing any
biological observations. Separate:
 * P(seed-predator presence | preregistered head),
 * mean damage fraction | seed-predator present,
 * head-weighted fraction damaged across enrolled heads,
 * achene-weighted fraction damaged across enrolled heads,
 * independently measured filled-achenes/head.

An unassessed uninfested head is NOT automatically assigned zero damage.
The audited fractions are as recorded; they do not independently establish
larval species, stigma pollen receipt, germination or causal true-spine action.

Related background:
Zhang et al. 2018 DOI 10.1086/698766 (34 Asteraceae, opposite trends
in incidence and conditional severity); Xi et al. 2021 DOI
10.1111/1365-2656.13480 (24 Asteraceae, pollination can outweigh seed loss).
These papers do not identify a Cirsium orientation x armature causal effect.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

KEYS=("population_id","individual_id","capitulum_id")
NEEDED=set(KEYS)|{
    "taxon","seed_predator_presence","seed_predator_damage_fraction",
    "total_achenes","filled_achenes","capitulum_diameter_mm",
    "spine_length_mm","natural_orientation_deg","orientation_treatment",
}
MISSING={"","NA","NULL","not_assessed","unknown"}
N_MIN_PLANTS=3

def read_file(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        return list(reader.fieldnames or []),list(reader)

def integer_count(value,name, positive=False):
    if value is None or value.strip() in MISSING:
        return None
    try:v=float(value)
    except ValueError:raise ValueError(f"{name}_NOT_A_COUNT") from None
    if not math.isfinite(v) or not v.is_integer() or v<0 or (positive and v==0):
        raise ValueError(f"{name}_INVALID_NONNEGATIVE_INTEGER")
    return int(v)

def measured_fraction(value,name):
    if value is None or value.strip() in MISSING:
        return None
    try:v=float(value)
    except ValueError:raise ValueError(f"{name}_NOT_NUMERIC") from None
    if not math.isfinite(v) or not 0<=v<=1:
        raise ValueError(f"{name}_OUTSIDE_UNIT_INTERVAL")
    return v

def measured_trait(value,name):
    if value is None or value.strip() in MISSING:
        return None
    try:v=float(value)
    except ValueError:raise ValueError(f"{name}_NOT_NUMERIC") from None
    if not math.isfinite(v) or v<0:
        raise ValueError(f"{name}_INVALID_VALUE")
    return v

def summarize_complete_cases(rows):
    """Identity checks hold at BOTH head and achene unit of aggregation."""
    n=len(rows)
    if n==0:
        return None
    infected=[r for r in rows if r["infested"]==1]
    p=sum(r["infested"] for r in rows)/n
    severity=sum(r["fraction"] for r in infected)/len(infected) if infected else None
    head_loss=sum(r["fraction"] for r in rows)/n
    total=sum(r["total"] for r in rows)
    selected_total=sum(r["total"] for r in infected)
    selected_damage=sum(r["total"]*r["fraction"] for r in infected)
    achene_loss=sum(r["total"]*r["fraction"] for r in rows)/total
    achene_p=selected_total/total
    achene_severity=(selected_damage/selected_total
                     if selected_total>0 else None)
    if infected:
        assert abs(p*severity-head_loss)<1e-10
        assert abs(achene_p*achene_severity-achene_loss)<1e-10
    else:
        assert p==head_loss==achene_loss==0
    n_plants=len({(r["population"],r["plant"]) for r in rows})
    has_controlled_orientation=sum(r["treated"] for r in rows)
    n_natural_traits=sum(all(r[v] is not None for v in
                              ("orientation","spine","head_diameter"))
                         for r in rows)
    return {
        "n_heads_complete":n,
        "n_independent_plants":n_plants,
        "n_seed_predator_present":len(infected),
        "n_seed_predator_absent":n-len(infected),
        "head_weighted_infested_fraction":p,
        "mean_fraction_damaged_if_infested":severity,
        "head_weighted_fraction_damaged_all_heads":head_loss,
        "achene_weighted_infested_head_exposure_fraction":achene_p,
        "achene_weighted_damage_if_infested":achene_severity,
        "achene_weighted_fraction_damaged_all_achenes":achene_loss,
        "head_versus_achene_weighted_loss_difference":achene_loss-head_loss,
        "mean_filled_achenes_per_head":sum(r["filled"] for r in rows)/n,
        "pooled_filled_achene_fraction":sum(r["filled"] for r in rows)/total,
        "n_heads_with_measured_orientation_spine_size":n_natural_traits,
        "n_heads_with_nonempty_orientation_treatment_flag":has_controlled_orientation,
        "head_incidence_x_conditional_severity_identity_confirmed":True,
        "achene_incidence_x_conditional_severity_identity_confirmed":True,
        "cannot_claim":"CORRELATION_OR_CAUSAL_FITNESS_OR_ADAPTATION_FROM_SUMMARY",
    }

def audit_rows(rows):
    if not rows:
        return {
            "status":"NO_REAL_PREREGISTERED_HEAD_DATA",
            "n_source_heads":0,
            "n_complete_head_incidence_severity":0,
            "n_missing_or_unassessed_head_outcomes":0,
            "summaries":[],
            "biological_result":"NOT_IDENTIFIABLE",
        }
    seen=set()
    subsets=defaultdict(list)
    omissions=defaultdict(int)
    for i,row in enumerate(rows,1):
        if not NEEDED.issubset(row):
            raise ValueError("ORIGINAL_HEAD_SCHEMA_MISSING_COLUMNS")
        key=tuple((row.get(k) or "").strip() for k in KEYS)
        if not all(key) or key in seen:
            raise ValueError("EMPTY_OR_DUPLICATE_REGISTERED_HEAD")
        seen.add(key)
        taxon=(row.get("taxon") or "").strip()
        if not taxon:
            raise ValueError("REGISTERED_HEAD_UNRESOLVED_TAXON")
        group=(key[0],taxon)
        flag=(row.get("seed_predator_presence") or "").strip()
        if flag not in {"0","1"} and flag not in MISSING:
            raise ValueError("SEED_PREDATOR_PRESENCE_INVALID")
        infested=int(flag) if flag in {"0","1"} else None
        fraction=measured_fraction(row.get("seed_predator_damage_fraction"),
                                  "SEED_PREDATOR_DAMAGE_FRACTION")
        total=integer_count(row.get("total_achenes"),"TOTAL_ACHENES",positive=True)
        filled=integer_count(row.get("filled_achenes"),"FILLED_ACHENES")
        if total is not None and filled is not None and filled>total:
            raise ValueError("FILLED_ACHENES_EXCEED_TOTAL")
        if infested==0 and fraction is not None and fraction>0:
            raise ValueError("DAMAGE_FRACTION_CONTRADICTS_NO_SEED_PREDATOR")
        orientation=measured_trait(row.get("natural_orientation_deg"),
                                   "ORIENTATION_DEG")
        spine=measured_trait(row.get("spine_length_mm"),"TRUE_SPINE_MM")
        size=measured_trait(row.get("capitulum_diameter_mm"),
                            "HEAD_DIAMETER_MM")
        treatment=(row.get("orientation_treatment") or "").strip()
        # Any filled/count value without the other counts is still a partial
        # record. Preserve it but do not create a complete loss/fitness rate.
        if any(v is None for v in (infested,fraction,total,filled)):
            omissions[group]+=1
            continue
        subsets[group].append({
            "population":key[0],"plant":key[1],"head":key[2],
            "infested":infested,"fraction":fraction,
            "total":total,"filled":filled,
            "orientation":orientation,"spine":spine,"head_diameter":size,
            "treated":treatment not in MISSING,
        })
    outputs=[]
    for group in sorted(set(subsets)|set(omissions)):
        items=subsets[group]
        summary=summarize_complete_cases(items)
        outputs.append({
            "population_id":group[0],"taxon":group[1],
            "n_missing_or_unassessed_head_outcomes":omissions[group],
            "complete_case_summary":summary,
            "status":("NO_COMPLETE_HEADS" if summary is None
                      else "DESCRIPTIVE_LOW_BIOLOGICAL_REPLICATION"
                      if summary["n_independent_plants"]<N_MIN_PLANTS
                      else "DESCRIPTIVE_COMPLETE_CASES_NOT_CAUSAL"),
            "selection_warning":(
                "Incomplete head endpoints cannot be recoded as intact or "
                "undamaged, and complete cases may be nonrandomly selected."
            ),
        })
    n_complete=sum(len(v) for v in subsets.values())
    n_missing=sum(omissions.values())
    return {
        "status":("HOLD_NO_COMPLETE_HEAD_INCIDENCE_SEVERITY_OUTCOMES"
                  if n_complete==0 else
                  "DESCRIPTIVE_COMPLETE_CASES_SELECTION_UNRESOLVED"),
        "n_source_heads":len(rows),
        "n_complete_head_incidence_severity":n_complete,
        "n_missing_or_unassessed_head_outcomes":n_missing,
        "n_population_by_taxon_strata":len(outputs),
        "summaries":outputs,
        "biological_result":"DESCRIPTIVE_HEAD_LEVEL_ONLY_NO_STRUCTURAL_SELECTION",
        "missingness":"SOURCE_NA_IS_NOT_ZERO",
        "no_post_infested_head_selection":True,
    }

def synthetic_tests():
    def head(n,presence,damage,filled,plant,diameter):
        return {
            "individual_id":plant,"population_id":"SYNT_POP",
            "capitulum_id":"H_"+plant,"taxon":"SYNT_CIRSIUM",
            "seed_predator_presence":presence,
            "seed_predator_damage_fraction":damage,
            "total_achenes":str(n),"filled_achenes":str(filled),
            "capitulum_diameter_mm":str(diameter),
            "spine_length_mm":"3","natural_orientation_deg":"70",
            "orientation_treatment":"",
        }
    h0=head(10,"0","0",8,"P0",9)
    h1=head(10,"1","0.8",2,"P1",9)
    h2=head(100,"1","0.1",60,"P2",19)
    r=audit_rows([h0,h1,h2])
    assert r["status"]=="DESCRIPTIVE_COMPLETE_CASES_SELECTION_UNRESOLVED"
    q=r["summaries"][0]["complete_case_summary"]
    assert abs(q["head_weighted_infested_fraction"]-2/3)<1e-12
    assert abs(q["mean_fraction_damaged_if_infested"]-.45)<1e-12
    assert abs(q["head_weighted_fraction_damaged_all_heads"]-.3)<1e-12
    assert abs(q["achene_weighted_infested_head_exposure_fraction"]-110/120)<1e-12
    assert abs(q["achene_weighted_damage_if_infested"]-18/110)<1e-12
    assert abs(q["achene_weighted_fraction_damaged_all_achenes"]-.15)<1e-12
    assert q["n_independent_plants"]==3
    assert audit_rows([])["status"]=="NO_REAL_PREREGISTERED_HEAD_DATA"
    unmeasured=head(20,"0","",15,"P3",9)
    o=audit_rows([h0,unmeasured])
    assert o["n_missing_or_unassessed_head_outcomes"]==1
    assert o["n_complete_head_incidence_severity"]==1
    assert abs(o["summaries"][0]["complete_case_summary"]["head_weighted_fraction_damaged_all_heads"])<1e-12
    assert o["summaries"][0]["status"]=="DESCRIPTIVE_LOW_BIOLOGICAL_REPLICATION"
    no_complete=audit_rows([unmeasured])
    assert no_complete["status"]=="HOLD_NO_COMPLETE_HEAD_INCIDENCE_SEVERITY_OUTCOMES"
    cases=[
        ([h0,h0],"EMPTY_OR_DUPLICATE_REGISTERED_HEAD"),
        ([{**h0,"seed_predator_presence":"1.0"}],"SEED_PREDATOR_PRESENCE_INVALID"),
        ([{**h0,"seed_predator_damage_fraction":"0.1"}],
         "DAMAGE_FRACTION_CONTRADICTS_NO_SEED_PREDATOR"),
        ([{**h0,"total_achenes":"0"}],"TOTAL_ACHENES_INVALID_NONNEGATIVE_INTEGER"),
        ([{**h0,"filled_achenes":"11"}],"FILLED_ACHENES_EXCEED_TOTAL"),
        ([{**h0,"seed_predator_damage_fraction":"nan"}],
         "SEED_PREDATOR_DAMAGE_FRACTION_OUTSIDE_UNIT_INTERVAL"),
    ]
    for input_rows,needle in cases:
        try:audit_rows(input_rows)
        except ValueError as exc:
            assert needle in str(exc),(needle,str(exc))
        else:raise AssertionError("INVALID_CIRSIUM_HEAD_LOSS_DATA_WAS_ACCEPTED")
    return {
        "status":"PASS_SYNTHETIC_HEAD_WEIGHTED_VS_ACHENE_WEIGHTED_HURDLE",
        "head_weighted_loss":q["head_weighted_fraction_damaged_all_heads"],
        "achene_weighted_loss":q["achene_weighted_fraction_damaged_all_achenes"],
        "n_negative_guards":len(cases),
        "n_actual_biological_records":0,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--heads",type=Path,required=True)
    ap.add_argument("--out",type=Path)
    args=ap.parse_args()
    header,rows=read_file(args.heads)
    if not NEEDED.issubset(header):
        raise ValueError("ORIGINAL_EAZAMI_HEAD_COLUMNS_MISSING")
    payload={
        "version":"capitulum_incidence_conditional_severity_v1",
        "date":"2026-10-13",
        "synthetic_verification":synthetic_tests(),
        "source_audit":audit_rows(rows),
        "data_definition":{
            "prevalence":"Fraction of all registered heads with explicitly scored infestation",
            "severity":"Mean fraction of damaged achenes only among infested heads",
            "head_weighted":"Each registered head contributes equal unit weight",
            "achene_weighted":"Each recorded achene contributes equal weight via observed total_achenes",
            "filled_achenes":"Actual reported filled-head proxy not independent germination",
            "missing":"Unassessed damage fraction never converted to zero",
        },
        "source":"Original EAzami Aim2 head field ledger, no new field collection",
        "inferential_ceiling":"DESCRIPTIVE_ONLY_NO_TRAIT_CAUSAL_OR_ADAPTIVE_INFERENCE",
    }
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",
                            encoding="utf-8")
    print(json.dumps(payload,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
