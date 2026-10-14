#!/usr/bin/env python3
"""Reproducible ingestion of REAL Cirsium canescens capitulum experiment data.

SOURCE: Russell, Taylor & Louda (2025) Dryad doi:10.5061/dryad.kh18932dc,
version 3, Plattethistleseedproduction.csv, SHA256 pinned below.
The repository does NOT contain raw Dryad files and does NOT possess an
API download token. Without an authenticated, hash-verified input file,
return SOURCE_NOT_MATERIALIZED, not an empty biological effect.

Design constraints:
* HEAD is nested in randomized PLANT; FnlTrt is normally recorded once
  per plant. Never propagate treatment across plant boundaries.
* HdDam 5/6 represent DIFFERENT damaged organs/entry paths, not
  quantitative severity 5/6; stage Stg recorded at a FINAL census.
* 2018 Tanglefoot killed heads unintentionally; EstSdProdTang is the
  author's corrected estimated seed count for 68 heads.
* No true spine length, calibrated gravitational inclination, filmed
  physical access or independently verified pollen deposition is present.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

DATASET_DOI="10.5061/dryad.kh18932dc"
RAW_FILENAME="Plattethistleseedproduction.csv"
EXPECTED_SHA256="d8c584b0c67a80f7e637e95df2262c7cb16788b0bd09d7e6e4ae702e00508b50"
EXPECTED_BYTES=348096
DRYAD_VERSION=3
NEEDED={
    "Project","Year","Site","Plant","FnlTrt","HdN","Stg","HdDam",
    "RcEg","FinPRGD","EstSdProd","EstSdProdTang","TPltSDs"
}
MISS={"","NA","N/A","NULL","null","None","not_assessed"}
HEAD_DAMAGE={
    0:"externally_unmarked",
    1:"ambiguous_external_damage",
    2:"phyllaries_lt_5_percent_damaged",
    3:"phyllaries_6_to_30_percent_damaged",
    4:"phyllaries_gt_30_percent_damaged",
    5:"stem_mining_near_capitulum",
    6:"chewed_hole_into_capitulum",
}

def missing(x):
    return x is None or str(x).strip() in MISS

def count(x,name,allow_missing=True):
    if missing(x):
        if allow_missing:
            return None
        raise ValueError("MISSING_"+name)
    try:f=float(str(x).strip())
    except ValueError:raise ValueError("NOT_NUMERIC_"+name) from None
    if not math.isfinite(f) or not f.is_integer() or f<0:
        raise ValueError("INVALID_NONNEGATIVE_INTEGER_"+name)
    return int(f)

def nonnegative(x,name):
    if missing(x):return None
    try:f=float(str(x).strip())
    except ValueError:raise ValueError("NOT_NUMERIC_"+name) from None
    if not math.isfinite(f) or f<0:
        raise ValueError("INVALID_NONNEGATIVE_"+name)
    return f

def fraction(x,name):
    f=nonnegative(x,name)
    if f is not None and f>1:
        raise ValueError("NOT_PROPORTION_"+name)
    return f

def head_position(hdn):
    """Architecture rank, NOT real orientation; preserve unknown code."""
    if hdn==100:return "main_apical"
    if hdn>=200 and hdn%100==0:return "lateral_terminal"
    if hdn>=200 and hdn%100>0:return "lateral_axillary_or_nonterminal"
    return "other_or_undocumented_head_number"

def read_verified_csv(path:Path):
    data=path.read_bytes()
    digest=hashlib.sha256(data).hexdigest()
    if len(data)!=EXPECTED_BYTES or digest!=EXPECTED_SHA256:
        raise ValueError("SOURCE_SHA256_OR_LENGTH_NOT_FROZEN_DRYAD_V3")
    try:
        decoded=data.decode("utf-8-sig")
    except UnicodeError:
        raise ValueError("SOURCE_UTF8_DECODING_REQUIRES_INSPECTION") from None
    r=csv.DictReader(io.StringIO(decoded))
    if not NEEDED.issubset(r.fieldnames or ()):
        raise ValueError("SOURCE_COLUMNS_NOT_DOCUMENTED_AS_PUBLISHED")
    return list(r),digest

def analyze(rows):
    if not rows:
        return {"status":"NO_ARCHIVE_HEAD_ROWS_NOT_A_NEGATIVE_ECOLOGICAL_RESULT",
                "n_heads":0,"n_plants":0,"source_results":[]}
    plant_buckets=defaultdict(list)
    seen_heads=set()
    for lineno,row in enumerate(rows,2):
        if not NEEDED.issubset(row):
            raise ValueError("ARCHIVE_MISSING_FIELD_ROW")
        plantkey=tuple((row.get(k) or "").strip() for k in
                       ("Project","Year","Site","Plant"))
        if not all(plantkey):
            raise ValueError("MISSING_PROJECT_YEAR_SITE_OR_PLANT_ID")
        if plantkey[1] not in {"2017","2018","2019"}:
            raise ValueError("UNEXPECTED_STUDY_YEAR")
        hdn=count(row["HdN"],"HDN",allow_missing=False)
        headkey=plantkey+(hdn,)
        if headkey in seen_heads:
            raise ValueError("DUPLICATED_CAPITULUM_ID_WITHIN_PLANT")
        seen_heads.add(headkey)
        plant_buckets[plantkey].append(row)

    plants=[]
    n_corrections=0
    n_corrections_non2018=0
    damaged_modes=Counter()
    stage_counts=Counter()
    for plantkey,head_rows in sorted(plant_buckets.items()):
        nonempty_treatment={
            r["FnlTrt"].strip() for r in head_rows if not missing(r["FnlTrt"])
        }
        if nonempty_treatment not in ({"I"},{"C"}):
            raise ValueError("MISSING_OR_CONFLICTING_PLANT_RANDOMIZED_TREATMENT")
        treatment=next(iter(nonempty_treatment))
        plant_seed_corrected=[]
        plant_seed_uncorrected=[]
        n_heads_missing_corrected=0
        n_eggs_known=n_egg_events=0
        n_ext_damage_known=n_explicit_damage=0
        n_stage_known=0
        known_floret_viability=0
        head_descriptives=[]
        original_totals=[]
        for row in head_rows:
            hdn=count(row["HdN"],"HDN",allow_missing=False)
            head_stage=count(row["Stg"],"FINAL_AUGUST_STAGE")
            if head_stage is not None and head_stage not in range(1,7):
                raise ValueError("INVALID_FINAL_AUGUST_HEAD_STAGE")
            if head_stage is not None:
                n_stage_known+=1
                stage_counts[head_stage]+=1
            damage_code=count(row["HdDam"],"HD_DAMAGE")
            if damage_code is not None and damage_code not in HEAD_DAMAGE:
                raise ValueError("UNSUPPORTED_NONORDINAL_EXTERNAL_HEAD_DAMAGE_CODE")
            if damage_code is not None:
                damaged_modes[HEAD_DAMAGE[damage_code]]+=1
                if damage_code!=1:
                    n_ext_damage_known+=1
                    n_explicit_damage+=int(damage_code>=2)
            eggs=count(row["RcEg"],"RHINOCYLLUS_OVIPOSITION_COUNT")
            if eggs is not None:
                n_eggs_known+=1
                n_egg_events+=eggs
            viable=fraction(row["FinPRGD"],"ESTIMATED_VIABLE_PROPORTION")
            known_floret_viability+=int(viable is not None)
            raw_seed=nonnegative(row["EstSdProd"],"UNCORRECTED_SEED_ESTIMATE")
            fixed_seed=nonnegative(row["EstSdProdTang"],"TANGLEFOOT_CORRECTED_SEED_ESTIMATE")
            if fixed_seed is None:
                n_heads_missing_corrected+=1
            else:
                plant_seed_corrected.append(fixed_seed)
            if raw_seed is not None:
                plant_seed_uncorrected.append(raw_seed)
            if raw_seed is not None and fixed_seed is not None and abs(raw_seed-fixed_seed)>1e-7:
                n_corrections+=1
                if plantkey[1]!="2018":
                    n_corrections_non2018+=1
            published_total=nonnegative(row["TPltSDs"],"PUBLISHED_PLANT_SEED_TOTAL")
            if published_total is not None:
                original_totals.append(published_total)
            head_descriptives.append({
                "position":head_position(hdn),
                "final_census_stage":head_stage,
                "external_damage_mechanism":HEAD_DAMAGE.get(damage_code),
                "rhinocyllus_eggs_assessed":eggs is not None,
                "estimated_viable_fraction_assessed":viable is not None,
                "corrected_seed_estimate_assessed":fixed_seed is not None,
            })
        if original_totals and max(original_totals)-min(original_totals)>1e-5:
            raise ValueError("CONFLICTING_PUBLISHED_PLANT_TOTAL_WITHIN_PLANT")
        full=(n_heads_missing_corrected==0)
        source_sum=sum(plant_seed_corrected) if full else None
        published=original_totals[0] if original_totals else None
        total_diff=(source_sum-published
                    if source_sum is not None and published is not None
                    else None)
        plants.append({
            "project":plantkey[0],"year":plantkey[1],"site":plantkey[2],
            "plant":plantkey[3],"treatment":treatment,
            "n_heads":len(head_rows),
            "n_missing_corrected_seed_estimates":n_heads_missing_corrected,
            "plant_seed_estimate_corrected_all_heads":source_sum,
            "published_plant_total_if_assessed":published,
            "plant_total_discrepancy":total_diff,
            "n_heads_with_external_damage_assessed":n_ext_damage_known,
            "n_heads_with_obvious_external_damage":n_explicit_damage,
            "n_heads_with_weevil_eggs_assessed":n_eggs_known,
            "number_recorded_weevil_ovipositions":n_egg_events,
            "n_heads_with_viability_fraction_assessed":known_floret_viability,
            "n_heads_with_final_stage_assessed":n_stage_known,
            "head_position_counts":dict(Counter(x["position"] for x in head_descriptives)),
        })

    blocks=defaultdict(lambda:{"I":[],"C":[]})
    for p in plants:
        block=(p["project"],p["year"],p["site"])
        blocks[block][p["treatment"]].append(p)

    contrasts=[]
    for (project,year,site),treats in sorted(blocks.items()):
        group={}
        for trt in ("I","C"):
            pp=treats[trt]
            complete=[p for p in pp if p["plant_seed_estimate_corrected_all_heads"] is not None]
            vals=[p["plant_seed_estimate_corrected_all_heads"] for p in complete]
            group[trt]={
                "n_randomized_plants":len(pp),
                "n_complete_plant_corrected_seed_totals":len(complete),
                "mean_corrected_plant_seed_output_complete_cases":(
                    sum(vals)/len(vals) if vals else None),
                "n_heads_with_external_damage_assessed":sum(
                    p["n_heads_with_external_damage_assessed"] for p in pp),
                "n_heads_with_external_damage":sum(
                    p["n_heads_with_obvious_external_damage"] for p in pp),
                "n_heads_with_weevil_ovipositions_assessed":sum(
                    p["n_heads_with_weevil_eggs_assessed"] for p in pp),
                "sum_weevil_ovipositions_if_assessed":sum(
                    p["number_recorded_weevil_ovipositions"] for p in pp),
            }
        a=group["I"]["mean_corrected_plant_seed_output_complete_cases"]
        b=group["C"]["mean_corrected_plant_seed_output_complete_cases"]
        contrasts.append({
            "project":project,"year":year,"site":site,
            "insecticide":group["I"],"control":group["C"],
            "descriptive_plant_difference_I_minus_C":(
                a-b if a is not None and b is not None else None),
            "NO_head_as_independent_treatment_replicate":True,
            "NO_spine_or_orientation_trait_effect_identified":True,
            "complete_case_selection_may_bias_comparison":True,
        })

    return {
        "status":"ARCHIVAL_HEAD_AND_PLANT_DATA_VALIDATED_DESCRIPTIVE_ONLY",
        "n_raw_head_rows":len(rows),
        "n_unique_randomized_plants":len(plants),
        "n_project_year_site_blocks":len(contrasts),
        "n_plants_with_complete_corrected_seed_estimates":sum(
            p["plant_seed_estimate_corrected_all_heads"] is not None for p in plants),
        "n_heads_with_2018_tanglefoot_corrections":n_corrections,
        "n_seed_estimate_corrections_outside_2018":n_corrections_non2018,
        "published_readme_says_68_corrections_check":(
            "MATCH_68" if n_corrections==68 and not n_corrections_non2018
            else "AUDIT_REQUIRED_NOT_EQUIVALENT_TO_68_README_HEADS"),
        "external_head_damage_categorical_counts":dict(damaged_modes),
        "final_census_stage_counts_not_baseline":dict(stage_counts),
        "descriptive_site_year_project_contrasts":contrasts,
        "observational_limit":[
            "Insecticide was randomized at plant level: do not pseudoreplicate heads.",
            "Corrected seed estimate differs in 2018 because Tanglefoot treatment damaged heads.",
            "Post-treatment head stage and actual damage are not baseline covariates.",
            "HdDam 5 and 6 are distinct physical damage mechanisms not larger ordinal scores.",
            "RcEg counts one invasive biocontrol weevil, not all seed-feeder guilds.",
            "HdN is branch/head position, not gravity-oriented capitulum angle.",
            "No original botanical spine length or involucre access-gaps.",
            "Raw archived species effect is not an orientation x spine adaptive effect.",
        ],
    }

def synthetic_tests():
    def head(plant,hdn,trt,damage,eggs,raw,corrected,year="2018"):
        return {
            "Project":"2","Year":year,"Site":"SITE_A","Plant":plant,
            "FnlTrt":trt,"HdN":str(hdn),"Stg":"5",
            "HdDam":str(damage),"RcEg":str(eggs),"FinPRGD":"0.5",
            "EstSdProd":str(raw),"EstSdProdTang":str(corrected),
            "TPltSDs":"25" if plant=="P1" and hdn==100 else
                       "20" if plant=="P2" and hdn==100 else "NA",
        }
    rows=[
        head("P1",100,"I",2,0,5,10),
        head("P1",200,"NA",6,1,15,15),
        head("P2",100,"C",0,0,10,10),
        head("P2",210,"NA",1,0,10,10)
    ]
    res=analyze(rows)
    assert res["n_unique_randomized_plants"]==2
    assert res["n_raw_head_rows"]==4
    assert res["n_heads_with_2018_tanglefoot_corrections"]==1
    assert res["external_head_damage_categorical_counts"]["chewed_hole_into_capitulum"]==1
    assert res["external_head_damage_categorical_counts"]["ambiguous_external_damage"]==1
    group=res["descriptive_site_year_project_contrasts"][0]
    assert group["descriptive_plant_difference_I_minus_C"]==5
    assert group["insecticide"]["n_randomized_plants"]==1
    assert group["control"]["n_randomized_plants"]==1
    assert head_position(100)=="main_apical"
    assert head_position(200)=="lateral_terminal"
    assert head_position(210)=="lateral_axillary_or_nonterminal"
    for corrupt,needle in [
        (rows+[rows[0]],"DUPLICATED_CAPITULUM_ID_WITHIN_PLANT"),
        ([{**rows[0],"FnlTrt":"NA"},rows[1],rows[2],rows[3]],
            "MISSING_OR_CONFLICTING_PLANT_RANDOMIZED_TREATMENT"),
        ([rows[0],{**rows[1],"FnlTrt":"C"},rows[2],rows[3]],
            "MISSING_OR_CONFLICTING_PLANT_RANDOMIZED_TREATMENT"),
        ([{**rows[0],"HdDam":"7"},rows[1],rows[2],rows[3]],
            "UNSUPPORTED_NONORDINAL_EXTERNAL_HEAD_DAMAGE_CODE"),
        ([{**rows[0],"RcEg":"-1"},rows[1],rows[2],rows[3]],
            "INVALID_NONNEGATIVE_INTEGER_RHINOCYLLUS_OVIPOSITION_COUNT"),
        ([{**rows[0],"FinPRGD":"1.3"},rows[1],rows[2],rows[3]],
            "NOT_PROPORTION_ESTIMATED_VIABLE_PROPORTION"),
    ]:
        try:analyze(corrupt)
        except ValueError as err:
            assert needle in str(err),(needle,str(err))
        else:
            raise AssertionError("Archive source corruption was accepted: "+needle)
    missing_rows=[rows[0],{**rows[1],"EstSdProdTang":"NA"},rows[2],rows[3]]
    hold=analyze(missing_rows)
    assert hold["n_plants_with_complete_corrected_seed_estimates"]==1
    assert hold["descriptive_site_year_project_contrasts"][0]["descriptive_plant_difference_I_minus_C"] is None
    assert analyze([])["status"]=="NO_ARCHIVE_HEAD_ROWS_NOT_A_NEGATIVE_ECOLOGICAL_RESULT"
    return {
        "status":"PASS_SYNTHETIC_CIRSIUM_ARCHIVE_DESIGN_AND_PROVENANCE_GATES",
        "n_synthetic_head_rows":len(rows),
        "n_negative_controls":6,
        "tanglefoot_correction_flagged":True,
        "head_position_codes_resolved":True,
        "treatment_propagated_only_within_plant":True,
        "real_archival_rows_analyzed":0,
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--source-csv",type=Path)
    a.add_argument("--out",type=Path)
    args=a.parse_args()
    tests=synthetic_tests()
    source={
        "dataset_doi":DATASET_DOI,
        "dryad_version":DRYAD_VERSION,
        "filename":RAW_FILENAME,
        "expected_sha256":EXPECTED_SHA256,
        "expected_bytes":EXPECTED_BYTES,
        "file_download_API_requires_authentication":True,
    }
    if args.source_csv is None:
        result={
            "status":"SOURCE_CSV_NOT_MATERIALIZED_NO_EMPIRICAL_RECOMPUTATION",
            "n_actual_archival_head_rows_analyzed":0,
            "authentic_archive_conclusion":"PUBLISHED_2025_STUDY_ONLY",
            "archive_fitness_data":None,
        }
    else:
        original,digest=read_verified_csv(args.source_csv)
        result={
            "status":"HASH_VERIFIED_ARCHIVAL_DESCRIPTIVE_RESULTS_NOT_SPINE_ADAPTATION",
            "n_actual_archival_head_rows_analyzed":len(original),
            "source_sha256":digest,
            "archive_fitness_data":analyze(original),
        }
    output={
        "version":"audit_real_cirsium_canescens_head_archive_v1",
        "as_of":"2026-10-14",
        "source_manifest":source,
        "synthetic_source_integrity_tests":tests,
        "actual_source_analysis":result,
    }
    if args.out is not None:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(output,ensure_ascii=False,indent=2)+"\n",
                            encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "synthetic_integrity_status":tests["status"],
        "real_Cirsium_head_rows_recomputed":result["n_actual_archival_head_rows_analyzed"],
    },ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
