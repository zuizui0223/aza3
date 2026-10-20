#!/usr/bin/env python3
"""Source-faithful Cirsium arvense 2022 four-arm dual outcome + insect-file audit.

Original seed heads are nested in site x experimental block x treatment plant.
Insects in the separately identified Herbivore.csv are NOT assumed to match
sampled heads by vial, even if site-block-treatment totals coincidentally match.
Emitted stage contrasts are descriptive of source treatment assignments, NOT
orientation x real-spine fitness, florivore feeding or pollination mediation.

Original Karpate et al 2025 Research Square preprint, 2022 field experiment:
DOI 10.21203/rs.3.rs-6634172/v1; author archive DOI 10.5281/zenodo.15364152.
"""
from __future__ import annotations
from collections import Counter, defaultdict
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

from audit_cirsium_2022_four_arm_head_outcomes_v1 import (
    SOURCE_REPO, SOURCE_COMMIT, SOURCE_GIT_BLOB, SITES, ARMS,
    parse_rows, prep, contrast, quantile, check, mean, sign_counts,
)

IDENTIFIED_BLOB="b0578e8c814b371e910da119c66e45575bc8b675"
IDENTIFIED_REQUIRED={"Site_Code","Block_Code","Treatment_Code","Vial_Code",
                     "Feeding_Guild","Number"}
GUILDS=("Herbivore","Parasitoid")
SEED_METRICS=("viable_fraction","viable_per_head","emerged_insects_per_head",
              "observed_damaged_fraction","inferred_eaten_signed_fraction")
PAIRS=(("Pollinator Exclusion","Control"),
       ("Predator Exclusion","Control"),
       ("Predator Model","Predator Exclusion"),
       ("Predator Model","Control"))
SEED=20261020


def git_blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()


def int_value(value,field):
    try:
        v=float(value)
    except (ValueError,TypeError):
        raise ValueError("INVALID_INSECT_NUMBER_"+field) from None
    if not math.isfinite(v) or not v.is_integer() or v<0:
        raise ValueError("INVALID_INSECT_NUMBER_"+field)
    return int(v)


def plant_key(row):
    """Join ONLY at source site×block×treatment; NEVER assume vial is a head id."""
    site=(row.get("Site_Code") or "").strip().lower()
    block=(row.get("Block_Code") or "").strip().upper()
    arm=(row.get("Treatment_Code") or "").strip()
    if not site or not block or not arm:
        raise ValueError("MISSING_PLANT_JOIN_KEY")
    return site,block,arm


def parse_identification_rows(raw: bytes,strict=True):
    sha=git_blob(raw)
    if strict:
        check(sha==IDENTIFIED_BLOB,"ORIGINAL_HERBIVORE_CSV_BLOB_CHANGED")
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    check(IDENTIFIED_REQUIRED.issubset(reader.fieldnames or ()),
          "IDENTIFICATION_SOURCE_SCHEMA_DRIFT")
    rows=list(reader)
    if strict:check(len(rows)==181,"IDENTIFIED_ROWS_NOT_181")
    for row in rows:
        row["_number"]=int_value(row.get("Number"),"COUNT")
        if row.get("Feeding_Guild") not in GUILDS:
            raise ValueError("UNEXPECTED_GUILD")
        plant_key(row)
    if strict:
        check(Counter(r["Feeding_Guild"] for r in rows)=={
              "Herbivore":168,"Parasitoid":13},"GUILD_ROW_CLASS_DRIFT")
        check(Counter({g:sum(r["_number"] for r in rows
                          if r["Feeding_Guild"]==g)
                       for g in GUILDS})==
              {"Herbivore":498,"Parasitoid":24},"GUILD_NUMBER_DRIFT")
    return rows,sha


def claim_limits():
    return {
      "confirmed_head_or_individual_insect_identity":False,
      "confirmed_seed_feeder_suppression":False,
      "no_herbivore_row_is_not_undamaged_head":True,
      "identifies_physical_spine_barrier":False,
      "identifies_pollination_resource_path":False,
      "source_herbivore_labels_are_verified_seed_consumers":False,
      "meshed_pollinator_exclusion_is_compound_intervention":True,
      "predator_exclusion_is_adhesive_plus_manual_not_true_spine":True,
      "predator_model_is_plastic_spider_plus_exclusion_not_predator_readdition":True,
      "independent_biological_field_samples_added":0,
    }


def linkage_summary(seed_rows, identification_rows, strict=True):
    eligible=[r for r in seed_rows
              if (r["Site_Code"] in SITES and r["Treatment_Code"] in ARMS)]
    heads=defaultdict(list)
    for r in eligible:
        heads[plant_key(r)].append(r)
    insects=defaultdict(list)
    for r in identification_rows:
        z=dict(r)
        z["_number"]=int_value(z.get("Number"),"COUNT")
        if z.get("Feeding_Guild") not in GUILDS:
            raise ValueError("UNEXPECTED_GUILD")
        insects[plant_key(z)].append(z)
    summary=[]
    for k,head_rows in sorted(heads.items()):
        ident=insects.get(k,[])
        seed_n=sum(int(float(z["Insects"])) for z in head_rows)
        labels_n=sum(z["_number"] for z in ident)
        summary.append({
            "key":k,
            "insect_count_on_seed_heads":seed_n,
            "identified_count_total":labels_n,
            "n_identified_source_rows":len(ident),
            "match":seed_n==labels_n
        })
    counts=Counter(k[2] for k in heads)
    byarm={}
    for arm in ARMS:
        group=[x for x in summary if x["key"][2]==arm]
        byarm[arm]={
            "n_proxy_plants":len(group),
            "n_matching_counts":sum(x["match"] for x in group),
            "n_mismatching_counts":sum(not x["match"] for x in group),
            "n_proxy_plants_without_identification_rows":
                sum(x["n_identified_source_rows"]==0 for x in group),
            "head_source_emerged_insects":sum(x["insect_count_on_seed_heads"] for x in group),
            "recorded_identification_individuals":sum(x["identified_count_total"] for x in group),
        }
    joined_keys=set(heads)
    outside=[r for r in identification_rows if plant_key(r) not in joined_keys]
    guild_counts={g:sum(int_value(r.get("Number"),"COUNT")
                        for r in identification_rows if r["Feeding_Guild"]==g)
                  for g in GUILDS}
    result={
      "n_original_seed_heads":len(seed_rows),
      "n_eligible_seed_heads":len(eligible),
      "n_proxy_plants":len(summary),
      "n_identification_rows":len(identification_rows),
      "n_identification_source_groups":len(insects),
      "n_matching_aggregate_counts":sum(x["match"] for x in summary),
      "n_mismatching_aggregate_counts":sum(not x["match"] for x in summary),
      "n_identified_greater_than_source_insects":
          sum(x["identified_count_total"]>x["insect_count_on_seed_heads"] for x in summary),
      "n_identified_lower_than_source_insects":
          sum(x["identified_count_total"]<x["insect_count_on_seed_heads"] for x in summary),
      "n_groups_with_no_identification_rows":
          sum(x["n_identified_source_rows"]==0 for x in summary),
      "n_groups_no_identification_but_positive_seed_count":
          sum(x["n_identified_source_rows"]==0 and x["insect_count_on_seed_heads"]>0
              for x in summary),
      "insect_count_from_all_eligible_seed_heads":
          sum(x["insect_count_on_seed_heads"] for x in summary),
      "identified_guild_count":guild_counts,
      "all_identification_count":sum(guild_counts.values()),
      "n_identification_rows_not_matched_to_seed_proxy":len(outside),
      "identification_count_not_matched_to_seed_proxy":
          sum(int_value(r.get("Number"),"COUNT") for r in outside),
      "by_treatment":byarm,
      "confirmed_head_or_individual_insect_identity":False,
      "confirmed_seed_feeder_suppression":False,
      "claim_limits":claim_limits(),
      "interpretation":"Count reconciliation at plant-proxy aggregation grain; equality does not authorize one-to-one insect-to-head matching or biological guild-specific fecundity mediation.",
    }
    if strict:
        check(len(seed_rows)==615 and len(eligible)==612 and len(summary)==164,
              "ORIGINAL_ELIGIBLE_SOURCE_DRIFT")
        check(len(identification_rows)==181 and len(insects)==93,
              "ORIGINAL_IDENTIFICATION_GROUP_DRIFT")
        check(result["n_matching_aggregate_counts"]==114 and
              result["n_mismatching_aggregate_counts"]==50,
              "ORIGINAL_114_MATCHED_50_MISMATCHED_GROUPS_DRIFT")
        check(result["n_groups_with_no_identification_rows"]==73 and
              result["n_groups_no_identification_but_positive_seed_count"]==4,
              "IDENTIFICATION_ABSENCE_DRIFT")
        check(result["n_identified_greater_than_source_insects"]==25 and
              result["n_identified_lower_than_source_insects"]==25,
              "DIRECTIONAL_MISMATCH_COUNTS_DRIFT")
        check(result["insect_count_from_all_eligible_seed_heads"]==497 and
              result["all_identification_count"]==522 and
              result["identification_count_not_matched_to_seed_proxy"]==2,
              "SAMPLE_UNIVERSE_COUNT_DRIFT")
        check(result["identified_guild_count"]=={"Herbivore":498,"Parasitoid":24},
              "SOURCE_GUILD_COUNT_DRIFT")
        check({a:z["n_mismatching_counts"] for a,z in byarm.items()}==
              {"Control":21,"Pollinator Exclusion":4,
               "Predator Exclusion":14,"Predator Model":11},
              "TREATMENT_SPECIFIC_SOURCE_MISMATCH_DRIFT")
    return result


def _bootstrap_rng(seed):
    state=seed & 0xffffffff
    if not state: raise ValueError("RNG_STATE_MUST_NONZERO")
    while True:
        state ^= (state << 13)&0xffffffff
        state ^= state >> 17
        state ^= (state << 5)&0xffffffff
        state &= 0xffffffff
        yield state


def dual_stage_summary(four_arm_blocks,reps=5000,seed=SEED):
    if reps<100 or reps>30000:raise ValueError("BOOTSTRAP_BUDGET_OUT_OF_RANGE")
    sites={s:sorted([z for z in four_arm_blocks if z["site"]==s],
                    key=lambda v:v["block"]) for s in SITES}
    if any(not group for group in sites.values()):
        raise ValueError("ALL_FOUR_SITES_REQUIRED_FOR_CONDITIONAL_SENSITIVITY")
    ordered=[x for site in SITES for x in sites[site]]
    n=len(ordered)
    expected=Counter(z["site"] for z in four_arm_blocks)
    source={}
    for a,b in PAIRS:
        name=a+" MINUS "+b
        d=[contrast(row,a,b) for row in ordered]
        fields={}
        for m in SEED_METRICS:
            values=[x[m] for x in d]
            fields[m]={
                "mean":mean(values),
                "sign":sign_counts(values),
                "site_means":{s:mean([x[m] for row,x in zip(ordered,d)
                                     if row["site"]==s]) for s in SITES},
                "site_equal_weight_mean":mean([
                    mean([x[m] for row,x in zip(ordered,d) if row["site"]==s])
                    for s in SITES]),
                "leave_one_site_out":{
                    s:mean([x[m] for row,x in zip(ordered,d) if row["site"]!=s])
                    for s in SITES
                },
            }
        source[name]=fields
    # The same bootstrap sampled blocks are used for ALL contrasts and both
    # seed and emergence metrics: independent RCT sampling is NOT fabricated.
    sequences={n:{m:[] for m in SEED_METRICS} for n in source}
    rng=_bootstrap_rng(seed)
    for _ in range(reps):
        selected=[]
        for s in SITES:
            zz=sites[s]
            for _j in range(len(zz)):
                idx=next(rng)%len(zz)
                selected.append(zz[idx])
        for a,b in PAIRS:
            name=a+" MINUS "+b
            for m in SEED_METRICS:
                sequences[name][m].append(mean([
                    row["arms"][a][m]-row["arms"][b][m]
                    for row in selected
                ]))
    for name in source:
        for m in SEED_METRICS:
            v=sequences[name][m]
            source[name][m]["conditional_within_site_block_bootstrap_ci95"]=[
                quantile(v,.025),quantile(v,.975)]
    for m in SEED_METRICS:
        check(math.isclose(source["Predator Model MINUS Control"][m]["mean"],
            source["Predator Exclusion MINUS Control"][m]["mean"]+
            source["Predator Model MINUS Predator Exclusion"][m]["mean"],
            rel_tol=0,abs_tol=1e-11),"NESTED_TREATMENT_DIFFERENCE_IDENTITY_BROKEN")
    return {
        "n_four_arm_blocks":n,
        "n_source_heads_treated_as_independent_replicates":0,
        "site_n":dict(expected),
        "bootstrap_draws_conditional_on_fixed_sites":reps,
        "bootstrap_seed":seed,
        "resampling_scope":"fixed four original 2022 local study sites; sample site x block jointly, retain all four treatment plants within each block; inference cannot generalize to another pollinator regime or prove traits",
        "contrasts":source,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--seeds",type=Path,required=True)
    p.add_argument("--identified-insects",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=5000)
    args=p.parse_args()
    rows,seed_sha=parse_rows(args.seeds.read_bytes(),strict=True)
    idrows,id_sha=parse_identification_rows(args.identified_insects.read_bytes(),strict=True)
    _keep,_proxies,_blocks,four=prep(rows,strict=True)
    result={
       "version":"cirsium_2022_guild_linkage_and_dual_outcome_audit_v2",
       "generated_at_date":"2026-10-20",
       "source":{
           "repo":SOURCE_REPO,"commit":SOURCE_COMMIT,
           "seed_heads_blob":seed_sha,
           "insect_identification_blob":id_sha,
           "original_head_rows":len(rows),
           "identified_insect_rows":len(idrows),
           "preprint_doi":"10.21203/rs.3.rs-6634172/v1",
           "year_of_field_experiment":2022,
       },
       "insect_file_reconciliation":linkage_summary(rows,idrows,strict=True),
       "dual_source_matched_block_outcomes":dual_stage_summary(
           four,reps=args.bootstrap,seed=SEED),
       "claim_ceiling":{
          **claim_limits(),
          "independent_guild_to_seed_fitness_path_identified":False,
          "no_full_factorial_four_treatment_interaction":True,
          "no_morphological_head_armature_intervention":True,
          "site_and_source_conditional_not_genus_general":True,
          "observation_scope":"One four-site experiment and two non-identically enumerated source data files; identification totals cannot substitute for seed-file insect counts.",
       },
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                                   allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
      "status":"SOURCE_PINNED_REAL_STUDY",
      "matched":result["insect_file_reconciliation"]["n_matching_aggregate_counts"],
      "mismatched":result["insect_file_reconciliation"]["n_mismatching_aggregate_counts"],
      "complete_blocks":result["dual_source_matched_block_outcomes"]["n_four_arm_blocks"],
      "contrasts":{k:{m:v[m]["mean"] for m in
          ("viable_fraction","emerged_insects_per_head")} for k,v in
          result["dual_source_matched_block_outcomes"]["contrasts"].items()}
    },indent=2))

if __name__=="__main__":
    main()
