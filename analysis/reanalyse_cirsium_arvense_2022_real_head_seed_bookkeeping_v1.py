#!/usr/bin/env python3
"""Source-verified RAW 2022 Cirsium arvense pollinator-exclusion seed audit.

Karpate et al 2025, Research Square DOI 10.21203/rs.3.rs-6634172/v1.
Open CC BY 4.0 author data/code: DOI 10.5281/zenodo.15364152;
Yogita25/Multitrophic_Interactions Git commit 8a151c62ed069f775a16b9bc5fbdce209cbb8691.

This is a REAL source re-analysis (not synthetic simulation), using 2022
heads nested in inferred SITE×BLOCK×TREATMENT mothers. Contrasts are
DESCRIPTIVE treatment (bag) differences and never direct causal effects of
pollination mediated by seeds, adult seed predation, or head-spine structure.

Critical arithmetic identity confirmed on every raw head:
  Flowers = Viable_Seeds + NViable_Seeds + Damaged_Seeds + Missing_Seeds
  Eaten_Seeds (author R) = Damaged_Seeds + Missing_Seeds
                         = Flowers - Viable_Seeds - NViable_Seeds
Therefore inferred predation (Eaten) is NOT an independent count of damaged
seeds; 'Missing' is a residual that can be negative in the source.
Only transparent endpoints and honest denominator sensitivity reported.

This independently quotes all original row numbers, retains unverifiable
head/vial pairings as unresolved, and computes 20k block-paired stratified
bootstrap percentile intervals for source-observed POLLINATOR EXCLUSION vs
CONTROL, with SITE and original BLOCK intact. No frame replication or
unapproved manipulations.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import json
import math
import random
from pathlib import Path

SOURCE_COMMIT="8a151c62ed069f775a16b9bc5fbdce209cbb8691"
SOURCE_REPO="Yogita25/Multitrophic_Interactions"
DOI_DATA="10.5281/zenodo.15364152"
DOI_PAPER="10.21203/rs.3.rs-6634172/v1"
GIT_BLOBS={
  "seeds.csv":"ca99c762162d15f006f0c73b4603ca71fd320a29",
  "Herbivore.csv":"b0578e8c814b371e910da119c66e45575bc8b675",
  "Pollinators.csv":"1723e96b8fc8ea755832daa786c387c5f3c188a8",
  "Herbivores.R":"6016c74ef3a8783a376f6a6f5e50a4c35f81a75a",
}
SITES={"ditch","construction","campus","rapseed"}
TREATMENTS=("Control","Pollinator Exclusion",
            "Predator Exclusion","Predator Model")
SEED_FIELDS={
  "Site_Code","Block_Code","Treatment_Code","Vial_Code",
  "Viable_Seeds","NViable_Seeds","Damaged_Seeds","Missing_Seeds",
  "Insects","Flowers"
}
MISSING={"","NA","N/A","NaN","null","None"}
BOOT_SEED=20261020
METRICS=(
    "viable_fraction_all",
    "nonviable_fraction_all",
    "directly_recorded_damaged_fraction_all",
    "inferred_eaten_signed_fraction_all",
    "missing_residual_signed_fraction_all",
    "emerged_insects_per_head_all",
    "inferred_eaten_fraction_positive_only",
    "inferred_eaten_fraction_nonnegative_missing_only",
)

def require(ok,message):
    if not ok:raise ValueError(message)

def verify_blob(path,expected):
    data=Path(path).read_bytes()
    git_sha=hashlib.sha1(
        b"blob "+str(len(data)).encode()+b"\0"+data
    ).hexdigest()
    require(git_sha==expected,"PINNED_AUTHOR_GIT_BLOB_CHECKSUM_MISMATCH:"+str(path))
    return data

def source_csv(path,blob):
    raw=verify_blob(path,GIT_BLOBS[blob])
    with Path(path).open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        return list(reader.fieldnames or []),list(reader)

def number(z,key):
    v=str(z.get(key,"")).strip()
    require(v not in MISSING,"MISSING_ACTUAL_NUMERIC_SEED_FIELD_"+key)
    try:f=float(v)
    except ValueError:
        raise ValueError("NONNUMERIC_SEED_FIELD_"+key) from None
    require(math.isfinite(f) and f.is_integer(),"NONINTEGER_OR_NONFINITE_SEED_FIELD_"+key)
    return int(f)

def head_key(r):
    """Vial_Code is NOT a unique head identifier in these raw records."""
    return (r["Site_Code"].strip().lower(),
            r["Block_Code"].strip(),r["Treatment_Code"].strip(),
            r["Vial_Code"].strip())

def parse_head_rows(source):
    validated=[]
    for i,r in enumerate(source,2):
        values={k:number(r,k) for k in
                ("Viable_Seeds","NViable_Seeds","Damaged_Seeds",
                 "Missing_Seeds","Insects","Flowers")}
        require(values["Flowers"]>0,"NONPOSITIVE_FLOWERS_AT_RAW_LINE_"+str(i))
        require(all(values[k]>=0 for k in (
            "Viable_Seeds","NViable_Seeds","Damaged_Seeds","Insects")),
            "OBSERVED_SEED_OR_EMERGENCE_COUNT_NEGATIVE_AT_"+str(i))
        require(sum(values[k] for k in (
            "Viable_Seeds","NViable_Seeds","Damaged_Seeds",
            "Missing_Seeds"))==values["Flowers"],
            "SEED_MASS_BALANCE_SOURCE_IDENTITY_FAILED_AT_"+str(i))
        inferred=values["Damaged_Seeds"]+values["Missing_Seeds"]
        require(inferred==values["Flowers"]-values["Viable_Seeds"]-
                values["NViable_Seeds"],
                "AUTHOR_EATEN_VARIABLE_IS_NOT_COMPONENT_RESIDUAL")
        row={
          "raw_line":i,
          "site":r["Site_Code"].strip().lower(),
          "block":r["Block_Code"].strip(),
          "treatment":r["Treatment_Code"].strip(),
          "vial":r["Vial_Code"].strip(),
          **values,
          "inferred_eaten":inferred,
          "raw_rhead_id":r.get("","").strip(),
        }
        validated.append(row)
    require(len(validated)==615,"PUBLIC_ARCHIVE_615_RAW_HEAD_ROWS_CHANGED")
    return validated

def fraction_sum(rows,numerator):
    if not rows:return None
    denom=sum(r["Flowers"] for r in rows)
    return sum(numerator(r) for r in rows)/denom if denom>0 else None

def plant_stats(rows):
    require(rows,"EMPTY_PLANT_WITH_NO_HEADS")
    selected_positive=[r for r in rows if r["inferred_eaten"]>0]
    selected_clean=[r for r in rows if r["Missing_Seeds"]>=0]
    sf=lambda key,arr=rows:fraction_sum(arr,lambda x:x[key])
    return {
      "n_heads":len(rows),
      "flowers":sum(x["Flowers"] for x in rows),
      "n_heads_with_negative_missing_residual":
          sum(x["Missing_Seeds"]<0 for x in rows),
      "n_heads_with_negative_inferred_eaten":
          sum(x["inferred_eaten"]<0 for x in rows),
      "n_heads_with_positive_inferred_eaten":len(selected_positive),
      "n_heads_with_nonnegative_missing":len(selected_clean),
      "viable_fraction_all":sf("Viable_Seeds"),
      "nonviable_fraction_all":sf("NViable_Seeds"),
      "directly_recorded_damaged_fraction_all":sf("Damaged_Seeds"),
      "inferred_eaten_signed_fraction_all":sf("inferred_eaten"),
      "missing_residual_signed_fraction_all":sf("Missing_Seeds"),
      "emerged_insects_per_head_all":sum(r["Insects"] for r in rows)/len(rows),
      "inferred_eaten_fraction_positive_only":(
          sf("inferred_eaten",selected_positive) if selected_positive else None),
      "inferred_eaten_fraction_nonnegative_missing_only":(
          sf("inferred_eaten",selected_clean) if selected_clean else None),
    }

def percentile(xs,p):
    s=sorted(xs)
    i=(len(s)-1)*p
    lo=int(math.floor(i));hi=int(math.ceil(i))
    return s[lo]*(hi-i)+s[hi]*(i-lo) if hi!=lo else s[lo]

def stratified_pair_bootstrap(paired,metric,n_boot,seed):
    eligible=[
        r for r in paired
        if r["control"][metric] is not None
        and r["pollinator_exclusion"][metric] is not None
    ]
    by_site=collections.defaultdict(list)
    for r in eligible:
        by_site[r["site"]].append(
            r["pollinator_exclusion"][metric]-r["control"][metric])
    n=sum(len(x) for x in by_site.values())
    if not n:return None
    obs=sum(sum(a) for a in by_site.values())/n
    rng=random.Random(seed)
    boots=[]
    sites=sorted(by_site)
    for _ in range(n_boot):
        s=0.
        for site in sites:
            vals=by_site[site]
            for _ in vals:
                s+=vals[rng.randrange(len(vals))]
        boots.append(s/n)
    signs=[v for arr in by_site.values() for v in arr]
    return {
       "n_paired_blocks":n,
       "n_paired_blocks_by_site":{x:len(by_site[x]) for x in sites},
       "mean_paired_plant_fraction_difference":obs,
       "stratified_block_bootstrap_percentile_CI95":[
           percentile(boots,.025),percentile(boots,.975)],
       "bootstrap_replicates":n_boot,
       "bootstrap_seed":seed,
       "n_paired_blocks_positive_difference":sum(v>0 for v in signs),
       "n_paired_blocks_negative_difference":sum(v<0 for v in signs),
       "n_paired_blocks_zero_difference":sum(v==0 for v in signs),
       "not_randomization_exact_or_multi_trial_causal_test":True,
    }

def summarize(tested,herbivore,pollinators,n_boot):
    clean=[
        x for x in tested
        if x["site"] in SITES and x["treatment"] in TREATMENTS
    ]
    excluded=[x for x in tested if x not in clean]
    require(len(clean)==612 and len(excluded)==3,
            "EXPECTED_612_OF_615_DEFINED_EXPERIMENTAL_ROWS_CHANGED")
    groups=collections.defaultdict(list)
    for r in clean:
        k=(r["site"],r["block"],r["treatment"])
        groups[k].append(r)
    require(len(groups)==164,
            "SITE_BLOCK_TREATMENT_INFERRED_PLANT_UNIT_CHANGED")
    blocks=collections.defaultdict(set)
    for site,block,t in groups:
        blocks[(site,block)].add(t)
    require(len(blocks)==45,"REAL_SITE_BLOCK_UNIVERSE_CHANGED")
    plants={
        k:plant_stats(rows) for k,rows in groups.items()
    }
    same_keys=collections.Counter(
        (r["site"],r["block"],r["treatment"],r["vial"]) for r in clean
    )
    duplicate_vial_keys=sum(v>1 for v in same_keys.values())
    inferred_pairs=[]
    for (site,block),ts in sorted(blocks.items()):
        ka=(site,block,"Control")
        kb=(site,block,"Pollinator Exclusion")
        if ka in plants and kb in plants:
            inferred_pairs.append({
                "site":site,"block":block,
                "control":plants[ka],
                "pollinator_exclusion":plants[kb],
            })
    require(len(inferred_pairs)==39,
            "MATCHED_CONTROL_VS_POLLINATOR_EXCLUSION_BLOCK_COUNT_CHANGED")
    comparisons={
        metric:stratified_pair_bootstrap(
            inferred_pairs,metric,n_boot,BOOT_SEED+i*33)
        for i,metric in enumerate(METRICS)
    }
    by_arm={}
    for trt in TREATMENTS:
        rows=[r for r in clean if r["treatment"]==trt]
        positive=[r for r in rows if r["inferred_eaten"]>0]
        nonnegative=[r for r in rows if r["Missing_Seeds"]>=0]
        by_arm[trt]={
          "n_rows":len(rows),
          "n_inferred_mother_plant_IDs":
              len({(r["site"],r["block"]) for r in rows}),
          "n_negative_missing":sum(x["Missing_Seeds"]<0 for x in rows),
          "n_negative_inferred_eaten":sum(x["inferred_eaten"]<0 for x in rows),
          "n_zero_inferred_eaten":sum(x["inferred_eaten"]==0 for x in rows),
          "n_inferred_eaten_strictly_positive":len(positive),
          "n_nonnegative_missing":len(nonnegative),
          "all_heads_achene_weighted_viable":fraction_sum(
              rows,lambda x:x["Viable_Seeds"]),
          "all_heads_achene_weighted_visible_damaged":fraction_sum(
              rows,lambda x:x["Damaged_Seeds"]),
          "all_heads_achene_weighted_signed_inferred_eaten":fraction_sum(
              rows,lambda x:x["inferred_eaten"]),
          "positive_only_achene_weighted_inferred_eaten":fraction_sum(
              positive,lambda x:x["inferred_eaten"]),
          "nonnegative_missing_achene_weighted_inferred_eaten":fraction_sum(
              nonnegative,lambda x:x["inferred_eaten"]),
          "n_total_emerged_insect_count_as_scored_in_seeds_file":
              sum(x["Insects"] for x in rows),
        }
    # Treat external insect emergence table as potentially ambiguous
    # because multiple original rows share one vial identifier.
    key_to_heads=collections.Counter(
      (r["site"],r["block"],r["treatment"],r["vial"]) for r in clean
    )
    link_status=collections.Counter()
    for e in herbivore:
        key=(e["Site_Code"].strip().lower(),
             e["Block_Code"].strip(),e["Treatment_Code"].strip(),
             e["Vial_Code"].strip())
        n=key_to_heads[key]
        link_status["SOURCE_NO_HEAD_KEY_MATCH" if n==0 else
                    "MULTIPLE_HEAD_ROWS_SHARE_VIAL_KEY" if n>1
                    else "ONE_SOURCE_HEAD_KEY_MATCH"]+=1
    return {
      "status":"REAL_PUBLIC_2022_CIRSIUM_ARVENSE_MULTITROPHIC_DATA_REANALYZED_SOURCE_BOUNDED",
      "original_repository":SOURCE_REPO,
      "original_git_commit":SOURCE_COMMIT,
      "original_zenodo_doi":DOI_DATA,
      "study_preprint_doi":DOI_PAPER,
      "archive_license":"CC BY 4.0",
      "audit":{
        "original_seed_rows":len(tested),
        "usable_site_treatment_seed_rows":len(clean),
        "rows_with_missing_site_or_treatment_excluded":len(excluded),
        "paper_reported_mature_head_count":609,
        "head_rows_minus_manuscript_stated_count":len(tested)-609,
        "original_head_Vial_Code_not_unique_within_mother_plant":True,
        "duplicated_site_block_treatment_vial_keys":duplicate_vial_keys,
        "inferred_site_block_treatment_mother_plant_units":len(groups),
        "distinct_site_block_sets":len(blocks),
        "complete_four_arm_site_blocks":sum(len(z)==4 for z in blocks.values()),
        "paired_control_vs_pollinator_exclusion_blocks":len(inferred_pairs),
        "all_seed_rows_have_component_mass_balance_identity":True,
        "all_seed_rows_Eaten_equals_Flowers_minus_viable_minus_nonviable":True,
        "all_raw_seed_rows_with_negative_Missing":sum(
            r["Missing_Seeds"]<0 for r in tested),
        "all_raw_seed_rows_with_negative_Eaten":sum(
            r["inferred_eaten"]<0 for r in tested),
        "all_raw_seed_rows_with_Eaten_equal_zero":sum(
            r["inferred_eaten"]==0 for r in tested),
        "all_raw_seed_rows_with_Eaten_positive":sum(
            r["inferred_eaten"]>0 for r in tested),
        "original_author_Herbivores_R_filters_Eaten_Seed_gt_zero":True,
        "original_author_Herbivore_csv_rows":len(herbivore),
        "original_author_Pollinators_csv_observation_windows":len(pollinators),
        "herbivore_table_individual_head_join_quality":dict(link_status),
      },
      "treatment_arm_head_reductions":by_arm,
      "matched_site_block_plant_unit_contrasts":{
          "description":"Pollinator exclusion organza BAGS minus unbagged control, plant aggregated head seed ratio, 39 originally paired site×block sets. 20k stratified block bootstrap percentile CIs. Not direct causal pollen mediation or anatomical armature effect.",
          "data":comparisons,
      },
      "hard_biological_limits":[
        "Missing_Seeds is the accounting residual of Flowers minus three counted categories, not independently scored physical predation.",
        "Eaten_Seeds=Damaged_Seeds+Missing_Seeds algebraically cancels the directly recorded damaged seed count.",
        "Negative Missing and Eaten are impossible biological seed counts: preserve as data quality flags; no silent zeroing.",
        "Confirmed visible Damaged_Seeds, inferred missing seed loss and reproductive seed failure are distinct.",
        "The author R predation models select Eaten_Seeds>0; this excludes unassessed or true-zero damage heads and can induce selection.",
        "Plant IDs are inferred from site × block × treatment following original analysis; original head vial ID is duplicated, so insect species to head identity sometimes ambiguous.",
        "Organza pollinator-exclusion bag could also restrict seed-feeder adult/ovipositor access and alter microclimate; no independent pollination-only intervention.",
        "Seed categories did not assay head spine length, orientation to gravity, phyllary gap or physical contact; no Q2 evolutionary trait payoff estimated.",
        "A head fraction uses Flowers as denominator, not number of actual fertile achenes or independently germination-tested viable seeds.",
        "Source article is a Research Square preprint; result replication in new analysis is descriptive and distinct from publication claims.",
        "Site-block paired bootstrap ignores uncertainty in inferred plant identification and is not a pre-registered causal effect estimator.",
        "Only one published Cirsium species and observational contrast of one experimental intervention, not a cross-lineage gain estimate or Q3 re-expression."
      ],
    }

def tests():
    fixture=[
      {"Site_Code":"ditch","Block_Code":"B1",
       "Treatment_Code":"Control","Vial_Code":"1",
       "Flowers":"10","Viable_Seeds":"5","NViable_Seeds":"3",
       "Damaged_Seeds":"1","Missing_Seeds":"1","Insects":"2","":"1"},
      {"Site_Code":"ditch","Block_Code":"B1",
       "Treatment_Code":"Pollinator Exclusion","Vial_Code":"1",
       "Flowers":"10","Viable_Seeds":"2","NViable_Seeds":"7",
       "Damaged_Seeds":"0","Missing_Seeds":"1","Insects":"0","":"2"},
      {"Site_Code":"ditch","Block_Code":"B1",
       "Treatment_Code":"Control","Vial_Code":"2",
       "Flowers":"10","Viable_Seeds":"8","NViable_Seeds":"4",
       "Damaged_Seeds":"0","Missing_Seeds":"-2","Insects":"0","":"3"},
    ]
    x=fixture[:]
    # The complete real-data contract requires 615 rows; for fixtures,
    # just test the underlying count identity as a local source invariant.
    p=[]
    for row in x:
        v={k:number(row,k) for k in ("Flowers","Viable_Seeds",
           "NViable_Seeds","Damaged_Seeds","Missing_Seeds","Insects")}
        assert v["Flowers"]==sum(v[k] for k in
               ("Viable_Seeds","NViable_Seeds","Damaged_Seeds","Missing_Seeds"))
        v["inferred_eaten"]=v["Damaged_Seeds"]+v["Missing_Seeds"]
        assert v["inferred_eaten"]==v["Flowers"]-v["Viable_Seeds"]-v["NViable_Seeds"]
        p.append(v)
    assert p[0]["inferred_eaten"]==2
    assert p[2]["inferred_eaten"]==-2
    assert sum(z["inferred_eaten"]>0 for z in p)==2
    c=plant_stats(p)
    assert c["n_heads_with_negative_missing_residual"]==1
    assert c["n_heads_with_positive_inferred_eaten"]==2
    assert abs(c["inferred_eaten_signed_fraction_all"]-1/30)<1e-10
    return {
      "status":"PASS_SYNTHETIC_SEED_BOOKKEEPING_AND_NEGATIVE_RESIDUAL_GATES",
      "n_synthetic_head_rows":len(fixture),
      "n_real_plant_analyses_simulated":0,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--seeds",type=Path,required=True)
    p.add_argument("--herbivore",type=Path,required=True)
    p.add_argument("--pollinators",type=Path,required=True)
    p.add_argument("--original-r-code",type=Path,required=True)
    p.add_argument("--boot",type=int,default=20000)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    require(100<=a.boot<=100000,"BOOTSTRAP_REPLICATE_LIMIT")
    s_head,seeds=source_csv(a.seeds,"seeds.csv")
    h_head,herb=source_csv(a.herbivore,"Herbivore.csv")
    p_head,pol=source_csv(a.pollinators,"Pollinators.csv")
    require(SEED_FIELDS.issubset(s_head),"ORIGINAL_SEED_DATA_HEADER_CHANGED")
    require({"Site_Code","Block_Code","Treatment_Code","Vial_Code",
             "Feeding_Guild","Number"}.issubset(h_head),
            "ORIGINAL_INSECT_HEADER_CHANGED")
    require({"Site_Code","Block_Code","Treatment_Code",
             "Number_openFlowerhead"}.issubset(p_head),
            "ORIGINAL_POLLINATOR_HEADER_CHANGED")
    old_R=verify_blob(a.original_r_code,GIT_BLOBS["Herbivores.R"]).decode("utf-8")
    require("seeds.eaten <- seeds[seeds$Eaten_Seeds>0,]" in old_R,
            "AUTHOR_PUBLISHED_POSITIVE_EATEN_SELECTION_CODE_NOT_CONFIRMED")
    require("seeds$Eaten_Seeds <- seeds$Damaged_Seeds + seeds$Missing_Seeds" in old_R,
            "AUTHOR_EATEN_ALGEBRA_SOURCE_HAS_CHANGED")
    heads=parse_head_rows(seeds)
    analysis=summarize(heads,herb,pol,a.boot)
    payload={
        "version":"cirsium_arvense_2022_original_head_seed_accounting_reanalysis_v1",
        "date":"2026-10-20",
        "source":{"repository":SOURCE_REPO,
                  "original_commit":SOURCE_COMMIT,
                  "Zenodo_DOI":DOI_DATA,
                  "preprint_DOI":DOI_PAPER,
                  "git_blobs":GIT_BLOBS,
                  "source_type":"actual published author field CSV and author R code, not synthetic experiment",
                  "license":"CC BY 4.0"},
        "synthetic_negative_controls":tests(),
        "real_Cirsium_source_reanalysis":analysis,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",
                     encoding="utf-8")
    q=analysis["matched_site_block_plant_unit_contrasts"]["data"]
    print(json.dumps({
        "status":analysis["status"],
        "original_rows":analysis["audit"]["original_seed_rows"],
        "valid_heads":analysis["audit"]["usable_site_treatment_seed_rows"],
        "negative_residual_heads":analysis["audit"]["all_raw_seed_rows_with_negative_Missing"],
        "matched_blocks":analysis["audit"]["paired_control_vs_pollinator_exclusion_blocks"],
        "paired_mean_differences":{
            k:{
              "diff":v["mean_paired_plant_fraction_difference"],
              "CI95":v["stratified_block_bootstrap_percentile_CI95"]
            } if v else None for k,v in q.items()
        },
        "data_quality":"INFERRED_EATEN_IS_COMPOSITIONAL_RESIDUAL_NOT_AN_INDEPENDENT_PHYSICAL_DAMAGE_COUNT",
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
