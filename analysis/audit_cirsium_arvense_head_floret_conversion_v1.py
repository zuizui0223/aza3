#!/usr/bin/env python3
"""Real Cirsium arvense head fecundity: number of florets vs good-seed conversion.

From ORIGINAL Karpate et al (2022 fieldwork; 2025 preprint) original seeds.csv
on author fixed Git commit. No field data, anatomy, genetics or new plants.

Within each inferred site×block×arm proxy:
  good seeds / head = (source Flowers / sampled heads) * (Viable_Seeds / Flowers).
For an arm A versus arm B, the exact symmetric two-factor decomposition is:
  difference in good seeds/head =
  (florets_A-florets_B)*(p_A+p_B)/2 +
  (p_A-p_B)*(florets_A+florets_B)/2.
This is a Shapley/Kitagawa ACCOUNTING IDENTITY, not a causal mediation model.
Source Flowers (rather than pre-treatment phyllary anatomy or floret buds)
is a possibly post-treatment denominator. Original signed Missing_Seeds is an
accounting residual, not observed feeding. Four-arm site-block 36 only;
treating 615 head records as 615 randomized replicates is invalid.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path

from audit_cirsium_2022_four_arm_head_outcomes_v1 import (
    parse_rows, prep, SITES, ARMS, SOURCE_GIT_BLOB, SOURCE_COMMIT,
    SOURCE_REPO, ARCHIVE_DOI, PAPER_DOI, mean, quantile, check,
)

METRICS=(
 "florets_per_head_delta", "good_fraction_delta", "good_seeds_per_head_delta",
 "floret_quantity_component", "good_fraction_component",
 "nonviable_fraction_delta", "direct_damage_fraction_delta",
 "signed_missing_fraction_delta",
)
COMPARISONS=(
 ("bag_minus_control","Pollinator Exclusion","Control"),
 ("predator_exclusion_minus_control","Predator Exclusion","Control"),
 ("predator_model_minus_control","Predator Model","Control"),
 ("predator_model_minus_predator_exclusion","Predator Model","Predator Exclusion"),
 ("bag_minus_predator_exclusion","Pollinator Exclusion","Predator Exclusion"),
)
SEED=20261023
N_BOOT=10000
MASK=0xffffffff

def rng32(seed=SEED):
    check(isinstance(seed,int) and 0<seed<2**32,"INVALID_RNG_SEED")
    s=seed
    while True:
        s^=(s<<13)&MASK
        s^=s>>17
        s^=(s<<5)&MASK
        s&=MASK
        yield s

def signs(v):
    return {"positive":sum(x>1e-10 for x in v),
            "negative":sum(x< -1e-10 for x in v),
            "zero":sum(abs(x)<=1e-10 for x in v)}

def plant_head_counts(heads):
    check(len(heads)>0,"UNOBSERVED_HEAD_GROUP")
    def summed(key):
        return sum(int(float(h[key])) for h in heads)
    f=summed("Flowers")
    v=summed("Viable_Seeds")
    n=summed("NViable_Seeds")
    d=summed("Damaged_Seeds")
    m=summed("Missing_Seeds")
    check(f>0 and f==v+n+d+m,"SOURCE_FLORET_CATEGORY_IDENTITY_BROKEN")
    total=len(heads)
    out={
      "n_heads":total,"reported_florets":f,
      "good_seeds":v,"nonviable_seeds":n,
      "damaged_seeds":d,"signed_missing":m,
      "florets_per_head":f/total,
      "good_seeds_per_head":v/total,
      "good_fraction":v/f,
      "nonviable_fraction":n/f,
      "direct_damage_fraction":d/f,
      "signed_missing_fraction":m/f,
    }
    check(math.isclose(out["florets_per_head"]*out["good_fraction"],
                       out["good_seeds_per_head"],rel_tol=0,abs_tol=1e-10),
          "INVALID_HEAD_RATE_PRODUCT")
    return out

def compare_one(a,b,site,block):
    x=a["florets_per_head"]; y=b["florets_per_head"]
    p=a["good_fraction"]; q=b["good_fraction"]
    dx=x-y
    dp=p-q
    dv=a["good_seeds_per_head"]-b["good_seeds_per_head"]
    qterm=dx*(p+q)/2
    rterm=dp*(x+y)/2
    dn=a["nonviable_fraction"]-b["nonviable_fraction"]
    dd=a["direct_damage_fraction"]-b["direct_damage_fraction"]
    dm=a["signed_missing_fraction"]-b["signed_missing_fraction"]
    check(math.isclose(qterm+rterm,dv,rel_tol=0,abs_tol=1e-9),
          "SHAPLEY_PRODUCT_IDENTITY_BROKEN")
    check(math.isclose(dp+dn+dd+dm,0,rel_tol=0,abs_tol=1e-9),
          "SIGNED_REMAINDER_SEED_CATEGORY_IDENTITY_BROKEN")
    return {
      "site":site,"block":block,
      "florets_per_head_delta":dx,"good_fraction_delta":dp,
      "good_seeds_per_head_delta":dv,
      "floret_quantity_component":qterm,"good_fraction_component":rterm,
      "nonviable_fraction_delta":dn,"direct_damage_fraction_delta":dd,
      "signed_missing_fraction_delta":dm,
    }

def summarize_blocks(blocks):
    check(len(blocks)>=2,"REQUIRE_AT_LEAST_TWO_ORIGINAL_BLOCKS")
    by_site={s:[v for v in blocks if v["site"]==s] for s in SITES}
    check(all(by_site.values()),"FOUR_ORIGINAL_SITES_REQUIRED")
    result={}
    for slug,arm_a,arm_b in COMPARISONS:
        paired=[]
        for z in blocks:
            paired.append(compare_one(z["arms"][arm_a],z["arms"][arm_b],
                                      z["site"],z["block"]))
        m={}
        for field in METRICS:
            data=[x[field] for x in paired]
            m[field]={
              "mean":mean(data),
              "sign":signs(data),
              "site_means":{
                s:mean([v[field] for v in paired if v["site"]==s])
                for s in SITES},
              "leave_one_site_out":{
                s:mean([v[field] for v in paired if v["site"]!=s])
                for s in SITES},
            }
        result[slug]={
          "arm_a":arm_a,"arm_b":arm_b,
          "n_blocks":len(paired),
          "measurements":m,"paired_block_rows":paired,
          "total_identity_mean_residual":m["good_seeds_per_head_delta"]["mean"]-(
            m["floret_quantity_component"]["mean"]+
            m["good_fraction_component"]["mean"]),
          "good_fraction_categories_identity_mean_residual":
            sum(m[k]["mean"] for k in (
              "good_fraction_delta","nonviable_fraction_delta",
              "direct_damage_fraction_delta","signed_missing_fraction_delta")),
        }
    return result

def original_blocks(raw):
    rows,sourcehash=parse_rows(raw,strict=True)
    eligible,proxy,blocks,complete=prep(rows,True)
    # Sorting by actual source site & block ID, not assumed temporal order.
    ordered=sorted(complete,key=lambda b:(b["site"],b["block"]))
    result=[]
    for z in ordered:
        s=z["site"];b=z["block"]
        result.append({"site":s,"block":b,
                       "arms":{a:plant_head_counts(proxy[(s,b,a)]) for a in ARMS}})
    counts={
      "raw_heads":len(rows),"eligible_heads":len(eligible),
      "plant_proxies":len(proxy),"site_blocks":len(blocks),
      "complete_blocks":len(result),
      "complete_sites":{s:sum(z["site"]==s for z in result) for s in SITES},
    }
    check(counts=={
      "raw_heads":615,"eligible_heads":612,"plant_proxies":164,
      "site_blocks":45,"complete_blocks":36,
      "complete_sites":{"campus":12,"construction":9,"ditch":8,"rapseed":7}},
      "AUTHENTIC_36_FOUR_ARM_SOURCE_DRIFT")
    return result,counts,sourcehash

def bootstrap(comp,blocks,reps=N_BOOT):
    check(isinstance(reps,int) and 100<=reps<=20000,
          "INVALID_BOOTSTRAP_DRAWS")
    pool={s:[i for i,b in enumerate(blocks) if b["site"]==s] for s in SITES}
    check(all(pool.values()),"CANNOT_DROP_ORIGINAL_STUDY_SITE")
    prng=rng32(SEED)
    draws={slug:{name:[] for name in METRICS} for slug,_,_ in COMPARISONS}
    for _ in range(reps):
        indices=[]
        for s in SITES:
            a=pool[s]
            for _j in range(len(a)):
                indices.append(a[next(prng)%len(a)])
        for slug,_,_ in COMPARISONS:
            rows=comp[slug]["paired_block_rows"]
            for metric in METRICS:
                draws[slug][metric].append(
                    mean([rows[i][metric] for i in indices]))
    for slug,_,_ in COMPARISONS:
        for metric in METRICS:
            vals=draws[slug][metric]
            comp[slug]["measurements"][metric]["bootstrap_ci95"]=[
                quantile(vals,.025),quantile(vals,.975)]
    return comp

def simulate_control():
    a={"florets_per_head":100,"good_fraction":.7,
       "good_seeds_per_head":70,"nonviable_fraction":.2,
       "direct_damage_fraction":.1,"signed_missing_fraction":0}
    b={"florets_per_head":80,"good_fraction":.5,
       "good_seeds_per_head":40,"nonviable_fraction":.3,
       "direct_damage_fraction":.2,"signed_missing_fraction":0}
    y=compare_one(a,b,"campus","synthetic")
    check(y["good_seeds_per_head_delta"]==30 and
          math.isclose(y["floret_quantity_component"],12,abs_tol=1e-12) and
          math.isclose(y["good_fraction_component"],18,abs_tol=1e-12),
          "SYNTHETIC_DECOMPOSITION_NOT_30_12_18")
    return "PASS_SOURCE_FLORET_PRODUCT_AND_TWO_CATEGORY_IDENTITIES"

def analyze(raw,reps=N_BOOT):
    blocks,counts,h=original_blocks(raw)
    comparisons=summarize_blocks(blocks)
    bootstrap(comparisons,blocks,reps)
    main=comparisons["bag_minus_control"]["measurements"]
    check(main["good_fraction_delta"]["sign"]==
          {"positive":0,"negative":36,"zero":0},
          "UNEXPECTED_ORIGINAL_BAGGED_GOOD_RATE_SIGN")
    check(main["good_fraction_component"]["mean"]< -40 and
          abs(main["floret_quantity_component"]["mean"])<1,
          "ORIGINAL_STAGE_COMPONENT_DRIFT")
    return {
      "version":"aza3_arvense_head_floret_conversion_original_Python_v1",
      "as_of":"2026-10-23",
      "source":{
         "original_author":"Yogita Karpate et al",
         "study_year":2022,"preprint_doi":PAPER_DOI,
         "data_doi":ARCHIVE_DOI,
         "repo":SOURCE_REPO,"commit":SOURCE_COMMIT,
         "seeds_csv_git_blob":h,
      },
      "source_counts":counts,
      "design":{
        "unit":"site-block x treatment plant proxy (not known biological plant identifier)",
        "method":"symmetric Shapley/Kitagawa two-factor accounting identity applied to all 36 original four-arm complete blocks",
        "product":"viable_seeds_per_head = (source Flowers/number sampled heads) * (source Viable_Seeds/source Flowers)",
        "split":"quantity = (Florets_A-Florets_B)*(p_A+p_B)/2; rate = (p_A-p_B)*(Florets_A+Florets_B)/2",
        "fraction_identity":"p_good + p_nonviable + p_direct_damaged + p_signed_missing == 1",
        "source_signed_missing_is_residual_not_confirmed_eaten":True,
      },
      "bootstrap":{
        "seed":SEED,"draws":reps,
        "unit":"complete original site block",
        "sites_fixed":True,"resample_within_each_original_site":True,
        "site_order":list(SITES),
        "sort_order":"site then original block lexicographic",
        "algorithm":"xorshift32",
        "conditional_descriptive_not_causal":True,
      },
      "comparisons":comparisons,
      "bag_fraction_component_share_of_mean_good_seed_difference":
        main["good_fraction_component"]["mean"]/
        main["good_seeds_per_head_delta"]["mean"],
      "limits":{
        "not_true_pre_treatment_floret_count":True,
        "flower_count_is_author_Flowers_not_calibrated_head_geometry":True,
        "ratio_of_components_not_separate_causal_path_effects":True,
        "bag_blocks_both_pollinators_and_some_antagonists":True,
        "signed_missing_cannot_be_called_direct_consumption":True,
        "zero_new_fieldwork_or_anatomical_measurements":0,
        "no_trait_partner_genetic_LEGO_causal_evidence":True,
      },
      "synthetic_control":simulate_control(),
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--seeds",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--bootstrap",type=int,default=N_BOOT)
    q=p.parse_args()
    z=analyze(q.seeds.read_bytes(),q.bootstrap)
    q.out.parent.mkdir(parents=True,exist_ok=True)
    q.out.write_text(json.dumps(z,indent=2,ensure_ascii=False,allow_nan=False)
                     +"\n",encoding="utf-8")
    key=z["comparisons"]["bag_minus_control"]["measurements"]
    print(json.dumps({"source_counts":z["source_counts"],
                      "bag_mean_components":{k:key[k] for k in (
                        "good_seeds_per_head_delta","floret_quantity_component",
                        "good_fraction_component","florets_per_head_delta",
                        "good_fraction_delta")},
                      "conversion_share":z[
                        "bag_fraction_component_share_of_mean_good_seed_difference"]},
                     ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
