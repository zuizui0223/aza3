#!/usr/bin/env python3
"""Original USDA Cirsium undulatum lateral-head insect damage and seed outcomes.

CC0 West & Louda (2021), DOI 10.1002/ajb2.1744, USDA DOI
10.15482/USDA.ADC/1522649; original 106 mother / 323 flower-head rows
already source-locked in aza3 PR #73.

Primary source grain is individual heads with an actual mother/year ID.
Treat two 2007 source-contradictory maternal treatment labels as ineligible
for the MAIN head-linked comparison (but keep them and the 14 orphan heads
in explicit QA). 297 heads remain, 87 mothers have >=1 lateral head,
82 mothers have paired scored apex and lateral injury. Original 0.88 is
the HIGHEST CODED head damage score, NOT independently measured 88% head loss.

Mother-cluster analysis selects plants with one or more recorded lateral heads;
it is DESCRIPTIVE conditional on having such heads and is not the
intention-to-treat effect among the originally allocated 106 mothers.
Consequently exact mother-level "any lateral good seed" tests below use ALL
original assigned mother rows, never observed-lateral-only rows.

Head damage is itself observed after treatment; adjusting for it can select
on the treatment outcome. Good seeds are author-classified intact-looking
seeds, not experimental germination. No direct pollinator service, head
spine manipulation, heritable LEGO architecture or trait re-expression.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

from audit_cirsium_undulatum_maternal_hurdle_size_v1 import (
    BLOBS, YEARS, ARMS, read_locked, source_groups, value,
    fisher_two_sided, bh_qvalues,
)

SOURCE_DOI="10.15482/USDA.ADC/1522649"
PAPER_DOI="10.1002/ajb2.1744"
METRICS=("mean_damage","fraction_max","fraction_positive","mean_goodseed")
BOOT_SEED=20261022
BOOT_REPS=5000
MASK=0xffffffff
SOURCE_CONFLICTS=("2007|323","2007|375")


def require(ok,code):
    if not ok:raise ValueError(code)


def identity(row):
    return row["year"]+"|"+row["plt"]


def mean(numbers):
    require(bool(numbers),"EMPTY_MOTHER_HEAD_GROUP")
    return sum(numbers)/len(numbers)


def source_audit(mothers,heads):
    ix={identity(m):m for m in mothers}
    require(len(ix)==106 and len(heads)==323,
            "ORIGINAL_MOTHER_HEAD_SAMPLE_COUNTS_CHANGED")
    matched=[h for h in heads if identity(h) in ix]
    orphan=[h for h in heads if identity(h) not in ix]
    conflicts=[h for h in matched if (
        h["dam"]!=ix[identity(h)]["dam"] or
        h["ins"]!=ix[identity(h)]["ins"])]
    clean=[h for h in matched if identity(h) not in {
        identity(x) for x in conflicts}]
    matched_ids=sorted(set(identity(h) for h in matched))
    unique_heads=[identity(h)+"|"+h["hd"] for h in heads]
    require(len(unique_heads)==len(set(unique_heads)),
            "DUPLICATE_ORIGINAL_HEAD_ID")
    require(len(matched)==309 and len(orphan)==14 and
            len(conflicts)==12 and len(clean)==297 and
            sorted({identity(z) for z in conflicts})==list(SOURCE_CONFLICTS),
            "ORIGINAL_MOTHER_HEAD_SOURCE_JOIN_CHANGED")
    require(all(h["dam"]==ix[identity(h)]["dam"] and
                h["ins"]==ix[identity(h)]["ins"] for h in clean),
            "SILENT_HEAD_TREATMENT_REASSIGNMENT")
    return ix,clean,{
        "n_source_mothers":len(mothers),
        "n_source_head_rows":len(heads),
        "n_head_with_mother":len(matched),
        "n_head_without_mother":len(orphan),
        "n_head_treatment_label_mismatch":len(conflicts),
        "n_head_source_aligned":len(clean),
        "n_mothers_with_scored_nonapical_heads":None,
        "n_matched_heads_with_nonmissing_damage_score":
            sum(value(h,"av.dam") is not None for h in clean),
        "n_apical_lateral_damage_pairs":None,
        "head_conflicting_source_mother_ids":
            sorted(set(identity(h) for h in conflicts)),
    }


def real_mother_nonapical(mothers,head_rows):
    heads=defaultdict(list)
    for h in head_rows:
        heads[identity(h)].append(h)
    results=[]
    for mother in mothers:
        x=heads.get(identity(mother),[])
        laterals=[h for h in x if h["hd"]!="100"]
        if not laterals:
            continue
        scored=[h for h in laterals if value(h,"av.dam") is not None]
        require(scored,"NO_SCORED_LATERAL_HEAD_ON_MOTHER_"+identity(mother))
        apical=next((h for h in x if h["hd"]=="100"),None)
        result={
            "id":identity(mother),
            "year":mother["year"],
            "arm":mother["trt"],
            "n_lateral":len(laterals),
            "n_lateral_scored":len(scored),
            "mean_damage":mean([value(h,"av.dam") for h in scored]),
            "fraction_max":mean([
                 float(value(h,"av.dam")==0.88) for h in scored]),
            "fraction_positive":mean([
                 float(value(h,"hd.gsd")>0) for h in laterals]),
            "mean_goodseed":mean([
                 value(h,"hd.gsd") for h in laterals]),
            "apical_damage":value(apical,"av.dam") if apical else None,
            "apical_seed":value(apical,"hd.gsd") if apical else None,
        }
        require(all(math.isfinite(result[k]) for k in METRICS),
                "NONFINITE_HEAD_QUALITY")
        results.append(result)
    return sorted(results,key=lambda z:z["id"])


def group_mother_stats(records):
    return {
        "n_mothers":len(records),
        "n_head_rows":sum(r["n_lateral"] for r in records),
        "mean_scored_damage_per_mother":mean([r["mean_damage"] for r in records]),
        "mean_fraction_max_damage_per_mother":mean([r["fraction_max"] for r in records]),
        "mean_fraction_seed_positive_heads_per_mother":
            mean([r["fraction_positive"] for r in records]),
        "mean_seed_count_per_lateral_head_with_mother_equal_weight":
            mean([r["mean_goodseed"] for r in records]),
        "n_lateral_heads_per_mother":mean([r["n_lateral"] for r in records]),
    }


def year_arm_records(records):
    out={
        y:{a:sorted([r for r in records if r["year"]==y and r["arm"]==a],
                    key=lambda z:z["id"]) for a in ARMS}
        for y in YEARS
    }
    require({y:{a:len(out[y][a]) for a in ARMS} for y in YEARS}=={
        "2007":{"ndni":7,"dni":8,"ndi":7,"di":10},
        "2008":{"ndni":11,"dni":12,"ndi":17,"di":15},
    },"LATERAL_HEAD_MOTHER_COHORT_CHANGED")
    return out


def xorshift32(seed=BOOT_SEED):
    require(isinstance(seed,int) and 0<seed<2**32,"INVALID_BOOTSTRAP_SEED")
    state=seed
    while True:
        state^=(state<<13)&MASK
        state^=state>>17
        state^=(state<<5)&MASK
        state&=MASK
        yield state


def quantile(data,p):
    x=sorted(data)
    require(bool(x) and 0<=p<=1,"INVALID_BOOTSTRAP_QUANTILE")
    t=(len(x)-1)*p
    lo=math.floor(t)
    frac=t-lo
    return x[lo]*(1-frac)+x[min(lo+1,len(x)-1)]*frac


def cluster_bootstrap(groups,reps=BOOT_REPS,seed=BOOT_SEED):
    require(100<=reps<=20000,"BOOTSTRAP_REPETITIONS_INVALID")
    labels=(("intact","ndi","ndni"),("damaged","di","dni"))
    samples={y:{name:{metric:[] for metric in METRICS}
                for name,_,_ in labels} for y in YEARS}
    rng=xorshift32(seed)
    for _ in range(reps):
        for y in YEARS:
            group={}
            for a in ARMS:
                original=groups[y][a]
                vals={k:0. for k in METRICS}
                for _j in range(len(original)):
                    r=original[next(rng)%len(original)]
                    for k in METRICS:
                        vals[k]+=r[k]
                group[a]={k:v/len(original) for k,v in vals.items()}
            for name,reduced,ambient in labels:
                for k in METRICS:
                    samples[y][name][k].append(
                        group[reduced][k]-group[ambient][k])
    return {
        y:{
          name:{
            metric:{
              "estimated_reduced_minus_ambient":
                 mean([r[metric] for r in groups[y][reduced]])-
                 mean([r[metric] for r in groups[y][ambient]]),
              "maternal_percentile_CI95":[
                  quantile(samples[y][name][metric],.025),
                  quantile(samples[y][name][metric],.975)],
            }
            for metric in METRICS}
          for name,reduced,ambient in labels}
        for y in YEARS}


def head_arm_summary(original_heads,mothers):
    ix={identity(m):m for m in mothers}
    out={y:{} for y in YEARS}
    for y in YEARS:
        for arm in ARMS:
            hs=[h for h in original_heads if h["year"]==y
                and h["hd"]!="100" and ix[identity(h)]["trt"]==arm]
            scored=[h for h in hs if value(h,"av.dam") is not None]
            require(scored,"EMPTY_HEAD_DAMAGE_GROUP")
            out[y][arm]={
               "n_lateral_head_records":len(hs),
               "n_lateral_damage_scored":len(scored),
               "n_highest_recorded_damage_score":
                   sum(value(h,"av.dam")==.88 for h in scored),
               "n_any_good_seed_heads":
                   sum(value(h,"hd.gsd")>0 for h in hs),
               "mean_lateral_damage_score_scored":
                   mean([value(h,"av.dam") for h in scored]),
               "high_score_is_max_observed_0p88_not_exact_fraction_destroyed":True,
            }
    require(sum(z["n_lateral_head_records"] for y in YEARS
                for z in out[y].values())==193,
            "SOURCE_CLEAN_LATERAL_HEAD_TOTAL_CHANGED")
    return out


def mother_any_lateral(plant_groups):
    output={y:{} for y in YEARS}
    for y in YEARS:
        for arm in ARMS:
            mothers=plant_groups[y][arm]
            output[y][arm]={
                "n_mothers":len(mothers),
                "n_positive_nonapical_seed_mothers":
                    sum(value(m,"o.hds")>0 for m in mothers),
            }
    return output


def exact_any_lateral(output):
    results={y:{} for y in YEARS}
    ps={}
    for y in YEARS:
        for name,reduced,ambient in (
            ("intact","ndi","ndni"),("damaged","di","dni")
        ):
            a=output[y][reduced]
            b=output[y][ambient]
            a1=a["n_positive_nonapical_seed_mothers"]
            b1=b["n_positive_nonapical_seed_mothers"]
            p=fisher_two_sided(a1,a["n_mothers"]-a1,b1,b["n_mothers"]-b1)
            ps[y+"_"+name]=p
            results[y][name]={
                "arms_reduced_ambient":[reduced,ambient],
                "n_mothers_reduced_ambient":[a["n_mothers"],b["n_mothers"]],
                "n_any_lateral_seeds_reduced_ambient":[a1,b1],
                "p_any_lateral_seed_difference":
                  a1/a["n_mothers"]-b1/b["n_mothers"],
                "fisher_p":p,
            }
    qs=bh_qvalues(ps)
    for y in YEARS:
        for name in ("intact","damaged"):
            results[y][name]["BH_q_four_year_damage_strata"]=qs[y+"_"+name]
    return results


def within_mother_apical_lateral(groups):
    output={y:{} for y in YEARS}
    for y in YEARS:
        for arm in ARMS:
            paired=[r for r in groups[y][arm] if r["apical_damage"] is not None]
            d=[r["mean_damage"]-r["apical_damage"] for r in paired]
            eps=1e-12
            output[y][arm]={
               "n_mothers_with_both_position_damage_scores":len(d),
               "nonapical_minus_apical_damage_mother_mean":mean(d),
               "n_lateral_more_damaged":sum(x>eps for x in d),
               "n_lateral_less_damaged":sum(x<-eps for x in d),
               "n_equal":sum(abs(x)<=eps for x in d),
            }
    return output


def analyze(mothers,heads,bootstrap_reps=BOOT_REPS):
    grouped=source_groups(mothers)
    idx,clean,quality=source_audit(mothers,heads)
    records=real_mother_nonapical(mothers,clean)
    by=year_arm_records(records)
    paired=within_mother_apical_lateral(by)
    quality["n_mothers_with_scored_nonapical_heads"]=len(records)
    quality["n_apical_lateral_damage_pairs"]=sum(
        paired[y][arm]["n_mothers_with_both_position_damage_scores"]
        for y in YEARS for arm in ARMS)
    require(len(records)==87 and quality["n_apical_lateral_damage_pairs"]==82
            and quality["n_matched_heads_with_nonmissing_damage_score"]==290,
            "SOURCE_OBSERVED_LATERAL_DAMAGE_GRAIN_CHANGED")
    results={
       "version":"aza3_20261022_undulatum_nonapical_head_2007_2008_source_JS_v1",
       "source":{
           "repo":"zuizui0223/aza3",
           "original_head_git_blob":BLOBS["flower_heads.csv"],
           "original_mother_git_blob":BLOBS["mother_plants.csv"],
           "original_dictionary_git_blob":BLOBS["source_dictionary.txt"],
           "dataset_doi":SOURCE_DOI,
           "published_paper_doi":PAPER_DOI,
       },
       "source_audit":quality,
       "head_record_group_by_year_treatment":head_arm_summary(clean,mothers),
       "mother_balanced_nonapical_head_groups":{
           y:{a:group_mother_stats(by[y][a]) for a in ARMS}
           for y in YEARS},
       "year_apical_strata_reduced_minus_ambient_mother_cluster_bootstrap":
           cluster_bootstrap(by,bootstrap_reps,BOOT_SEED),
       "bootstrap":{
           "n_draws":bootstrap_reps,
           "seed":BOOT_SEED,
           "PRNG":"xorshift32 with unsigned modulo",
           "sample_order":"replicate > year ascending > source treatment arm ndni,dni,ndi,di > mother IDs sorted lexicographically",
           "percentile_interpolation":"linear over 0-based sorted rank (N-1)*p",
           "n_biological_sample_unit_mothers":106,
       },
       "original_mothers_any_lateral_good_seed_by_year_arm":
           mother_any_lateral(grouped),
       "exact_maternal_any_lateral_seed_comparisons":
           exact_any_lateral(mother_any_lateral(grouped)),
       "within_mother_apical_vs_lateral_damage":paired,
       "inference_limits":{
           "individual_source_head_rows_are_not_independent_mother_replicates":True,
           "head_damage_score_0p88_is_max_coded_score_not_direct_proportion_of_florets_lost":True,
           "head_damage_score_is_measured_after_access_treatment_and_not_causal_independent_covariate":True,
           "original_head_labels_disagree_for_two_2007_mothers_not_repaired":True,
           "original_mothers_with_no_lateral_heads_not_imputed_damage_zero":True,
           "high_lateral_goodseed_success_under_insecticide_not_genetic_lego_effect":True,
           "head_rank_not_orientation_or_phyllary_spine_geometry":True,
           "randomized_main_design_not_verified_here_against_original_assignment_logs":True,
           "no_florivore_oviposition_pollinator_delivery_or_seed_germination_at_same_head":True,
           "new_biological_field_samples_added":0,
       },
    }
    return results


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mothers",required=True,type=Path)
    p.add_argument("--heads",required=True,type=Path)
    p.add_argument("--dictionary",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--bootstrap",type=int,default=BOOT_REPS)
    q=p.parse_args()
    mothers=read_locked(q.mothers,"mother_plants.csv",
        ("plt","year","dam","ins","trt","o.hds","plt.gsd"))
    heads=read_locked(q.heads,"flower_heads.csv",
        ("plt","year","dam","ins","hd","av.dam","hd.gsd"))
    dict_text=read_locked(q.dictionary,"source_dictionary.txt")
    for span in ("Apical Damage Treatment","Flowerhead position on plant",
                 "Average proportion flowerhead area damaged"):
        require(span in dict_text,"USDA_HEAD_SOURCE_DICTIONARY_CHANGED")
    result=analyze(mothers,heads,q.bootstrap)
    q.out.parent.mkdir(parents=True,exist_ok=True)
    q.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                               allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
       "status":"ACTUAL_297_ALIGNED_HEAD_ROWS_87_MATERNAL_LATERAL_CLUSTERS",
       "source_audit":result["source_audit"],
       "year2008_mother_cluster":result["mother_balanced_nonapical_head_groups"]["2008"],
       "year2008_lateral_damage_and_seed_cluster_bootstrap":
          result["year_apical_strata_reduced_minus_ambient_mother_cluster_bootstrap"]["2008"],
       "year2008_exact_any_lateral_seed_test":
          result["exact_maternal_any_lateral_seed_comparisons"]["2008"],
    },indent=2))


if __name__=="__main__":
    main()
