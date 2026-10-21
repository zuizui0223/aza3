#!/usr/bin/env python3
"""Actually observed, mother-linked capitulum compensation in Cirsium undulatum.

Source: West & Louda 2021, American Journal of Botany, 10.1002/ajb2.1744,
public US Department of Agriculture 10.15482/USDA.ADC/1522649 CC0.
2007 and 2008 original full source files: 106 maternal plants, 323 heads.

Treatments: dam = experimental apical damage d/nd; ins = later head
insecticide i vs water ni. Do NOT equate treatment d to physical loss of
every terminal head. Plant seed measure is visually 'good' seeds, not
germinated viable offspring. Experimental potential compensation in nonapical
heads != return to equal total fecundity. Pooled lateral seed share increases
mechanically when apical seeds are lost even without absolute lateral gain.

A four-arm randomized-looking source design is not enough to infer causal
cross-year effects, which are two different reproductive cohorts. Site,
maternal size and head development are possible confounders. All tests are
descriptive stratified bootstrap (plants, NEVER heads as independent n).
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from statistics import median

BASE="data/external/usda_wavyleaf_2007_2008/"
SOURCE_BLOBS={
 "mother_plants.csv":"458800272c3d0a03c7653702027a8f8578d9e88c",
 "flower_heads.csv":"637a4dcd5194d0bf5b26e7f9cf0dde56a50d7616",
 "source_dictionary.txt":"5e20c9d4609abc982287bd12f739c76d7c6d1d2f",
}
REPO="zuizui0223/aza3"
SOURCE_DOI="10.15482/USDA.ADC/1522649"
PAPER_DOI="10.1002/ajb2.1744"
YEARS=("2007","2008")
ARMS=("ndni","dni","ndi","di")
METRICS=("total","lateral","apical","flowered","has_seeds")
CONTRAST_NAMES=(
 "damage_under_reduced_insects",
 "damage_under_ambient_insects",
 "insecticide_under_apical_damage",
 "insecticide_under_protected_apex",
 "damage_insecticide_interaction",
)
MOTHER_HEADERS=("plt","year","av.dam","no.flw","n.fl.sds","plt.gsd",
 "dam","ins","trt","ap.hd","o.hds","n.hds","r.diam",
 "n.br","plt.ht","pnt.flw","pnt.ohds","sds.flhd","sd.fh.sd")
HEAD_HEADERS=("plt","year","dam","ins","hd","av.dam","hd.gsd")
SENSITIVITY_EXCLUSIONS=("2007|323","2007|375","2007|327")
BOOT_SEED=20261021
BOOT_REPS=10000
MASK=(1<<32)-1


def require(ok,code):
    if not ok: raise ValueError(code)


def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def load(path,expected_name,required_fields=None):
    raw=Path(path).read_bytes()
    require(git_blob(raw)==SOURCE_BLOBS[expected_name],
            "GIT_BLOB_SOURCE_DRIFT_"+expected_name)
    if required_fields is None:
        return raw.decode("utf-8-sig")
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    require(tuple(reader.fieldnames or ())==required_fields,
            "SOURCE_CSV_HEADER_CHANGED_"+expected_name)
    rows=list(reader)
    require(all(None not in r and all(v is not None for v in r.values())
                for r in rows),"MALFORMED_SOURCE_CSV_ROW")
    return rows


def plant_id(r):
    return str(r["year"])+"|"+str(r["plt"])


def head_id(r):
    return plant_id(r)+"|"+str(r["hd"])


def number(r,key):
    x=r.get(key)
    require(x is not None and x not in ("", "NA", "na", "N/A", "."),
            "MISSING_NUMERIC_PLANT_OUTCOME_"+key)
    try:
        y=float(x)
    except (TypeError, ValueError):
        raise ValueError("INVALID_NUMERIC_"+key) from None
    require(math.isfinite(y) and y>=0,"INVALID_NUMERIC_"+key)
    return y


def metric(r,k):
    if k=="total":return number(r,"plt.gsd")
    if k=="lateral":return number(r,"o.hds")
    if k=="apical":return number(r,"ap.hd")
    if k=="flowered":return number(r,"no.flw")
    if k=="has_seeds":return float(number(r,"plt.gsd")>0)
    raise ValueError("UNKNOWN_OUTCOME")


def mean(vals):
    require(len(vals)>0,"EMPTY_GROUP")
    return sum(vals)/len(vals)


def group_stats(rows):
    require(rows,"EMPTY_TREATMENT_ARM")
    result={"n":len(rows)}
    for k in METRICS:
        y=[metric(r,k) for r in rows]
        result[k]=mean(y)
        if k=="total":
            result["total_median"]=float(median(y))
            result["total_sum"]=sum(y)
            result["zero_seed_mothers"]=sum(z==0 for z in y)
    require(result["total_sum"]>0,"UNDEFINED_POOLED_LATERAL_SHARE")
    result["lateral_seed_share_of_pooled_good_seed_total"]=(
       sum(metric(r,"lateral") for r in rows)/result["total_sum"])
    return result


def contrast(groupmeans,metric_name):
    x={a:groupmeans[a][metric_name] for a in ARMS}
    return {
        "damage_under_reduced_insects":x["di"]-x["ndi"],
        "damage_under_ambient_insects":x["dni"]-x["ndni"],
        "insecticide_under_apical_damage":x["di"]-x["dni"],
        "insecticide_under_protected_apex":x["ndi"]-x["ndni"],
        "damage_insecticide_interaction":
            (x["di"]-x["ndi"])-(x["dni"]-x["ndni"]),
    }


def xs32(seed):
    require(isinstance(seed,int) and 0<seed<=MASK,"RNG_SEED_INVALID")
    state=seed
    while True:
        state ^= (state<<13)&MASK
        state ^= state>>17
        state ^= (state<<5)&MASK
        state &=MASK
        yield state


def percentile(values,q):
    require(values and 0<=q<=1,"INVALID_PERCENTILE")
    v=sorted(values)
    t=(len(v)-1)*q
    lo=math.floor(t)
    frac=t-lo
    return v[lo]*(1-frac)+v[min(len(v)-1,lo+1)]*frac


def year_groups(mothers):
    groups={year:{arm:[r for r in mothers if r["year"]==year
                            and r["trt"]==arm] for arm in ARMS}
            for year in YEARS}
    require({y:{a:len(groups[y][a]) for a in ARMS} for y in YEARS}=={
             "2007":{"ndni":9,"dni":9,"ndi":10,"di":11},
             "2008":{"ndni":15,"dni":12,"ndi":21,"di":19}},
            "ORIGINAL_FACTORIAL_ARM_SAMPLE_COUNTS_CHANGED")
    return groups


def bootstrap_contrasts(groups,reps=BOOT_REPS,seed=BOOT_SEED):
    require(100<=reps<=30000,"BOOTSTRAP_BUDGET_INVALID")
    point={y:{m:contrast({a:group_stats(groups[y][a]) for a in ARMS},m)
              for m in METRICS} for y in YEARS}
    samplers={
        y:{m:{c:[] for c in CONTRAST_NAMES} for m in METRICS}
        for y in YEARS}
    rng=xs32(seed)
    for _ in range(reps):
        for y in YEARS:
            group={ }
            for arm in ARMS:
                original=groups[y][arm]
                totals={m:0. for m in METRICS}
                for draw in range(len(original)):
                    source=original[next(rng)%len(original)]
                    for m in METRICS:
                        totals[m]+=metric(source,m)
                group[arm]={m:totals[m]/len(original) for m in METRICS}
            for m in METRICS:
                for name,val in contrast(group,m).items():
                    samplers[y][m][name].append(val)
    return {
        y:{m:{name:{
             "estimate":point[y][m][name],
             "bootstrap_percentile_CI95":[
                 percentile(samplers[y][m][name],.025),
                 percentile(samplers[y][m][name],.975)
             ],
        } for name in CONTRAST_NAMES} for m in METRICS}
        for y in YEARS
    }


def validate_mothers(mothers):
    require(len(mothers)==106,"USDA_ORIGINAL_MOTHER_N_CHANGED")
    ids=[plant_id(x) for x in mothers]
    require(len(ids)==len(set(ids)),"DUPLICATE_MOTHER_YEAR_PLANT_ID")
    for r in mothers:
        require(r["year"] in YEARS and r["dam"] in ("d","nd")
                and r["ins"] in ("i","ni")
                and r["trt"]==r["dam"]+r["ins"],
                "MOTHER_TREATMENT_CODING_CHANGED")
        for key in ("plt.gsd","o.hds","ap.hd","no.flw","n.fl.sds"):
            number(r,key)
    return mothers


def quality(mothers,heads):
    require(len(heads)==323,"USDA_ORIGINAL_FLOWERHEAD_N_CHANGED")
    x={plant_id(r):r for r in mothers}
    hid=[head_id(h) for h in heads]
    require(len(set(hid))==len(hid),"DUPLICATE_SOURCE_HEAD_ID")
    groups={}
    for head in heads:
        require(head["year"] in YEARS and
                head["dam"] in ("d","nd") and head["ins"] in ("i","ni"),
                "HEAD_TREATMENT_CATEGORY_CHANGED")
        number(head,"hd.gsd")
        groups.setdefault(plant_id(head),[]).append(head)
    matches=[h for h in heads if plant_id(h) in x]
    orphan=[h for h in heads if plant_id(h) not in x]
    conflicts=[h for h in matches if
               h["dam"]!=x[plant_id(h)]["dam"] or
               h["ins"]!=x[plant_id(h)]["ins"]]
    mids=[r for r in mothers if plant_id(r) not in groups]
    full_bad=[r for r in mothers if
              plant_id(r) in groups and
              sum(number(h,"hd.gsd") for h in groups[plant_id(r)])!=
              number(r,"plt.gsd")]
    part_bad=[r for r in mothers if
              number(r,"ap.hd")+number(r,"o.hds")!=
              number(r,"plt.gsd")]
    out={
      "n_distinct_mothers":len(mothers),
      "n_source_head_records":len(heads),
      "n_unique_head_ids":len(hid),
      "n_heads_with_matching_mother":len(matches),
      "n_head_records_without_mother":len(orphan),
      "n_head_only_parent_ids":len({plant_id(h) for h in orphan}),
      "head_only_parent_ids":sorted({plant_id(h) for h in orphan}),
      "mother_ids_with_no_head_records":sorted(plant_id(r) for r in mids),
      "n_mothers_with_matching_head_seed_sum_except_absent_head_records":
        len(mothers)-len(mids)-len(full_bad),
      "n_matched_parent_head_sum_mismatches":len(full_bad),
      "n_head_rows_treatment_label_conflicts":len(conflicts),
      "treatment_label_conflict_mother_ids":
          sorted({plant_id(h) for h in conflicts}),
      "n_mothers_with_apical_plus_lateral_accounting_mismatch":
          len(part_bad),
      "component_mismatch_mother_ids":sorted(plant_id(r) for r in part_bad),
      "plant_vs_head_flowering_count_not_equivalent":
        "Head file includes flowerhead records not necessarily equal to the plant-level number reaching anthesis; do not equate row count to no.flw",
      "head_level_conflict_no_auto_reassignment":True,
    }
    require(len(matches)==309 and len(orphan)==14 and
            len(out["head_only_parent_ids"])==6 and
            out["mother_ids_with_no_head_records"]==
               ["2007|329","2007|367"] and
            len(full_bad)==0 and len(conflicts)==12 and
            out["treatment_label_conflict_mother_ids"]==
              ["2007|323","2007|375"] and
            out["component_mismatch_mother_ids"]==["2007|327"],
            "USDA_EXTERNAL_HEAD_JOIN_OR_ACCOUNTING_DRIFT")
    return out


def analyze(mothers,heads,bootstrap_reps=BOOT_REPS):
    validate_mothers(mothers)
    qa=quality(mothers,heads)
    groups=year_groups(mothers)
    stats={y:{arm:group_stats(groups[y][arm]) for arm in ARMS} for y in YEARS}
    outcomes=bootstrap_contrasts(groups,bootstrap_reps,BOOT_SEED)
    excluded=set(SENSITIVITY_EXCLUSIONS)
    sensitivity={}
    for y in YEARS:
        sub={a:group_stats([r for r in groups[y][a]
                            if plant_id(r) not in excluded])
             for a in ARMS}
        sensitivity[y]=sub
        sensitivity[y]["damage_effects"]={
            m:{
              "reduced":sub["di"][m]-sub["ndi"][m],
              "ambient":sub["dni"][m]-sub["ndni"][m],
            } for m in ("total","lateral","apical","flowered")
        }
    return {
      "version":"cirsium_undulatum_2007_2008_original_mother_head_python_analysis_v1",
      "source":{
        "doi":SOURCE_DOI,
        "paper_doi":PAPER_DOI,
        "license":"CC0",
        "git_blob_source_locks":SOURCE_BLOBS,
      },
      "cohort":{
        "n_mothers":len(mothers),"n_heads":len(heads),
        "year_n":{y:sum(r["year"]==y for r in mothers) for y in YEARS},
        "year_treatment_mothers":{
            y:{a:len(groups[y][a]) for a in ARMS} for y in YEARS},
        "head_rows_are_not_independent_maternal_replicates":True,
      },
      "measured_group_outcomes":stats,
      "year_specific_contrasts_plant_bootstrap":outcomes,
      "bootstrap_contract":{
         "n_draws":bootstrap_reps,"seed":BOOT_SEED,
         "prng":"xorshift32 unsigned mask then modulo group size",
         "unit":"MOTHER_PLANT_WITHIN_YEAR_TREATMENT",
         "resample_group_order":list(ARMS),
         "resample_year_order":list(YEARS),
         "outcome_order":list(METRICS),
         "CI_is_descriptive_stratum_bootstrap_not_randomization_based":True,
      },
      "source_reconciliation":qa,
      "sensitivity_drop_two_treatment_mismatch_and_one_component_mismatch":sensitivity,
      "claims":{
        "two_years_are_distinct_individual_cohorts":True,
        "good_seeds_are_not_seed_germination_assay":True,
        "apical_damage_treatment_does_not_guarantee_all_apical_seed_loss":True,
        "high_pooled_lateral_share_does_not_prove_absolute_lateral_increase":True,
        "source_treatment_conflict_ids_are_not_silently_recoded":True,
        "this_reanalysis_not_a_discovery_of_compensation_first_reported_2021":True,
        "same_head_pollinator_oviposition_and_seed_fitness_measured":False,
        "genetically_selected_lego_or_ancestral_memory_identified":False,
        "external_independent_plant_samples_newly_collected":0,
      },
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mothers",type=Path,required=True)
    ap.add_argument("--heads",type=Path,required=True)
    ap.add_argument("--dictionary",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--bootstrap",type=int,default=BOOT_REPS)
    args=ap.parse_args()
    moms=load(args.mothers,"mother_plants.csv",MOTHER_HEADERS)
    heads=load(args.heads,"flower_heads.csv",HEAD_HEADERS)
    dictionary=load(args.dictionary,"source_dictionary.txt")
    for t in ("Apical Damage Treatment","Insecticide Treatment",
              "Flowerhead position on plant","per plant"):
        require(t in dictionary,"PUBLIC_SOURCE_DICTIONARY_CONTENT_DRIFT")
    z=analyze(moms,heads,args.bootstrap)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+"\n",
                        encoding="utf-8")
    print(json.dumps({
     "status":"ACTUAL_CC0_MOTHER_AND_HEAD_OBSERVATIONS_REPLAYED",
     "n_mothers":z["cohort"]["n_mothers"],"n_heads":z["cohort"]["n_heads"],
     "source_join":z["source_reconciliation"],
     "contrasts":{
         y:{m:{k:round(d["estimate"],4) for k,d in
                   z["year_specific_contrasts_plant_bootstrap"][y][m].items()}
            for m in ("total","lateral","has_seeds")}
         for y in YEARS},
    },indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
