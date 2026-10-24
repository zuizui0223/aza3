#!/usr/bin/env python3
"""Independent genuine Ipomopsis outgroup: pollination stage vs seed-predator escape.

Reanalysis of Campbell et al. (2022), American Naturalist 199:808-823
DOI 10.1086/716740, open Zenodo DOI 10.7280/D1KM49,
file Selection_by_pollinators_and_predators.csv, MD5 as published by the archive:
12ed3dbe6d269c5f808f9c6b230d468c.

KEY SCIENTIFIC BOUNDARIES:
  - The 437 original rows combine plant species/hybrids, 2010–18 gardens
    and year. No original plant identifier is supplied. Treat source rows,
    not proven independent individuals or independent evolutionary origins.
  - Original author "Type" is A=I. aggregata, H=hybrid, T=I. tenuituba.
    "Site" source hyb & ten are parts of the original upper site; both the
    original 3-label coding and author-merged upper site are sensitivities.
  - "width" is corolla width and "sepal" is sepal width, NOT Cirsium spines.
  - seedsinit is the number of seeds initiated (author estimate); flyescape
    is proportion of FRUITS escaping fly predation, NOT a per-seed guarantee.
    Do not impose seeds == seedsinit * flyescape: >initiated can occur for
    the separately measured author's final seed estimate.
  - Associations in strata Site×Year×Type are not intervention-specific
    selection gradients, genetic LEGO, phylogenetic transitions or causal
    trait-mediated fitness. Original Campbell paper used richer multi-trait
    selection models; these descriptive analyses do NOT refute them.
  - Shuffle original widths/sepals within source labels as a CONDITIONAL,
    exploratory association benchmark, not random allocation of plant sex,
    pollinators or predation. Use 4-test BH correction as guard, not proof.
  - Reproduces independent JS receipt from original source downloaded live.
"""
from __future__ import annotations

import argparse
from collections import defaultdict, Counter
import csv
import hashlib
import io
import json
import math
from pathlib import Path

EXPECTED_MD5="12ed3dbe6d269c5f808f9c6b230d468c"
DATA_DOI="10.7280/D1KM49"
PAPER_DOI="10.1086/716740"
SOURCE_URL="https://zenodo.org/records/5083399"
NAME="Selection_by_pollinators_and_predators.csv"
READ_ME="ReadMe_Selection_by_pollinators_and_predators.csv"
COLS=["Site","Year","Type","width","color","nectar","sepal",
      "totpinene","indoleN","seeds","seedsinit","flyescape"]
SITES=["agg","hyb","ten"]
YEARS=list(range(2010,2019))
TYPES=["A","H","T"]
NUMERIC=("width","color","nectar","sepal","totpinene",
         "indoleN","seeds","seedsinit","flyescape")
SEED_WIDTH=20261024
SEED_SEPAL=20261025
N_PERM=10000
N_BOOT=3000
MASK=0xffffffff

def check(ok,error):
    if not ok:
        raise ValueError(error)

def load(raw,verify_hash=True):
    if verify_hash:
        check(hashlib.md5(raw).hexdigest()==EXPECTED_MD5,
              "ORIGINAL_ZENODO_MD5_MISMATCH")
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig"),newline=""))
    check(list(reader.fieldnames or [])==COLS,
          "ORIGINAL_ZENODO_12_COLUMNS_CHANGED")
    rows=[]
    source_strings=[]
    for i,r in enumerate(reader,2):
        check(None not in r and all(x is not None for x in r.values()),
              "ORIGINAL_ZENODO_SOURCE_WIDTH_CHANGED")
        source_strings.append(tuple(r[k] for k in COLS))
        o={"Site":r["Site"],"Year":int(r["Year"]), "Type":r["Type"],
           "original_file_row":i}
        check(o["Site"] in SITES and o["Year"] in YEARS and
              o["Type"] in TYPES,"SOURCE_CONTEXT_INVALID")
        for k in NUMERIC:
            if r[k]=="." or r[k]=="":
                o[k]=None
            else:
                try:
                    val=float(r[k])
                except ValueError:
                    raise ValueError("SOURCE_NUMERIC_INVALID_"+k) from None
                check(math.isfinite(val),"SOURCE_NONFINITE_"+k)
                o[k]=val
        check(o["flyescape"] is None or 0<=o["flyescape"]<=1,
              "ORIGINAL_FRUIT_FLY_ESCAPE_OUT_OF_BOUNDS")
        check(o["seedsinit"] is None or o["seedsinit"]>=0,
              "ORIGINAL_SEEDS_INITIATED_NEGATIVE")
        check(o["seeds"] is None or o["seeds"]>=0,
              "ORIGINAL_SEEDS_FINAL_NEGATIVE")
        rows.append(o)
    check(len(rows)==437,"ORIGINAL_437_SOURCE_ROWS_CHANGED")
    check(len(source_strings)==len(set(source_strings)),
          "ORIGINAL_EXACT_DUPLICATE_SOURCE_ROW_CHANGED")
    return rows

def mean(arr):
    check(len(arr)>0,"MEAN_EMPTY")
    return sum(arr)/len(arr)

def pearson(x,y):
    check(len(x)==len(y),"PEARSON_SHAPE_MISMATCH")
    if len(x)<4:
        return None
    ax=mean(x);ay=mean(y)
    c=xx=yy=0.0
    for a,b in zip(x,y):
        a-=ax;b-=ay
        c+=a*b;xx+=a*a;yy+=b*b
    if not xx or not yy:
        return None
    return c/math.sqrt(xx*yy)

def quantile(arr,q):
    check(arr and 0<=q<=1,"INVALID_BOOTSTRAP_QUANTILE")
    seq=sorted(arr)
    ix=(len(seq)-1)*q
    i=math.floor(ix);f=ix-i
    return seq[i]*(1-f)+seq[min(i+1,len(seq)-1)]*f

def xorshift32(seed):
    check(type(seed) is int and 0<seed<2**32,
          "INVALID_OR_ZERO_XORSHIFT_SEED")
    state=seed
    while True:
        state^=(state<<13)&MASK
        state^=state>>17
        state^=(state<<5)&MASK
        state&=MASK
        yield state

def source_groups(rows,trait,sitecode="original",drop_year=None):
    check(trait in ("width","sepal"),"ONLY_MEASURED_HEAD_ANALOGUES_TESTED")
    check(sitecode in ("original","merged"),"INVALID_ORIGINAL_SITE_CODE")
    raw=[r for r in rows if r[trait] is not None and
         r["seedsinit"] is not None and r["flyescape"] is not None and
         (drop_year is None or r["Year"]!=drop_year)]
    groups={}
    for r in raw:
        site=(r["Site"] if sitecode=="original" else
              ("agg" if r["Site"]=="agg" else "upper"))
        key=f'{site}|{r["Year"]}|{r["Type"]}'
        groups.setdefault(key,[]).append(r)
    eligible=[]
    for key,rr in groups.items():
        if len(rr)<2 or len({z[trait] for z in rr})<=1:
            continue
        eligible.append({
            "key":key,
            "rows":[{
                "file_row":z["original_file_row"],
                "x":z[trait],
                "init":math.log1p(z["seedsinit"]),
                "escape":z["flyescape"],
            } for z in rr],
        })
    return raw,eligible

def source_summary(rows,trait,sitecode="original",drop_year=None):
    raw,groups=source_groups(rows,trait,sitecode,drop_year)
    xx=[r[trait] for r in raw]
    yy=[math.log1p(r["seedsinit"]) for r in raw]
    zz=[r["flyescape"] for r in raw]
    centered=[]
    sizes=[]
    for g in groups:
        rr=g["rows"]
        mx=mean([r["x"] for r in rr])
        my=mean([r["init"] for r in rr])
        mz=mean([r["escape"] for r in rr])
        sizes.append(len(rr))
        centered.extend({
            "x":r["x"]-mx,
            "init":r["init"]-my,
            "escape":r["escape"]-mz,
        } for r in rr)
    return {
        "n_source_complete_records":len(raw),
        "n_within_group_comparable_records":len(centered),
        "n_groups":len(groups),
        "source_group_size_min":min(sizes),
        "source_group_size_max":max(sizes),
        "raw_correlations":{
            "width_or_sepal_vs_log1p_initiated":pearson(xx,yy),
            "width_or_sepal_vs_fly_escape":pearson(xx,zz),
        },
        "within_group_correlations":{
            "trait_vs_log1p_initiated":pearson(
                [r["x"] for r in centered],
                [r["init"] for r in centered]),
            "trait_vs_fly_escape":pearson(
                [r["x"] for r in centered],
                [r["escape"] for r in centered]),
        },
    }

def conditional_inference(rows,trait,seed,perms=N_PERM,boots=N_BOOT):
    check(perms==N_PERM and boots==N_BOOT,
          "FREEZE_NUMERICAL_INFERENCE_BUDGET")
    raw,source=source_groups(rows,trait)
    config=[]
    for g in source:
        rr=g["rows"]
        mu=mean([r["x"] for r in rr])
        my=mean([r["init"] for r in rr])
        mz=mean([r["escape"] for r in rr])
        config.append({
            "key":g["key"],"n":len(rr),
            "xs":[r["x"]-mu for r in rr],
            "yi":[r["init"]-my for r in rr],
            "ye":[r["escape"]-mz for r in rr],
            "raw":rr,
        })
    svx=svy=sve=cov1=cov2=0.0
    for g in config:
        for x,y,e in zip(g["xs"],g["yi"],g["ye"]):
            svx+=x*x;svy+=y*y;sve+=e*e
            cov1+=x*y;cov2+=x*e
    check(svx>0 and svy>0 and sve>0,
          "NO_WITHIN_SOURCE_STUDY_VARIABILITY")
    r1=cov1/math.sqrt(svx*svy)
    r2=cov2/math.sqrt(svx*sve)
    delta=r1-r2
    rng=xorshift32(seed)
    count1=count2=count_difference=0
    for _ in range(perms):
        cv1=cv2=0.0
        for g in config:
            ids=list(range(g["n"]))
            for j in range(len(ids)-1,0,-1):
                k=next(rng)%(j+1)
                ids[j],ids[k]=ids[k],ids[j]
            for j in range(g["n"]):
                cv1+=g["xs"][ids[j]]*g["yi"][j]
                cv2+=g["xs"][ids[j]]*g["ye"][j]
        v1=cv1/math.sqrt(svx*svy)
        v2=cv2/math.sqrt(svx*sve)
        count1+=abs(v1)>=abs(r1)-1e-12
        count2+=abs(v2)>=abs(r2)-1e-12
        count_difference+=abs(v1-v2)>=abs(delta)-1e-12
    drawn1=[];drawn2=[];drawn_diff=[]
    for _ in range(boots):
        vxx=vyy=vee=vxy=vxe=0.0
        for g in config:
            sampled=[g["raw"][next(rng)%g["n"]] for j in range(g["n"])]
            mx=mean([r["x"] for r in sampled])
            my=mean([r["init"] for r in sampled])
            me=mean([r["escape"] for r in sampled])
            for r in sampled:
                dx=r["x"]-mx;dy=r["init"]-my;de=r["escape"]-me
                vxx+=dx*dx;vyy+=dy*dy;vee+=de*de
                vxy+=dx*dy;vxe+=dx*de
        check(vxx>0 and vyy>0 and vee>0,"BOOTSTRAP_SOURCE_DEGENERATE")
        z1=vxy/math.sqrt(vxx*vyy)
        z2=vxe/math.sqrt(vxx*vee)
        drawn1.append(z1);drawn2.append(z2);drawn_diff.append(z1-z2)
    return {
        "design":{
            "source_eligible_strata":[{
                "key":g["key"],
                "n":g["n"],
                "source_file_rows":[r["file_row"] for r in g["raw"]],
            } for g in config],
            "permutations":perms,"bootstrap_draws":boots,"seed":seed,
            "prng":"xorshift32",
            "stratum_field_order":"Site|Year|Type",
            "individual_record_id_not_supplied_by_author":True,
            "unknown_repeat_individual_possible":True,
            "outcome_transform":"log1p(seedsinit) and raw flyescape fraction",
            "sitecode":(
                "original three source site labels (original README notes hyb "
                "and ten same upper site)"),
        },
        "summary":source_summary(rows,trait),
        "r_init":r1,"r_escape":r2,
        "r_init_minus_r_escape":delta,
        "permutation_p_two_sided":{
            "init":(count1+1)/(perms+1),
            "escape":(count2+1)/(perms+1),
            "stage_difference":(count_difference+1)/(perms+1),
        },
        "within_stratum_row_resample_ci95":{
            "init":[quantile(drawn1,.025),quantile(drawn1,.975)],
            "escape":[quantile(drawn2,.025),quantile(drawn2,.975)],
            "stage_difference":[quantile(drawn_diff,.025),
                                quantile(drawn_diff,.975)],
        },
    }

def analyze(raw):
    rows=load(raw)
    width=conditional_inference(rows,'width',SEED_WIDTH)
    sepal=conditional_inference(rows,'sepal',SEED_SEPAL)
    ps=[
        ("width_init",width["permutation_p_two_sided"]["init"]),
        ("width_escape",width["permutation_p_two_sided"]["escape"]),
        ("sepal_init",sepal["permutation_p_two_sided"]["init"]),
        ("sepal_escape",sepal["permutation_p_two_sided"]["escape"]),
    ]
    ps.sort(key=lambda v:v[1])
    qs={}
    latest=1.0
    for j in range(len(ps)-1,-1,-1):
        key,p=ps[j]
        latest=min(latest,p*len(ps)/(j+1))
        qs[key]=latest
    years={str(year):source_summary(rows,'width','original',year)
           for year in YEARS}
    audit={
        "n_original_data_records":len(rows),
        "n_original_fields":len(COLS),
        "source_header":list(COLS),
        "missing_count_by_column":{
            k:sum(r[k] is None for r in rows) for k in COLS
        },
        "site_counts":dict(Counter(r["Site"] for r in rows)),
        "year_counts":dict(Counter(str(r["Year"]) for r in rows)),
        "type_counts":dict(Counter(r["Type"] for r in rows)),
        "seeds_greater_than_seedsinit":sum(
            r["seeds"] is not None and r["seedsinit"] is not None
            and r["seeds"]>r["seedsinit"]+1e-8 for r in rows),
        "source_plant_identifier_field_exists":False,
        "source_exact_duplicate_full_rows":0,
    }
    check(audit["site_counts"]=={"agg":221,"hyb":144,"ten":72}
          and audit["type_counts"]=={"H":246,"A":175,"T":16}
          and audit["missing_count_by_column"]["sepal"]==329
          and audit["seeds_greater_than_seedsinit"]==41,
          "IPOMOPSIS_MISSINGNESS_STAGE_DEFINITION_DRIFT")
    check(width["summary"]["n_within_group_comparable_records"]==354
          and sepal["summary"]["n_within_group_comparable_records"]==90,
          "IPOMOPSIS_WITHIN_STRATUM_SAMPLE_DRIFT")
    return {
        "schema":"aza3_general_cross_stage_original_ipomopsis_source_v1",
        "as_of":"2026-10-24",
        "source":{
            "publication":"Campbell et al (2022), American Naturalist",
            "doi":PAPER_DOI,"open_data_doi":DATA_DOI,
            "zenodo_record":SOURCE_URL,
            "author_csv":NAME,
            "metadata_author_md5":EXPECTED_MD5,
            "author_readme":READ_ME,
            "published_author_data_not_an_aza3_field_study":True,
        },
        "audit":audit,
        "analyses":{
            "width":width,"sepal":sepal,
            "four_analyses_multiplicity_BH_q":qs,
            "merged_actual_upper_site_sensitivity":{
                "width":source_summary(rows,'width','merged'),
                "sepal":source_summary(rows,'sepal','merged'),
            },
            "leave_one_year_out_width":years,
        },
        "limits":{
            "Ipomopsis_not_Cirsium":True,
            "traits_are_corolla_and_sepal_width_not_Cirsium_spines":True,
            "observational_within_source_strata_not_true_trait_randomized_intervention":True,
            "independent_pollinator_visitation_and_oviposition_video_not_in_source":True,
            "source_type_A_H_T_is_taxon_class_not_three_independent_lineages":True,
            "original_source_has_no_unique_plant_identifier":True,
            "staged_outcomes_related_not_seeds_equals_initiated_times_flyescape":True,
            "flyescape_measures_fraction_of_fruits_avoiding_predation":True,
            "selection_model_after_original_article_not_reconstructed":True,
            "permutation_exchangeability_not_proven":True,
            "bootstrap_within_site_year_type_descriptive_only":True,
            "analyses_exploratory_after_source_selection":True,
            "new_own_field_observations":0,
        },
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    result=analyze(a.csv.read_bytes())
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                                allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "source":result["audit"],
        "width":{
            "n":result["analyses"]["width"]["summary"]["n_within_group_comparable_records"],
            "r_initiation":result["analyses"]["width"]["r_init"],
            "r_escape":result["analyses"]["width"]["r_escape"],
            "conditional_p":result["analyses"]["width"]["permutation_p_two_sided"],
        },
        "sepal":{
            "n":result["analyses"]["sepal"]["summary"]["n_within_group_comparable_records"],
            "r_escape":result["analyses"]["sepal"]["r_escape"],
            "conditional_p":result["analyses"]["sepal"]["permutation_p_two_sided"],
        },
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
