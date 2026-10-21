#!/usr/bin/env python3
"""Separate Cirsium undulatum (2007/08) seed-bearing probability from yield.

Actual original USDA CC0 mother plant CSV 106 records, West & Louda (2021)
10.1002/ajb2.1744, original data 10.15482/USDA.ADC/1522649.
Original source text and cohort from aza3 PR #73, not newly sampled plants.

Whole-plant "good seed" is author-classified seeds, NOT germination success.
A 2-part arithmetic identity E[Y] = P(Y>0) * E[Y|Y>0] is descriptive,
NOT a statistical mediation decomposition of an assigned treatment.
Fixed contrast exact Fisher uses plants as n and 8 BH-adjusted comparisons.
Log1p count OLS+HC3 on 56 (2008) or 36 (2007) mothers with measured
rosette diameter is SENSITIVITY: timing of diameter measurement and the
missingness process are not verified. Do not infer size-adjusted causation.
Zero-inflated 2008 cells (insecticide positive 40/40) preclude naive logit.

The official USDA dictionary says no.flw / n.fl.sds count flowered / seed-
bearing heads, but 10 / 17 mothers respectively show positive *apical-only*
seeds despite zero counts. These fields may exclude apical heads in practice.
We flag this as SEMANTIC AMBIGUITY, not fabrication or 10 errors.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import math
from pathlib import Path

DOI="10.15482/USDA.ADC/1522649"
PAPER_DOI="10.1002/ajb2.1744"
BLOBS={
    "mother_plants.csv":"458800272c3d0a03c7653702027a8f8578d9e88c",
    "flower_heads.csv":"637a4dcd5194d0bf5b26e7f9cf0dde56a50d7616",
    "source_dictionary.txt":"5e20c9d4609abc982287bd12f739c76d7c6d1d2f",
}
YEARS=("2007","2008")
ARMS=("ndni","dni","ndi","di")
CONTRAST={
    "insecticide_intact":("ndi","ndni"),
    "insecticide_damaged":("di","dni"),
    "damage_ambient":("dni","ndni"),
    "damage_reduced":("di","ndi"),
}
LINEAR_CONTRAST={
    "insecticide_intact":(0,0,1,0),
    "insecticide_damaged":(0,0,1,1),
    "damage_ambient":(0,1,0,0),
    "damage_reduced":(0,1,0,1),
}
POSITIVE_ZERO_FLOWERING_ID_COUNT=10
POSITIVE_ZERO_SEED_HEAD_ID_COUNT=17


def require(ok,reason):
    if not ok:
        raise ValueError(reason)


def git_blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()


def read_locked(path,name,required=None):
    raw=Path(path).read_bytes()
    require(git_blob(raw)==BLOBS[name],"ORIGINAL_USDA_GIT_BLOB_DRIFT_"+name)
    if required is None:
        return raw.decode("utf-8-sig")
    rd=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    rows=list(rd)
    require(set(required).issubset(rd.fieldnames or ()),
            "ORIGINAL_SOURCE_COLUMNS_MISSING_"+name)
    require(all(None not in row and all(v is not None for v in row.values())
                for row in rows),"ORIGINAL_SOURCE_CSV_MALFORMED")
    return rows


def value(r,k):
    s=r.get(k)
    if s is None or str(s).strip().upper() in ("","NA","N/A","."):
        return None
    try:
        x=float(s)
    except (TypeError,ValueError):
        raise ValueError("INVALID_SOURCE_NUMBER_"+k) from None
    require(math.isfinite(x),"NONFINITE_SOURCE_NUMBER_"+k)
    return x


def mean(xs):
    require(bool(xs),"MEAN_EMPTY")
    return sum(xs)/len(xs)


def group_stats(rows):
    require(bool(rows),"EMPTY_SOURCE_GROUP")
    outcomes=[value(r,"plt.gsd") for r in rows]
    lateral=[value(r,"o.hds") for r in rows]
    require(all(x is not None and x>=0 for x in outcomes+lateral),
            "MISSING_OR_NEGATIVE_MATERNAL_SEEDS")
    positive=[x for x in outcomes if x>0]
    npos=len(positive)
    return {
        "n_mothers":len(rows),
        "n_seed_bearing_mothers":npos,
        "n_seed_free_mothers":len(rows)-npos,
        "p_any_good_seeds":npos/len(rows),
        "mean_good_seeds_all_mothers":mean(outcomes),
        "mean_good_seeds_given_positive":mean(positive) if positive else None,
        "mean_lateral_good_seeds_all_mothers":mean(lateral),
        "p_any_lateral_good_seeds":
            sum(x>0 for x in lateral)/len(lateral),
        "n_rosette_measurements":
            sum(value(r,"r.diam") is not None for r in rows),
        "mean_rosette_cm_if_measured":
            mean([value(r,"r.diam") for r in rows
                  if value(r,"r.diam") is not None])
            if any(value(r,"r.diam") is not None for r in rows) else None,
        "mean_good_seeds_if_size_known":
            mean([value(r,"plt.gsd") for r in rows
                  if value(r,"r.diam") is not None])
            if any(value(r,"r.diam") is not None for r in rows) else None,
    }


def logchoose(n,k):
    if k<0 or k>n:return -math.inf
    return (math.lgamma(n+1)-math.lgamma(k+1)
            -math.lgamma(n-k+1))


def fisher_two_sided(success1, fail1, success2, fail2):
    counts=(success1,fail1,success2,fail2)
    require(all(type(v)==int and v>=0 for v in counts)
            and success1+fail1>0 and success2+fail2>0,
            "INVALID_FISHER_GROUP_COUNTS")
    n=sum(counts)
    sample=success1+fail1
    successes=success1+success2
    def pmf(k):
        return math.exp(logchoose(successes,k)+
                        logchoose(n-successes,sample-k)-
                        logchoose(n,sample))
    observed=pmf(success1)
    p=sum(pmf(k) for k in range(max(0,sample-(n-successes)),
                                min(sample,successes)+1)
          if pmf(k) <= observed*(1+1e-9))
    return min(1.0,p)


def bh_qvalues(ps):
    names=list(ps)
    ordered=sorted(names,key=lambda name:ps[name])
    k=len(names)
    out={}
    tail=1.0
    for rank in range(k,0,-1):
        name=ordered[rank-1]
        p=ps[name]
        require(type(p) in (float,int) and
                math.isfinite(p) and 0<=p<=1,
                "INVALID_P_FOR_BH")
        tail=min(tail,k*p/rank)
        out[name]=tail
    return out


def exact_seed_presence_tests(groups):
    raw={}
    for y in YEARS:
        for name,(arm1,arm0) in CONTRAST.items():
            a=group_stats(groups[y][arm1])
            b=group_stats(groups[y][arm0])
            tab=[
                [a["n_seed_bearing_mothers"],a["n_seed_free_mothers"]],
                [b["n_seed_bearing_mothers"],b["n_seed_free_mothers"]]
            ]
            p=fisher_two_sided(*tab[0],*tab[1])
            raw[y+"_"+name]=p
    q=bh_qvalues(raw)
    return {y:{
        name:{
           "tested_arms":[arm1,arm0],
           "mothers_positive_negative":[
                [group_stats(groups[y][arm])["n_seed_bearing_mothers"],
                 group_stats(groups[y][arm])["n_seed_free_mothers"]]
                for arm in (arm1,arm0)],
           "difference_in_seed_bearing_probability":
                group_stats(groups[y][arm1])["p_any_good_seeds"]-
                group_stats(groups[y][arm0])["p_any_good_seeds"],
           "fisher_exact_two_sided_p":raw[y+"_"+name],
           "BH_q_across_eight_predefined_contrasts":q[y+"_"+name],
        }
        for name,(arm1,arm0) in CONTRAST.items()
    } for y in YEARS}


def hurdle_identity(groups):
    out={}
    for y in YEARS:
        out[y]={}
        for name,(arm1,arm0) in CONTRAST.items():
            a=group_stats(groups[y][arm1])
            b=group_stats(groups[y][arm0])
            # Comparisons here include zeros by design: if one group has no
            # positives, conditional mean decomposition is undefined, not 0.
            if a["n_seed_bearing_mothers"] and b["n_seed_bearing_mothers"]:
                prob_ratio=a["p_any_good_seeds"]/b["p_any_good_seeds"]
                positive_ratio=(
                   a["mean_good_seeds_given_positive"]/
                   b["mean_good_seeds_given_positive"])
                total_ratio=(a["mean_good_seeds_all_mothers"]/
                             b["mean_good_seeds_all_mothers"])
                require(math.isclose(prob_ratio*positive_ratio,total_ratio,
                                     rel_tol=0,abs_tol=1e-9),
                        "HURDLE_ARITHMETIC_IDENTITY_FAILED")
                out[y][name]={
                    "observed_mean_good_seed_ratio":total_ratio,
                    "probability_any_seed_ratio":prob_ratio,
                    "positive_mother_conditional_mean_seed_ratio":positive_ratio,
                    "product_identity_verified":True,
                    "reproductive_processes_are_causally_decomposed":False,
                }
            else:
                out[y][name]={
                    "observed_mean_good_seed_ratio":None,
                    "probability_any_seed_ratio":None,
                    "positive_mother_conditional_mean_seed_ratio":None,
                    "product_identity_verified":False,
                    "reproductive_processes_are_causally_decomposed":False,
                }
    return out


def solve_inverse(a):
    n=len(a)
    v=[[float(z) for z in row]+[
        1. if i==j else 0. for j in range(n)]
        for i,row in enumerate(a)]
    for j in range(n):
        pos=max(range(j,n),key=lambda i:abs(v[i][j]))
        require(abs(v[pos][j])>1e-11,"SINGULAR_DESIGN_MATRIX")
        v[j],v[pos]=v[pos],v[j]
        t=v[j][j]
        v[j]=[z/t for z in v[j]]
        for i in range(n):
            if i!=j:
                fac=v[i][j]
                v[i]=[v[i][k]-fac*v[j][k] for k in range(2*n)]
    return [row[n:] for row in v]


def dot(a,b):
    return sum(u*v for u,v in zip(a,b))


def fit_log1p_hc3(rows,adjust_for_rosette):
    """Sensitivity linear contrast, not a log-count GLM causal effect."""
    n=len(rows)
    npar=5 if adjust_for_rosette else 4
    require(n>=3*npar,"TOO_FEW_FOR_DESCRIPTIVE_REGRESSION")
    X=[]
    Y=[]
    for row in rows:
        d=int(row["dam"]=="d")
        ins=int(row["ins"]=="i")
        xi=[1.,float(d),float(ins),float(d*ins)]
        if adjust_for_rosette:
            rd=value(row,"r.diam")
            require(rd is not None and rd>0,"SIZE_MISSING_IN_ADJUSTED_ANALYSIS")
            xi.append(math.log(rd))
        X.append(xi)
        Y.append(math.log1p(value(row,"plt.gsd")))
    xx=[[sum(z[i]*z[j] for z in X) for j in range(npar)]
        for i in range(npar)]
    inv=solve_inverse(xx)
    xy=[sum(z[i]*y for z,y in zip(X,Y)) for i in range(npar)]
    beta=[dot(z,xy) for z in inv]
    resid=[y-dot(x,beta) for x,y in zip(X,Y)]
    meat=[[0.]*npar for _ in range(npar)]
    for x,e in zip(X,resid):
        h=dot(x,[dot(z,x) for z in inv])
        require(h<1-1e-8,"HC3_SATURATED_LEVERAGE")
        factor=e*e/(1-h)**2
        for i in range(npar):
            for j in range(npar):
                meat[i][j]+=x[i]*x[j]*factor
    tmp=[[sum(inv[i][k]*meat[k][j] for k in range(npar))
         for j in range(npar)] for i in range(npar)]
    cov=[[sum(tmp[i][k]*inv[j][k] for k in range(npar))
         for j in range(npar)] for i in range(npar)]
    results={}
    for name,c in LINEAR_CONTRAST.items():
        cc=list(c)+([0.] if adjust_for_rosette else [])
        est=dot(cc,beta)
        err=math.sqrt(max(0.,dot(cc,[dot(z,cc) for z in cov])))
        results[name]={
            "difference_log1p_good_seeds":est,
            "se_HC3":err,
            "normal_approximate_95_interval":[est-1.96*err,est+1.96*err],
            "exp_difference_is_ratio_of_geometric_1plus_counts":math.exp(est),
        }
    return {
        "n_mothers":n,
        "n_parameters":npar,
        "design":"intercept + damage + insecticide + interaction "+
                 ("+log(rosette)" if adjust_for_rosette else ""),
        "coefficient_order":[
            "intercept","apical_damage","insecticide","damage_by_insecticide"
        ]+(["log_rosette_cm"] if adjust_for_rosette else []),
        "betas":beta,
        "HC3_log1p_contrasts":results,
        "claim":"DESCRIPTIVE_COMPLETE_CASE_COVARIATE_SENSITIVITY_NOT_CAUSAL_ADJUSTED_MATERNAL_SEED_COUNTS",
    }


def audit_head_semantics(mothers,heads):
    m={r["year"]+"|"+r["plt"]:r for r in mothers}
    head_matching=[x for x in heads if
                   x["year"]+"|"+x["plt"] in m]
    mismatches=[x for x in head_matching if
                x["dam"]!=m[x["year"]+"|"+x["plt"]]["dam"] or
                x["ins"]!=m[x["year"]+"|"+x["plt"]]["ins"]]
    zero_flowering_positive=[r for r in mothers if
            value(r,"plt.gsd")>0 and value(r,"no.flw")==0]
    zero_seedhead_positive=[r for r in mothers if
            value(r,"plt.gsd")>0 and value(r,"n.fl.sds")==0]
    def apical_only(r):
        return value(r,"o.hds")==0 and math.isclose(
            value(r,"ap.hd"),value(r,"plt.gsd"),abs_tol=1e-9)
    require(len(zero_flowering_positive)==10 and
            len(zero_seedhead_positive)==17 and
            all(apical_only(r) for r in
                zero_flowering_positive+zero_seedhead_positive),
            "APICAL_ONLY_HEAD_COUNT_SEMANTIC_SOURCE_CHANGED")
    require(len(head_matching)==309 and
            len(mismatches)==12 and
            sorted({z["year"]+"|"+z["plt"] for z in mismatches})==
            ["2007|323","2007|375"],
            "SOURCE_HEAD_PARENT_TREATMENT_JOIN_CHANGED")
    return {
       "n_source_heads":len(heads),
       "n_head_rows_matching_original_mothers":len(head_matching),
       "n_head_rows_with_treatment_disagreement":len(mismatches),
       "n_positive_mothers_no_reported_flowered_heads":
           len(zero_flowering_positive),
       "n_positive_mothers_zero_reported_seed_bearing_heads":
           len(zero_seedhead_positive),
       "zero_reported_flowered_but_positive_seed_mother_ids":[
           r["year"]+"|"+r["plt"] for r in zero_flowering_positive],
       "zero_reported_seed_head_but_positive_seed_mother_ids":[
           r["year"]+"|"+r["plt"] for r in zero_seedhead_positive],
       "all_such_positive_seed_cases_are_apical_only":True,
       "no_flw_and_n_fl_sds_count_full_flowering_capitula_confirmed":False,
       "these_rows_are_proven_data_errors":False,
       "documented_variable_dictionary_may_use_nonapical_stage_count":True,
       "original_head_treatment_discrepancies_are_not_autocorrected":True,
    }


def source_groups(mothers):
    require(len(mothers)==106,"ORIGINAL_106_MOTHERS_CHANGED")
    ids=[r["year"]+"|"+r["plt"] for r in mothers]
    require(len(set(ids))==106,"DUPLICATE_SOURCE_MOTHER_ID")
    groups={y:{a:[r for r in mothers
                   if r["year"]==y and r["trt"]==a] for a in ARMS}
            for y in YEARS}
    require({y:{a:len(groups[y][a]) for a in ARMS}
             for y in YEARS}=={
       "2007":{"ndni":9,"dni":9,"ndi":10,"di":11},
       "2008":{"ndni":15,"dni":12,"ndi":21,"di":19}},
       "USDA_ORIGINAL_TWO_YEAR_ARM_COUNTS_CHANGED")
    for r in mothers:
        require(r["trt"]==r["dam"]+r["ins"]
                and r["trt"] in ARMS,"SOURCE_TREATMENT_NOT_VERIFIED")
    return groups


def analyze(mothers,heads,dictionary):
    require("no.flw" in dictionary and "n.fl.sds" in dictionary
            and "per plant" in dictionary,"MISSING_USDA_AUTHOR_DICTIONARY")
    groups=source_groups(mothers)
    audit=audit_head_semantics(mothers,heads)
    group={y:{a:group_stats(groups[y][a]) for a in ARMS} for y in YEARS}
    fits={}
    for y in YEARS:
        yr=[r for a in ARMS for r in groups[y][a]]
        keep=[r for r in yr if value(r,"r.diam") is not None
              and value(r,"r.diam")>0]
        fits[y]={
          "all_mothers_no_size_adjustment":fit_log1p_hc3(yr,False),
          "complete_cases_no_size_adjustment":fit_log1p_hc3(keep,False),
          "complete_cases_with_rosette_covariate":fit_log1p_hc3(keep,True),
          "n_missing_or_invalid_rosette":len(yr)-len(keep),
        }
    return {
       "version":"aza3_cirsium_undulatum_original_two_part_hurdle_exact_size_v1",
       "source":{
          "dataset_doi":DOI,"paper_doi":PAPER_DOI,
          "source_git_blobs":BLOBS,
          "original_mother_n":len(mothers),
          "original_head_n":len(heads),
          "real_original_sources_not_simulated":True,
       },
       "maternal_year_treatment_hurdle":group,
       "maternal_year_probability_exact_tests":exact_seed_presence_tests(groups),
       "maternal_group_mean_identity_factorization":hurdle_identity(groups),
       "rosette_complete_case_log1p_HC3_sensitivity":fits,
       "source_head_count_semantics":audit,
       "scientific_limits":{
          "reproductive_success_means_at_least_one_author_good_seed_not_germination":True,
          "mean_output_prob_times_positive_conditional_mean_is_identity_not_causal_mediation":True,
          "fisher_exact_tests_use_mother_not_head_independent_n":True,
          "year_group_treatment_populations_are_not_matched_repeated_mothers":True,
          "normal_HC3_interval_does_not_prove_randomization_causal_effect":True,
          "rosette_measurement_timing_is_not_verified_pre_treatment":True,
          "eleven_mothers_in_2008_have_missing_rosette":True,
          "head_count_zero_with_apical_seeds_may_reflect_variable_ontology":True,
          "insecticide_is_not_selective_to_flowerhead_seed_predators_alone":True,
          "independent_C_canescens_original_134_mother_rows_unavailable":True,
          "genetic_module_evolution_and_trait_reexpression_unmeasured":True,
          "new_field_plants_or_heads_collected":0,
       },
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mothers",type=Path,required=True)
    ap.add_argument("--heads",type=Path,required=True)
    ap.add_argument("--dictionary",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    mothers=read_locked(args.mothers,"mother_plants.csv",
                        ("plt","year","trt","dam","ins","plt.gsd",
                         "o.hds","ap.hd","r.diam","no.flw","n.fl.sds"))
    heads=read_locked(args.heads,"flower_heads.csv",
                      ("plt","year","dam","ins","hd","hd.gsd"))
    d=read_locked(args.dictionary,"source_dictionary.txt")
    result=analyze(mothers,heads,d)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                                     allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
       "source":result["source"],
       "year2008_hurdle":result["maternal_year_treatment_hurdle"]["2008"],
       "year2008_exact":result["maternal_year_probability_exact_tests"]["2008"],
       "year2008_model_adjusted":result["rosette_complete_case_log1p_HC3_sensitivity"]["2008"]["complete_cases_with_rosette_covariate"]["HC3_log1p_contrasts"],
       "head_semantic_issue":{
          "positive_but_flowered_zero":result["source_head_count_semantics"]["n_positive_mothers_no_reported_flowered_heads"],
          "positive_but_seed_head_zero":result["source_head_count_semantics"]["n_positive_mothers_zero_reported_seed_bearing_heads"]
       }
    },indent=2))


if __name__=="__main__":
    main()
