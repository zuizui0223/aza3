#!/usr/bin/env python3
"""Cirsium undulatum: reproductive *route* co-occurrence vs independent marginals.

Original USDA CC0 maternal data, West & Louda (2021) DOI 10.1002/ajb2.1744,
data DOI 10.15482/USDA.ADC/1522649. Source from merged aza3 PR #73,
original apical/lateral classification from merged PR #76, injury stages PR #77.

Actual maternal plant ID year|plt is the observational unit. Exclude ONLY
2007|327 (0 apical + 1 lateral but maternal total 20) from JOINT states;
it remains in the published 106-mother marginal outcome analyses.

Null comparison: within each of the eight fixed year x original factorial
treatment strata, apical and lateral seed-success indicators are independent,
conditional on their two observed marginal positive counts. Under the
conditional hypergeometric law, exact convolution gives the number of
mother plants with BOTH active routes (and, at fixed margins, with NEITHER).
This is not a randomized intervention or a phylogenetically independent
evolutionary comparison. A pooled 2x2 test mixes very different treatments.

Positive 'good seed' = author-scored seed, NOT proven viability/germination.
Organ allocation, size and pollinators may jointly determine these outcomes;
negative residual within-stratum co-occurrence is NOT proof of adaptive
risk spreading, resource trade-offs or genetic independent LEGO modules.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

DATA_DOI="10.15482/USDA.ADC/1522649"
PAPER_DOI="10.1002/ajb2.1744"
SOURCE_GIT_BLOB="458800272c3d0a03c7653702027a8f8578d9e88c"
YEARS=("2007","2008")
ARMS=("ndni","dni","ndi","di")
EXPECTED_STATES={
  "2007":{
    "ndni":(0,2,2,5),"dni":(2,0,6,0),
    "ndi":(0,3,2,5),"di":(0,1,6,4)
  },
  "2008":{
    "ndni":(12,1,1,1),"dni":(5,0,7,0),
    "ndi":(0,7,1,13),"di":(0,4,8,7)
  }
}
STATE_ORDER=("neither","apical_only","lateral_only","both")
INCONSISTENT_ORIGINAL_ID="2007|327"


def require(ok,code):
    if not ok:raise ValueError(code)


def git_blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()


def load_source(path):
    data=Path(path).read_bytes()
    require(git_blob(data)==SOURCE_GIT_BLOB,"ORIGINAL_USDA_MOTHER_BLOB_CHANGED")
    rows=list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))
    require(len(rows)==106,"ORIGINAL_MOTHER_TOTAL_CHANGED")
    ids=[r["year"]+"|"+r["plt"] for r in rows]
    require(len(ids)==len(set(ids)),"SOURCE_MOTHER_DUPLICATION")
    for r in rows:
        require(r["year"] in YEARS and r["trt"] in ARMS and
                r["trt"]==r["dam"]+r["ins"],
                "ORIGINAL_FACTORIAL_TREATMENT_CHANGED")
    return rows


def classify_original(rows):
    result={}
    excluded=[]
    for r in rows:
        name=r["year"]+"|"+r["plt"]
        try:
            ap=float(r["ap.hd"])
            side=float(r["o.hds"])
            total=float(r["plt.gsd"])
        except (KeyError,TypeError,ValueError):
            raise ValueError("MOTHER_ROUTE_NUMERIC_SOURCE_INVALID") from None
        require(all(math.isfinite(v) and v>=0 for v in (ap,side,total)),
                "MOTHER_ROUTE_NUMERIC_SOURCE_INVALID")
        if not math.isclose(ap+side,total,abs_tol=1e-9,rel_tol=0):
            excluded.append(name)
            continue
        state=("both" if ap>0 and side>0 else
               "apical_only" if ap>0 else
               "lateral_only" if side>0 else "neither")
        result[name]={
            "year":r["year"],"arm":r["trt"],"state":state,
            "apical":int(ap>0),"lateral":int(side>0),
        }
    require(excluded==[INCONSISTENT_ORIGINAL_ID] and len(result)==105,
            "UNRESOLVED_SOURCE_MOTHER_COMPONENT_DISCREPANCY")
    return result,excluded


def distribution(n,a,b):
    """P(K=k | n, A-positive=a, L-positive=b), K intersection.

    K is hypergeometric: C(a,k) C(n-a,b-k) / C(n,b).
    This conditions on source-observed route margins, not natural history
    family allocation randomization or group size equalities.
    """
    require(all(type(v)==int for v in (n,a,b)) and
            n>=1 and 0<=a<=n and 0<=b<=n,
            "INVALID_HYPERGEOMETRIC_MARGINAL")
    lower=max(0,a+b-n)
    upper=min(a,b)
    den=math.comb(n,b)
    probs={
        k:(math.comb(a,k)*math.comb(n-a,b-k))/den
        for k in range(lower,upper+1)
    }
    require(math.isclose(sum(probs.values()),1.,rel_tol=0,abs_tol=1e-11),
            "HYPERGEOMETRIC_PROBABILITY_MASS_DRIFT")
    return probs


def convolve(probs,dist):
    out={}
    for a,pa in probs.items():
        for b,pb in dist.items():
            out[a+b]=out.get(a+b,0.0)+pa*pb
    require(math.isclose(sum(out.values()),1.,rel_tol=0,abs_tol=1e-9),
            "NULL_CONVOLUTION_MASS_DRIFT")
    return out


def summarized_stratum(rows,year,arm):
    rs=[r for r in rows.values() if r["year"]==year and r["arm"]==arm]
    require(rs,"NO_ORIGINAL_MOTHERS_IN_STRATUM")
    n=len(rs)
    state_counts={k:sum(r["state"]==k for r in rs) for k in STATE_ORDER}
    expect=EXPECTED_STATES[year][arm]
    require(tuple(state_counts.values())==expect,
            "ORIGINAL_FOUR_STATE_MOTHER_CELLS_CHANGED")
    a=state_counts["apical_only"]+state_counts["both"]
    b=state_counts["lateral_only"]+state_counts["both"]
    both=state_counts["both"]
    none=state_counts["neither"]
    require(n-a-b+both==none,"MOTHER_ROUTE_IDENTITY_FAIL")
    h=distribution(n,a,b)
    return {
        "year":year,"arm":arm,"n_mothers":n,
        "route_counts":state_counts,
        "apical_positive":a,"lateral_positive":b,
        "both_positive_observed":both,"neither_positive_observed":none,
        "independent_expected_both_conditioning_on_source_margins":a*b/n,
        "independent_expected_neither_conditioning_on_source_margins":
            (n-a)*(n-b)/n,
        "overlap_sharp_support":[min(h),max(h)],
        "overlap_null_distribution":{str(k):v for k,v in h.items()},
        "exact_conditioning_on_original_margins":True,
        "is_a_randomized_head_module_independence_test":False,
    }


def exact_route_analysis(rows):
    strata=[summarized_stratum(rows,year,arm)
            for year in YEARS for arm in ARMS]
    dp={0:1.}
    for z in strata:
        d={int(k):v for k,v in z["overlap_null_distribution"].items()}
        dp=convolve(dp,d)
    n=sum(z["n_mothers"] for z in strata)
    ap=sum(z["apical_positive"] for z in strata)
    lat=sum(z["lateral_positive"] for z in strata)
    both=sum(z["both_positive_observed"] for z in strata)
    none=sum(z["neither_positive_observed"] for z in strata)
    expect=sum(z["independent_expected_both_conditioning_on_source_margins"]
               for z in strata)
    expect_none=sum(z["independent_expected_neither_conditioning_on_source_margins"]
                    for z in strata)
    require((n,ap,lat,both,none)==(105,53,68,35,19),
            "ORIGINAL_105_MOTHER_JOINT_COUNTS_CHANGED")
    lower=sum(v for k,v in dp.items() if k<=both)
    upper=sum(v for k,v in dp.items() if k>=both)
    primary=min(1.0,2*min(lower,upper))
    observed_mass=dp[both]
    probability_ordered=sum(v for v in dp.values()
                            if v<=observed_mass+1e-14)
    expected_pooled=ap*lat/n
    observed_or=(both*none)/((ap-both)*(lat-both))
    numerator=sum(z["route_counts"]["both"]*
                  z["route_counts"]["neither"]/z["n_mothers"]
                  for z in strata)
    denominator=sum(z["route_counts"]["apical_only"]*
                    z["route_counts"]["lateral_only"]/z["n_mothers"]
                    for z in strata)
    require(denominator>0,"STRATIFIED_ODDS_REFERENCE_UNDEFINED")
    years={}
    for year in YEARS:
        r=[s for s in strata if s["year"]==year]
        d={0:1.}
        for x in r:d=convolve(d,{int(k):v for k,v in x["overlap_null_distribution"].items()})
        c=sum(x["both_positive_observed"] for x in r)
        lower_y=sum(v for k,v in d.items() if k<=c)
        upper_y=sum(v for k,v in d.items() if k>=c)
        years[year]={
            "n_mothers":sum(x["n_mothers"] for x in r),
            "both_observed":c,
            "both_expected_under_conditional_independence":sum(
                x["independent_expected_both_conditioning_on_source_margins"]
                for x in r),
            "exact_two_inclusive_tails_doubled_min":
                min(1.0,2*min(lower_y,upper_y)),
            "null_overlap_distribution":{str(k):v for k,v in sorted(d.items())},
        }
    return {
        "version":"cirsium_undulatum_two_reproductive_routes_exact_null_v1",
        "source":{
            "paper_doi":PAPER_DOI,"dataset_doi":DATA_DOI,
            "mother_csv_git_blob":SOURCE_GIT_BLOB,
            "original_mother_n":106,"joint_state_coherent_n":105,
            "excluded_mother_ids":[INCONSISTENT_ORIGINAL_ID],
        },
        "eight_original_year_by_treatment_strata":strata,
        "pooled_ignoring_year_and_treatment":{
            "n_mothers":n,"apical_positive":ap,
            "lateral_positive":lat,"both_positive":both,
            "neither_positive":none,
            "apical_only":ap-both,"lateral_only":lat-both,
            "observed_minus_independence_expected_both":both-expected_pooled,
            "independent_expected_both":expected_pooled,
            "pooled_odds_ratio_both_vs_mutual_failure":observed_or,
            "positive_residual_relative_to_pooled_null":both>expected_pooled,
            "pooled_expectation_not_valid_under_stratified_treatment":True,
        },
        "stratified_conditional_marginal_independence":{
            "observed_both":both,"independent_expected_both":expect,
            "observed_minus_expected_both":both-expect,
            "observed_neither":none,"independent_expected_neither":expect_none,
            "observed_minus_expected_neither":none-expect_none,
            "prob_neither_observed":none/n,
            "prob_neither_if_each_stratum_marginals_independent":expect_none/n,
            "observed_conditional_MH_odds_ratio_numerator_over_denominator":
                numerator/denominator,
            "exact_conditional_lower_tail_inclusive":lower,
            "exact_conditional_upper_tail_inclusive":upper,
            "exact_two_inclusive_tails_doubled_min":primary,
            "probability_ordered_exact_two_sided_sensitivity":probability_ordered,
            "null_total_overlap_distribution":{
                str(k):v for k,v in sorted(dp.items())},
            "independence_is_statistical_conditional_reference_not_measured":
                True,
            "evidence_supports_nonzero_negative_genetic_covariance":False,
        },
        "year_conditioned_sensitivity":years,
        "scientific_claim_limits":{
            "source_seed_routes_are_author_scored_good_seeds_not_confirmed_seedlings":
                True,
            "conditioned_on_year_and_herbivory_intervention_margins":True,
            "two_observed_routes_are_not_lost_ancestral_genes":True,
            "route_correlation_can_reflect_maternal_quality_resource_allocation_or_sampling":
                True,
            "no_independent_oviposition_or_effective_pollination_at_same_head":True,
            "head_orientation_true_spine_or_phyllary_shape_not_measured":True,
            "post_treatment_selection_cannot_identify_optimal_resource_allocation":True,
            "negative_covariance_direction_is_not_significant_at_doubled_tail_0p05":
                primary>=.05,
            "head_pair_null_does_not_establish_selected_bet_hedging":True,
            "does_not_infer_species_or_lifehistory_wide_conclusion":True,
            "newly_observed_individual_plants_or_heads":0,
        }
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mothers",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    q=p.parse_args()
    rows=load_source(q.mothers)
    joint,bad=classify_original(rows)
    result=exact_route_analysis(joint)
    q.out.parent.mkdir(parents=True,exist_ok=True)
    q.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                                allow_nan=False)+"\n",encoding="utf-8")
    x=result["stratified_conditional_marginal_independence"]
    print(json.dumps({
        "actual_coherent_mother_n":len(joint),
        "incoherent_original_row":bad,
        "pooled_positive_overlap_contrast":
             result["pooled_ignoring_year_and_treatment"][
                  "observed_minus_independence_expected_both"],
        "stratified_negative_overlap_contrast":
             x["observed_minus_expected_both"],
        "stratified_doubled_min_tail_p":
             x["exact_two_inclusive_tails_doubled_min"],
        "claim":"DESCRIPTIVE_SOURCE_STRATIFICATION_SIGN_REVERSAL_NOT_ADAPTIVE_RISK_INSURANCE",
    },indent=2))


if __name__=="__main__":
    main()
