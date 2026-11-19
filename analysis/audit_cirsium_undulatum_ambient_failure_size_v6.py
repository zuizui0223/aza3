#!/usr/bin/env python3
"""Measured maternal seed-zero vs size and scored-head abundance, real USDA study.

Only 2008 ambient-insect Cirsium undulatum maternal plants (27) are used.
West & Louda (2021) DOI 10.1002/ajb2.1744; original USDA CC0
10.15482/USDA.ADC/1522649. The original 106 mothers and 323 head
records are source-Git-SHA checked by the parent established audit.

This is a POST-HOC, post-treatment COVARIATE association screen, not a
selection or causal mediation model. Conditional label permutations hold
observed seed-success counts fixed within apical damage levels. They are a
conditional descriptive null, not randomized outcome assignments or a
realistic genetic/environmental process. Source-ID prefixes are NOT
verified independent field sites. Multiple testing is BH over exactly six
exploratory variables. Zero/positive 'good seeds' is author-scored, NOT
germination, descendant fitness or ANTO native evolutionary regain.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
from functools import reduce
from itertools import combinations, product
import json
import math
from pathlib import Path

from audit_cirsium_undulatum_scored_head_failure_v4 import (
    source_locked_tables, head_mother_crosswalk, good_seed_int
)
from audit_cirsium_undulatum_lateral_head_positions_v5 import source_prefix_group

YEAR="2008"
ARMS=("ndni","dni")
PREFIXES=("ES","SD","WG","WR")
SCALE=10**6
FEATURES={
    "plant_height_cm":"plt.ht",
    "basal_rosette_diameter_cm":"r.diam",
    "flower_bud_count":"n.hds",
    "branch_count":"n.br",
    "mean_head_damage_fraction":"av.dam",
    "scored_lateral_head_count":"HEAD_SOURCE"
}


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def exact_decimal_units(raw):
    """Use exact fixed six-decimal units, never float bin comparisons."""
    if raw is None or str(raw).strip().upper() in ("NA","","N/A","."):
        return None
    try:
        x=Decimal(str(raw).strip())
    except (ValueError,InvalidOperation):
        raise ValueError("INVALID_VALUE_NONNUMERIC") from None
    require(x.is_finite() and x>=0, "INVALID_VALUE_NEGATIVE_OR_NONFINITE")
    scaled=x*SCALE
    require(scaled==scaled.to_integral_value(),
            "INVALID_VALUE_EXCEEDS_SIX_DECIMAL_PLACES")
    return int(scaled)


def exact_stratified_binary(rows,strata_keys):
    """Finite exact two-sided centered positive-sum with fixed stratum K.

    Let stratum g have n measured plants, k positive seed mothers,
    x_ig measured covariate. Statistic:
          T = sum_g [ sum_(positive i) x_ig - k_g mean_g(x) ]
    (all in original measurement units).
    Exhaustively enumerate all combinations of k positive labels within
    each eligible nondegenerate stratum, conditional on that stratum's
    measured plants and observed k. A degenerate all-zero/all-one stratum
    contributes no variation and is not used for evidence.

    This tests only a conditional associational null. There is no
    identification of temporal causal exposure or random genotype assignment.
    """
    require(isinstance(strata_keys,(list,tuple)) and strata_keys,
            "INVALID_STRATUM_KEYS")
    require(isinstance(rows,(tuple,list)) and rows,
            "INVALID_STRATUM_EMPTY_ROWS")
    valid=[]
    missing=0
    for row in rows:
        require(type(row.get("positive")) is int and row["positive"] in (0,1),
                "INVALID_POSITIVE_BINARY_OUTCOME")
        require(all(key in row and str(row[key]) for key in strata_keys),
                "INVALID_STRATUM_MISSING_KEY")
        num=exact_decimal_units(row.get("value"))
        if num is None:
            missing+=1
            continue
        valid.append({
            "positive":row["positive"],
            "values_scaled":num,
            "group":tuple(str(row[k]) for k in strata_keys),
        })
    require(valid,"INVALID_STRATUM_NO_MEASURED_PLANTS")
    groups={}
    for row in valid:
        groups.setdefault(row["group"],[]).append(row)
    active=[]
    uninformative=0
    for name in sorted(groups):
        g=groups[name]
        n=len(g)
        k=sum(r["positive"] for r in g)
        if k==0 or k==n:
            uninformative+=n
            continue
        vals=[r["values_scaled"] for r in g]
        active.append({
            "keys":list(name), "n":n, "positive":k,
            "values":vals, "total":sum(vals),
            "observed_positive_sum":sum(r["values_scaled"] for r in g
                                        if r["positive"]),
        })
    require(active, "INVALID_STRATUM_NO_WITHIN_GROUP_OUTCOME_VARIATION")

    denom=math.lcm(*(group["n"] for group in active))
    obs=0
    distributions=[]
    for g in active:
        factor=denom//g["n"]
        shift=g["positive"]*g["total"]*factor
        obs+=denom*g["observed_positive_sum"]-shift
        distributions.append([
            denom*sum(chosen)-shift
            for chosen in combinations(g["values"],g["positive"])
        ])
    nways=math.prod(len(d) for d in distributions)
    require(nways>0 and nways<=10**7,
            "INVALID_STRATUM_TOO_MANY_EXACT_ASSIGNMENTS")
    tail=0
    # Event abs(T_perm) >= abs(T_observed). Exact INTEGER arithmetic.
    # Python loop is bounded; no naive sampling or asymptotic approximation.
    target=abs(obs)
    for realization in product(*distributions):
        if abs(sum(realization))>=target:
            tail+=1
    return {
        "n_measured_mothers":len(valid),
        "n_missing":missing,
        "n_informative_strata":len(active),
        "uninformative_stratum_mothers":uninformative,
        "total_assignments":nways,
        "tail_assignments":tail,
        "p_two_sided":tail/nways,
        "observed_centered_positive_sum":obs/(denom*SCALE),
        "fixed_positive_counts_by_informative_stratum":[
            {"stratum":g["keys"],"n":g["n"],"n_positive":g["positive"],
             "n_combinations":math.comb(g["n"],g["positive"])}
            for g in active
        ],
        "outcome_not_randomly_permuted_by_experiment":True,
    }


def bh_qvalues_six(pvalues):
    require(isinstance(pvalues,dict) and len(pvalues)==6,
            "SIX_SCREENED_VARIABLES_REQUIRED")
    require(all(type(p) in (float,int) and math.isfinite(p)
                and 0<=p<=1 for p in pvalues.values()),
            "INVALID_MULTIPLICITY_P")
    ranked=sorted(pvalues,key=lambda key:pvalues[key])
    n=len(ranked)
    out={}
    q=1.
    for i in range(n-1,-1,-1):
        key=ranked[i]
        q=min(q,n*pvalues[key]/(i+1))
        out[key]=q
    return out


def source_ambient_mothers(mothers=None,heads=None,strict=True):
    if mothers is None and heads is None:
        mothers,heads=source_locked_tables()
    require(mothers is not None and heads is not None,
            "NEED_BOTH_REAL_SOURCE_TABLES")
    verified,source_2007=head_mother_crosswalk(
        mothers,heads,strict2008=strict)
    ambient=[]
    for row in verified:
        mom=row["mother"]
        if mom["ins"]!="ni":
            continue
        arm=mom["trt"]
        require(arm in ARMS,"ORIGINAL_TREATMENT_AMBIENT_NOT_VERIFIED")
        entry={
            "id":str(mom["plt"]),
            "prefix":source_prefix_group(mom["plt"]),
            "arm":arm,
            "positive":int(good_seed_int(mom,"plt.gsd")>0),
            "values":{
                var:str(len(row["lateral"])) if src=="HEAD_SOURCE"
                    else str(mom[src])
                for var,src in FEATURES.items()
            },
        }
        ambient.append(entry)
    ambient.sort(key=lambda x:x["id"])
    require(len(ambient)==27 and
            sum(x["positive"] for x in ambient)==10,
            "ORIGINAL_2008_AMBIENT_27_MOTHERS_CHANGED")
    require({a:sum(x["arm"]==a for x in ambient) for a in ARMS}
            =={"ndni":15,"dni":12},"ORIGINAL_2008_AMBIENT_ARMS_CHANGED")
    require(set(x["prefix"] for x in ambient)==set(PREFIXES),
            "ORIGINAL_SOURCE_ID_PREFIX_SET_CHANGED")
    return ambient


def _rows_for(ambient,var):
    return [
        {"value":item["values"][var],"positive":item["positive"],
         "arm":item["arm"],"prefix":item["prefix"]}
        for item in ambient
    ]


def _variable_screen(ambient,var):
    original=_rows_for(ambient,var)
    measured=[r for r in original if exact_decimal_units(r["value"]) is not None]
    success=[float(r["value"]) for r in measured if r["positive"]==1]
    failure=[float(r["value"]) for r in measured if r["positive"]==0]
    require(success and failure,"ORIGINAL_VARIABLE_ONLY_ONE_OUTCOME")
    p=exact_stratified_binary(original,["arm"])
    require(p["n_measured_mothers"]==len(measured),
            "COVARIATE_MISSINGNESS_DRIFT")
    return {
        "n_measured_mothers":len(measured),
        "seed_positive":len(success),
        "seed_zero":len(failure),
        "n_missing":len(original)-len(measured),
        "mean_in_zero_seed_mothers":sum(failure)/len(failure),
        "mean_in_seed_positive_mothers":sum(success)/len(success),
        "mean_difference_positive_minus_zero":
            sum(success)/len(success)-sum(failure)/len(failure),
        "exact_stratified_permutation":p,
        "q_BH_six_exploratory_covariates":None,
    }


def audit_ambient_size(mothers=None, heads=None, strict=True):
    ambient=source_ambient_mothers(mothers,heads,strict=strict)
    screened={var:_variable_screen(ambient,var) for var in FEATURES}
    qs=bh_qvalues_six({
        var:screened[var]["exact_stratified_permutation"]["p_two_sided"]
        for var in FEATURES
    })
    for var in FEATURES:
        screened[var]["q_BH_six_exploratory_covariates"]=qs[var]

    height_rows=_rows_for(ambient,"plant_height_cm")
    leave_prefix={}
    for prefix in PREFIXES:
        cut=[r for r in height_rows if r["prefix"]!=prefix]
        r=exact_stratified_binary(cut,["arm"])
        leave_prefix[prefix]={
            "removed_prefix":prefix,
            "n_remaining_mothers":len(cut),
            "p_two_sided":r["p_two_sided"],
            "n_assignments":r["total_assignments"],
        }
    tallest=max(
        (v for v in ambient if v["positive"]==1),
        key=lambda v:exact_decimal_units(v["values"]["plant_height_cm"])
    )
    without=[r for v,r in zip(ambient,height_rows) if v["id"]!=tallest["id"]]
    result_no_tallest=exact_stratified_binary(without,["arm"])
    within_prefix=exact_stratified_binary(
        height_rows,["arm","prefix"])
    output={
        "schema":"aza3_Cirsium_undulatum_2008_ambient_mother_source_size_screen_v6",
        "source":{
            "dataset_doi":"10.15482/USDA.ADC/1522649",
            "article_doi":"10.1002/ajb2.1744",
            "original_mothers_all_years":106,
            "original_heads_all_years":323,
            "complete_2008_mothers":67,
            "ambient_insect_2008_mothers":len(ambient),
            "ambient_apical_intact":sum(r["arm"]=="ndni" for r in ambient),
            "ambient_apical_damaged":sum(r["arm"]=="dni" for r in ambient),
            "any_seed_positive":sum(r["positive"] for r in ambient),
            "both_seed_routes_failed":sum(not r["positive"] for r in ambient),
            "original_source_id_prefixes":list(PREFIXES),
            "no_flw_not_in_six_covariates":
                "ORIGINAL_ZERO_FLOWERED_BUT_POSITIVE_APICAL_SEED_CONFLICT",
            "not_seven_posthoc_favorable_covariates":True,
            "raw_mother_and_head_Git_SHA_checked":True,
            "same_2008_population_as_prior_aza3_source_audits":True,
        },
        "six_exploratory_covariates":screened,
        "height_sensitivity":{
            "leave_out_source_prefix":leave_prefix,
            "leave_out_highest_successful_height":{
                "source_plant_identifier":tallest["id"],
                "height_cm":float(tallest["values"]["plant_height_cm"]),
                "n_mothers":len(without),
                "p_two_sided":result_no_tallest["p_two_sided"],
                "n_assignments":result_no_tallest["total_assignments"],
            },
            "source_prefix_and_apical_treatment_conditional":
                within_prefix,
        },
        "interpretive_ceiling":{
            "source_prefixes_are_verified_field_sites":False,
            "flowering_fields_are_complete_verified_anthesis":False,
            "post_treatment_plant_height_is_pre_exposure_causal_covariate":False,
            "rosette_size_measurements_are_complete_and_pretreatment":False,
            "number_of_heads_is_non_endogenous_unconfounded_treatment":False,
            "observed_damage_score_identifies_pollinator_vs_seed_predator":False,
            "height_p_012_is_confirmatory_after_six_BH":False,
            "joint_failure_due_to_one_arthropod_established":False,
            "head_number_independent_of_maternal_condition":False,
            "SLK_candidate_architecture_relative_fitness_identified":False,
            "ANTO_ancestral_genetic_re_expression_or_heritability_estimated":False,
            "parent_good_seed_positive_means_seedling_established":False,
            "independent_new_plant_experiment_done_by_aza3":False,
        },
        "scientific_conclusion":
            "AMBIENT_MATERNAL_SEED_FAILURE_NOT_DETERMINISTICALLY_EXPLAINED_BY_SIZE_OR_SCORED_DAMAGE_NO_CAUSE_IDENTIFIED",
    }
    if strict:
        h=screened["plant_height_cm"]
        require(h["n_measured_mothers"]==27 and
                screened["basal_rosette_diameter_cm"]["n_measured_mothers"]==24,
                "SOURCE_SIZE_MISSINGNESS_CHANGED")
        require(math.isclose(h["exact_stratified_permutation"]["p_two_sided"],
                             .012218337218337219,rel_tol=0,abs_tol=1e-14),
                "SOURCE_HEIGHT_EXACT_PERMUTATION_CHANGED")
        require(math.isclose(h["q_BH_six_exploratory_covariates"],
                             .07331002331002331,rel_tol=0,abs_tol=1e-14),
                "SOURCE_HEIGHT_MULTIPLICITY_CHANGED")
        require(math.isclose(leave_prefix["WR"]["p_two_sided"],
                             .2945165945165945,rel_tol=0,abs_tol=1e-14),
                "SOURCE_WR_LEAVE_PREFIX_CHANGED")
        require(math.isclose(result_no_tallest["p_two_sided"],
                             .07473082473082474,rel_tol=0,abs_tol=1e-14),
                "SOURCE_TALLEST_LEAVE_ONE_CHANGED")
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    r=audit_ambient_size(strict=True)
    payload=json.dumps(r,sort_keys=True,indent=2,allow_nan=False)+"\n"
    if args.output is None:
        print(payload)
    else:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(payload,encoding="utf-8")


if __name__=="__main__":
    main()
