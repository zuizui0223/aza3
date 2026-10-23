#!/usr/bin/env python3
"""Real maternal-plant apical/lateral reproductive routes in Cirsium undulatum.

West & Louda 2021 DOI 10.1002/ajb2.1744, genuine USDA-ADC source
DOI 10.15482/USDA.ADC/1522649, CC0 (source locked on aza3 main).
This is NOT an invented split of unknown heads: ap.hd (apical seed)
and o.hds (sum nonapical good seeds) coexist in each original maternal row.

A plant's apical vs lateral seed-production routes are NOT ancestral trait
regain, developmental gene module retention or independent anatomical
fitness effects. Apical "damage" is not the same as complete organ removal.
Source 2007 plant 327 has 0 apical + 1 nonapical vs maternal total 20:
exclude ONE row from four-state route claims without correcting it.
The original 106 plants remain available in parent hurdle analysis.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import json
import math

from audit_cirsium_undulatum_maternal_hurdle_size_v1 import (
    YEARS, ARMS, CONTRAST, read_locked, source_groups,
    value, fisher_two_sided, bh_qvalues,
)

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data/external/usda_wavyleaf_2007_2008"
ROUTES=("no_good_seed","apical_only","lateral_only","both_sites")


def require(ok,msg):
    if not ok:
        raise ValueError(msg)


def classify_route(row):
    """Return anatomical seed outcome or None for known inconsistent totals."""
    try:
        apical=value(row,"ap.hd")
        lateral=value(row,"o.hds")
        whole=value(row,"plt.gsd")
    except (ValueError,TypeError):
        raise ValueError("MOTHER_ROUTE_INVALID_SEED_COMPONENT") from None
    require(all(x is not None and math.isfinite(x) and x>=0
                for x in (apical,lateral,whole)),
            "MOTHER_ROUTE_NEGATIVE_OR_MISSING_SEED_COMPONENT")
    # The source's one inconsistent row is a missing JOINT observation.
    # It remains in the baseline 106-mother marginal seed dataset.
    if not math.isclose(apical+lateral,whole,rel_tol=0,abs_tol=1e-8):
        return None
    if apical and lateral:
        return "both_sites"
    if apical:
        return "apical_only"
    if lateral:
        return "lateral_only"
    return "no_good_seed"


def summarize_year_arm(rows):
    """One independently verified maternal plant = one multinomial unit."""
    route_counts={r:0 for r in ROUTES}
    bad_ids=[]
    all_apical=all_lateral=all_whole=0
    for r in rows:
        state=classify_route(r)
        if state is None:
            bad_ids.append(f"{r['year']}|{r['plt']}")
            continue
        route_counts[state]+=1
        all_apical+=value(r,"ap.hd")
        all_lateral+=value(r,"o.hds")
        all_whole+=value(r,"plt.gsd")
    n=sum(route_counts.values())
    require(n>0,"NO_CONSISTENT_MOTHER_ROUTE_CASES")
    n_ap=route_counts["apical_only"]+route_counts["both_sites"]
    n_lat=route_counts["lateral_only"]+route_counts["both_sites"]
    nseed=n-route_counts["no_good_seed"]
    require(math.isclose(all_apical+all_lateral,all_whole,
                         abs_tol=1e-8),"MOTHER_ROUTE_SUM_NO_LONGER_CONSISTENT")
    return {
        "n":n,
        "original_n_mothers":len(rows),
        "excluded_incoherent_mothers":len(bad_ids),
        "states":route_counts,
        "apical_positive_mothers":n_ap,
        "lateral_positive_mothers":n_lat,
        "seed_positive_mothers":nseed,
        "apical_success_fraction":n_ap/n,
        "lateral_success_fraction":n_lat/n,
        "p_lateral_only":route_counts["lateral_only"]/n,
        "p_no_good_seed":route_counts["no_good_seed"]/n,
        "apical_mean_good_seed_count":all_apical/n,
        "lateral_mean_good_seed_count":all_lateral/n,
        "maternal_mean_good_seed_count":all_whole/n,
        "among_apical_negative_lateral_success_fraction":
            route_counts["lateral_only"]/
            (route_counts["lateral_only"]+route_counts["no_good_seed"])
            if route_counts["lateral_only"]+route_counts["no_good_seed"] else None,
        "pollen_deposition_not_measured":True,
        "apical_and_lateral_genotype_not_two_independent_mothers":True,
    }


def fixed_fisher_eight_contrasts(year_states):
    """Fisher-BH on EIGHT descriptive parent lateral-presence contrasts.

    All four named contrasts per study year are included; a nominal p<0.05
    is NOT automatically evidential after multiple contrasts. This is a
    post-hoc source audit, not a preregistered trial significance claim.
    """
    result={}
    pvals={}
    for year in YEARS:
        for name,(treat,control) in CONTRAST.items():
            a=year_states[year][treat]
            b=year_states[year][control]
            tab=[
                [a["lateral_positive_mothers"],a["n"]-a["lateral_positive_mothers"]],
                [b["lateral_positive_mothers"],b["n"]-b["lateral_positive_mothers"]]
            ]
            exact=fisher_two_sided(*tab[0],*tab[1])
            key=f"{year}_{name}"
            pvals[key]=exact
            result[key]={
                "compared_arms":[treat,control],
                "lateral_positive_negative_counts":tab,
                "difference_lateral_success_probability":
                    a["lateral_success_fraction"]-b["lateral_success_fraction"],
                "p_exact":exact,
                "association_not_measured_causal_pollen_herbivory_decomposition":True
            }
    q=bh_qvalues(pvals)
    for key in result:
        result[key]["q_BH_eight"]=q[key]
    return result


def maternal_routes(mothers,strict=True):
    """Re-use original aza3 source bounds rather than importing from models."""
    source_groups(mothers)  # checks original 106, 2×2 assignments, year arms
    rows=source_groups(mothers)
    head=read_locked(BASE/"flower_heads.csv","flower_heads.csv",
                     ("plt","year","dam","ins","hd","hd.gsd"))
    dictionary=read_locked(BASE/"source_dictionary.txt",
                           "source_dictionary.txt")
    require("per plant" in dictionary and "o.hds" in dictionary,
            "ORIGINAL_USDA_DICTIONARY_NOT_AUTHENTICATED")
    require(len(head)==323,"ORIGINAL_USDA_HEAD_ROWS_CHANGED")
    observed={}
    invalid=[]
    for year in YEARS:
        observed[year]={}
        for arm in ARMS:
            records=rows[year][arm]
            for row in records:
                if classify_route(row) is None:
                    invalid.append(f"{row['year']}|{row['plt']}")
            observed[year][arm]=summarize_year_arm(records)
    valid=sum(observed[y][a]["n"] for y in YEARS for a in ARMS)
    require(valid+len(invalid)==len(mothers),
            "MOTHER_ROUTE_BAD_ROW_ACCOUNTING")
    if strict:
        require(valid==105 and invalid==["2007|327"],
                "ORIGINAL_APICAL_LATERAL_MASS_BALANCE_DRIFT")
        require(observed["2008"]["dni"]["states"]=={
            "no_good_seed":5,"apical_only":0,
            "lateral_only":7,"both_sites":0},
            "FIELD_2008_DAMAGE_AMBIENT_STATES_CHANGED")
        require(observed["2008"]["di"]["states"]=={
            "no_good_seed":0,"apical_only":4,
            "lateral_only":8,"both_sites":7},
            "FIELD_2008_DAMAGE_REDUCED_STATES_CHANGED")
    return {
        "version":"aza3_usda_Cirsium_undulatum_real_maternal_reproductive_routes_v2",
        "source":{
            "dataset_doi":"10.15482/USDA.ADC/1522649",
            "original_paper_doi":"10.1002/ajb2.1744",
            "published_experiment_ownership":"West and Louda (2021)",
            "original_mother_rows":len(mothers),
            "coherent_mother_rows":valid,
            "inconsistent_mother_identifiers":invalid,
            "original_heads_in_related_USDA_source":len(head),
            "source_verification":"existing_parent_USDA_git_sha_and_exact_row_count_contract",
            "missing_fraction_rejected_not_imputed":True,
        },
        "maternal_route_states":observed,
        "lateral_presence_fisher_tests":fixed_fisher_eight_contrasts(observed),
        "inference":{
            "observed_entity":"maternal_good_seed_anatomical_route",
            "biological_conclusion":
                "OBSERVED_REAL_SPATIAL_REPRODUCTIVE_FALLBACK_NOT_ADAPTIVE_SELECTION_OR_GENETIC_REGAIN",
            "original_author_damage_compensation_already_prior_art":True,
            "2008_ambient_apical_damaged_lateral_only_mothers":7,
            "2008_ambient_apical_damaged_zero_good_seed_mothers":5,
            "2008_insecticide_apical_damaged_apical_positive_mothers":11,
            "exclude_year2007_bad_sum_but_do_not_erase_from_parent_study":True,
        },
        "scientific_limits":{
            "apical_damage_means_apical_function_completely_removed":False,
            "same_mother_lateral_positive_means_extra_compensation_caused_by_damage":False,
            "new_treatment_effect_is_randomized_between_plant_genotypes_proven":False,
            "posthoc_fisher_nominal_p_is_robust_BH_significance":False,
            "insecticide_selectively_eliminated_seed_predator_without_pollinators":False,
            "maternal_seed_presence_means_established_seedling":False,
            "documented_heritable_anatomical_reexpression_or_gene_toolkit":False,
            "SLK_continuous_architecture_payoff_or_rare_invader_estimated":False,
            "ANTO_lost_ancestral_trait_regained_by_native_allele":False,
            "different_years_are_repeated_same_maternal_genotypes":False,
            "public_original_plant_data_collected_by_this_repository":False,
        }
    }


def main():
    root=ROOT
    parser=argparse.ArgumentParser()
    parser.add_argument("--mothers",type=Path,
                        default=BASE/"mother_plants.csv")
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    r=read_locked(args.mothers,"mother_plants.csv",
                  ("plt","year","trt","dam","ins","plt.gsd","ap.hd","o.hds"))
    result=maternal_routes(r,strict=True)
    payload=json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(payload,encoding="utf-8")
    else:
        print(payload)


if __name__=="__main__":
    main()
