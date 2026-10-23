#!/usr/bin/env python3
"""Observed, SAME-mother coverage of apical and lateral Cirsium reproduction.

USDA CC0 original 2007-08 Cirsium undulatum 106 mother rows,
West & Louda (2021), DOI 10.1002/ajb2.1744, USDA DOI
10.15482/USDA.ADC/1522649. Source-git-blob verification and invalid
plant 2007|327 exclusion are inherited from original aza3 v2.

P(A or L)-max(P(A),P(L)) is the observed extra UNION coverage beyond
the better observed route. It is an algebraic identity for success
events on the same individual, NEVER fitness gained by removing the
other head, selection for head division of labour, or recovered ancestry.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import json
import math

from audit_cirsium_undulatum_maternal_reproductive_routes_v2 import (
    maternal_routes, ROUTES, BASE
)
from audit_cirsium_undulatum_maternal_hurdle_size_v1 import (
    read_locked, YEARS, ARMS
)

ORIGINAL_DOI = "10.15482/USDA.ADC/1522649"
ORIGINAL_PAPER = "10.1002/ajb2.1744"


def require(test, reason):
    if not test:
        raise ValueError(reason)


def exact_discordance_test(apical_only, lateral_only):
    """Two-sided exact McNemar via Binomial(n discordant, p=.5).

    This asks whether within-mother route-only frequencies differ,
    not whether either route *caused* another or genetically reappeared.
    """
    a,b=apical_only,lateral_only
    require(type(a) is int and type(b) is int and a>=0 and b>=0,
            "INVALID_PAIRED_DISCORDANCE_COUNTS")
    n=a+b
    if not n:
        return 1.
    low=min(a,b)
    return min(1., 2*sum(math.comb(n,k) for k in range(low+1))/(2**n))


def bh_correction(ps):
    require(bool(ps) and all(
        type(v) in (float,int) and math.isfinite(v) and 0<=v<=1
        for v in ps.values()),"INVALID_MULTIPLICITY_P_VALUES")
    order=sorted(ps,key=ps.get)
    n=len(ps)
    q={}
    running=1.
    for i in range(n-1,-1,-1):
        name=order[i]
        running=min(running, n*ps[name]/(i+1))
        q[name]=running
    return q


def coverage_from_counts(counts):
    """The four mutually exclusive maternal outcomes for ONE arm/year."""
    require(isinstance(counts,dict) and set(counts)==set(ROUTES),
            "INVALID_FOUR_ANATOMICAL_ROUTES")
    require(all(type(v) is int and v>=0 for v in counts.values()),
            "INVALID_ROUTE_COUNT")
    n=sum(counts.values())
    require(n>0,"EMPTY_MATERNAL_ARM")
    zero=counts["no_good_seed"]
    ap=counts["apical_only"]
    lat=counts["lateral_only"]
    both=counts["both_sites"]
    pa=(ap+both)/n
    pl=(lat+both)/n
    union=(ap+lat+both)/n
    pn=both/n
    require(abs(pa+pl-pn-union)<1e-13,
            "MATERNAL_JOINT_PROBABILITY_IDENTITY")
    s=ap+lat
    no_ap=zero+lat
    no_lat=zero+ap
    return {
        "n_mothers":n,
        "four_anatomical_routes":dict(counts),
        "P_apical_good_seed":pa,
        "P_lateral_good_seed":pl,
        "P_both_apical_and_lateral":pn,
        "P_joint_failure_both_routes":zero/n,
        "P_any_maternal_good_seed":union,
        "observed_union_surplus_over_apical_route":lat/n,
        "observed_union_surplus_over_lateral_route":ap/n,
        "observed_union_surplus_over_better_single_route":
            min(ap,lat)/n,
        "fraction_positive_mothers_using_exactly_one_route":
            s/(n-zero) if n!=zero else None,
        "P_lateral_positive_given_apical_negative":
            lat/no_ap if no_ap else None,
        "P_apical_positive_given_lateral_negative":
            ap/no_lat if no_lat else None,
        "conditional_independence_benchmark_P_union":
            pa+pl-pa*pl,
        "observed_minus_independence_benchmark_P_union":
            union-(pa+pl-pa*pl),
        "independence_benchmark_is_only_product_of_marginals":True,
        "conditional_independence_is_a_realistic_biological_null":False,
        "observed_correlated_routes_are_molecularly_independent":False,
        "organ_removal_causal_selection_estimated":False,
        "vegetative_clonal_or_future_generation_success_not_measured":True,
    }


def audit_maternal_route_coverage(mothers=None,strict=True):
    if mothers is None:
        mothers=read_locked(BASE/"mother_plants.csv","mother_plants.csv",
                            ("plt","year","trt","dam","ins",
                             "plt.gsd","ap.hd","o.hds"))
    old=maternal_routes(mothers,strict=strict)
    route={year:{
        arm:coverage_from_counts(
            old["maternal_route_states"][year][arm]["states"])
        for arm in ARMS
    } for year in YEARS}
    # Do not invent an independently measured treatment-specific route
    # model. Existing v2 exact Fisher tests are reused with their
    # eight-group multiplicity adjustment.
    paired={}
    p={}
    for year in YEARS:
        for arm in ARMS:
            x=route[year][arm]
            a=x["four_anatomical_routes"]["apical_only"]
            b=x["four_anatomical_routes"]["lateral_only"]
            name=year+"_"+arm
            pv=exact_discordance_test(a,b)
            p[name]=pv
            paired[name]={
                "apical_only_vs_lateral_only":[a,b],
                "number_of_discordant_mothers":a+b,
                "p_two_sided":pv,
                "paired_unit":"maternal_individual_not_two_independent_heads",
                "route_outcome_not_natural_gene_reexpression":True,
            }
    q=bh_correction(p)
    for name in paired:
        paired[name]["q_BH_eight"]=q[name]
    src=old["source"]
    require(src["dataset_doi"]==ORIGINAL_DOI and
            src["original_paper_doi"]==ORIGINAL_PAPER,
            "USDA_SOURCE_PROVENANCE_CHANGED")
    if strict:
        require(src["coherent_mother_rows"]==105 and
                src["inconsistent_mother_identifiers"]==["2007|327"],
                "SOURCE_JOINT_ACCOUNTING_CHANGED")
        require(route["2008"]["di"]["four_anatomical_routes"]=={
            "no_good_seed":0,"apical_only":4,"lateral_only":8,"both_sites":7},
                "SOURCE_2008_PROTECTED_DAMAGED_ROUTE_CHANGED")
        require(route["2008"]["dni"]["four_anatomical_routes"]=={
            "no_good_seed":5,"apical_only":0,"lateral_only":7,"both_sites":0},
                "SOURCE_2008_AMBIENT_DAMAGED_ROUTE_CHANGED")
    return {
        "schema":"aza3_original_Cirsium_maternal_joint_route_coverage_v3",
        "source":{
            "dataset_doi":ORIGINAL_DOI,
            "original_study_DOI":ORIGINAL_PAPER,
            "source_owner":"West & Louda (2021) USDA CC0 original experiment",
            "genuine_maternal_rows":src["original_mother_rows"],
            "jointly_consistent_mothers":src["coherent_mother_rows"],
            "original_individual_head_rows":
                src["original_heads_in_related_USDA_source"],
            "bad_joint_maternal_route_IDs":
                src["inconsistent_mother_identifiers"],
            "USDA_original_git_blob_mothers":
                "458800272c3d0a03c7653702027a8f8578d9e88c",
            "not_a_new_independent_Cirsium_sample":True,
        },
        "year_arm_observed_coverage":route,
        "paired_route_discordance_exact_tests":paired,
        "original_alternative_route_presence_tests":
            old["lateral_presence_fisher_tests"],
        "biology":{
            "sibling_heads_are_present_before_this_season_not_evolutionarily_resurrected":True,
            "same_mother_but_not_same_genetic_program_assayed":True,
            "joint_structural_seed_states_not_causal_architecture_manipulations":True,
            "apical_damage_not_equal_complete_anatomical_apex_removal":True,
            "normal_insects_and_reduced_insects_change_multiple_interactors":True,
            "unadjusted_conditional_outcome_not_independent_genetic_trait":True,
        },
        "biological_claim_limits":{
            "common_antagonist_proven_to_cause_both_failures":False,
            "cause_of_lateral_seed_success_separately_identified":False,
            "literal_removal_of_apical_route_experiment_conducted":False,
            "alternative_heads_reconstruct_lost_ancestral_stamen_gene":False,
            "structural_division_of_labour_fitness_quantified":False,
            "ancestral_genetic_retention_or_regain_measured":False,
            "established_descendant_recruitment_recorded":False,
            "treatment_is_selective_to_antagonists_excluding_pollinators":False,
            "two_years_are_same_individual_longitudinal_records":False,
            "source_floral_seed_stage_equals_pollinator_pollen_deposition":False,
        },
        "conclusion":
            "OBSERVED_CONDITIONAL_MATERNAL_ROUTE_COVERAGE_NOT_CAUSAL_ARCHITECTURE_FITNESS_OR_NATIVE_REGAIN",
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=audit_maternal_route_coverage(strict=True)
    data=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(data,encoding="utf-8")
    else:
        print(data)


if __name__=="__main__":
    main()
