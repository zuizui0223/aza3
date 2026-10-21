#!/usr/bin/env python3
"""What can original Cirsium pitcheri rank-specific damage MARGINS identify?

Source 2020 Gijsman, Havens & Vitt, GEC DOI 10.1016/j.gecco.2020.e00945,
is represented only by transparently transcribed published summaries.
Original authors state individual mother/head data are available on request,
not freely downloaded. The 2021 Dryad workbook is NOT accessible in the
current project through public API (HTTP 401); no individual raw file input.

The 2016 WFD proportions 27/50/80% are head-weighted rank marginal field
summaries. Fréchet-Hoeffding calculations conditional on treating them as
marginals of ONE hypothetical head per rank are mathematically SHARP.
They are NOT necessarily bounds for the actual source plant-weighted
within-parent co-infestation because plants have unequal numbers of heads.
Even not-infested does NOT mean seed maturation or seedling recruitment.

The equal-effort, rank-invariant .60 injury and whole-head-loss scenarios
below are explicitly COUNTERFACTUAL assumptions, not fitted field effects.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


SOURCE_DOI = "10.1016/j.gecco.2020.e00945"
EXPECTED_P = {"terminal": .27, "secondary": .50, "tertiary": .80}
SOURCE_2020_ACCESS = "AUTHOR_REQUEST_REQUIRED_FOR_INDIVIDUAL_HEAD_OR_PLANT_DATA"
SCENARIO_RATIOS = (.25, .50, .75, 1.0)


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def _unit(p):
    return type(p) in (int, float) and math.isfinite(p) and 0 <= p <= 1


def hypothetical_one_per_rank(p_list):
    """Sharp bounds conditional on exact marginals for co-present ranks.

    Fréchet joint upper = min(p), lower = max(0,sum(p)-n+1).
    At least one uninfested = 1 - P(all infested).
    The product is a diagnostic independence REFERENCE, not an estimate.
    """
    require(len(p_list) >= 2, "NEED_AT_LEAST_TWO_HEAD_RANKS")
    require(all(_unit(p) for p in p_list), "INVALID_MARGINAL")
    n=len(p_list)
    lo=max(0.0, sum(p_list)-(n-1))
    hi=min(p_list)
    independent_all=math.prod(p_list)
    unaffected=[1-hi, 1-lo]
    extra_over_first=[unaffected[0]-(1-p_list[0]),
                      unaffected[1]-(1-p_list[0])]
    # avoid subtraction-induced negative zero in report
    extra_over_first=[0.0 if abs(z)<1e-14 else z for z in extra_over_first]
    return {
        "n_hypothetical_rank_positions":n,
        "input_published_like_rank_marginals":[float(x) for x in p_list],
        "joint_all_infested_sharp_bounds":[lo,hi],
        "at_least_one_free_sharp_bounds":unaffected,
        "independent_reference_all_infested":independent_all,
        "independent_reference_at_least_one_free":1-independent_all,
        "independent_reference_extra_over_terminal_alone":
            p_list[0]-independent_all,
        "sharp_extra_chance_over_terminal_alone":extra_over_first,
        "assumption_exact_same_plant_one_per_rank_marginals_verified":False,
        "independence_of_head_infestation_observed":False,
        "not_infested_equals_reproductive_success":False,
        "not_inferred_offspring_fitness":True,
    }


def equal_effort_mechanism_scenarios(pt, ps, pr, loss=.60):
    """Mean allocation break-even under additional explicitly false/unverified assumptions.

    One reproductive-effort unit if invested terminal has yield:
       (1-d)*(1-L*pt)*fT ; secondary: (1-L*ps)*fS
    where d is INDEPENDENT catastrophic terminal loss, L is a UNIFORM
    loss per infested head, and fS/fT is intrinsic productivity per effort.
    Equating gives d* = 1 - (fS/fT)(1-L ps)/(1-L pt).
    Neither d nor this rank-uniform L is independently observed in source.
    This is a prospective mechanism discriminator, not a source estimate.
    """
    require(all(_unit(x) for x in (pt, ps, pr, loss)),
            "SCENARIO_PARAMETER_OUTSIDE_PROBABILITY_DOMAIN")
    f={name:1-loss*p for name,p in
       (("terminal",pt),("secondary",ps),("tertiary",pr))}
    require(f["terminal"]>0,"TERMINAL_BASELINE_EXPECTED_YIELD_ZERO")
    return {
        "scenario": "EQUAL_REPRODUCTIVE_EFFORT__RANK_INVARIANT_WEevil_DAMAGE__INDEPENDENT_TERMINAL_CATASTROPHE",
        "published_pooled_infested_head_mature_seed_count_deficit_used_only_as_hypothetical_uniform_loss":loss,
        "remaining_mature_seed_proxy_if_weevil_loss_rate_rank_invariant":f,
        "whole_terminal_loss_fraction_break_even_for_secondary_equal_cost_resource_ratios":{
            str(r):max(0.,1.-r*f["secondary"]/f["terminal"])
            for r in SCENARIO_RATIOS},
        "assumptions_confirmed_in_original_study":False,
        "rank_invariant_weevil_effect_is_not_empirically_established":True,
        "whole_terminal_loss_independent_of_weevil_risk_is_not_empirically_established":True,
        "secondary_intrinsic_fecundity_per_cost_not_measured":True,
        "terminal_absence_2017_can_be_plugged_into_2016_weevil_scenario_as_estimate":False,
        "arithmetic_mean_fitness_benefit_of_diversification_is_not_observed":True,
        "no_assumption_that_secondary_or_tertiary_heads_are_costless":True,
    }


def partial_recovery_bound(missing_terminal_fraction):
    """With no observed P(viable secondary | missing terminal), only [0,pA] is sharp."""
    require(_unit(missing_terminal_fraction),"MISSING_TERMINAL_FRACTION_INVALID")
    return {
        "missing_terminal_fraction":missing_terminal_fraction,
        "joint_missing_terminal_and_successful_secondary_sharp_bound_without_secondary_marginal":
            [0.,missing_terminal_fraction],
        "possible_successful_secondary_given_terminal_missing":[0.,1.],
        "has_original_plant_level_joint_seed_counts":False,
        "observed_deer_cause_for_each_missing_terminal":False,
        "scope":"Bound on population share with potentially compensatory secondary heads; not a causal contribution or a probability estimate of viable seedlings.",
    }


def analyze(source):
    s=source["source"]
    require(s["doi"]==SOURCE_DOI and
            s["primary_published_article_transcription_not_raw_head_rows"] is True,
            "SOURCE_DOI_OR_GRAIN_DRIFT")
    w=source["rank_weevil_2016_WFD"]
    p=w["infestation_evidence_fractions_approx"]
    require(set(p)==set(EXPECTED_P) and all(
        abs(p[k]-EXPECTED_P[k])<1e-12 for k in EXPECTED_P),
        "ORIGINAL_SOURCE_WEevil_RATE_CHANGED")
    require(w["rank_causal_weevil_assignment"] is False,
            "NO_RANDOMIZED_WEevil_ACCESS")
    require(w["mature_seed_infested_to_uninfested_approx_ratio"]==.40,
            "ORIGINAL_SEED_COUNT_DAMAGE_RATIO_CHANGED")
    m=source["terminal_head_absence_by_site_year"]
    require(m["specific_herbivore_causation_identified"] is False,
            "DEER_CAUSE_UNVERIFIED")
    require(m["2017"]=={"SC":.27,"WFD":.33}
            and m["2016"]=={"SC":0,"WFD":.089},
            "ORIGINAL_MISSING_TERMINAL_FRACTIONS_CHANGED")
    y=source["field_seed_output_by_rank_year"]
    require(y["tertiary_2017_mature_seeds_both_sites"]==0
            and abs(y["secondary_year_difference_p"]-.24)<1e-12,
            "SOURCE_OFFSPRING_STAGE_EVIDENCE_DRIFT")
    limits=source["material_limit"]
    require(limits["original_plant_and_head_rows_not_downloaded"] is True,
            "ORIGINAL_PLANT_HEAD_DATA_NOT_HELD")
    require(limits["2020_and_2021_seeds_not_joined_by_original_head_id"] is True,
            "RAW_MOTHER_SEED_LINK_INVENTED")
    require(source["dryad_blocked_separately"]["no_verified_raw_xlsx"] is True,
            "UNVERIFIED_DRYAD_DATA_ACCESS_MISLABELED")
    return {
        "version":"cirsium_pitcheri_2020_original_marginal_risk_partial_identification_v1",
        "source_doi":SOURCE_DOI,
        "source_numerical_status":"PUBLISHED_ROUNDED_HEAD_RANK_MARGINS_NOT_ORIGINAL_WITHIN_MOTHER_DATA",
        "source_2020_individual_head_data_access":SOURCE_2020_ACCESS,
        "published_input":{
            "original_WFD_2016_per_head_weevil_rate":p,
            "source_sampled_head_rank_not_fixed_per_plant":True,
            "observed_2017_terminal_absence_fraction":m["2017"],
            "2017_no_mature_tertiary_seed_in_sampled_heads":True,
        },
        "hypothetical_if_exact_same_plant_one_head_per_rank_marginals":{
            "terminal_secondary":hypothetical_one_per_rank(
                [p["terminal"],p["secondary"]]),
            "terminal_tertiary":hypothetical_one_per_rank(
                [p["terminal"],p["tertiary"]]),
            "secondary_tertiary":hypothetical_one_per_rank(
                [p["secondary"],p["tertiary"]]),
            "terminal_secondary_tertiary":hypothetical_one_per_rank(
                [p["terminal"],p["secondary"],p["tertiary"]]),
        },
        "observed_terminal_absence_and_untested_recovery":{
            site:partial_recovery_bound(m["2017"][site]) for site in ("SC","WFD")
        },
        "equal_cost_mechanism_only_scenarios":equal_effort_mechanism_scenarios(
            p["terminal"],p["secondary"],p["tertiary"],
            1.-w["mature_seed_infested_to_uninfested_approx_ratio"]),
        "inference_limit":{
            "no_plant_specific_compensation_claim":True,
            "not_infested_does_not_imply_filled_achene_output":True,
            "rank_is_not_orientation_or_phyllary_geometry":True,
            "adaptive_rank_insurance_established":False,
            "within_parent_infestation_correlation_measured":False,
            "same_year_joint_weevil_deer_exposure_identified":False,
            "missing_terminal_mechanism_directly_identified":False,
            "source_head_rank_margins_not_guaranteed_as_plant_equal_weighted_margins":True,
            "pooled_head_infestation_cost_not_rank_specific":True,
            "seed_germination_of_same_attacked_head_unavailable":True,
            "new_head_access_field_record_count":0,
        },
        "n_new_independent_biological_field_observations":0,
        "next_primary_data_needed":{
            "source":"2020 ORIGINAL AUTHOR HEAD/MOTHER TABLE",
            "required":["parent_id","site","year","capitulum_id",
                        "rank","head_present","head_age","weevil_damage",
                        "mature_achene_count","filled_achene_mass",
                        "reason_for_missing_head_if_verified",
                        "parent_total_mature_achenes"],
            "needed_to_discriminate":{
                "within_parent_weevil_risk_dependence":
                    "joint_infested_indicator_over_multiple_ranks",
                "actual_categorical_compensation":
                    "P(secondary_mature_seeds>0 | terminal_head_missing)",
                "whole_plant_fitness":
                    "sum mature seeds across ALL pre-enrolled heads per maternal plant",
                "enemy_cause":
                    "verified browse evidence separate from unexplained absent heads",
            },
            "original_source_contact_needed":True,
            "do_not_use_unverified_Dryad_as_matching_2020_head_table":True,
        },
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    q=ap.parse_args()
    src=json.loads(q.source.read_text(encoding="utf-8"))
    result=analyze(src)
    q.out.parent.mkdir(parents=True,exist_ok=True)
    q.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                                allow_nan=False)+"\n",encoding="utf-8")
    p=result["hypothetical_if_exact_same_plant_one_head_per_rank_marginals"]["terminal_secondary"]
    print(json.dumps({"scientific_status":"2020_RAW_MOTHER_LEVEL_RISK_COMPENSATION_NOT_IDENTIFIABLE",
                      "conditional_pair_frechet_bounds":p["at_least_one_free_sharp_bounds"],
                      "conditional_independent_reference":p["independent_reference_at_least_one_free"],
                      "no_actual_plant_pair_risk_measured":True},indent=2))


if __name__=="__main__":
    main()
