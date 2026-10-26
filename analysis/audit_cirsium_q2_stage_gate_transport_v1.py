#!/usr/bin/env python3
"""Actual original Cirsium exclosure -> *conditional* later-stage payoff transport.

Real source: Klinkhamer et al. 1988 Journal of Ecology DOI 10.2307/2260602.
The input is the publisher-paper source Table 2 manual 10-row transcript
(1981 Sep-1982 Feb, five dates x two habitat groups x cage/net/open), NOT
the missing raw individual seed-container histories. Source consistency
is checked by the independent, already merged original-table audit.

What new biology is (and isn't) identifiable:
- At fixed head-level seed supply H_z and fixed common postdispersal survival
  q shared between head configurations, E[retained seeds_z] = H_z*q,
  so postdispersal loss attenuates a difference but CANNOT reverse its sign.
- If q depends on z via seed traits, actual dispersal microhabitat or predator
  access, E[retained seeds_z] = H_z*q_z and a reversal requires
  H_1/H_0 > q_0/q_1. Conditional algebra, not an evolutionary theorem claim.
- Published cage vs open effects identify *postdispersal physical access*,
  NOT how head orientation, phyllary posture or spines change q(z).
- Every hypothetical H_ratio is labelled as a math scenario. The authors'
  real 0.6 cm cage itself is no genotype or capitulum-trait manipulation.
- The C. arvense 2022 flowerhead data are a DIFFERENT taxon, year and cohort,
  with a reported proportion not H. The two sources may NOT be multiplied to
  manufacture lifetime plant fitness or a quantitative SLK R-K margin.
- Even a common upstream stem-survival factor cannot change sign, but source
  1982 "98% stems attacked" is not "98% died before making heads".
- Seed eaten/retained is NOT viability, germination, recruitment or fitness.

The real empirical contribution is a fully enumerated 10-cell
predation-retention tradeoff and a prospective necessary-measurement gate,
not simulated original plant records or a measured heritable trait value.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from audit_cirsium_vulgare_1988_multistage_herbivore_gate_v1 import (
    audit_input as verify_original_1988_paper_table
)

SCENARIOS = (0.5,0.75,0.8,0.9,1.0)
YEARS = (1981,1982)
SOURCE_DOI = "10.2307/2260602"
CROSS_SOURCE_DOI = "10.21203/rs.3.rs-6634172/v1"
NUM_TOL=1e-12

def require(statement,reason):
    if not statement:
        raise ValueError(reason)

def mean(numbers):
    require(len(numbers)>0,"EMPTY_SOURCE_GROUP")
    return sum(numbers)/len(numbers)

def pos_unit(value,nonzero=False):
    return (type(value) in (int,float) and math.isfinite(value)
            and 0<=value<=1 and (not nonzero or value>0))

def constant_survival_difference(H1,H0,q):
    """A non-mechanistic identity if every pre-dispersal cost is included in H."""
    require(all(type(x) in (int,float) and math.isfinite(x) for x in (H1,H0)),
            "HEAD_OUTPUT_UNDEFINED")
    require(H0>=0 and H1>=0 and pos_unit(q),
            "INVALID_LATE_STAGE_SURVIVAL")
    return (H1-H0)*q

def conditional_survival_difference(H1,H0,q1,q0):
    """End-point difference; does not identify whether H or q is due to a trait."""
    require(all(type(x) in (int,float) and math.isfinite(x) and x>=0
                for x in (H1,H0)) and
            all(pos_unit(q) for q in (q1,q0)),
            "INVALID_JOINT_STAGE_PRODUCT")
    return H1*q1-H0*q0

def rank_threshold_for_different_survival(q0,q1):
    require(pos_unit(q0) and pos_unit(q1,True),
            "CANNOT_IDENTIFY_SURVIVAL_RATIO")
    return q0/q1

def assert_no_cross_study_lifetime_merge(original,unlinked):
    a=original["source"]
    b=unlinked["source"]
    require(a["doi"]==SOURCE_DOI and
            b["doi"]==CROSS_SOURCE_DOI and
            b["study_year"]==2022 and
            unlinked["source_counts"]["heads_all"]==615 and
            unlinked["source_counts"]["complete_four_arm_blocks"]==36,
            "TRUE_SOURCE_COHORT_INTEGRITY_CHANGED")
    require(b["doi"]!=a["doi"] and
            b["study_year"] not in YEARS,
            "FALSE_SAME_PLANT_REPRODUCTIVE_COHORT")
    require(unlinked["limits"]["new_biological_field_observations"]==0,
            "FALSE_STUDY_REPLICATION")
    measured=unlinked["comparisons"]["bag_minus_control"]["measures"]
    require(measured["viable_fraction"]["mean"]< -.5,
            "ORIGINAL_ARVENSE_DIFFERENT_SOURCE_DRIFT")
    return {
        "source":(
            "Cirsium arvense original 2022 4-arm intervention, not "
            "Cirsium vulgare 1981/1982 exclosure"),
        "paper_doi":b["doi"],
        "original_head_rows":unlinked["source_counts"]["heads_all"],
        "matched_four_arm_blocks":
            unlinked["source_counts"]["complete_four_arm_blocks"],
        "author_defined_bag_minus_control_good_seed_fraction":
            measured["viable_fraction"]["mean"],
        "not_a_head_module_treatment":True,
        "not_a_count_H_or_genotype_matched":True,
        "no_original_individual_joint_source_link_to_1988":True,
    }

def joint_alignment_witness():
    """Constructive same-marginal opposite-joint-alignment counterexample.

    All values are SYNTHETIC. This demonstrates why product of marginal
    head fecundity means and mean per-head seed survival is generally not
    the expected surviving-seed output without source head-batch linkage.
    """
    H=[8,2]
    q_aligned=[.9,.1]
    q_opposed=[.1,.9]
    aH=mean(H)
    aq=mean(q_aligned)
    y_aligned=mean([h*q for h,q in zip(H,q_aligned)])
    y_opposed=mean([h*q for h,q in zip(H,q_opposed)])
    require(math.isclose(y_aligned,3.7) and
            math.isclose(y_opposed,1.3) and
            math.isclose(aH*aq,2.5),
            "SYNTHETIC_JOINT_MARGINAL_LOGIC_CHANGED")
    return {
        "status":"CONSTRUCTED_MATHEMATICAL_COUNTEREXAMPLE_NOT_ORIGINAL_SOURCE_RECORDS",
        "source_field_rows_used":0,
        "heterogeneous_expected_value_identity":(
            "E[H_i*q_i]=E[H_i]*E[q_i]+Cov(H_i,q_i). Even when head output "
            "and per-head retention have IDENTICAL MARGINAL distributions, "
            "expected surviving output can differ if their within-head "
            "pairing differs."),
        "both_scenarios_same_original_head_seed_counts_H":H,
        "both_scenarios_same_unweighted_retention_marginal_sorted":
            sorted(q_aligned),
        "mean_H":aH,
        "mean_unweighted_q":aq,
        "naive_marginal_product":aH*aq,
        "aligned_H_and_q_scenario":{
            "q_assigned_to_H_in_original_H_order":q_aligned,
            "mean_H_times_q":y_aligned,
            "covariance_H_q":y_aligned-aH*aq,
            "head_seed_weighted_retention":y_aligned/aH,
        },
        "opposed_H_and_q_scenario":{
            "q_assigned_to_H_in_original_H_order":q_opposed,
            "mean_H_times_q":y_opposed,
            "covariance_H_q":y_opposed-aH*aq,
            "head_seed_weighted_retention":y_opposed/aH,
        },
        "difference_in_expected_retained_seeds_between_scenarios":
            y_aligned-y_opposed,
        "joint_seed_batch_trace_required_to_identify_which_alignment_applies":
            True,
        "stage_1988_table_marginal_lacks_head_link":True,
    }


def analyze(source,unlinked):
    verify_original_1988_paper_table(source)
    other=assert_no_cross_study_lifetime_merge(source,unlinked)
    t=source["postdispersal_exclosure"]
    original=t["rows"]
    require(len(original)==10 and len({r["date"] for r in original})==5
            and len({r["habitat"] for r in original})==2,
            "WRONG_EXCLOSURE_DATE_HABITAT_GRAIN")
    obs=[]
    for r in original:
        retention={
            k:1-r[k]/100 for k in ("cage","net","open")
        }
        require(all(pos_unit(v,nonzero=True) for v in retention.values()),
                "SOURCE_RETENTION_NOT_IN_UNIT_INTERVAL")
        obs.append({
            "date":r["date"],
            "habitat":r["habitat"],
            "source_seed_predation_percentage":{
                k:r[k] for k in ("cage","net","open")},
            "measured_seed_retention_per_exposure_treatment":retention,
            "source_recorded_cage_minus_open_retention":
                retention["cage"]-retention["open"],
            "source_recorded_net_minus_open_retention":
                retention["net"]-retention["open"],
            "head_seed_yield_ratio_threshold_for_hypothetical_cage_route_to_exceed_open":
                rank_threshold_for_different_survival(
                    retention["open"],retention["cage"]),
            "cage_over_open_per_seed_retention_ratio":
                retention["cage"]/retention["open"],
        })
    qs={k:[r["measured_seed_retention_per_exposure_treatment"][k]
             for r in obs] for k in ("cage","net","open")}
    pred={k:[r["source_seed_predation_percentage"][k] for r in obs]
          for k in ("cage","net","open")}
    avgret={k:mean(v) for k,v in qs.items()}
    threshold=[r["head_seed_yield_ratio_threshold_for_hypothetical_cage_route_to_exceed_open"]
               for r in obs]
    summary={
        "n_habitat_date_groups":len(obs),
        "n_dates":len({r["date"] for r in obs}),
        "n_habitats":len({r["habitat"] for r in obs}),
        "mean_predation_percentage":{
            k:mean(pred[k]) for k in ("cage","net","open")},
        "mean_seed_retention_percentage":{
            k:100*avgret[k] for k in ("cage","net","open")},
        "within_habitat_date_cage_retention_exceeds_open":sum(
            r["source_recorded_cage_minus_open_retention"]>0 for r in obs),
        "within_habitat_date_cage_retention_exceeds_net":sum(
            r["measured_seed_retention_per_exposure_treatment"]["cage"]>
            r["measured_seed_retention_per_exposure_treatment"]["net"]
            for r in obs),
        "mean_matched_cage_minus_open_survival_percentage_points":
            100*mean([r["source_recorded_cage_minus_open_retention"]
                      for r in obs]),
        "ratio_of_group_mean_retention_cage_divided_by_open":
            avgret["cage"]/avgret["open"],
        "mean_of_cellwise_ratios_cage_divided_by_open":
            mean([r["cage_over_open_per_seed_retention_ratio"] for r in obs]),
        "source_break_even_head_yield_ratio_min":min(threshold),
        "source_break_even_head_yield_ratio_max":max(threshold),
        "ratio_of_group_mean_retention_open_divided_by_cage":
            avgret["open"]/avgret["cage"],
        "n_source_cells_where_bird_net_retention_is_below_open":sum(
            r["measured_seed_retention_per_exposure_treatment"]["net"]<
            r["measured_seed_retention_per_exposure_treatment"]["open"]
            for r in obs),
    }
    scen=[]
    for ratio in SCENARIOS:
        check=[
            {
                "cell":r["date"]+"|"+r["habitat"],
                "normalized_survivor_advantage_if_hypothetical_head_ratio_and_cage_vs_open":
                    conditional_survival_difference(
                        ratio,1,
                        r["measured_seed_retention_per_exposure_treatment"]["cage"],
                        r["measured_seed_retention_per_exposure_treatment"]["open"]),
            }
            for r in obs
        ]
        values=[x["normalized_survivor_advantage_if_hypothetical_head_ratio_and_cage_vs_open"]
                for x in check]
        scen.append({
            "hypothetical_predispersal_head_seed_output_ratio_protected_over_open":
                ratio,
            "n_habitat_dates_where_protected_head_route_would_yield_more_surviving_exposed_seeds":
                sum(x>NUM_TOL for x in values),
            "n_habitat_dates_where_open_route_would_yield_more":
                sum(x< -NUM_TOL for x in values),
            "n_exact_ties":sum(abs(x)<=NUM_TOL for x in values),
            "normalized_differences_by_source_cell":check,
            "head_output_ratio_is_observed":False,
            "trait_changes_seed_caging_or_deposition_observed":False,
        })
    result={
        "version":"aza3_q2_cirsium_stage_transport_conditional_access_original_JS_v1",
        "as_of":"2026-10-26",
        "source":{
            "source_doi":source["source"]["doi"],
            "paper_year":1988,
            "source_type":"published original Table 1/2 re-transcription only, not individual raw data",
            "exclosure_original_years":[1981,1982],
            "head_stem_years":[1982,1984],
            "site":"Meijendel, Netherlands",
            "original_number_of_exclosure_habitat_date_cells":len(original),
            "other_head_source":other,
        },
        "definitions":{
            "units_of_yield":"mature seeds per original head entering original postdispersal exposure; expected retained seeds, NOT viable seedlings",
            "q_definition":"1 - table2 treatment-specific source percentage seeds consumed/100; conditional on reported exposure experiment",
            "H_definition":("Expected mature head output including pre-dispersal "
                "costs; heterogeneous cohorts need joint H_i × q_i "
                "measurements, not product of separate unweighted group means"),
            "no_same_seed_ID_link_across_1988_cohorts":True,
            "independent_cross_source_numeric_product_forbidden":True,
            "use_of_cage_as_if_genetic_module_could_change_postdispersal_route_is_scenario_not_empirical_trait_effect":True,
        },
        "source_original_exclosure_cell_replay":obs,
        "source_aggregates":summary,
        "counterfactual_algebra_only":{
            "constant_downstream_route":{
                "expected_surviving_seed_proxy":(
                    "E[Y_z|matched_background]=q * H_z, with H_z measured "
                    "mature seed output inclusive of all pretreatment head "
                    "costs; q common to the two head phenotype states; no "
                    "genotype-dependent postdispersal selection"),
                "difference":"W_1-W_0=q*(H_1-H_0)",
                "condition":"q>0 implies sign(W_1-W_0)=sign(H_1-H_0); q=0 implies both W=0",
                "no_seed_independence_across_individuals_needed_for_expectation":
                    True,
                "only_same_late_s_and_same_background_stage_cost_scale":True,
                "does_not_identify_head_effect_from_original_table":True,
            },
            "variable_downstream_route":{
                "expected_surviving_seed_proxy":(
                    "W_z=q_z*H_z only if q_z denotes the SEED-WEIGHTED "
                    "retention factor for the same z-defined source seed "
                    "cohort (E[H_i*q_i]/E[H_i] when cohort has "
                    "heterogeneous head counts), or every head has common "
                    "retention; multiplying separate unweighted cohort "
                    "means requires zero relevant covariance."),
                "break_even":(
                    "If H_0>0 and q_1>0, W_1>W_0 iff H_1/H_0>q_0/q_1. "
                    "Actual random cage/open treatment is not a genetically "
                    "controlled Cirsium capitulum trait."),
                "reverse_possible_only_if_correlated_trait_or_environment_conditional_postdispersal_route":
                    True,
                "matched_outcome_trial_required":True,
            },
            "stem_stage":{
                "head_eligible_survival":(
                    "W_z=a_z*q_z*H_z with a_z probability of making "
                    "seed-bearing heads; a_z can be correlated with prior "
                    "stem attack/plant size; not equal to 1-% any stem attack"),
                "if_a_z_equal_and_positive":(
                    "The common upstream factor a cannot reverse a head-only "
                    "ranking, but attenuates its absolute reproductive increment"),
                "original_1982_any_stem_attack_fraction_is_not_stem_death_probability":
                    True,
            },
            "costs_and_limits":{
                "head_output_H_must_include_all_pre_dispersal_head_costs":True,
                "offspring_establishment_or_late_architecture_costs_need_separate_term":
                    True,
                "no_claim_of_architecture_global_advantage_slk_Phi_or_SCH_conflict":
                    True,
                "1988_cage_effect_cannot_be_inferred_for_head_morphology":True,
            },
            "declared_illustrative_Hratios":scen,
            "joint_alignment_nonidentifiability_witness":
                joint_alignment_witness(),
        },
        "measurement_gate":{
            "same_mother_and_genotyped_original_head_records_observed":False,
            "head_module_measurements_orientation_bracts_true_spines_stickiness":
                False,
            "matched_pollinator_deposition_and_oviposition_observed":False,
            "original_heads_mature_seed_output_measured_but_not_joint_trait_treatment":
                True,
            "same_head_seed_batch_split_to_randomized_postdispersal_access":False,
            "actual_seed_borne_genotype_specific_retention_estimated":False,
            "proposed_actual_field_design":(
                "prospectively register original head IDs and phenotype; "
                "measure seed output under valid original P x E sham "
                "interventions, split real ID-tied mature seeds among matched "
                "microhabitat predator-exclusion/open blocks with permissions, "
                "measure retention and germination without multiplying "
                "unlinked years or species; preserve same-source head-seed "
                "alignment or estimate seed-weighted retention to avoid "
                "E[H]E[q] error"),
            "status":"SOURCE_RETENTION_QUANTIFIED_BUT_HEAD_MODULE_TRANSMISSION_NOT_IDENTIFIABLE",
            "new_field_observations":0,
        },
    }
    require(summary["n_habitat_date_groups"]==10 and
            summary["within_habitat_date_cage_retention_exceeds_open"]==10 and
            abs(summary["mean_matched_cage_minus_open_survival_percentage_points"]-57.3)<NUM_TOL,
            "REAL_CIRSIUM_EXCLOSURE_STAGE_GATE_DRIFT")
    require(scen[2]["n_habitat_dates_where_protected_head_route_would_yield_more_surviving_exposed_seeds"]==8,
            "HYPOTHETICAL_0P8_RATIO_CELL_GATE_DRIFT")
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-1988",type=Path,required=True)
    p.add_argument("--separate-source-2022",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    source=json.loads(args.source_1988.read_text(encoding="utf-8"))
    head=json.loads(args.separate_source_2022.read_text(encoding="utf-8"))
    receipt=analyze(source,head)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(receipt,ensure_ascii=False,indent=2,
                                   allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "source":receipt["source"]["source_doi"],
        "aggregates":receipt["source_aggregates"],
        "stage_gate":receipt["measurement_gate"]["status"],
        "scenario_counts":[{
            "ratio":x["hypothetical_predispersal_head_seed_output_ratio_protected_over_open"],
            "positive":x["n_habitat_dates_where_protected_head_route_would_yield_more_surviving_exposed_seeds"]}
            for x in receipt["counterfactual_algebra_only"]["declared_illustrative_Hratios"]],
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
