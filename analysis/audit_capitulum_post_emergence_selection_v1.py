#!/usr/bin/env python3
"""Adversarial post-treatment emergence-selection audit for head-level seed damage.

Xi et al. (2015), DOI 10.1111/1365-2656.12361:
- Cirsium setosum is part of five-species observational comparisons.
- Saussurea nigrescens experimental arms introduced flies +/- parasitoids.
- The paper selected fly-only heads with emerging flies and fly+wasp
  heads with emerging wasps. This is conditioned on POST treatment emergence.
No access/spine/orientation/viability effects on Cirsium are identified.

The deterministic examples prove that DIFFERENT post-treatment inclusion
mechanisms can create a selected-head difference when the TRUE unselected
treatment effect is exactly zero. They are NOT re-analyses of Xi's raw data.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

def group_damage(high_count=50,low_count=50,high_damage=20.0,low_damage=5.0):
    if min(high_count,low_count)<0 or not high_count+low_count:
        raise ValueError("Nonempty nonnegative head group needed")
    return (high_count*high_damage+low_count*low_damage)/(high_count+low_count)

def selected_mean(n_high, n_low, high_damage=20.0, low_damage=5.0):
    if not all(isinstance(x,int) and x>=0 for x in (n_high,n_low)):
        raise ValueError("Counts must be nonnegative integers")
    if n_high+n_low==0:
        return None
    return (n_high*high_damage+n_low*low_damage)/(n_high+n_low)

def pure_selection_counterexample():
    # Both arms start with the SAME 50:50 low/high head damage composition.
    # Every head's potential outcome Y(0) = Y(1), so intervention causal
    # effect is identically zero. Only the post-treatment inclusion changes.
    all_control=group_damage()
    all_wasp=group_damage()
    fly_only=selected_mean(10,40)   # 10 high / 40 low heads selected
    fly_and_wasp=selected_mean(33,17)  # 33 high / 17 low selected
    assert all_control==all_wasp==12.5
    assert abs(fly_only-8)<1e-12
    assert abs(fly_and_wasp-14.9)<1e-12
    assert abs(fly_and_wasp/fly_only-1.8625)<1e-12
    assert fly_and_wasp>fly_only and all_control==all_wasp
    return {
        "status":"PASS_DETERMINISTIC_POST_TREATMENT_SELECTION_COUNTEREXAMPLE",
        "full_randomized_arm_damage_means":{"flies":all_control,
                                           "flies_plus_wasps":all_wasp},
        "unconditional_true_causal_difference":all_wasp-all_control,
        "conditional_head_means_after_differential_emergence_selection":{
            "fly_emerged_selection":fly_only,
            "wasp_emerged_selection":fly_and_wasp,
        },
        "observed_selected_ratio":fly_and_wasp/fly_only,
        "heads_selected":{"fly_only":50,"flies_and_wasps":50},
        "differential_selection_can_mimic_positive_damage_effect":True,
        "biological_observations":0,
        "not_an_estimate_of_bias_in_Xi_2015":True,
    }

def publishable_gate(data):
    # Real head-level causal contrast requires intention-to-observe cohorts
    # with original treatment assignment and all outcomes, rather than only
    # emergent-insect conditional samples.
    necessary={
        "randomized_head_or_enclosure_assignment",
        "all_heads_before_outcome",
        "outcome_or_documented_reason_for_missing",
        "enclosure_id_and_biological_replication",
        "host_stage_and_prior_infestation",
        "emergence_as_post_treatment_variable_only",
        "viable_filled_achenes_or_seed_damage",
        "true_anatomical_armature_orientation",
    }
    if not isinstance(data,dict):
        raise ValueError("Audit gates must be named")
    met={x for x in necessary if data.get(x) is True}
    status="PRE_SPECIFIED_ESTIMAND_DESIGN_READY_NOT_CAUSALLY_FITTED" if len(met)==len(necessary) else "HOLD_SELECTION_OR_TRAIT_IDENTIFIABILITY"
    return {"status":status,"n_requirements":len(necessary),
            "met":sorted(met),"missing":sorted(necessary-met)}

def tests():
    sim=pure_selection_counterexample()
    none=publishable_gate({})
    all_ok=publishable_gate({
        x:True for x in (
          "randomized_head_or_enclosure_assignment","all_heads_before_outcome",
          "outcome_or_documented_reason_for_missing",
          "enclosure_id_and_biological_replication",
          "host_stage_and_prior_infestation",
          "emergence_as_post_treatment_variable_only",
          "viable_filled_achenes_or_seed_damage",
          "true_anatomical_armature_orientation"
        )
    })
    assert none["status"]=="HOLD_SELECTION_OR_TRAIT_IDENTIFIABILITY"
    assert len(none["missing"])==8
    assert all_ok["status"]=="PRE_SPECIFIED_ESTIMAND_DESIGN_READY_NOT_CAUSALLY_FITTED"
    assert not all_ok["missing"]
    assert selected_mean(0,0) is None
    assert selected_mean(5,5)==12.5
    return {
        "status":"PASS_EMERGENCE_SELECTION_AND_GO_NO_GO_TESTS",
        "selection_counterexample":sim,
        "real_Cirsium_head_effect":"NOT_IDENTIFIABLE",
        "prospective_requirements":none["missing"],
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    report={
        "version":"capitulum_emergence_conditioned_selection_v1",
        "date":"2026-10-08",
        "empirical_source":{
            "paper":"Xi, Eisenhauer & Sun 2015",
            "doi":"10.1111/1365-2656.12361",
            "focal_experiment_plant":"Saussurea nigrescens",
            "Cirsium_setosum_part":"observational only",
            "selection_rule":"Fly-only treatment retained heads with emerging flies; fly+parasitoid treatment retained heads with emerging wasps. Heads in other arms were randomly sampled.",
            "reported_per_capitulum_damaged_seed_means":{"flies":8.04,
                "flies_plus_wasps":14.83},
            "paper_conclusion":"Corroborating microcosm shows prolonged development; selected head comparison itself is not a full intention-to-treat estimate.",
        },
        "synthetic_test":tests(),
        "design_status":publishable_gate({}),
        "prohibited_interpretations":[
            "Selected-head association identifies population-average treatment effect",
            "The illustrative selection mechanism explains the observed Xi difference",
            "Saussurea estimates are Cirsium true-spine causal effects",
            "Host emergence proves early seed rescue",
            "Number of emerged wasps measures plant fitness",
        ],
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",
                        encoding="utf-8")
    print(json.dumps({
      "status":report["synthetic_test"]["status"],
      "unconditional_causal_difference":0,
      "illustrative_selected_damage_ratio":
       report["synthetic_test"]["selection_counterexample"]["observed_selected_ratio"],
      "biological_observations":0,
    },indent=2))

if __name__=="__main__":
    main()
