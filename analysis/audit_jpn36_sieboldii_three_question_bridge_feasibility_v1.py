#!/usr/bin/env python3
"""Cirsium sieboldii JPN36 bridge: real published history and frozen FIELD plans.

A plant with phenotype and fitness data is not automatically a genotyped
LEGO case. A JPN36 manipulative feasibility head is not automatically one
of the 40 planned C. sieboldii Phase-A genetic individuals. No sampling,
field intervention, DNA sequencing or source observation is authorized here.

Source:
 EAzami pinned c9870e91604ed1514344ddc8f765eae5237dbf27
 data/evidence/jpn36_phyllary_access_pilot_{contract,validation}_v1.json
 sampling/jpn36_phyllary_access_{eligible_heads,assignment,observations}_v1.csv
 aza3 planning/genomic panel 298 individuals, SIE target n=40.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PANEL=ROOT/"data/planning/aza3_nature_genomic_sampling_panel_v1.csv"
BRIDGE=ROOT/"data/contracts/aza3_anto_empirical_lego_bridge_v1.json"
LEGO=ROOT/"data/contracts/aza3_capitulum_higher_order_lego_v2.json"
REF="c9870e91604ed1514344ddc8f765eae5237dbf27"
COMMON=("population_id","individual_id","capitulum_id")
EXPECTED_RECORDS={
 "eligible":set(COMMON)|{
     "site_id","phenological_stage","pair_id","minimum_access_gap_mm",
     "natural_head_orientation_deg","permit_record_id","eligible"},
 "assignment":set(COMMON)|{
     "pair_id","assignment","tissue_injury","treatment_integrity"},
 "observation":set(COMMON)|{
     "pair_id","assignment","observation_bout_id",
     "legitimate_visitor_approaches","enemy_approaches",
     "filled_achenes","total_achenes","visitor_observation_minutes"},
}
FIELDS=("eligible","assignment","observation")

def require(x,reason):
    if not x:raise ValueError(reason)

def csv_rows(path,needed,label):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        r=csv.DictReader(f)
        col=set(r.fieldnames or [])
        require(needed.issubset(col),"PINNED_"+label+"_SCHEMA_CHANGED")
        return list(r)

def source_data(meta,validation,eligible,assigned,observed,panel,bridge,lego):
    require(meta["contract_version"]=="jpn36_phyllary_access_pilot_v1",
            "JPN36_CONTRACT_VERSION_CHANGED")
    require(meta["taxon"]["paper_japan_member_id"]=="JPN_36"
            and meta["taxon"]["accepted_name"]=="Cirsium sieboldii",
            "SOURCE_PHYLOGENETIC_TAXON_NOT_RESOLVED")
    require(validation["fdt8_alignment"]["taxon"]=="Cirsium sieboldii"
            and validation["fdt8_alignment"]["field_first_target"]=="JPN_36",
            "FROZEN_EAZAMI_STATE_NOT_MATCHED")
    require(meta["scope"]["current_field_execution_authorized"] is False
            and validation["field_execution_authorized"] is False,
            "FIELD_AUTHORIZATION_STATE_HAS_CHANGED_REAUDIT_BEFORE_EXECUTION")
    require(meta["scope"]["empirical_result_available"] is False
            and validation["empirical_result_available"] is False,
            "REPORTED_PILOT_RESULTS_NEED_SEPARATE_VERIFICATION")

    rights=meta["rights_and_conservation_gate"]
    reference=meta["natural_reference_gate"]
    device=meta["device_qualification_gate"]
    assignment=validation["sample_allocation"]
    fraction=rights["maximum_fraction_of_counted_flowering_individuals"]
    cap=rights["maximum_experimental_heads"]
    n=assignment["individuals"]
    pairs=assignment["pairs"]
    require(isinstance(fraction,(int,float)) and 0<fraction<=1,
            "SOURCE_FRACTIONAL_CONSERVATION_CAP_INVALID")
    require(n==24 and cap==24 and pairs==12
            and assignment["one_head_per_individual"] is True
            and rights["one_experimental_head_per_individual"] is True,
            "ORIGINAL_JPN36_PAIRED_EXPERIMENT_CONTRACT_CHANGED")
    require(sorted(assignment["arms"])==["access_proxy","sham"],
            "JPN36_ACTIVE_VS_SHAM_TREATMENTS_CHANGED")
    require(reference["minimum_observation_only_individuals"]==30
            and reference["minimum_observation_only_heads"]==30,
            "REFERENCE_COHORT_DENOMINATOR_CHANGED")
    require(device["minimum_heads_per_arm"]==6
            and device["duration_hours"]==24,
            "DEVICE_QUALIFICATION_COHORT_CHANGED")
    required=set(rights["required_before_randomization"])
    require({"land_manager_authorization_approved",
             "site_specific_manipulation_authorization_approved",
             "conservation_review_approved",
             "mature_head_or_achene_collection_authorization_approved",
             "viability_assay_authorization_approved",
             "population_census_complete",
             "device_qualification_passed"}.issubset(required),
            "SOURCE_CONSERVATION_AUTHORIZATION_REQUIREMENTS_LOST")
    require(reference["same_population_and_phenophase_as_pilot"] is True,
            "REFERENCE_POPULATION_OR_STAGE_COMPARABILITY_LOST")

    for label,rows in zip(FIELDS,(eligible,assigned,observed)):
        unique=set()
        for row in rows:
            k=tuple((row.get(c) or "").strip() for c in COMMON)
            require(all(k),"SOURCE_"+label+"_EMPTY_BIOLOGICAL_ID")
            if label!="observation":
                require(k not in unique,"SOURCE_"+label+"_DUPLICATE_BIOLOGICAL_HEAD")
            unique.add(k)
    e_ids={tuple(r[c] for c in COMMON) for r in eligible}
    a_ids={tuple(r[c] for c in COMMON) for r in assigned}
    o_ids={tuple(r[c] for c in COMMON) for r in observed}
    require(a_ids.issubset(e_ids),"ASSIGNED_PILOT_HEAD_NEVER_QUALIFIED")
    require(o_ids.issubset(a_ids),"OBSERVATION_OF_NONASSIGNED_PILOT_HEAD")

    planned=[p for p in panel if p["panel_id"]=="SIE"
             and p["phase"]=="A" and p["taxon"]=="Cirsium sieboldii"]
    require(len(planned)==1 and planned[0]["recommended_n"]=="40",
            "SIEBOLDII_PHASE_A_GENOMIC_TARGET_CHANGED")
    planned_bridge=[p for p in bridge["samplings"]["by_taxon_target"]
                    if p["taxon"]=="Cirsium sieboldii"]
    require(len(planned_bridge)==1 and planned_bridge[0]["planned_n"]==40
            and lego["study_design"]["planned_phase_A_total"]==40,
            "THREE_INDEPENDENT_SOURCE_TARGETS_DISAGREE")
    require(bridge["samplings"]["phase_A_collected_n_currently_admitted"]==0
            and lego["explicit_claim_ceilings"]["actual_collected_individuals_in_this_new_protocol"]==0,
            "PROSPECTIVE_NOT_COLLECTED_STATUS_CHANGED")

    # Necessary aggregate bound for exactly 24 heads at most 2% of
    # counted flowering individuals: do NOT interpret it as permission;
    # site-specific quotas and actual flowering census are still required.
    necessary_census=math.ceil(n/fraction)
    no_actual_overlap_count=0 if not eligible and not assigned and not observed else None
    return {
        "status":"JPN36_SOURCE_PLAN_CONCORDANT_NO_FIELD_EXECUTION_AUTHORIZED",
        "actual_JPN36_source_rows":{
            "eligible_heads":len(eligible),
            "randomized_head_assignments":len(assigned),
            "observational_bouts":len(observed),
            "actual_sample_identity_crosswalk_to_genomic_panel":no_actual_overlap_count,
        },
        "taxon":"Cirsium sieboldii",
        "historical_EAzami":{
            "tip_label":"JPN_36",
            "posture":"appressed",
            "sticky_state":"nonsticky_or_nearly_nonsticky",
            "phylogenetic_terminal_state_is_not_a_causal_functional_test":True,
            "ufboot_terminal_forced_fraction":validation["fdt8_alignment"]["ufboot_terminal_forced_fraction"],
            "n_UFBoot_trees":1000,
        },
        "current_EAzami_pilot":{
            "nonmanipulated_same_population_reference_individuals":reference["minimum_observation_only_individuals"],
            "nonmanipulated_reference_heads":reference["minimum_observation_only_heads"],
            "minimum_device_qualification_heads_per_arm":device["minimum_heads_per_arm"],
            "device_qualification_duration_hours":device["duration_hours"],
            "randomized_independent_plants":n,
            "randomized_pairs":pairs,
            "heads_per_randomized_plant":1,
            "maximum_allowed_experimental_heads":cap,
            "fraction_of_counted_flowering_population_cap":fraction,
            "aggregate_necessary_counted_flowering_individuals_for_full_24_head_pilot":necessary_census,
            "necessary_census_not_sufficient_approval":True,
            "independent_local_site_quota_and_permission_must_be_verified":True,
            "field_execution_authorized":False,
            "all_7_authorization_gates_required":sorted(required),
            "not_an_adaptation_or_effect_size_study":True,
        },
        "prospective_aza3_molecular_design":{
            "separate_Phase_A_SIE_target_individuals":40,
            "actually_admitted_new_phase_A_genotyped_individuals":0,
            "not_proved_24_are_subset_of_40":True,
            "not_proved_30_reference_are_subset_of_40":True,
            "do_not_add_24_and_40_as_independent_unique_individuals":True,
            "same_specimen_ID_must_be_explicitly_crosswalked":True,
        },
        "three_question_readiness":{
            "Q1_lineage_memory":"REAL_EAZAMI_JPN36_SINGLE_TIP_HISTORY_NOT_CROSS_LINEAGE_PREDICTION",
            "Q2_partner_payoff":"REAL_PHYL_POSTURE_WITH_PLANNED_RANDOMIZED_SHAM_BUT_NO_DATA_OR_PERMISSION",
            "Q3_hidden_capacity":"PLANNED_40_SAME_SPECIES_RESOURCE_INVENTORY_WITH_ZERO_FUNCTIONALLY_GENOTYPED_INDIVIDUALS",
            "Q2_x_Q3_same_plant_bridge":"NONE_VERIFIED_24_AND_40_ARE_DISTINCT_PLANNING_DENOMINATORS",
        },
        "first_actionable_gate":"OBTAIN_SITE_APPROVED_OBSERVATION_ONLY_REFERENCE_OF_30_HEADS_FROM_30_DISTINCT_PLANTS_BEFORE_ANY_24_HEAD_MANIPULATION",
        "new_field_sampling_authorized":False,
        "can_generate_phylogenetic_memory_slope_from_one_jpn36_tip":False,
        "citation":"EAzami pinned "+REF+", original JPN36 protocol/validation; aza3 main SIE Phase A",
    }

def synthetic_checks(meta,validation,panel,bridge,lego):
    # Positive in the sense of the *source contract*, NOT experiment outcome.
    s=source_data(meta,validation,[],[],[],panel,bridge,lego)
    require(s["current_EAzami_pilot"]["aggregate_necessary_counted_flowering_individuals_for_full_24_head_pilot"]==1200,
            "PROTOCOL_POPULATION_QUOTA_ARITHMETIC")
    require(s["actual_JPN36_source_rows"]["eligible_heads"]==0,
            "ACTUAL_SOURCE_NOT_EMPTY")
    scenarios=[]
    for flowering in (0,200,1000,1199,1200,1400):
        eligible=min(24,math.floor(.02*flowering))
        scenarios.append({"hypothetical_flowering_individuals":flowering,
                          "mathematical_maximum_heads_under_two_percent_and_24_cap":eligible,
                          "can_fit_all_24_by_number_alone":eligible>=24,
                          "hypothetical_not_a_real_site_census":True,
                          "even_if_true_authorization_still_required":True})
    require([v["mathematical_maximum_heads_under_two_percent_and_24_cap"]
             for v in scenarios]==[0,4,20,23,24,24],
            "CONSERVATION_QUOTA_LIMITS_INCORRECT")
    # Bad protocol data cannot be rescued by the fact that no observations
    # have yet been entered.
    import copy
    for mk,needle in (
        (lambda v:v["scope"].update(current_field_execution_authorized=True),
         "FIELD_AUTHORIZATION_STATE_HAS_CHANGED_REAUDIT_BEFORE_EXECUTION"),
        (lambda v:v["rights_and_conservation_gate"].update(
            maximum_experimental_heads=25),
         "ORIGINAL_JPN36_PAIRED_EXPERIMENT_CONTRACT_CHANGED"),
        (lambda v:v["natural_reference_gate"].update(
            minimum_observation_only_individuals=0),
         "REFERENCE_COHORT_DENOMINATOR_CHANGED"),
    ):
        altered=copy.deepcopy(meta)
        mk(altered)
        try:source_data(altered,validation,[],[],[],panel,bridge,lego)
        except ValueError as e:
            require(needle in str(e),
                    "WRONG_ADVERSARIAL_CONSERVATION_TEST_ERROR:"+str(e))
        else:raise AssertionError("INVALID_EAZAMI_CONSERVATION_PLAN_ACCEPTED")
    return {"status":"PASS_FROZEN_JPN36_CONSERVATION_AND_GENOMIC_COHORT_SEPARATION",
            "hypothetical_census_sensitivity":scenarios,
            "negative_controls":3,
            "new_real_plant_observations":0}

def main():
    p=argparse.ArgumentParser()
    for arg in ("jpn36-contract","jpn36-validation",
                "jpn36-eligible","jpn36-assignment","jpn36-observations"):
        p.add_argument("--"+arg,type=Path,required=True)
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    meta=json.loads(a.jpn36_contract.read_text(encoding="utf-8"))
    validation=json.loads(a.jpn36_validation.read_text(encoding="utf-8"))
    field={
        "eligible":csv_rows(a.jpn36_eligible,EXPECTED_RECORDS["eligible"],"ELIGIBLE"),
        "assignment":csv_rows(a.jpn36_assignment,EXPECTED_RECORDS["assignment"],"ASSIGNMENT"),
        "observation":csv_rows(a.jpn36_observations,EXPECTED_RECORDS["observation"],"OBSERVATION"),
    }
    planned=csv_rows(PANEL,{"panel_id","phase","taxon","recommended_n"},"PANEL")
    b=json.loads(BRIDGE.read_text(encoding="utf-8"))
    lego=json.loads(LEGO.read_text(encoding="utf-8"))
    result={
        "version":"aza3_jpn36_three_question_bridge_feasibility_v1",
        "as_of":"2026-10-19",
        "original_EAzami_source_ref":REF,
        "source_integrity_tests":synthetic_checks(meta,validation,planned,b,lego),
        "real_pinned_source_reconciliation":source_data(
            meta,validation,field["eligible"],field["assignment"],
            field["observation"],planned,b,lego),
    }
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                         encoding="utf-8")
    r=result["real_pinned_source_reconciliation"]
    print(json.dumps({"status":r["status"],
        "planned_independent_field_plants":r["current_EAzami_pilot"]["randomized_independent_plants"],
        "minimum_counted_flowering_individuals_if_all_24_heads_used":r["current_EAzami_pilot"]["aggregate_necessary_counted_flowering_individuals_for_full_24_head_pilot"],
        "planned_independent_phase_A_genomic_individuals":r["prospective_aza3_molecular_design"]["separate_Phase_A_SIE_target_individuals"],
        "actually_enrolled_field_heads":r["actual_JPN36_source_rows"]["eligible_heads"],
        "authorization":r["current_EAzami_pilot"]["field_execution_authorized"]
    },ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
