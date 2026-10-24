#!/usr/bin/env python3
"""Q2 original STUDIES (not made-up individuals): when is sticky head defence paid off?

Source papers:
- Willson, Anderson & Thomas 1983 AMN DOI 10.2307/2425228;
  source abstract excerpts reused in scholarship review; 1980 experiment
  C. discolor reports increased damage after sticky exudate occlusion,
  C. flodmanii reported no detectable effect in 1980 and 1982.
- Thomas 2003 Great Lakes Entomol DOI 10.22543/0090-0222.2085;
  publisher original abstract, reports failure to support increased seed
  predators and decreased seed production under exudate neutralization
  during study of C. discolor + P. digitalis in 1982–1987.
- Thomas 2007 DOI 10.22543/0090-0222.2188; original journal abstract
  documents pollinator NONcontact, some seed feeder successful bypass,
  ant/aphid deterrence, and prey gleaning arthropods on C. discolor.
- Independent externally-shaped comparison: Sun & Huang 2015 Pedicularis
  rex DOI 10.1093/aobpla/plv019; water drainage via scissor cut,
  significant seed damage and final seed change, and no detected visit
  change in just one visitor population. This is NOT thistle adhesion.

This is exact publication/cohort evidence AUDIT, not reanalysis of
individual head outcomes. No raw sample means/n/p-values supplied for
1983/2003; they are NULL/UNKNOWN, not zeros. A reported nonsignificant
treatment difference is not equivalence and is not a sign reversal from
a significant positive effect. Study environment, year, insect community
and exudate neutralization procedure are not controlled BETWEEN papers.
No trait adaptation, heritable independence, physical head angle or
restored evolutionary developmental pathway is concluded.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

SOURCES={
    "WAT1983":"10.2307/2425228",
    "T2003":"10.22543/0090-0222.2085",
    "T2007":"10.22543/0090-0222.2188",
}
COHORT_OUTCOMES={
    "WAT1983_discolor_Illinois_1980":
      "GREATER_SEED_PREDATION_WHEN_STICKY_EXUDATE_OCCLUDED",
    "WAT1983_flodmani_Minnesota_1980":"NO_EFFECT_DETECTED",
    "WAT1983_flodmani_Minnesota_1982":"NO_EFFECT_DETECTED",
    "T2003_discolor_Illinois_1982_1987":
      "NO_SUPPORT_FOR_GREATER_SEED_PREDATORS_OR_LOWER_SEED_PRODUCTION_WHEN_NEUTRALIZED",
}
ORIGINAL_ROLES={
    "reported_pollinators":"DO_NOT_CONTACT_STICKY_TRAPS",
    "specialist_seed_consumers":"SOME_SEED_PREDATORS_BYPASS_STICKY_TRAPS",
    "small_prey_gleaning_arthropod_predators":
        "MINUTE_PIRATE_BUG_AND_SMALL_SALTICID_SPIDER_SEEM_TO_GLEAN_TRAPPED_PREY",
    "ants_and_aphids":"REPORTED_DETERRED_BY_STICKY_TRAPS",
}

def check(statement,error):
    if not statement:
        raise ValueError(error)

def audit(registry,comparator):
    check(registry["version"]==
          "aza3_real_cirsium_sticky_exudate_state_dependent_function_original_study_registry_v1",
          "REGISTRY_VERSION_MISMATCH")
    studies=registry["primary_sources"]
    check({x["source_id"] for x in studies}==set(SOURCES)
          and len(studies)==3,
          "SOURCE_STUDIES_NOT_THREE_WITH_UNIQUE_IDS")
    for s in studies:
        check(s["doi"]==SOURCES[s["source_id"]],
              "STICKY_STUDY_PRIMARY_DOI_DRIFT")
        check(s["source_grade"].endswith("NOT_FULL_PDF_REVERIFIED"),
              "FALSE_CLAIM_FULL_PRIMARY_PDF_AUTHENTICATED")
        check(s.get("access_urls"),
              "UNATTRIBUTED_STICKY_SOURCE_FINDING")
    byid={s["source_id"]:s for s in studies}
    # Explicitly preserve original STUDY vs species/year COHORT grain.
    cohorts=[(s["source_id"],x) for s in studies
             for x in s.get("organismal_samples",[])]
    check(len(cohorts)==4,
          "SOURCE_TRIAL_COHORT_NOT_FOUR")
    got={x["cohort_id"]:x["authors_reported_response"] for sid,x in cohorts}
    check(got==COHORT_OUTCOMES,
          "PUBLISHED_STICKY_SEED_OUTCOME_CODE_CHANGED")
    for sid,x in cohorts:
        check(x["source_p_value"] is None and
              x["source_mean_seed_loss"] is None and
              x["source_individual_head_records_available"] is False,
              "UNVERIFIED_STICKY_NUMERIC_EFFECT_OR_HEAD_RECORD_FABRICATED")
        check(x["experimental_treatment"] in (
                "original sticky bract exudates experimentally occluded",
                "sticky inflorescence exudates neutralized"),
              "INCORRECT_STICKINESS_INTERVENTION")
        check(x["species"] in ("Cirsium discolor","Cirsium flodmanii"),
              "UNVERIFIED_SOURCE_CIRSIUM_TAXON")
    discolor=[x for _,x in cohorts if x["species"]=="Cirsium discolor"]
    check(len(discolor)==2 and
          discolor[0]["field_years"]==[1980] and
          discolor[1]["field_years"]==[1982,1983,1984,1985,1986,1987],
          "FALSE_SAME_YEAR_DISCOLOR_REPLICATE_OR_COMPLETE_EVERY_YEAR")
    check(discolor[1]["source_result_reported_jointly_across_plant_species"] is True,
          "FALSE_2003_DISCOLOR_EXCLUSIVE_ESTIMAND")
    check(cohorts[0][0]=="WAT1983" and cohorts[-1][0]=="T2003",
          "CHRONOLOGY_CONTRAST_CHANGED")
    insects=byid["T2007"]["observed_interfaces"]
    check(len(insects)==4 and
          {x["role"]:x["original_observation"] for x in insects}==ORIGINAL_ROLES,
          "PUBLISHED_ARTHROPOD_INTERFACE_ROLES_CHANGED")
    check(not byid["T2007"]["direct_original_individual_event_ids_available"],
          "FALSE_STICKY_SURFACE_INSECT_CONTACT_EVENTS")
    # Compare the published Pedicularis water defence experimental result
    # from a previously audited aza3 source, without copying or pretending
    # the traits/organisms or treatments are homologous.
    cmp=registry["other_genus_original_comparator"]
    p=comparator["sources"]["drainage_2015"]
    check(cmp["doi"]==p["doi"]=="10.1093/aobpla/plv019",
          "WRONG_ORIGINAL_PEDICULARIS_CONTROL")
    check(cmp["distinct_from_sticky_exudate_mechanism"] is True and
          p["seed_populations"]==6 and
          p["pollinator_video_populations"]==1 and
          p["conservative_BH_significant_sites"]==4,
          "PEDICULARIS_SOURCE_NOT_DISTINCT_OR_INCORRECT_GRAIN")
    check(cmp["source_population_tests_nominal_significant"]==5 and
          cmp["new_aza3_six_test_BH_conservative_significant"]==4 and
          cmp["source_single_population_pollinator_visitation_p"]==.958 and
          cmp["source_initial_seed_set_across_six_populations_p"]==.906,
          "INCORRECT_POSITIVE_EXTERNAL_FUNCTION_COMPARATOR")
    evidence=registry["biology_contrast"]
    check(evidence["study_year_x_method_x_arthropod_community_confounding_unresolved"] is True
          and evidence["zero_effect_statistical_equivalence_not_tested"] is True
          and evidence["source_study_environment_and_treatment_comparable"] is False
          and evidence["no_meta_effect_size_without_original_group_means_and_unit_counts"] is True,
          "CAUSAL_TEMPORAL_REGIME_SHIFT_OR_META_EFFECT_FALSELY_IDENTIFIED")
    need=registry["provenance_caveats"]
    check(all(need[k] for k in (
          "source_1983_full_original_pdf_verified",
          "source_2003_full_original_pdf_verified",
          "source_2007_full_original_pdf_verified")) is False
          and all(not need[k] for k in (
          "source_1983_full_original_pdf_verified",
          "source_2003_full_original_pdf_verified",
          "source_2007_full_original_pdf_verified")),
          "FALSE_FULL_TEXT_SOURCE_VERIFICATION")
    check(evidence["no_new_own_field_observations"]==0,
          "FALSE_NEW_FIELD_WORK")
    # Positive/null SUPPORT counts only. These are not fixed effects,
    # empirical proportions of all thistles, or independent origin counts.
    seed_classes=Counter(x["report_significance"] for _,x in cohorts)
    study_id=sorted({sid for sid,_ in cohorts})
    review={
        "version":"cirsium_sticky_head_access_payoff_regime_scientific_gate_v1",
        "provenance":{
            "primary_publications":len(studies),
            "original_source_organism_cohorts":len(cohorts),
            "source_publication_dois":SOURCES,
            "original_reports_not_raw_observation_matrix":True,
            "1980_discolor_and_1982_1987_discolor_are_not_same_field_cohort":True,
            "C_fodmani_1980_and_1982_are_not_two_independent_evolutionary_origins":True,
        },
        "original_seed_predation_support":{
            "Cirsium discolor 1980":"DETECTED_GREATER_PREDATION_WHEN_STICKINESS_BLOCKED",
            "Cirsium discolor 1982-1987":"NO_SUPPORT_FOR_SAME_EFFECT_IN_ORIGINAL_ABSTRACT",
            "Cirsium flodmanii 1980":"NO_DETECTED_EFFECT",
            "Cirsium flodmanii 1982":"NO_DETECTED_EFFECT",
            "significant_source_reports":seed_classes["REPORTED_STATISTICALLY_SIGNIFICANT"],
            "reported_nondetections_not_null_effect_size":sum(
                x["report_significance"]!="REPORTED_STATISTICALLY_SIGNIFICANT"
                for _,x in cohorts),
            "pooled_experimental_effect_size":None,
            "confidence_interval_for_sticky_defence_benefit":None,
            "sign_reversal_demonstrated":False,
        },
        "independent_observed_visitor_access_2007":{
            x["role"]:x["original_observation"] for x in insects
        },
        "positive_non_Cirsium_comparator":{
            "organism":"Pedicularis rex",
            "mechanism":"cupular bract water-retention + physical hole",
            "source_populations_seed_outcome":p["seed_populations"],
            "source_populations_visitor_frequency":p["pollinator_video_populations"],
            "site_specific_predation_significant_after_conservative_BH":
                p["conservative_BH_significant_sites"],
            "not_homologous_to_sticky_Cirsium_head_intervention":True,
        },
        "biological_inference":{
            "same_sticky_structure_yields_different_published_detection_categories":True,
            "sticky_exudate_protective_net_seed_fitness_in_all_environments_identified":False,
            "universal_functional_trait_payoff_transfers_to_other_Cirsium_lineages":False,
            "trait_by_true_insect_contact_measured_in_same_head_to_fitness":False,
            "pollinator_noncontact_not_zero_pollen_service":True,
            "no_true_spines_capitulum_gravity_angle_measured":True,
            "no_ancestral_regain_or_genetic_LEGO":True,
            "new_original_organismal_field_rows":0,
        },
    }
    check(review["original_seed_predation_support"]["significant_source_reports"]==1 and
          review["original_seed_predation_support"]
                ["reported_nondetections_not_null_effect_size"]==3,
          "PUBLISHED_DISCLOSURE_PATTERN_CHANGED")
    return review

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--pedicularis-comparator",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    data=json.loads(a.source.read_text(encoding="utf-8"))
    comp=json.loads(a.pedicularis_comparator.read_text(encoding="utf-8"))
    report=audit(data,comp)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,ensure_ascii=False,
                                allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
