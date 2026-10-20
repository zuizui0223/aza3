#!/usr/bin/env python3
"""Actual original-literature Cirsium sticky-bract payoff context audit.

Source archive:
- Willson, Anderson & Thomas 1983, American Midland Naturalist 110:212-214,
  DOI 10.2307/2425228 (published abstract accessible; full numeric not yet extracted)
- Thomas 2003, Great Lakes Entomologist DOI 10.22543/0090-0222.2085
  (original publisher abstract: neutralization did not improve seed predation,
   predators or seeds, in C. discolor 1982–1987)
- Thomas 2007, Great Lakes Entomologist DOI 10.22543/0090-0222.2188
  (original publisher abstract: pollinators avoid sticky, several seed
   predators bypass, ants/aphids deterred, predatory arthropods sometimes glean)

EAzami's frozen 15-row source seed included the 2003 NULL but NOT the 1983
positive/no-detected two-species manipulation. Do not change upstream EAzami.
This script returns real published *categorical effect status*, not invented
effect estimates, statistical heterogeneity, cross-species treatment effects,
or guild-conditioned adaptive selection.

No original field data are added. No manipulation is authorized.
"""
from __future__ import annotations
import argparse
import collections
import copy
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
NEW=ROOT/"data/evidence/cirsium_stickiness_real_primary_context_evidence_v1.csv"

NEEDED={
 "source_id","publication_year","experiment_years","species","location",
 "study_design","manipulated_head_module","intervention","endpoint",
 "reported_direction","source_access","doi","source_url",
 "quant_effect_estimate","independent_study_key","interpretation_ceiling"
}
UPSTREAM_NEEDED={"evidence_id","study_id","taxon","capitulum_module",
                 "design","fitness_endpoint","direction","doi"}
BOTANICAL_NEEDED={"paper_japan_member_id","paper_taxon_concept",
                  "orientation_state","phyllary_posture","stickiness_state",
                  "authority_match_status","source_scope"}
FROZEN={
 "WILLSON_1983_DISCOLOR":{
   "publication_year":"1983","experiment_years":"1980",
   "species":"Cirsium discolor",
   "reported_direction":"significant_increase_after_occlusion",
   "independent_study_key":"Willson_Anderson_Thomas_1983",
   "doi":"10.2307/2425228",
 },
 "WILLSON_1983_FLODMANI":{
   "publication_year":"1983","experiment_years":"1980_and_1982",
   "species":"Cirsium flodmani",
   "reported_direction":"no_effect_detected",
   "independent_study_key":"Willson_Anderson_Thomas_1983",
   "doi":"10.2307/2425228",
 },
 "THOMAS_2003_DISCOLOR":{
   "publication_year":"2003","experiment_years":"1982_to_1987",
   "species":"Cirsium discolor",
   "reported_direction":"no_effect_detected",
   "independent_study_key":"Thomas_2003",
   "doi":"10.22543/0090-0222.2085",
 },
 "THOMAS_2007_DISCOLOR":{
   "publication_year":"2007",
   "species":"Cirsium discolor",
   "reported_direction":"pollinators_avoid_sticky_surface_some_seed_predators_bypass_ants_and_aphids_deterred",
   "independent_study_key":"Thomas_2007",
   "doi":"10.22543/0090-0222.2188",
 },
}
CATEGORIES={
 "pollinators":{
   "role":"legitimate_flower_visitors",
   "observed_contact":"avoid_sticky_surfaces",
   "defence_outcome":"NOT_ASSESSED",
 },
 "some_seed_predators":{
   "role":"pre_dispersal_insect_seed_predators",
   "observed_contact":"bypass_or_avoid_sticky_surfaces",
   "defence_outcome":"NOT_ASSESSED_BY_THIS_BEHAVIORAL_STUDY",
 },
 "ants_and_aphids":{
   "role":"other_arthropods",
   "observed_contact":"deterred",
   "defence_outcome":"DOES_NOT_ESTABLISH_SEED_SAVING",
 },
 "predatory_arthropods":{
   "role":"predatory_bugs_and_spiders",
   "observed_contact":"some_appear_to_glean_from_sticky_traps",
   "defence_outcome":"CONCURRENT_GLEANING_NOT_CAUSAL_ENEMY_RESCUE",
 },
}

def require(ok,why):
    if not ok:raise ValueError(why)

def csvrows(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as stream:
        r=csv.DictReader(stream)
        return set(r.fieldnames or []),list(r)

def audit(rows,upstream):
    require(len(rows)==4,
            "FULL_FROZEN_FOUR_CIRSIUM_INTERACTION_SOURCE_COMPONENTS_REQUIRED")
    ids={}
    for row in rows:
        require(NEEDED.issubset(row),
                "STICKY_EVIDENCE_SCHEMA_CHANGED")
        identifier=row["source_id"]
        require(identifier in FROZEN and identifier not in ids,
                "STICKY_EVIDENCE_MISSING_DUPLICATE_OR_UNFROZEN_CASE")
        expected=FROZEN[identifier]
        for k,v in expected.items():
            require(row[k]==v,"ORIGINAL_PUBLISHED_SOURCE_DRIFT_"+k)
        require(row["manipulated_head_module"]=="bracteal_stickiness",
                "SUBSTITUTED_SPINES_FOR_EXUDATES")
        require(row["quant_effect_estimate"]=="NOT_EXTRACTED",
                "FIELD_EFFECT_SIZE_FABRICATED_FROM_ABSTRACT")
        require(row["source_url"].startswith("https://"),
                "PUBLIC_SOURCE_CITATION_MISSING")
        if row["study_design"]=="field_manipulation":
            require(identifier!="THOMAS_2007_DISCOLOR",
                    "CONTACT_OBSERVATION_WAS_PRESENTED_AS_EXPERIMENT")
        elif identifier!="THOMAS_2007_DISCOLOR":
            raise ValueError("TREATMENT_FALSIFIED_AS_NONEXPERIMENTAL")
        require(row["source_access"] in {
            "abstract_report_secondary_no_full_numeric_table",
            "primary_publisher_abstract"},
            "CLAIMED_FULL_TEXT_OR_NUMERIC_EXTRACT_UNVERIFIED")
        ids[identifier]=row
    require(set(ids)==set(FROZEN),"KNOWN_SOURCE_MISSING")
    old={x["evidence_id"]:x for x in upstream}
    require(len(upstream)==15 and len(old)==15,
            "ORIGINAL_EAZAMI_CIRSIUM_EVIDENCE_SEED_NOT_15_ROWS")
    require("INT006" in old and
            old["INT006"]["doi"]=="10.22543/0090-0222.2085" and
            old["INT006"]["taxon"]=="Cirsium discolor" and
            old["INT006"]["direction"]=="null",
            "ORIGINAL_EAZAMI_DIRECT_NULL_SOURCE_MISSING_OR_DRIFTED")
    require(not any(row["doi"]=="10.2307/2425228" for row in upstream),
            "ORIGINAL_EAZAMI_SEED_ALREADY_CONTAINS_1983_SOURCE")
    require(sum(x["capitulum_module"]=="stickiness" for x in upstream)==1,
            "ORIGINAL_EAZAMI_STICKY_FOCAL_SEED_COUNT_CHANGED")

    manip=[x for x in rows if x["study_design"]=="field_manipulation"]
    by_species=collections.defaultdict(list)
    for x in manip:
        by_species[x["species"]].append({
           "source_id":x["source_id"],
           "experimental_period":x["experiment_years"],
           "reported_direction":x["reported_direction"],
           "study_key":x["independent_study_key"],
           "location":x["location"],
        })
    discolor=by_species["Cirsium discolor"]
    require(len(discolor)==2 and
            {x["reported_direction"] for x in discolor}==
            {"significant_increase_after_occlusion","no_effect_detected"},
            "SAME_SPECIES_DIFFERENT_EXPERIMENTAL_CONTEXT_WAS_ERASED")
    require(len(by_species["Cirsium flodmani"])==1,
            "FLODMANI_NULL_CONTEXT_IGNORED")
    require(len({x["independent_study_key"] for x in manip})==2,
            "SPLIT_1983_SPECIES_ROWS_FALSELY_COUNTED_AS_3_STUDIES")
    assert all(x["doi"]==FROZEN[x["source_id"]]["doi"] for x in rows)
    return {
      "status":"REAL_CIRSIUM_SAME_SPECIES_STICKINESS_BENEFIT_DETECTION_HETEROGENEITY_FOUND",
      "original_EAzami_study_seed_rows":len(upstream),
      "original_EAzami_direct_stickiness_rows":1,
      "previously_missing_1983_study_identified":True,
      "original_manipulative_1983_and_2003_publications":2,
      "actual_original_manipulative_species_by_study_cases":3,
      "two_year_period_discordance_on_same_C_discolor":discolor,
      "second_C_flodmani_null":[*by_species["Cirsium flodmani"]],
      "observational_2007_real_insect_guild_surface_interactions":CATEGORIES,
      "original_manipulative_mean_effect_sizes_extracted":0,
      "has_statistically_demonstrated_sign_reversal":False,
      "has_shown_partner_identity_causes_difference":False,
      "has_measured_fitness_for_multiple_partner_guilds_under_same_spine_or_stickiness_state":False,
      "source_limitation":{
        "1983":"Published abstract/status accessible; original full numerical tables not independently extracted",
        "2003":"Original publisher abstract does not include standardized effect and confidence interval",
        "2007":"Original publisher abstract gives behavior by arthropod guild; not parallel manipulations or seed fitness",
        "comparison":"Site, year, treatment details and insect community not aligned; positive versus nondetection is NOT positive versus negative effect",
      },
      "scientific_Q2_decision":{
        "primary_rebuttal_to_universal_stickiness":
           "Same C. discolor exudate treatment showed detectable seed-predation deterrence in 1980, but no detected effects under the later 1982–1987 test.",
        "biological_interpretation":
           "Sticky surface can be bypassed by some seed feeders; pollinators avoid it; ants and aphids appear deterred. Partner contact is a plausible modifier but was NOT isolated across site-years.",
        "Q2_falsifiable_next_effect":
           "Head contact opportunity x sticky exudate manipulation x independently observed guild identity predicting actual oviposition and filled achenes, with original intervention randomized and shams",
        "do_not_claim":
           "Original 1983 and 2003 results constitute a same-study guild treatment factorial or a statistical interaction",
      },
    }

def source_sister_pair_audit(traits,bridge):
    names={}
    for r in traits:
        key=r.get("paper_japan_member_id")
        require(key not in names and bool(key),
                "DUPLICATE_CIRSIUM_BOTANICAL_TAXON")
        names[key]=r
    expected={
      "JPN_06":("Cirsium dipsacolepis","upward_or_erect",
                "nonsticky_or_nearly_nonsticky","detail_page"),
      "JPN_15":("Cirsium lineare","upward_or_erect",
                "sticky","index_catchphrase"),
      "JPN_36":("Cirsium sieboldii","downward_or_nodding",
                "nonsticky_or_nearly_nonsticky","detail_page"),
    }
    for key,vals in expected.items():
        require(key in names,"REQUIRED_CIRSIUM_HEAD_MORPH_MISSING")
        r=names[key]
        require((
          r["paper_taxon_concept"],r["orientation_state"],
          r["stickiness_state"],r["source_scope"]
        )==vals,"SOURCE_CIRSIUM_SISTER_HEAD_TRAIT_STATUS_DRIFT")
        require(r["authority_match_status"]=="exact_authority_concept_match",
                "NO_DIRECT_CIRSIUM_AUTHORITY_TAXON_CONCEPT")
    correction=bridge["state_direction_correction"]
    require(correction["JPN06_state"]=="nonsticky_or_nearly_nonsticky"
            and correction["JPN15_state"]=="sticky"
            and correction["JPN06_JPN15_sister_support"]==
            "100/100_on_canonical_ML_tree"
            and bridge["module_decisions"]["stickiness"][
                "best_ancestry_matched_pair"]==["JPN_06","JPN_15"],
            "AUTHORITATIVE_HEAD_STICKINESS_SISTER_SUPPORT_DRIFT")
    require(names["JPN_06"]["orientation_state"]==
            names["JPN_15"]["orientation_state"]=="upward_or_erect",
            "NUCLEAR_SISTER_PAIR_ORIENTATION_NO_LONGER_MATCHES")
    return {
      "status":"REAL_SISTER_CONCEPT_PAIR_HEAD_TRAITS_DECOUPLED_NOT_CAUSAL",
      "source":"EAzami frozen Japanese NMNS botanical states plus frozen canonical ML sister-pair summary",
      "JPN06":{"taxon":names["JPN_06"]["paper_taxon_concept"],
               "head_orientation":"upward_or_erect",
               "exudate_code":"nonsticky_or_nearly_nonsticky",
               "authority_source":"detail_page"},
      "JPN15":{"taxon":names["JPN_15"]["paper_taxon_concept"],
               "head_orientation":"upward_or_erect",
               "exudate_code":"sticky",
               "authority_source":"index_catchphrase",
               "field_stickiness_and_exudate_composition_measured":False},
      "candidate_phylogeny":"100/100_on_canonical_ML_tree",
      "independent_species_pair_not_two_evolutionary_transition_replications":True,
      "shared_orientation_but_different_stickiness":True,
      "no_partner_identified_or_head_reproductive_trait_causal_effect_estimated":True,
      "first_required_observation":"Verify actual live secretion in JPN15 and JPN06; track pre-entry insect identity and whether the animal contacts the true bractal exudate; source record JPN15 sticky is from index not a measured plant.",
      "proposed_manipulation_status":"NOT_AUTHORIZED",
    }

def corrupt_test(rows,upstream):
    def fail(bad,needle):
        try:audit(bad,upstream)
        except ValueError as er:
            require(needle in str(er),
                    "NEGATIVE_FIXTURE_ERROR_WRONG:"+str(er))
            return
        raise AssertionError("BAD_SOURCE_CIRSIUM_CONTEXT_TREATED_AS_REAL")
    z=copy.deepcopy(rows)
    z[0]["reported_direction"]="no_effect_detected"
    fail(z,"ORIGINAL_PUBLISHED_SOURCE_DRIFT_reported_direction")
    z=copy.deepcopy(rows)
    z[1]["species"]="Cirsium discolor"
    fail(z,"ORIGINAL_PUBLISHED_SOURCE_DRIFT_species")
    z=copy.deepcopy(rows)
    z[2]["quant_effect_estimate"]="0.0"
    fail(z,"FIELD_EFFECT_SIZE_FABRICATED_FROM_ABSTRACT")
    z=copy.deepcopy(rows)
    z[0]["manipulated_head_module"]="involucral_spines"
    fail(z,"SUBSTITUTED_SPINES_FOR_EXUDATES")
    z=copy.deepcopy(rows)
    z[3]["study_design"]="field_manipulation"
    fail(z,"CONTACT_OBSERVATION_WAS_PRESENTED_AS_EXPERIMENT")
    z=copy.deepcopy(rows)
    z[0]["independent_study_key"]="Thomas_2003"
    fail(z,"ORIGINAL_PUBLISHED_SOURCE_DRIFT_independent_study_key")
    return {
      "status":"PASS_SYNTHETIC_SOURCE_CONTEXT_AND_NO_FALSE_GUILD_SELECTION_GATES",
      "source_corruptions_rejected":6,
      "real_raw_insect_or_seed_records_created":0,
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--new-evidence",type=Path,default=NEW)
    a.add_argument("--original-eazami-seed",type=Path,required=True)
    a.add_argument("--original-japan-head-states",type=Path,required=True)
    a.add_argument("--original-module-function-bridge",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    arg=a.parse_args()
    hdr,rows=csvrows(arg.new_evidence)
    original_head,old=csvrows(arg.original_eazami_seed)
    head_header,japan=csvrows(arg.original_japan_head_states)
    bridge=json.loads(arg.original_module_function_bridge.read_text(encoding="utf-8"))
    require(NEEDED.issubset(hdr) and UPSTREAM_NEEDED.issubset(original_head)
            and BOTANICAL_NEEDED.issubset(head_header),
            "FROZEN_SOURCE_CSV_HEADERS_MISSING")
    live_sister_source=source_sister_pair_audit(japan,bridge)
    report={
      "analysis_id":"cirsium_real_stickiness_context_pr50_Q2_v1",
      "as_of":"2026-10-19",
      "sources":"Original publisher abstract 2003/2007; original 1983 abstract and DOI metadata",
      "EAzami_source_commit":"c9870e91604ed1514344ddc8f765eae5237dbf27",
      "synthetic_corruption_tests":corrupt_test(rows,old),
      "actual_published_evidence_audit":audit(rows,old),
      "actual_source_Cirsium_JPN06_JPN15_sister_trait_contrast":live_sister_source,
      "no_new_Cirsium_experiment_performed":True,
    }
    arg.out.parent.mkdir(parents=True,exist_ok=True)
    arg.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",
                       encoding="utf-8")
    q=report["actual_published_evidence_audit"]
    print(json.dumps({
      "result":q["status"],
      "actual_independent_manipulation_studies":q["original_manipulative_1983_and_2003_publications"],
      "study_by_species_experimental_cases":q["actual_original_manipulative_species_by_study_cases"],
      "C_discolor_effect_detection_across_periods":[x["reported_direction"] for x in q["two_year_period_discordance_on_same_C_discolor"]],
      "partner_dependent_fitness_effect_measured":q["has_measured_fitness_for_multiple_partner_guilds_under_same_spine_or_stickiness_state"],
      "Japanese_sister_pair_trait_dissociation":live_sister_source["shared_orientation_but_different_stickiness"],
      "source_negative_controls":report["synthetic_corruption_tests"]["status"],
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
