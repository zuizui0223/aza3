#!/usr/bin/env python3
"""Validate observational event-level Cirsium head access ontology.

This cannot estimate effects without observations, denominators and viable
achenes. It guards ecological semantics: contact != pollination; egg-laying
!= seed damage; parasitoid searching != seed rescue. No actual synthetic
test episode is committed as biological data.
"""
from __future__ import annotations
import argparse
import csv
import io
import json
import math
from pathlib import Path

BIN = {"0","1","NA"}
BIN_COLS=[
 "approached","landed_or_hover_access","reproductive_zone_reached",
 "anther_stigma_contact_observed","oviposition_confirmed",
 "predator_or_parasitoid_search_observed","host_attacked_confirmed",
 "host_killed_before_seed_damage","seed_damage_timing_assessed",
 "seeds_saved_inferred","scorer_blinded_to_head_traits"
]
ID_COLS=[
 "record_id","individual_id","population_id","capitulum_id",
 "observation_bout_id","approach_episode_id","video_clip_id",
 "visitor_taxon_label","behavioral_role","role_evidence",
 "pre_entry_guild","pre_entry_guild_evidence",
 "phenological_stage","entry_route_world","entry_route_head"
]


def header(path):
    with open(path,encoding="utf-8-sig",newline="") as f:
        return next(csv.reader(f))


def validate(rows, contract, field_header, bout_header, plant_header):
    required=contract["required_fields"]
    assert len(required)==len(set(required))
    field_needed=set(contract["source_ledgers"][1]["key"])
    bout_needed=set(contract["source_ledgers"][0]["key"])
    plant_needed=set(contract["source_ledgers"][2]["key"])
    if not field_needed.issubset(set(field_header)):
        raise ValueError("FIELD_LEDGER_SCHEMA_MISMATCH")
    if not bout_needed.issubset(set(bout_header)):
        raise ValueError("BOUT_LEDGER_SCHEMA_MISMATCH")
    if not plant_needed.issubset(set(plant_header)):
        raise ValueError("PLANT_CENSUS_SCHEMA_MISMATCH")
    if not set(contract["head_traits_from_existing_field_ledger"]).issubset(
         set(field_header)):
        raise ValueError("HEAD_MORPHOLOGY_LEDGER_COLUMNS_MISSING")
    key_seen=set()
    observed_bouts=set()
    role_counts={}
    unblinded=0
    reviewed=0
    for i, row in enumerate(rows,1):
        for k in required:
            if k not in row:
                raise ValueError(f"ROW_{i}_FIELD_MISSING:{k}")
        for k in ID_COLS:
            if not str(row.get(k) or "").strip():
                raise ValueError(f"ROW_{i}_MISSING_IDENTIFIER:{k}")
        for label, options in contract["allowed"].items():
            if not isinstance(options,list):
                continue
            if label in row and label != "binary":
                if row[label] not in options:
                    raise ValueError(f"ROW_{i}_INVALID_{label}:{row[label]}")
        if row["pollen_deposition_assay"] not in contract["allowed"]["pollen_deposition_assay"]:
            raise ValueError(f"ROW_{i}_INVALID_POLLEN_ASSAY")
        for k in BIN_COLS:
            if row[k] not in BIN:
                raise ValueError(f"ROW_{i}_INVALID_BOOLEAN:{k}")
        if row["approached"]!="1":
            raise ValueError(f"ROW_{i}_APPROACH_NOT_DEFINED")
        if row["pre_entry_guild"]=="unresolved":
            if row["pre_entry_guild_evidence"]!="unresolved":
                raise ValueError(f"ROW_{i}_PRE_ENTRY_GUILD_WITHOUT_VALID_ASSIGNMENT")
        elif (row["pre_entry_guild_evidence"]=="unresolved" or
              not row["evidence_uri"].strip()):
            raise ValueError(f"ROW_{i}_POST_OUTCOME_ROLE_CANNOT_DEFINE_ARRIVAL_GUILD")
        if row["seeds_saved_inferred"]=="1":
            raise ValueError(f"ROW_{i}_OBSERVATION_CANNOT_CERTIFY_SEED_RESCUE")
        key=(row["individual_id"],row["population_id"],row["capitulum_id"],
             row["observation_bout_id"],row["approach_episode_id"])
        if key in key_seen: raise ValueError(f"ROW_{i}_REPEATED_APPROACH_ID")
        key_seen.add(key)
        observed_bouts.add(key[:4])
        try:
            secs=float(row["time_from_bout_start_s"])
        except (ValueError,TypeError):
            raise ValueError(f"ROW_{i}_INVALID_TIME")
        if not math.isfinite(secs) or secs<0:
            raise ValueError(f"ROW_{i}_INVALID_TIME")
        if row["reproductive_zone_reached"]=="1" and row["landed_or_hover_access"]!="1":
            raise ValueError(f"ROW_{i}_ZONE_WITHOUT_ENTRY")
        if row["anther_stigma_contact_observed"]=="1" and row["reproductive_zone_reached"]!="1":
            raise ValueError(f"ROW_{i}_CONTACT_WITHOUT_ZONE")
        if row["pollen_deposition_assay"]=="positive" and row["anther_stigma_contact_observed"]!="1":
            raise ValueError(f"ROW_{i}_POLLEN_ASSAY_NO_CONTACT")
        if row["oviposition_confirmed"]=="1" and row["behavioral_role"]!="ovipositing_seed_feeder":
            raise ValueError(f"ROW_{i}_OVIPOSITION_BUT_NO_VERIFIED_GUILD")
        if (row["oviposition_confirmed"]=="1" or
            row["anther_stigma_contact_observed"]=="1" or
            row["host_attacked_confirmed"]=="1") and row["role_evidence"]=="not_assessed":
            raise ValueError(f"ROW_{i}_BIOLOGICAL_EFFECT_WITHOUT_EVIDENCE")
        if row["predator_or_parasitoid_search_observed"]=="1" and row["behavioral_role"] not in ("enemy_parasitoid","enemy_predator"):
            raise ValueError(f"ROW_{i}_ENEMY_SEARCH_WITHOUT_ENEMY_IDENTITY")
        if row["host_attacked_confirmed"]=="1" and (
            row["predator_or_parasitoid_search_observed"]!="1" or
            not row["linked_host_id"].strip()):
            raise ValueError(f"ROW_{i}_HOST_ATTACK_WITHOUT_LINK")
        if row["host_killed_before_seed_damage"]=="1" and (
            row["host_attacked_confirmed"]!="1" or
            row["seed_damage_timing_assessed"]!="1" or
            not row["linked_host_id"].strip()):
            raise ValueError(f"ROW_{i}_SEED_RESCUE_TIMING_UNRESOLVED")
        if row["host_killed_before_seed_damage"]=="0" and row["seed_damage_timing_assessed"]!="1":
            raise ValueError(f"ROW_{i}_NEGATIVE_EARLY_KILL_REQUIRES_ASSESSED_TIMING")
        if row["role_evidence"] != "not_assessed":
            if not row["evidence_uri"].strip():
                raise ValueError(f"ROW_{i}_CLAIM_WITHOUT_VIDEO_OR_ASSAY")
            reviewed+=1
        if row["scorer_blinded_to_head_traits"]!="1": unblinded+=1
        role_counts[row["behavioral_role"]]=role_counts.get(row["behavioral_role"],0)+1
    return {
       "status":"PASS_STRUCTURAL_VALIDATION_ZERO_FITNESS_INFERENCE",
       "n_approach_episodes":len(rows),
       "n_observed_bouts_with_episodes":len(observed_bouts),
       "role_counts":role_counts,
       "n_event_rows_with_evidence":reviewed,
       "n_unblinded_or_unknown_head_trait_scorers":unblinded,
       "cannot_estimate_zero_visits_or_arrival_rate_from_event_table_alone":True,
       "requires_denominator_bout_effort_and_undetected_event_recall":True,
       "direct_trait_impact_on_seed_fitness":"NOT_IDENTIFIED",
       "claim_ceiling":"SCHEMA_AND_SEMANTIC_VERIFICATION_ONLY_NOT_REAL_ECOLOGICAL_RESULT"
    }


def test_fixture(contract):
    valid={
      k:"" for k in contract["required_fields"]
    }
    valid.update({
      "record_id":"TEST_ROW_1","individual_id":"TEST_PLANT",
      "population_id":"TEST_POP","capitulum_id":"TEST_HEAD",
      "observation_bout_id":"TEST_BOUT","approach_episode_id":"TEST_APPROACH_1",
      "video_clip_id":"test/synthetic.mp4","time_from_bout_start_s":"5",
      "phenological_stage":"full_anthesis","reproductive_sex_state":"hermaphroditic",
      "visitor_taxon_label":"unidentified insect",
      "pre_entry_guild":"unresolved","pre_entry_guild_evidence":"unresolved",
      "behavioral_role":"floral_forager","role_evidence":"video_behavior_confirmed",
      "entry_route_world":"side","entry_route_head":"disc_facing",
      "approached":"1","landed_or_hover_access":"1","reproductive_zone_reached":"1",
      "anther_stigma_contact_observed":"1","pollen_deposition_assay":"not_assessed",
      "oviposition_confirmed":"0","predator_or_parasitoid_search_observed":"0",
      "host_attacked_confirmed":"0","host_killed_before_seed_damage":"NA",
      "seed_damage_timing_assessed":"0","seeds_saved_inferred":"NA",
      "evidence_uri":"test/synthetic.mp4","scorer_blinded_to_head_traits":"1",
    })
    bout=["individual_id","population_id","capitulum_id","observation_bout_id",
          "observation_minutes","pollinator_visit_count","effective_contact_count"]
    field=["individual_id","population_id","capitulum_id",
           *contract["head_traits_from_existing_field_ledger"]]
    plant=["individual_id","population_id","phenology_census_id"]
    outcome=validate([valid],contract,field,bout,plant)
    assert outcome["n_approach_episodes"]==1
    invalid_cases=[
      ("REJECT_SEED_RESCUE_FROM_EVENT",
       {"seeds_saved_inferred":"1"}, "OBSERVATION_CANNOT_CERTIFY_SEED_RESCUE"),
      ("REJECT_CONTACT_IS_POLLEN_DEPOSITION",
       {"anther_stigma_contact_observed":"0","pollen_deposition_assay":"positive"},
        "POLLEN_ASSAY_NO_CONTACT"),
      ("REJECT_OVIPOSITION_AS_VISIT",
       {"oviposition_confirmed":"1","behavioral_role":"floral_forager"},
        "OVIPOSITION_BUT_NO_VERIFIED_GUILD"),
      ("REJECT_FALSE_EARLY_PARASITOID_BENEFIT",
       {"host_killed_before_seed_damage":"1","linked_host_id":"",
        "host_attacked_confirmed":"0","seed_damage_timing_assessed":"0"},
        "SEED_RESCUE_TIMING_UNRESOLVED"),
      ("REJECT_FALSE_NO_DAMAGE_AS_EARLY_KILL",
       {"host_killed_before_seed_damage":"0","seed_damage_timing_assessed":"0"},
        "NEGATIVE_EARLY_KILL_REQUIRES_ASSESSED_TIMING"),
      ("REJECT_NOT_PRECLASSIFIED_SEED_FEEDER",
       {"pre_entry_guild":"seed_feeder_candidate",
        "pre_entry_guild_evidence":"unresolved"},
       "POST_OUTCOME_ROLE_CANNOT_DEFINE_ARRIVAL_GUILD"),
      ("REJECT_PSEUDO_PRE_ENTRY_ASSIGNMENT",
       {"pre_entry_guild":"unresolved",
        "pre_entry_guild_evidence":"pre_entry_taxon_and_local_life_history_verified"},
       "PRE_ENTRY_GUILD_WITHOUT_VALID_ASSIGNMENT"),
      ("REJECT_NO_APPROACH",
       {"approached":"0"},"APPROACH_NOT_DEFINED"),
    ]
    for label,mutation,error in invalid_cases:
        current={**valid,**mutation}
        try:
            validate([current],contract,field,bout,plant)
        except ValueError as ex:
            assert error in str(ex),(label,str(ex))
        else:
            raise AssertionError(f"{label}: false biological claim admitted")
    return {
       "n_pass_cases":1, "n_expected_failure_checks":len(invalid_cases),
       "status":"PASS_SYNTHETIC_SEMANTIC_FIREWALL_NOT_FIELD_DATA"
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--contract",type=Path,required=True)
    parser.add_argument("--events",type=Path,required=True)
    parser.add_argument("--bout-header",type=Path)
    parser.add_argument("--field-header",type=Path)
    parser.add_argument("--plant-header",type=Path)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    contract=json.loads(args.contract.read_text(encoding="utf-8"))
    if header(args.events)!=contract["required_fields"]:
        raise SystemExit("ERROR_EVENT_HEADER_CONTRACT_DRIFT")
    with args.events.open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    bout=header(args.bout_header) if args.bout_header else [
        "individual_id","population_id","capitulum_id","observation_bout_id"]
    field=header(args.field_header) if args.field_header else [
        "individual_id","population_id","capitulum_id",
        *contract["head_traits_from_existing_field_ledger"]]
    plant=header(args.plant_header) if args.plant_header else [
        "individual_id","population_id","phenology_census_id"]
    out=validate(rows,contract,field,bout,plant)
    if args.self_test:
        out["synthetic_tests"]=test_fixture(contract)
    print(json.dumps(out,ensure_ascii=False,indent=2))
if __name__=="__main__":
    main()
