#!/usr/bin/env python3
"""Fail-closed, real-anatomy head contact ontology for Cirsium access.

Only a verified insect touching actual phyllary/spine tissue AND verified
attempt-level physical blockage can be called a physical barrier. Candidate
seed feeders and pollinators are independently classified before their
approach outcomes. Repeated attempt rows do not create independent heads.

No field events have been committed; tests are explicitly SYNTHETIC.
"""
from __future__ import annotations
import argparse
import copy
import csv
import json
import math
from pathlib import Path
from validate_capitulum_video_effort_denominator_v1 import (
    validate as validate_video_effort, read_csv as read_effort_csv
)

COLUMNS=[
 "record_id","individual_id","population_id","capitulum_id",
 "observation_bout_id","approach_episode_id","attempt_index",
 "video_clip_id","head_stage","contact_surface","contact_resolved",
 "spine_touch_verified","physical_blockage_verified",
 "entry_not_achieved_at_attempt","actual_spine_length_mm",
 "minimum_phyllary_gap_mm","anatomy_scale_evidence_uri",
 "video_contact_evidence_uri","pre_entry_taxon_evidence_uri",
 "scorer_blinded_to_oviposition_result"
]
KEYS=("individual_id","population_id","capitulum_id",
      "observation_bout_id","approach_episode_id")
ATTEMPT_KEY=KEYS+("attempt_index",)
SURFACES={
 "spine_tip","spine_shaft","phyllary_lamina","phyllary_gap",
 "floret_disc","stem","no_surface_contact","unresolved"
}
ARMATURE_SURFACES={"spine_tip","spine_shaft","phyllary_lamina","phyllary_gap"}
VERIFIED_GUILDS={
 "pre_entry_taxon_and_local_life_history_verified",
 "independently_assigned_before_access"
}
FLAGS={"0","1","NA"}
PARENT_ZONE_ALLOWED={
    "involucre_outer":{"spine_tip","spine_shaft","phyllary_lamina",
                       "no_surface_contact","unresolved"},
    "involucre_gap":{"spine_tip","spine_shaft","phyllary_lamina",
                     "phyllary_gap","no_surface_contact","unresolved"},
    "floret_disc":{"floret_disc","no_surface_contact","unresolved"},
    "receptacle_base":{"phyllary_lamina","no_surface_contact","unresolved"},
    "stem_entry":{"stem","no_surface_contact","unresolved"},
    "undetermined":SURFACES,
}

def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        return list(reader.fieldnames or []),list(reader)

def validate(rows, attempts, events):
    at_index={}
    for a in attempts:
        key=tuple(a.get(k,"") for k in ATTEMPT_KEY)
        if not all(key) or key in at_index:
            raise ValueError("DUPLICATE_OR_MISSING_PARENT_CONTACT_ATTEMPT")
        at_index[key]=a
    ev_index={}
    for ev in events:
        key=tuple(ev.get(k,"") for k in KEYS)
        if not all(key) or key in ev_index:
            raise ValueError("DUPLICATE_OR_MISSING_PARENT_APPROACH_EPISODE")
        ev_index[key]=ev
    expected=set()
    parent_guild_unresolved=0
    parent_coverage_incomplete=0
    for key,a in at_index.items():
        episode=ev_index.get(key[:len(KEYS)])
        if episode is None:
            raise ValueError("PARENT_ATTEMPT_WITHOUT_APPROACH_EPISODE")
        coverage=a.get("continuous_coverage")
        if coverage not in {"0","1"}:
            raise ValueError("PARENT_ATTEMPT_COVERAGE_INVALID")
        if coverage=="0":
            parent_coverage_incomplete+=1
            continue
        if (episode.get("pre_entry_guild")=="unresolved" or
            episode.get("pre_entry_guild_evidence") not in VERIFIED_GUILDS):
            parent_guild_unresolved+=1
            continue
        expected.add(key)
    if not rows and not attempts:
        return {"status":("HOLD_EPISODES_EXIST_WITHOUT_ATTEMPT_SEQUENCE" if events
                          else "NO_REAL_TRUE_ARMATURE_CONTACT_ANNOTATIONS"),
                "n_parent_attempts":0,"n_eligible_parent_attempts":0,
                "n_contact_rows":0, "n_verified_spine_touches":None,
                "n_verified_physical_blockages":None,
                "n_independent_heads":0,
                "empirical_inference":"NOT_IDENTIFIABLE"}
    if not expected and not rows:
        return {"status":"HOLD_NO_ELIGIBLE_PARENT_ATTEMPTS",
                "n_parent_attempts":len(attempts),
                "n_eligible_parent_attempts":0,
                "n_parent_attempts_unresolved_guild":parent_guild_unresolved,
                "n_parent_attempts_unscorable_video":parent_coverage_incomplete,
                "n_contact_rows":0, "n_verified_spine_touches":None,
                "n_verified_physical_blockages":None,
                "n_independent_heads":0,
                "empirical_inference":"NOT_IDENTIFIABLE"}
    counts={"touches":0,"blockages":0,"unknown":0}
    seen=set()
    heads=set()
    for i,row in enumerate(rows,1):
        if set(row)!=set(COLUMNS):
            raise ValueError("ANNOTATION_HEADER_OR_COLUMN_DRIFT")
        key=tuple(row[k].strip() for k in ATTEMPT_KEY)
        if not all(key) or key in seen:
            raise ValueError("ANNOTATION_EMPTY_OR_DUPLICATE_ATTEMPT")
        seen.add(key)
        a=at_index.get(key)
        e=ev_index.get(key[:len(KEYS)])
        if a is None or e is None:
            raise ValueError("ANNOTATION_MUST_MATCH_PARENT_ATTEMPT_AND_EPISODE")
        if row["video_clip_id"]!=a["video_clip_id"] or row["head_stage"]!=a["head_stage"]:
            raise ValueError("ANNOTATION_VIDEO_OR_STAGE_MISMATCH")
        if row["head_stage"]!=e["phenological_stage"]:
            raise ValueError("ANNOTATION_EPISODE_STAGE_MISMATCH")
        if (e["pre_entry_guild"]=="unresolved" or
            e["pre_entry_guild_evidence"] not in VERIFIED_GUILDS or
            not row["pre_entry_taxon_evidence_uri"].strip()):
            raise ValueError("INSECT_ROLE_CANNOT_BE_ASSIGNED_AFTER_EGG_OR_POLLEN_SUCCESS")
        if a["continuous_coverage"]!="1":
            raise ValueError("UNCOVERED_CONTACT_SEQUENCE_CANNOT_PROVE_BARRIER")
        if row["contact_surface"] not in SURFACES:
            raise ValueError("UNRECOGNIZED_TRUE_BOTANICAL_CONTACT_SURFACE")
        parent_zone=a.get("contact_zone","")
        if parent_zone not in PARENT_ZONE_ALLOWED:
            raise ValueError("PARENT_CONTACT_ZONE_MISSING_OR_INVALID")
        if row["contact_surface"] not in PARENT_ZONE_ALLOWED[parent_zone]:
            raise ValueError("TRUE_CONTACT_SURFACE_CONTRADICTS_PARENT_PORTAL")
        for col in ("contact_resolved","spine_touch_verified",
                    "physical_blockage_verified","entry_not_achieved_at_attempt"):
            if row[col] not in FLAGS:
                raise ValueError("INVALID_TRI_STATE_CONTACT_FLAG")
        if row["scorer_blinded_to_oviposition_result"] not in {"0","1"}:
            raise ValueError("REVIEW_BLINDING_FLAG_MUST_BE_EXPLICIT")
        if not row["video_contact_evidence_uri"].strip():
            raise ValueError("VIDEO_CITATION_REQUIRED_FOR_EVERY_REAL_ANNOTATION")
        surface=row["contact_surface"]
        spine=row["spine_touch_verified"]
        if surface=="unresolved":
            if row["contact_resolved"]!="0" or spine!="NA":
                raise ValueError("UNRESOLVED_SURFACE_CANNOT_BE_NO_SPINE_TOUCH")
            counts["unknown"]+=1
        else:
            if row["contact_resolved"]!="1":
                raise ValueError("CLASSIFIED_SURFACE_NEEDS_REVIEWED_CONTACT")
            if surface in {"spine_tip","spine_shaft"} and spine!="1":
                raise ValueError("TRUE_SPINE_SURFACE_REQUIRES_VERIFIED_TOUCH")
            if surface not in {"spine_tip","spine_shaft"} and spine!="0":
                raise ValueError("NONSPINE_OR_NONE_CANNOT_BE_SPINE_TOUCH")
        if spine=="1":
            counts["touches"]+=1
        if row["entry_not_achieved_at_attempt"]=="1" and a["access_observed"]!="0":
            raise ValueError("FAILED_ACCESS_CONTRADICTS_PARENT_ATTEMPT")
        if (row["physical_blockage_verified"]=="1" or
            row["entry_not_achieved_at_attempt"]=="1"):
            if row["contact_resolved"]!="1":
                raise ValueError("UNASSESSED_CONTACT_CANNOT_SUPPORT_CONFIRMED_FAILURE")
        if row["physical_blockage_verified"]=="1":
            if (surface not in ARMATURE_SURFACES or
                row["entry_not_achieved_at_attempt"]!="1" or
                a["obstruction_observed"]!="1" or
                a["access_observed"]!="0"):
                raise ValueError("PHYSICAL_BLOCKAGE_REQUIRES_REAL_BRACT_TOUCH_AND_PARENT_OBSTRUCTION")
            if not row["anatomy_scale_evidence_uri"].strip():
                raise ValueError("UNSCALED_ARMATURE_CANNOT_IDENTIFY_MECHANICAL_GATE")
            counts["blockages"]+=1
        for col in ("actual_spine_length_mm","minimum_phyllary_gap_mm"):
            value=row[col]
            if value not in ("NA",""):
                try:x=float(value)
                except ValueError:raise ValueError("ANATOMICAL_SCALE_NOT_NUMERIC") from None
                if not math.isfinite(x) or x<0:
                    raise ValueError("ANATOMICAL_DIMENSION_MUST_BE_NONNEGATIVE")
        if row["physical_blockage_verified"]=="1":
            required_dimension=("actual_spine_length_mm" if surface in
                                {"spine_tip","spine_shaft"} else
                                "minimum_phyllary_gap_mm")
            if row[required_dimension] in ("","NA"):
                raise ValueError("CONFIRMED_MECHANICAL_OBSTRUCTION_LACKS_RELEVANT_MEASURED_ANATOMY")
            if surface in {"spine_tip","spine_shaft"} and float(row[required_dimension])<=0:
                raise ValueError("POSITIVE_TRUE_SPINE_TOUCH_REQUIRES_NONZERO_SPINE_LENGTH")
        heads.add(tuple(row[k] for k in
                        ("population_id","individual_id","capitulum_id")))
    missing=expected-seen
    if missing:
        return {
            "status":"HOLD_INCOMPLETE_ELIGIBLE_ATTEMPT_ANNOTATION",
            "n_contact_rows":len(rows),
            "n_parent_attempts":len(attempts),
            "n_eligible_parent_attempts":len(expected),
            "n_unannotated_eligible_parent_attempts":len(missing),
            "n_parent_attempts_unresolved_guild":parent_guild_unresolved,
            "n_parent_attempts_unscorable_video":parent_coverage_incomplete,
            "n_verified_spine_touches":None,
            "n_verified_physical_blockages":None,
            "empirical_inference":"NOT_IDENTIFIABLE_DENOMINATOR_INCOMPLETE",
        }
    # A mechanically blocked FIRST route is not equivalent to preventing
    # reproductive-zone access. One candidate insect may switch portals
    # and succeed later during the SAME biological approach episode.
    from collections import defaultdict
    by_episode=defaultdict(list)
    for row in rows:
        key=tuple(row[k].strip() for k in ATTEMPT_KEY)
        try:
            idx=int(row["attempt_index"])
        except ValueError:
            raise ValueError("ATTEMPT_INDEX_INVALID_FOR_ROUTE_SEQUENCE") from None
        if idx<=0:
            raise ValueError("ATTEMPT_INDEX_INVALID_FOR_ROUTE_SEQUENCE")
        by_episode[key[:len(KEYS)]].append((
            idx, row["physical_blockage_verified"]=="1",
            at_index[key].get("access_observed")=="1",
        ))
    blocked_then_later_access=0
    episodes_with_any_verified_block=0
    for seq in by_episode.values():
        seq.sort(key=lambda x:x[0])
        if len({x[0] for x in seq})!=len(seq):
            raise ValueError("REPEATED_ATTEMPT_SEQUENCE_INDEX")
        blocked_indexes=[x[0] for x in seq if x[1]]
        if blocked_indexes:
            episodes_with_any_verified_block+=1
            if any(x[2] and x[0]>min(blocked_indexes) for x in seq):
                blocked_then_later_access+=1
    # Do not silently replace the candidate ARRIVAL denominator with
    # the selected subset of insects having photographed contact attempts.
    # Episodes with no sequenced attempt are UNKNOWN as to contact, not
    # demonstrated to have avoided the armature.
    eligible_event_index={
        key:ev for key,ev in ev_index.items()
        if ev.get("pre_entry_guild")!="unresolved"
        and ev.get("pre_entry_guild_evidence") in VERIFIED_GUILDS
    }
    no_documented_attempt=set(eligible_event_index)-set(by_episode)
    documented_final_service={
        "seed_feeder_candidate":{
            "confirmed_oviposition":0,"oviposition_not_confirmed_in_reviewed_episode":0,
            "oviposition_not_assessed":0
        },
        "legitimate_pollinator_candidate":{
            "positive_pollen_deposition_assay":0,
            "negative_pollen_deposition_assay":0,
            "pollen_deposition_not_assessed":0
        }
    }
    for key,ev in eligible_event_index.items():
        guild=ev.get("pre_entry_guild")
        if guild=="seed_feeder_candidate":
            state=ev.get("oviposition_confirmed","NA")
            if state not in {"0","1","NA"}:
                raise ValueError("OVIPOSITION_TRISTATE_INVALID")
            if state=="1" and (
                ev.get("behavioral_role")!="ovipositing_seed_feeder"
                or ev.get("role_evidence") in (None,"","not_assessed")
            ):
                raise ValueError("CONFIRMED_EGG_WITHOUT_INDEPENDENT_EVIDENCE")
            key_name={
                "1":"confirmed_oviposition",
                "0":"oviposition_not_confirmed_in_reviewed_episode",
                "NA":"oviposition_not_assessed"
            }[state]
            documented_final_service[guild][key_name]+=1
        elif guild=="legitimate_pollinator_candidate":
            state=ev.get("pollen_deposition_assay","not_assessed")
            if state not in {"positive","negative","not_assessed"}:
                raise ValueError("POLLEN_ASSAY_STATE_INVALID")
            # A claimed positive assay is invalid if the audited video
            # includes only a blocked outer-head attempt, without any
            # subsequent recorded successful floral-disc access. Likewise,
            # an event-level zone flag must agree with the contact sequence.
            event_key=tuple(ev.get(k,"") for k in KEYS)
            has_floral_access=any(
                attempt_key[:len(KEYS)]==event_key
                and attempt.get("contact_zone")=="floret_disc"
                and attempt.get("access_observed")=="1"
                for attempt_key,attempt in at_index.items()
            )
            if state=="positive" and (
                ev.get("anther_stigma_contact_observed")!="1"
                or ev.get("reproductive_zone_reached")!="1"
                or ev.get("role_evidence")!="pollen_deposition_measured"
                or not ev.get("evidence_uri","").strip()
                or not has_floral_access
            ):
                raise ValueError("POSITIVE_POLLEN_NEEDS_INDEPENDENT_ASSAY_AND_FLORAL_ACCESS")
            key_name={
                "positive":"positive_pollen_deposition_assay",
                "negative":"negative_pollen_deposition_assay",
                "not_assessed":"pollen_deposition_not_assessed",
            }[state]
            documented_final_service[guild][key_name]+=1
    return {
        "status":"DESCRIPTIVE_TRUE_CONTACT_GATES_ONLY",
        "n_independently_identified_candidate_approach_episodes":len(eligible_event_index),
        "n_candidate_episodes_without_recorded_attempt_sequence":len(no_documented_attempt),
        "no_recorded_attempt_is_NOT_verified_no_contact":True,
        "independent_functional_evidence_counts":documented_final_service,
        "functional_service_counts_are_NOT_viable_achene_effects":True,
        "n_observed_approach_episodes":len(by_episode),
        "n_episodes_with_any_verified_armature_block":episodes_with_any_verified_block,
        "n_episodes_with_verified_block_then_later_access":blocked_then_later_access,
        "blocking_without_later_access_is_NOT_confirmed_abandonment":True,
        "n_parent_attempts":len(attempts),
        "n_eligible_parent_attempts":len(expected),
        "n_unannotated_eligible_parent_attempts":0,
        "n_parent_attempts_unresolved_guild":parent_guild_unresolved,
        "n_parent_attempts_unscorable_video":parent_coverage_incomplete,
        "n_contact_rows":len(rows),
        "n_verified_spine_touches":counts["touches"],
        "n_verified_physical_blockages":counts["blockages"],
        "n_unresolved_contacts":counts["unknown"],
        "n_independent_heads":len(heads),
        "empirical_inference":"CONTACT_BEHAVIOUR_ONLY_NOT_CAUSAL_SELECTION_OR_FITNESS",
        "cannot_infer":["treatment effects", "seed set", "genotypic fitness",
                        "selection on trait covariance", "historical coevolution"],
    }

def synthetic_tests():
    # Independent pre-entry candidate guild, parent attempt recording
    # confirmed involucre obstruction with complete coverage.
    a={
        "individual_id":"I1","population_id":"POP1","capitulum_id":"H1",
        "observation_bout_id":"B1","approach_episode_id":"E1",
        "attempt_index":"1","video_clip_id":"synthetic://clip",
        "head_stage":"bud","contact_zone":"involucre_outer",
        "obstruction_observed":"1",
        "access_observed":"0","continuous_coverage":"1",
    }
    e={
        **{k:a[k] for k in KEYS},
        "pre_entry_guild":"seed_feeder_candidate",
        "pre_entry_guild_evidence":"independently_assigned_before_access",
        "phenological_stage":"bud",
    }
    x=dict(zip(COLUMNS,[
        "R1",a["individual_id"],a["population_id"],a["capitulum_id"],
        a["observation_bout_id"],a["approach_episode_id"],"1",
        a["video_clip_id"],"bud","spine_tip","1","1","1","1",
        "3.4","NA","synthetic://caliper","synthetic://clip",
        "synthetic://taxon","1"
    ]))
    def check(row=x, attempt=a, episode=e):
        return validate([copy.deepcopy(row)],[copy.deepcopy(attempt)],
                        [copy.deepcopy(episode)])
    assert check()["n_verified_physical_blockages"]==1
    assert check()["n_verified_spine_touches"]==1
    assert validate([],[],[])["n_verified_spine_touches"] is None
    assert validate([] ,[],[copy.deepcopy(e)])["status"]=="HOLD_EPISODES_EXIST_WITHOUT_ATTEMPT_SEQUENCE"
    unresolved=validate([],[copy.deepcopy(a)],[{**e,"pre_entry_guild":"unresolved",
                    "pre_entry_guild_evidence":"unresolved"}])
    assert unresolved["status"]=="HOLD_NO_ELIGIBLE_PARENT_ATTEMPTS"
    assert unresolved["n_verified_spine_touches"] is None
    # An annotator who only codes blocked approaches falsely inflates the
    # physical-barrier fraction. Missing successful attempts must block
    # ALL mechanistic summaries, even if positive contact is valid.
    a2={**a,"attempt_index":"2","contact_zone":"floret_disc",
        "obstruction_observed":"0","access_observed":"1"}
    x2={**x,"record_id":"R2","attempt_index":"2",
        "contact_surface":"floret_disc","spine_touch_verified":"0",
        "physical_blockage_verified":"0","entry_not_achieved_at_attempt":"0",
        "actual_spine_length_mm":"NA","minimum_phyllary_gap_mm":"NA"}
    missing=validate([copy.deepcopy(x)],[copy.deepcopy(a),a2],[copy.deepcopy(e)])
    assert missing["status"]=="HOLD_INCOMPLETE_ELIGIBLE_ATTEMPT_ANNOTATION"
    assert missing["n_unannotated_eligible_parent_attempts"]==1
    assert missing["n_verified_physical_blockages"] is None
    empty_annotation=validate([],[copy.deepcopy(a)],[copy.deepcopy(e)])
    assert empty_annotation["status"]=="HOLD_INCOMPLETE_ELIGIBLE_ATTEMPT_ANNOTATION"
    full=validate([copy.deepcopy(x),x2],[copy.deepcopy(a),a2],[copy.deepcopy(e)])
    assert full["n_contact_rows"]==full["n_eligible_parent_attempts"]==2
    assert full["n_verified_physical_blockages"]==1
    assert full["n_independent_heads"]==1
    assert full["n_observed_approach_episodes"]==1
    assert full["n_episodes_with_any_verified_armature_block"]==1
    assert full["n_episodes_with_verified_block_then_later_access"]==1
    assert full["n_independently_identified_candidate_approach_episodes"]==1
    assert full["n_candidate_episodes_without_recorded_attempt_sequence"]==0
    assert full["independent_functional_evidence_counts"]["seed_feeder_candidate"]["oviposition_not_assessed"]==1
    # Other reviewed candidate approaches can have NO annotated contact
    # episode. They remain an unknown contact state in the arrival
    # denominator, never classified as "avoided spine" or "no visit".
    e_no_contact={**e,"approach_episode_id":"E_NO_ATTEMPT"}
    extended=validate([copy.deepcopy(x),copy.deepcopy(x2)],
                      [copy.deepcopy(a),copy.deepcopy(a2)],
                      [copy.deepcopy(e),e_no_contact])
    assert extended["n_independently_identified_candidate_approach_episodes"]==2
    assert extended["n_candidate_episodes_without_recorded_attempt_sequence"]==1
    assert extended["n_verified_physical_blockages"]==1
    # An independently verified egg-laying end point must NOT be inferred
    # from a temporarily successful floret access attempt.
    assert full["independent_functional_evidence_counts"]["seed_feeder_candidate"]["confirmed_oviposition"]==0
    true_egg={**e,"oviposition_confirmed":"1","behavioral_role":"ovipositing_seed_feeder",
              "role_evidence":"egg_scar_and_video"}
    egg=validate([copy.deepcopy(x),copy.deepcopy(x2)],
                 [copy.deepcopy(a),copy.deepcopy(a2)],[true_egg])
    assert egg["independent_functional_evidence_counts"]["seed_feeder_candidate"]["confirmed_oviposition"]==1
    false_egg={**e,"oviposition_confirmed":"1"}
    try:
        validate([copy.deepcopy(x),copy.deepcopy(x2)],
                 [copy.deepcopy(a),copy.deepcopy(a2)],[false_egg])
    except ValueError as ex:
        assert "CONFIRMED_EGG_WITHOUT_INDEPENDENT_EVIDENCE" in str(ex)
    else:
        raise AssertionError("Successful arrival was mislabelled confirmed oviposition")
    assert check()["n_episodes_with_verified_block_then_later_access"]==0
    # A post-contact positive pollen assay cannot be fabricated from a
    # successful landing/attempt or an insect-order-only classification.
    pollinator={**e,"pre_entry_guild":"legitimate_pollinator_candidate",
        "pollen_deposition_assay":"positive"}
    bad_pollen=validate([copy.deepcopy(x)],[copy.deepcopy(a)],
                        [{**pollinator,"pollen_deposition_assay":"not_assessed"}])
    assert bad_pollen["independent_functional_evidence_counts"]["legitimate_pollinator_candidate"]["pollen_deposition_not_assessed"]==1
    try:
        validate([copy.deepcopy(x)],[copy.deepcopy(a)],[pollinator])
    except ValueError as ex:
        assert "POSITIVE_POLLEN_NEEDS_INDEPENDENT_ASSAY_AND_FLORAL_ACCESS" in str(ex)
    else:
        raise AssertionError("Pollen delivery was fabricated from mere access")
    # Even an independently labelled positive assay cannot be joined to a
    # video episode documenting ONLY a blocked outer spine attempt.
    false_arrival={**pollinator,"anther_stigma_contact_observed":"1",
        "reproductive_zone_reached":"1",
        "role_evidence":"pollen_deposition_measured",
        "evidence_uri":"synthetic://assay"}
    try:
        validate([copy.deepcopy(x)],[copy.deepcopy(a)],[false_arrival])
    except ValueError as ex:
        assert "POSITIVE_POLLEN_NEEDS_INDEPENDENT_ASSAY_AND_FLORAL_ACCESS" in str(ex)
    else:
        raise AssertionError("Blocked-only video incorrectly confirmed flower access")
    # A truly positive synthetic floral assay also needs a documented
    # later successful floret-disc access attempt on an open head.
    pollen_attempt1={**a,"head_stage":"full_anthesis"}
    pollen_attempt2={**a2,"head_stage":"full_anthesis"}
    pollen_row1={**x,"head_stage":"full_anthesis"}
    pollen_row2={**x2,"head_stage":"full_anthesis"}
    good_pollen={**false_arrival,"phenological_stage":"full_anthesis"}
    yes=validate([pollen_row1,pollen_row2],
                 [pollen_attempt1,pollen_attempt2],[good_pollen])
    assert yes["independent_functional_evidence_counts"]["legitimate_pollinator_candidate"]["positive_pollen_deposition_assay"]==1
    bad=[
        ({**x,"spine_touch_verified":"0"},a,e,"TRUE_SPINE_SURFACE"),
        ({**x,"anatomy_scale_evidence_uri":""},a,e,"UNSCALED_ARMATURE"),
        ({**x,"actual_spine_length_mm":"NA"},a,e,"LACKS_RELEVANT_MEASURED"),
        ({**x,"entry_not_achieved_at_attempt":"0"},a,e,"PHYSICAL_BLOCKAGE"),
        ({**x,"video_contact_evidence_uri":""},a,e,"VIDEO_CITATION"),
        (x,{**a,"continuous_coverage":"0"},e,"UNCOVERED_CONTACT"),
        (x,a,{**e,"pre_entry_guild_evidence":"unresolved"},"INSECT_ROLE"),
        ({**x,"contact_surface":"unresolved"},a,e,"UNRESOLVED_SURFACE"),
        ({**x,"approach_episode_id":"E_DOES_NOT_EXIST"},a,e,"ANNOTATION_MUST_MATCH"),
        ({**x,"actual_spine_length_mm":"-2"},a,e,"ANATOMICAL_DIMENSION_MUST_BE_NONNEGATIVE"),
        (x,{**a,"contact_zone":"floret_disc"},e,"TRUE_CONTACT_SURFACE_CONTRADICTS_PARENT_PORTAL"),
    ]
    for rr,aa,ee,needle in bad:
        try:check(rr,aa,ee)
        except ValueError as er:
            assert needle in str(er),(needle,str(er))
        else:
            raise AssertionError(f"INVALID_REAL_CONTACT_INFERRED:{needle}")
    return {
        "status":"PASS_SYNTHETIC_REAL_ARMATURE_AND_PARENT_LINK_GATES",
        "positive_fixtures":1,"rejected_bad_scientific_claims":len(bad),
        "real_organismal_field_records":0,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--contact",type=Path,default=Path(
        "data/intake/capitulum_true_armature_contact_v1.csv"))
    p.add_argument("--attempts",type=Path,default=Path(
        "data/intake/capitulum_contact_attempt_sequence_v1.csv"))
    p.add_argument("--events",type=Path,default=Path(
        "data/intake/capitulum_guild_access_event_ledger_v1.csv"))
    p.add_argument("--effort",type=Path,default=Path(
        "data/intake/capitulum_video_effort_denominator_v1.csv"))
    p.add_argument("--effort-contract",type=Path,default=Path(
        "data/contracts/capitulum_video_effort_denominator_v1.json"))
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    ch,contacts=read_csv(a.contact)
    ah,attempts=read_csv(a.attempts)
    eh,events=read_csv(a.events)
    if ch!=COLUMNS or not set(ATTEMPT_KEY).issubset(ah) or not set(KEYS).issubset(eh):
        raise ValueError("CONTACT_ATTEMPT_OR_EVENT_CSV_CONTRACT_DRIFT")
    effort_spec=json.loads(a.effort_contract.read_text(encoding="utf-8"))
    effort_header,efforts=read_effort_csv(a.effort)
    if effort_header!=effort_spec["columns"]:
        raise ValueError("AUTHORITY_VIDEO_EFFORT_SCHEMA_DRIFT")
    observed_effort=validate_video_effort(efforts,events,effort_spec)
    result={
        "verified_observation_effort":observed_effort,
        "version":"capitulum_true_armature_contact_gate_v1",
        "date":"2026-10-08",
        "synthetic_tests":synthetic_tests(),
        "actual_intake":validate(contacts,attempts,events),
        "science":"True structural barrier = independently identified before access, documented anatomical touch, verified blockage, unsuccessful attempt. No phylogenetic or plant-fitness effect from this alone.",
    }
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",
                         encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
