#!/usr/bin/env python3
"""Partial ID for Cirsium nippoense JPN19 adult phyllary posture in Japan.

NMNS authority records upright-or-erect orientation and an outer-phyllary
range from ascending to short recurving. The original EAzami strict state
is 'ascending_or_recurved': NEVER convert to observed binary 1 or binary 0.

Compute separate mathematically possible outcomes for the already verified
Japan dominant-radiation six-tip four-cell table. The outputs are BOUNDS
under alternative source-state interpretations, not new measured plants.
Within-individual plasticity, within-population polymorphism, pollination,
genetic evolution and population prevalence are NOT identifiable here.

Primary botanical source: NMNS C. nippoense detail no. 110 (Kahaku)
https://www.kahaku.go.jp/research/db/botany/azami/detail.html?no=110
Related original botanical description Kadota 2004 BNSM Ser.B 30(2) 63-69.
Frozen EAzami taxon-source data @ c9870e91604ed1514344ddc8f765eae5237dbf27.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

NMNS_BLOB="c4da867491e913119d78e4b07df800965a0f021d"
SOURCE_EAZAMI_COMMIT="c9870e91604ed1514344ddc8f765eae5237dbf27"
JPN19="JPN_19"
EXPECTED_BASE={
    "erect_not_outward":3,"erect_outward":0,
    "nodding_not_outward":2,"nodding_outward":1,
}
FLORA_URL="https://www.kahaku.go.jp/research/db/botany/azami/detail.html?no=110"


def require(test,msg):
    if not test:raise ValueError(msg)


def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def load_original_nmns(path):
    raw=Path(path).read_bytes()
    require(git_blob(raw)==NMNS_BLOB,"NMNS_ORIGINAL_GIT_BLOB_CHANGED")
    rows=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
    require(len(rows)==22,"NMNS_ORIGINAL_22_SOURCE_ROWS_CHANGED")
    match=[x for x in rows if x["paper_japan_member_id"]==JPN19]
    require(len(match)==1,"JPN19_UNIQUE_SOURCE_RECORD_MISSING")
    return match[0]


def claim_limit():
    return {
        "jpn19_fixed_outward_morph_assigned":False,
        "within_plant_reversible_phyllary_change_measured":False,
        "actual_adaptive_selection_or_pollination_measured":False,
        "ancestor_state_or_allele_reexpression_inferred":False,
        "independent_phylogenetic_events_added":False,
        "confirmed_fourth_state_in_Japan_dominant_radiation":False,
        "n_source_head_individuals_newly_measured":0,
    }


def outcome_envelope(row,baseline):
    """Both endpoints of a source category spanning either direction."""
    require(dict(baseline)==EXPECTED_BASE,"PRIOR_SOURCE_COHORT_CHANGED")
    require(row.get("paper_japan_member_id")==JPN19 and
            row.get("paper_taxon_concept")=="Cirsium nippoense" and
            row.get("authority_match_status")=="exact_authority_concept_match" and
            row.get("source_url")==FLORA_URL,
            "ORIGINAL_JPN19_IDENTITY_OR_PROVENANCE_CHANGED")
    require(row.get("orientation_state")=="upward_or_erect",
            "UNRESOLVED_ERECT_SOURCE_STATE")
    require(row.get("phyllary_posture")=="ascending_or_recurved",
            "ORIGINAL_AMBIGUOUS_SOURCE_CHANGED")
    require(row.get("stickiness_state")=="sticky",
            "JPN19_SOURCE_STICKINESS_CHANGED")
    a=dict(baseline)
    b=dict(baseline)
    a["erect_not_outward"]+=1
    b["erect_outward"]+=1
    def summary(tab):
        return {"n_taxon_concepts_if_JPN19_is_assigned":sum(tab.values()),
                "four_cell_counts":tab,
                "n_observed_combinations":sum(v>0 for v in tab.values()),
                "NOT_an_observed_resolution":True}
    return {
       "authority_taxon":"Cirsium nippoense",
       "authority_record_id":JPN19,
       "original_orientation":"upward_or_erect",
       "original_outer_phyllary_descriptor":"ascending_or_recurved",
       "sticky_original_botanical_description":True,
       "source_authority_url":FLORA_URL,
       "n_existing_confirmed_dominant_taxa":sum(baseline.values()),
       "source_mixed_descriptor_options":["ascending", "shortly_recurved"],
       "forced_new_source_tip_state":None,
       "if_ascending":summary(a),
       "if_recurved":summary(b),
       "n_observed_combinations_bound":[
          min(sum(x>0 for x in tab.values()) for tab in (a,b)),
          max(sum(x>0 for x in tab.values()) for tab in (a,b)),
       ],
       "within_single_plant_polymorphism_demonstrated":False,
       "within_single_population_frequency_or_dev_stage":None,
       "distinguishing_observation_needed":
          "Voucher-linked botanical closeup/field measurements of phyllary posture and capitulum orientation for C. nippoense, including plant, head, stage and site; do not disturb rare populations without authorization.",
       "interpretation":"A single source taxon describes a direction RANGE. Under 0 or 1 assignment Japan dominant state repertoire may be 3 or 4; neither is known to hold for a fixed individual or taxon phenotype.",
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--nmns",type=Path,required=True)
    ap.add_argument("--existing-receipt",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    opts=ap.parse_args()
    nmns=load_original_nmns(opts.nmns)
    lock=json.loads(opts.existing_receipt.read_text(encoding="utf-8"))
    source=lock["source"]
    aa=next(x for x in source if x["name"]=="jp_base")
    require(aa["source_sha"]==NMNS_BLOB and aa["git_commit"]==SOURCE_EAZAMI_COMMIT,
            "PREVIOUS_JP_BOTANICAL_SOURCE_NOT_FROZEN")
    prev=lock["observed_matrices"]["jp_main6"]
    require(prev["n"]==6,"PRIOR_JP6_SIZE_CHANGED")
    result={
       "version":"aza3_20261021_jpn19_mixed_phyllary_source_range_bounds_v1",
       "source":{"repo":"zuizui0223/EAzami","commit":SOURCE_EAZAMI_COMMIT,
                 "original_nmns_git_blob":NMNS_BLOB,
                 "NMNS_url":FLORA_URL,
                 "biological_data_added":0},
       "actual_source_state":outcome_envelope(nmns,prev["four_cell_counts"]),
       "claim_ceiling":claim_limit(),
    }
    opts.out.parent.mkdir(parents=True,exist_ok=True)
    opts.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
       "status":"ORIGINAL_BOTANICAL_STATE_RANGE_NOT_FIXED_TAXON_TRAIT",
       "new_species":"Cirsium nippoense",
       "strict_JP_existing":result["actual_source_state"]["n_existing_confirmed_dominant_taxa"],
       "possible_head_combination_count_range":result["actual_source_state"]["n_observed_combinations_bound"],
       "direct_new_fieldhead_measurements":0
    },indent=2))


if __name__=="__main__":
    main()
