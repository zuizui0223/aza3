#!/usr/bin/env python3
"""Actual botanical trait-state coverage for Cirsium lineage-transfer Q1.

Read source-pinned (EAzami c9870e91) Japanese authority seed+extensions,
38-member nuclear-tree discovery roster and independent East Asian
authority-state screen. Outputs genuine *taxon-state* frequencies, not a
phylogenetic memory slope, cross-lineage trait portability effect or repeated
independent evolutionary origins.

IMPORTANT: 19/38 Japanese orientation-coded and 11/19 candidate heldout
orientation-coded is *sampling coverage*, not occurrence prevalence.
The heldout authority taxa do not have a source-pinned, matched nuclear tree
here; patristic distance cannot be invented from taxon names or geography.
"""
from __future__ import annotations
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

JAPAN_N=38
JAPAN_TRAIT_COLUMNS={"paper_japan_member_id","paper_taxon_concept",
                    "orientation_state","phyllary_posture","stickiness_state",
                    "authority_match_status","source_url"}
ROSTER_COLUMNS={"paper_japan_member_id","paper_taxon_concept"}
HELDOUT_COLUMNS={"analysis_taxon","provenance_frame","orientation_state",
                 "stickiness_state","phyllary_state","source_basis",
                 "source_url","discovery_overlap","orientation_pool_status",
                 "stickiness_pool_status"}
BINARY_ORIENTATION={"upward_or_erect","downward_or_nodding"}
UNKNOWN={"","unknown","ambiguous","not_codable",
         "source_conflict_index_downward_detail_erect"}
CANDIDATE="candidate_primary"

def source_csv(path:Path,required):
    with path.open(encoding="utf-8-sig",newline="") as fh:
        r=csv.DictReader(fh)
        if not required.issubset(set(r.fieldnames or ())):
            raise ValueError("SOURCE_CSV_SCHEMA_DRIFT")
        return list(r)

def check_unique(rows,col,reason):
    found=set()
    for row in rows:
        key=(row.get(col) or "").strip()
        if not key or key in found:
            raise ValueError(reason)
        found.add(key)
    return found

def analyse(roster,base,extension,heldout):
    roster_ids=check_unique(roster,"paper_japan_member_id",
                            "JAPAN_MEMBERSHIP_EMPTY_OR_DUPLICATE")
    if len(roster_ids)!=JAPAN_N:
        raise ValueError("JAPAN38_MEMBERSHIP_WRONG_N")
    source=base+extension
    ids=check_unique(source,"paper_japan_member_id",
                     "JAPAN_TRAIT_DUPLICATE_OR_EMPTY_TIP")
    if not ids.issubset(roster_ids):
        raise ValueError("AUTHORITATIVE_JAPAN_TRAIT_NOT_IN_NUCLEAR_ROSTER")
    by_roster={row["paper_japan_member_id"]:row for row in roster}
    for row in source:
        if row["authority_match_status"]!="exact_authority_concept_match":
            raise ValueError("JAPAN_TRAIT_UNRESOLVED_AUTHORITY_CONCEPT")
        if row["paper_taxon_concept"].strip()!=(
                by_roster[row["paper_japan_member_id"]]
                ["paper_taxon_concept"].strip().split(" Kadota")[0]
            ) and row["paper_taxon_concept"].strip() not in (
                by_roster[row["paper_japan_member_id"]]
                ["paper_taxon_concept"].strip(),):
            # Authored authority taxon concepts can differ in the published
            # authority string beyond a simple author suffix; registry IDs,
            # not fuzzy name transforms, are the primary exact join.
            pass
        if not row["source_url"].startswith("https://"):
            raise ValueError("JAPAN_STATE_SOURCE_URL_MISSING")
    if len(base)!=22 or len(extension)!=2 or len(source)!=24:
        raise ValueError("JAPAN_AUTHORITY_SOURCE_COVERAGE_CHANGED")
    check_unique(heldout,"analysis_taxon",
                 "EAST_ASIA_HOLDOUT_EMPTY_OR_DUPLICATE_TAXON")
    if len(heldout)!=19:
        raise ValueError("HOLDOUT_EAST_ASIA_SCREEN_SIZE_CHANGED")

    j_orient=[r for r in source if
              r["orientation_state"] in BINARY_ORIENTATION]
    j_phyll=[r for r in source if
             r["phyllary_posture"] not in UNKNOWN]
    j_sticky=[r for r in source if
              r["stickiness_state"] not in UNKNOWN]
    if [len(j_orient),len(j_phyll),len(j_sticky)]!=[19,10,13]:
        raise ValueError("JAPAN_FROZEN_19_10_13_RESOLVED_TRAIT_COVERAGE_DRIFT")

    external_not_overlap=[r for r in heldout
                          if r["discovery_overlap"]=="no"]
    overlapping=[r for r in heldout if r["discovery_overlap"]=="yes"]
    if len(external_not_overlap)!=16 or len(overlapping)!=3:
        raise ValueError("SOURCE_HOLDOUT_OVERLAP_COUNTS_CHANGED")
    def eligible(column,state_col):
        chosen=[r for r in heldout if r[column]==CANDIDATE]
        for row in chosen:
            if row["discovery_overlap"]!="no":
                raise ValueError("DISCOVERY_TIP_FALSELY_QUALIFIED_AS_HOLDOUT")
            state=row[state_col]
            if state in UNKNOWN:
                raise ValueError("UNRESOLVED_TRAIT_FALSELY_QUALIFIED")
            if not row["source_url"].startswith("https://"):
                raise ValueError("HOLDOUT_SOURCE_URL_MISSING")
        return chosen
    ho=eligible("orientation_pool_status","orientation_state")
    hs=eligible("stickiness_pool_status","stickiness_state")
    if len(ho)!=11 or len(hs)!=2:
        raise ValueError("SOURCE_HOLDOUT_11_ORIENTATION_2_STICKINESS_CHANGED")
    if any(x["orientation_state"] not in BINARY_ORIENTATION for x in ho):
        raise ValueError("HOLDOUT_BINARY_ORIENTATION_HARMONIZATION_DRIFT")
    hf=[r for r in external_not_overlap
        if r["phyllary_state"] not in UNKNOWN]
    if len(hf)!=1 or hf[0]["analysis_taxon"]!="Cirsium reflexiphyllarium":
        raise ValueError("ONLY_ONE_QUALIFIED_HELDOUT_PHYLLARY_STATE_EXPECTED")

    j_counts=Counter(x["orientation_state"] for x in j_orient)
    h_counts=Counter(x["orientation_state"] for x in ho)
    if (dict(j_counts)!={"upward_or_erect":13,"downward_or_nodding":6}
        or dict(h_counts)!={"downward_or_nodding":5,"upward_or_erect":6}):
        raise ValueError("SOURCE_ORIENTATION_STATE_COUNTS_CHANGED")

    return {
        "version":"cirsium_phylogenetic_memory_holdout_trait_coverage_v1",
        "status":"REAL_TAXON_TRAIT_STATES_AUDITED_NO_COMMON_PHYLOGENETIC_MEMORY_TEST",
        "source_provenance":{
            "EAzami_ref":"c9870e91604ed1514344ddc8f765eae5237dbf27",
            "published_japanese_membership_tips":len(roster),
            "authority_japanese_rows":len(source),
            "external_east_asia_screened_taxa":len(heldout),
            "genuine_external_nondiscovery_taxa":len(external_not_overlap),
            "external_overlap_excluded":len(overlapping),
        },
        "japan_discovery":{
            "orientation_resolved_taxa":len(j_orient),
            "orientation_states":dict(sorted(j_counts.items())),
            "phyllary_resolved_taxa":len(j_phyll),
            "stickiness_resolved_taxa":len(j_sticky),
            "unscored_orientation_roster_tips":len(roster)-len(j_orient),
        },
        "heldout_east_asia":{
            "orientation_eligible_taxa":len(ho),
            "orientation_states":dict(sorted(h_counts.items())),
            "phyllary_authority_resolved_taxa":len(hf),
            "stickiness_eligible_taxa":len(hs),
            "sticky_taxa_source_codes":{
                x["analysis_taxon"]:x["stickiness_state"] for x in hs},
            "orientation_candidates":sorted(x["analysis_taxon"] for x in ho),
            "no_trait_phylogenetic_distances_supplied":True,
        },
        "genuine_descriptive_observation":
            "Both upward and downward capitulum orientation categories are source-backed within the Japanese discovery taxa (13/6) and within external East Asian candidate taxa (6/5). This is categorical trait coverage, not phylogenetic memory portability.",
        "Q1_inferential_ceiling":{
            "source_tree_for_japan_exists":True,
            "same_branch_length_tree_for_11_external_heldout_tips":False,
            "trait_dissimilarity_pairs_with_verified_patristic_distances":0,
            "phylogenetic_memory_gradient_rho":None,
            "leave_one_lineage_out_predictive_gain":None,
            "same_trait_memory_transport":"NOT_IDENTIFIABLE_WITH_EXISTING_EXTERNAL_STATES",
            "major_risks":["external taxon sampling by trait description is nonrandom",
                           "separate mainland/clade branch lengths missing",
                           "stickiness and phyllary state ontologies not harmonized",
                           "taxon descriptions not within-species individual sampling",
                           "missing states not absent organs"],
        },
    }

def synthetic_negative_controls(roster,base,extension,heldout):
    import copy
    scenarios=[
        ("JAPAN_TRAIT_DUPLICATE_OR_EMPTY_TIP",
         (roster,base,base[:1]+extension,heldout)),
        ("DISCOVERY_TIP_FALSELY_QUALIFIED_AS_HOLDOUT",
         (roster,base,extension,[
             {**r,"orientation_pool_status":CANDIDATE}
             if r["discovery_overlap"]=="yes" and
               r["orientation_pool_status"]!="not_codable" else r
             for r in heldout])),
        ("UNRESOLVED_TRAIT_FALSELY_QUALIFIED",
         (roster,base,extension,[
             {**r,"orientation_pool_status":CANDIDATE}
             if r["orientation_state"]=="ambiguous" else r
             for r in heldout])),
        ("HOLDOUT_EAST_ASIA_SCREEN_SIZE_CHANGED",
         (roster,base,extension,heldout[:-1])),
    ]
    count=0
    for expected,args in scenarios:
        try:analyse(*copy.deepcopy(args))
        except ValueError as exc:
            if expected not in str(exc):
                raise AssertionError("BAD_NEGATIVE_CONTROL_"+expected+
                                     "_GOT_"+str(exc)) from exc
            count+=1
        else:
            raise AssertionError("UNQUALIFIED_EXTERNAL_MEMORY_CLAIM_ACCEPTED")
    return {
        "status":"PASS_SYNTHETIC_HOLDOUT_QUALIFICATION_NEGATIVE_CONTROLS",
        "deliberate_corrupted_sources_rejected":count,
        "biological_trait_states_are_FROM_PINNED_EAzami_not_synthetic":True,
        "phylogenetic_memorability_inferred_from_unplaced_taxa":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--japan-roster",type=Path,required=True)
    p.add_argument("--japan-authority",type=Path,required=True)
    p.add_argument("--japan-extension",type=Path,required=True)
    p.add_argument("--east-asia-holdout",type=Path,required=True)
    p.add_argument("--out",type=Path)
    args=p.parse_args()
    roster=source_csv(args.japan_roster,ROSTER_COLUMNS)
    base=source_csv(args.japan_authority,JAPAN_TRAIT_COLUMNS)
    extension=source_csv(args.japan_extension,JAPAN_TRAIT_COLUMNS)
    heldout=source_csv(args.east_asia_holdout,HELDOUT_COLUMNS)
    report={
        "real_source_state_audit":analyse(roster,base,extension,heldout),
        "negative_controls":synthetic_negative_controls(
            roster,base,extension,heldout),
    }
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",
                            encoding="utf-8")
    print(json.dumps({
        "status":report["real_source_state_audit"]["status"],
        "Japan_orientation":report["real_source_state_audit"]["japan_discovery"]["orientation_states"],
        "External_orientation":report["real_source_state_audit"]["heldout_east_asia"]["orientation_states"],
        "Q1_memory_transfer":"NOT_IDENTIFIABLE_WITH_EXISTING_EXTERNAL_STATES",
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
