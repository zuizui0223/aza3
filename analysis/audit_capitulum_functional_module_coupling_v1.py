#!/usr/bin/env python3
"""An audit of *capitulum-only* coupling, not a selection test.

Use immutable published file blobs:
- Azami complete 9-construct / 36-pair image RV table
- EAzami Japanese-radiation orientation/phyllary/stickiness history
- EAzami authority specimen/taxon states
Never treat trait RV or transition overlap as evidence of selection or causal
genetic covariation. Fetch pinned ref + verify GitHub blob SHA before reading.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.request import Request, urlopen

SOURCES = {
 "azami_rv": {
   "url":"https://raw.githubusercontent.com/zuizui0223/azami/b463dbc938e454baa54b2e480f7f1d8badd7795d/reproducibility/current_reference/upgrade/complete18_construct_pairwise.csv",
   "git_blob_sha":"8be80b4a42fd0265ab95b8209871db0001d0352f",
 },
 "eazami_history": {
   "url":"https://raw.githubusercontent.com/zuizui0223/EAzami/c9870e91604ed1514344ddc8f765eae5237dbf27/data/evidence/japan38_multitrait_history_summary_v1.json",
   "git_blob_sha":"ddfbf0692c71ed83ac0a0d4eaa9404c338a74b57",
 },
 "eazami_traits": {
   "url":"https://raw.githubusercontent.com/zuizui0223/EAzami/c9870e91604ed1514344ddc8f765eae5237dbf27/data/evidence/japan38_nmns_capitulum_trait_seed_v1.csv",
   "git_blob_sha":"c4da867491e913119d78e4b07df800965a0f021d",
 },
}
EXPLICIT_CANDIDATE_PAIRS = {
 "orientation_x_head_elongation":("presentation_angle","head_elongation"),
 "orientation_x_projection":("presentation_angle","projection_prominence"),
 "orientation_x_projection_pattern":("presentation_angle","projection_pattern"),
 "floral_chroma_x_projection":("floral_chroma","projection_prominence"),
 "projection_strength_x_pattern":("projection_prominence","projection_pattern"),
 "orientation_x_involucre":("presentation_angle","involucre_form"),
 "lightness_x_hue":("floral_lightness","floral_hue"),
 "shape_x_involucre":("head_elongation","involucre_form"),
}
ALLOWED_CONSTRUCTS = {
 "presentation_angle","floral_lightness","floral_chroma","floral_hue",
 "head_elongation","head_compactness","involucre_form",
 "projection_prominence","projection_pattern",
}

def sha_git_blob(b):
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def fetch_sources():
    raw = {}
    for label, source in SOURCES.items():
        with urlopen(Request(source["url"], headers={"User-Agent":"aza3-capitulum-audit-v1"}),timeout=40) as r:
            binary=r.read()
        if sha_git_blob(binary)!=source["git_blob_sha"]:
            raise RuntimeError(f"ERROR_PINNED_UPSTREAM_BLOB_MISMATCH {label}")
        raw[label]=binary
    return raw

def run(raw):
    for label,source in SOURCES.items():
        if sha_git_blob(raw[label]) != source["git_blob_sha"]:
            raise ValueError(f"ERROR_SOURCE_DRIFT {label}")
    pairs=list(csv.DictReader(io.StringIO(raw["azami_rv"].decode("utf-8-sig"))))
    assert len(pairs)==36
    lookup={}
    for r in pairs:
        k=frozenset((r["construct_left"],r["construct_right"]))
        if len(k)!=2 or k in lookup: raise ValueError("MISSING_OR_DUPLICATE_RELATION")
        if not k.issubset(ALLOWED_CONSTRUCTS): raise ValueError("UNKNOWN_CONSTRUCT")
        w,a=float(r["within_taxon_rv"]),float(r["among_taxon_rv"])
        if not(0<=w<=1 and 0<=a<=1): raise ValueError("INVALID_RV_RANGE")
        lookup[k]={"within_taxon_rv":w,"among_taxon_rv":a,
                   "difference_among_minus_within":a-w}
    if len({t for e in lookup for t in e})!=9:
        raise ValueError("CONSTRUCT_COVERAGE_FAILURE")
    selected={}
    for label,(left,right) in EXPLICIT_CANDIDATE_PAIRS.items():
        entry=lookup[frozenset((left,right))]
        rv=entry["among_taxon_rv"]
        rank=1+sum(r["among_taxon_rv"]>rv for r in lookup.values())
        selected[label]=dict(construct_left=left,construct_right=right,
                             rank_of_36_among_taxon=rank,**entry,
                             statistical_significance="NOT_ESTIMATED_PER_PAIR",
                             causal_interpretation="NOT_IDENTIFIED")
    hist=json.loads(raw["eazami_history"])
    states=list(csv.DictReader(io.StringIO(raw["eazami_traits"].decode("utf-8-sig"))))
    if len(states)!=22: raise ValueError("UNEXPECTED_AUTHORITY_SAMPLE_SIZE")
    joint=[r for r in states
           if(r["orientation_state"].startswith(("upward","downward"))
              and r["phyllary_posture"]!="unknown")]
    examples=[]
    for r in joint:
        examples.append({"taxon":r["paper_taxon_concept"],
            "orientation":r["orientation_state"],
            "phyllary_posture":r["phyllary_posture"],
            "stickiness":r["stickiness_state"],
            "source_url":r["source_url"]})
    if len(joint)!=8: raise ValueError("JOINT_STATE_COVERAGE_DRIFT")
    overlaps=hist["common_lability_diagnostic"]
    bn=overlaps["branch_length_aware_ml_excess_transition_overlap"]
    tn=overlaps["topology_only_equal_branch_ufboot1000"]
    crit={}
    for label,source in [
        ("orientation_x_phyllary",("orientation_x_phyllary",)),
        ("orientation_x_stickiness",("orientation_x_stickiness",)),
        ("phyllary_x_stickiness",("phyllary_x_stickiness",)),
    ]:
        b=bn[label];t=tn[label]
        positive_across_layers=bool(
            b["rho"]>0
            and b["one_sided_stratified_p"]<0.05
            and t["q05"]>0
            and t["positive_fraction"]>=0.95
        )
        crit[label]={
            "branch_length_aware_rho":b["rho"],
            "branch_length_aware_stratified_one_sided_p":b["one_sided_stratified_p"],
            "topology_only_median_rho":t["median_rho"],
            "topology_only_q05":t["q05"],
            "topology_only_positive_fraction":t["positive_fraction"],
            "robustly_positive_both_layers":positive_across_layers,
        }
    # Refuse to report the fragile orientation x stickiness branch-length p
    # as an unqualified coevolution discovery.
    if any(v["robustly_positive_both_layers"] for v in crit.values()):
        raise ValueError("CRITICAL_SENSITIVITY_DECISION_DRIFT")
    if "No module pair is consistently positive" not in overlaps["decision"]:
        raise ValueError("EAZAMI_CONCLUSION_DRIFT")
    result={
        "version":"capitulum_functional_module_empirical_boundary_v1",
        "status_date":"2026-10-08",
        "inputs":SOURCES,
        "phenotypic_scale":"Azami: 1,734 observations/42 taxa; image-RV construct associations",
        "azami_all_36_pairs":{
            "n_constructs":9,"n_pairs":len(lookup),
            "n_among_stronger_than_within":
                sum(v["among_taxon_rv"]>v["within_taxon_rv"] for v in lookup.values()),
            "n_within_stronger_than_among":
                sum(v["within_taxon_rv"]>v["among_taxon_rv"] for v in lookup.values()),
            "candidate_subset":selected,
        },
        "history":{
            "minimum_changes":{
               "orientation":hist["minimum_change_history"]["orientation"],
               "phyllary_posture":hist["minimum_change_history"]["phyllary_posture"],
               "stickiness":hist["minimum_change_history"]["stickiness"],
            },
            "trait_pairs_sensitive_to_branch_model":crit,
            "robust_shared_transition_localization_pairs":0,
            "tested_pairs":3,
            "conclusion":"NO_ROBUST_SYNCHRONOUS_TRANSITION_LOCALIZATION_IN_THREE_SPARSE_MORPHOLOGICAL_STATES; a null result is NOT evidence of adaptive asynchronous change",
        },
        "authority_state_support":{
            "source_taxon_concepts":len(states),
            "joint_orientation_phyllary_concepts":len(joint),
            "examples":examples,
            "warning":"Phyllary posture not spine length; incomplete species states and species/taxon geography, not within-population repeated events"
        },
        "claims":{
            "current_visible_phenotypic_coupling":"SUPPORTED_SELECTED_PAIRS_BUT_NOT_EVERYWHERE",
            "real_spine_or_bract_mechanical_function":"NOT_DIRECTLY_MEASURED_IN_AZAMI",
            "genetic_trait_correlation":"NOT_IDENTIFIED",
            "historical_synchronous_evolution":"NOT_ROBUSTLY_SUPPORTED",
            "historical_adaptive_nonsynchrony":"NOT_IDENTIFIED",
            "ecological_trait_interaction_on_viable_seed":"NOT_IDENTIFIED",
            "guild_dependent_synchronous_vs_compensatory_hypothesis":"PROSPECTIVE_ONLY",
        },
        "scientific_boundaries":[
            "Within- and among-taxon RV is unsigned and not a transition-rate or fitness interaction.",
            "Only preselected illustrative edge RV magnitudes reported, not their per-edge p-values or phylogenetic significance.",
            "The 0/3 robust shared localization may reflect low state resolution, missingness, and uncertainty; do not promote absence of coevolution.",
            "Projection prominence is 2D involucre outline, not spine length/force; image angle need not be gravity-referenced.",
            "Historical module names reflect categorical phyllary posture, stickiness and orientation, not all nine Azami constructs.",
            "Direct biotic role and adaptive benefit must be demonstrated by separate guild, time and viable-achene contrasts."
        ]
    }
    if result["azami_all_36_pairs"]["n_among_stronger_than_within"]!=33:
        raise ValueError("AZAMI_RELATION_COUNT_DRIFT")
    return result

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    result=run(fetch_sources())
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
       "status":"PASS_36_PAIR_AUDIT",
       "n_pairs":result["azami_all_36_pairs"]["n_pairs"],
       "n_among_stronger":result["azami_all_36_pairs"]["n_among_stronger_than_within"],
       "illustrative_pairs":result["azami_all_36_pairs"]["candidate_subset"],
       "history":result["history"]["trait_pairs_sensitive_to_branch_model"],
       "claim_boundary":"NO_FITNESS_OR_ADAPTATION_IDENTIFIED",
    },indent=2))

if __name__=="__main__":main()
