#!/usr/bin/env python3
"""Source-pinned Cirsium botanical *state-label* transfer audit, not memory.

The Japanese Japan38 NMNS authority states are discovery; published East Asian
species with independent botanical state evidence are provisional held-out
*taxon labels*. They are NOT verified phylogenetically independent clades.

The Brier negative control asks ONLY whether the discovery marginal fraction
of upward heads predicts the held-out categorical labels better than 0.5.
It does NOT estimate patristic-memory rho, evolutionary rates, homoplasy,
selection or conservation of a genetic/developmental module.

Original EAzami SHA c9870e91604ed1514344ddc8f765eae5237dbf27.
No new Cirsium individual phenotypes, gene expression or seed fitness data.
"""
from __future__ import annotations
import argparse
import collections
import csv
import json
import math
from fractions import Fraction
from pathlib import Path

STATES=("upward_or_erect","downward_or_nodding")
SRC_COLUMNS={"paper_japan_member_id","paper_taxon_concept",
    "authority_match_status","orientation_state","phyllary_posture",
    "stickiness_state","source_url"}
HOLD_COLUMNS={"analysis_taxon","provenance_frame","orientation_state",
    "stickiness_state","phyllary_state","source_url","discovery_overlap",
    "orientation_pool_status","stickiness_pool_status"}
MEMBER_COLUMNS={"paper_japan_member_id","paper_taxon_concept","tree_codes"}
FROZEN_EAZAMI_SHA="c9870e91604ed1514344ddc8f765eae5237dbf27"

def read(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as fp:
        reader=csv.DictReader(fp)
        return set(reader.fieldnames or ()), list(reader)

def demand(ok,code):
    if not ok: raise ValueError(code)

def unique_rows(rows,field,code):
    data={}
    for row in rows:
        name=(row.get(field) or "").strip()
        demand(bool(name),code+"_EMPTY_NAME")
        demand(name not in data,code+"_DUPLICATE_NAME")
        data[name]=row
    return data

def canonical_concept(x):
    """Species-token guard only; never collapse species into varieties."""
    tokens=x.strip().split()
    demand(len(tokens)>=2 and tokens[0]=="Cirsium",
           "NOT_A_CIRSIUM_AUTHORITY_TAXON")
    if len(tokens)>=4 and tokens[2]=="var.":
        return " ".join(tokens[:4])
    return " ".join(tokens[:2])

def brier(y,p):
    demand(bool(y) and math.isfinite(float(p)) and 0<=p<=1,
           "INVALID_BRIER_INPUT")
    return sum((Fraction(int(v))-p)**2 for v in y)/len(y)

def audit(member,base,extra,held):
    members=unique_rows(member,"paper_japan_member_id","JAPAN38")
    demand(len(members)==38 and
           all(v.startswith("JPN_") for v in members),
           "SOURCE_JAPAN38_COHORT_IS_NOT_38")
    # One extra botanical source per separate paper taxon.
    source=unique_rows(base+extra,"paper_japan_member_id","JAPAN38_NMNS")
    demand(len(base)==22 and len(extra)==2 and len(source)==24,
           "SOURCE_NMNS_24_CONCEPT_COVERAGE_CHANGED")
    japonica=collections.Counter()
    known={z:0 for z in ("orientation","phyllary","stickiness")}
    source_species=set()
    for cid,r in source.items():
        demand(cid in members,"SOURCE_TAXON_NOT_PART_OF_38_PAPER_COHORT")
        demand(r["authority_match_status"]=="exact_authority_concept_match",
               "AUTHORITY_NAME_NOT_EXACT")
        demand(r["source_url"].startswith("https://"),
               "MISSING_BOTANICAL_PRIMARY_CITATION")
        # The paper membership contains taxonomic authorship; the authority
        # fields carry exact species/variety concepts, not specimen colour.
        sp=canonical_concept(r["paper_taxon_concept"])
        source_species.add(sp)
        demand(canonical_concept(members[cid]["paper_taxon_concept"])==sp,
               "JAPAN38_PAPER_AND_BOTANICAL_CONCEPT_DISAGREE")
        state=r["orientation_state"]
        demand(state in STATES+("unknown",),
               "UNRESOLVED_JAPANESE_ORIENTATION_WAS_FORCE_CODED")
        if state in STATES:
            japonica[state]+=1
            known["orientation"]+=1
        for field,category in [("phyllary_posture","phyllary"),
                               ("stickiness_state","stickiness")]:
            if r[field]!="unknown":
                known[category]+=1
    demand(known=={"orientation":19,"phyllary":10,"stickiness":13}
           and dict(japonica)=={
               "upward_or_erect":13,
               "downward_or_nodding":6},
           "ORIGINAL_JAPANESE_TRAIT_TALLY_CHANGED")

    held_index=unique_rows(held,"analysis_taxon","EAST_ASIAN_HELDOUT")
    demand(len(held_index)==19,"EAST_ASIAN_AUTHORITY_ROWS_NOT_19")
    pending=collections.Counter()
    eligible=[]
    nonoverlap=0
    for name,h in held_index.items():
        concept=canonical_concept(name)
        demand(h["source_url"].startswith("https://"),
               "HELDOUT_WITHOUT_BOTANICAL_CITATION")
        is_disjoint=h["discovery_overlap"]=="no"
        if is_disjoint:
            nonoverlap+=1
            demand(concept not in source_species,
                   "SAME_DISCOVERY_TAXON_LEAKED_TO_HELDOUT")
        else:
            demand(h["discovery_overlap"]=="yes",
                   "HELDOUT_OVERLAP_UNKNOWN_NOT_ASSUMED_INDEPENDENT")
        o=h["orientation_pool_status"]
        s=h["orientation_state"]
        pending[o]+=1
        if o=="candidate_primary":
            demand(is_disjoint,
                   "DISCOVERY_OVERLAPPING_TAXON_USED_TO_SCORE_PORTABILITY")
            demand(s in STATES,
                   "AMBIGUOUS_ORIENTATION_FORCED_TO_KNOWN_STATE")
            eligible.append(h)
        else:
            demand(o in {
                "not_codable","exclude_ambiguous",
                "exclude_overlap_and_ambiguity",
                "exclude_discovery_overlap",
            },"UNKNOWN_HELDOUT_SELECTION_STATUS")
    demand(nonoverlap==16 and len(eligible)==11,
           "HELDOUT_ORIENTATION_ELIGIBILITY_CHANGED")
    selected=collections.Counter(x["orientation_state"] for x in eligible)
    demand(dict(selected)=={
        "upward_or_erect":6,"downward_or_nodding":5},
        "HELDOUt_ORIENTATION_LABEL_MIX_CHANGED")
    sticky=[h for h in held if h["stickiness_pool_status"]=="candidate_primary"]
    demand(len(sticky)==2 and all(h["discovery_overlap"]=="no" for h in sticky),
           "STICKINESS_HOLDOUT_COUNT_OR_LEAKAGE_CHANGED")
    # Note: external nonsticky does NOT exactly equal the Japanese
    # merged category "nonsticky_or_nearly_nonsticky"; do not harmonize.
    phyllary=[h for h in held
              if h["discovery_overlap"]=="no" and
              h["phyllary_state"] not in {"unknown","ambiguous",""}]
    demand(len(phyllary)==1,
           "HELDOUT_PHYLLARY_ELIGIBILITY_CHANGED")

    p_j=Fraction(japonica["upward_or_erect"],known["orientation"])
    y=[h["orientation_state"]=="upward_or_erect" for h in eligible]
    out_of_group=brier(y,p_j)
    uninformative=brier(y,Fraction(1,2))
    # An in-sample heldout climatology is a descriptive oracle ONLY,
    # explicitly not an honest predictive benchmark.
    heldout_p=Fraction(selected["upward_or_erect"],len(y))
    oracle=brier(y,heldout_p)
    return {
       "status":"REAL_TAXON_STATE_NEGATIVE_CONTROL_COMPLETE_PHYLOGENETIC_PORTABILITY_HOLD",
       "source_snapshot":FROZEN_EAZAMI_SHA,
       "source_units":"taxon-level NMNS authority botanical labels, not individuals",
       "discovery":{
         "paper_cohort_taxa":len(members),
         "botanical_coded_taxon_concepts":len(source),
         "orientation_scored_taxa":known["orientation"],
         "phyllary_scored_taxa":known["phyllary"],
         "stickiness_scored_taxa":known["stickiness"],
         "orientation_states":dict(japonica),
       },
       "external_candidate_pool":{
         "source_backed_taxon_rows":len(held),
         "no_discovery_overlap_labels":nonoverlap,
         "orientation_predeclared_eligible_taxa":len(eligible),
         "orientation_state_counts":dict(selected),
         "orientation_status_including_exclusions":dict(sorted(pending.items())),
         "stickiness_predeclared_eligible_taxa":len(sticky),
         "phyllary_codable_nonoverlap_taxa":len(phyllary),
         "external_taxon_concepts_are_not_independent_phylogenetic_events":True,
         "phyllary_and_stickiness_codes_not_assumed_comparable":True,
       },
       "taxon_label_only_sanity_check":{
         "task":"predict each heldout up/down taxon label using Japanese marginal P(up)",
         "japan_upward_rate":float(p_j),
         "heldout_upward_rate":float(heldout_p),
         "trained_majority_prediction":"upward_or_erect",
         "trained_majority_heldout_correct_taxa":selected["upward_or_erect"],
         "trained_majority_heldout_total_taxa":len(eligible),
         "trained_majority_heldout_accuracy":float(heldout_p),
         "train_Japan_upward_probability_heldout_Brier":float(out_of_group),
         "predefined_half_probability_heldout_Brier":float(uninformative),
         "Japanese_marginal_improves_over_half":out_of_group<uninformative,
         "heldout_marginal_oracle_Brier_not_valid_forecast":float(oracle),
         "interpretation":"Source-selection-sensitive categorical occurrence only; NOT trait phylogenetic memory transfer, statistical test or historical common cause",
       },
       "Q1_phylogenetic_memory":{
         "inferential_status":"NOT_IDENTIFIABLE_MISSING_MATCHED_PHYLOGENETIC_DISTANCE_ENSEMBLE",
         "empirical_patristic_pairs_for_east_asian_holdout":0,
         "reason":"Public trait labels for 11 distinct external species do not supply independently validated nuclear tree placement and comparable branch lengths for those taxa",
         "two_Ryukyu_Arenicola_tip_pair_insufficient_for_memory_correlation":True,
         "required_to_test":[
            "shared time/branch-length phylogeny with source/holdout taxa and topology uncertainty",
            "homology- and stage-matched trait scores, uncertainty and sampling fraction by lineage",
            "family/clade holdout split defined without seeing outcomes; avoid dependent tip pairs as independent replicates",
            "compare full predictive loss on heldout phylogenetic distances against trait-name-only and within-clade baselines"
         ],
         "do_not_conclude":[
            "same evolutionary memory is portable or nonportable",
            "patristic memory rho from state frequency Brier loss",
            "11 taxa are 11 independent evolutionary transitions",
            "geography automatically establishes independent lineage history",
         ],
       },
    }

def synthetic_negative_controls(members,base,extra,held):
    import copy
    def rejects(ms,bs,es,hs,needle):
        try: audit(ms,bs,es,hs)
        except ValueError as err:
            demand(needle in str(err),
                   "WRONG_NEGATIVE_CONTROL_ERROR:"+str(err))
            return
        raise AssertionError("INVALID_EVOLUTIONARY_EVIDENCE_ACCEPTED:"+needle)
    observed=audit(members,base,extra,held)
    demand(observed["external_candidate_pool"]["orientation_predeclared_eligible_taxa"]==11,
           "ORIGINAL_REAL_TEST_NEGATIVE_CONTROL")
    z=copy.deepcopy(base);z.append(copy.deepcopy(z[0]))
    rejects(members,z,extra,held,"JAPAN38_NMNS_DUPLICATE_NAME")
    z=copy.deepcopy(extra);z[0]["paper_japan_member_id"]="JPN_99"
    rejects(members,base,z,held,"SOURCE_TAXON_NOT_PART_OF_38_PAPER_COHORT")
    z=copy.deepcopy(base);z[0]["orientation_state"]="ambiguous"
    rejects(members,z,extra,held,"UNRESOLVED_JAPANESE_ORIENTATION_WAS_FORCE_CODED")
    z=copy.deepcopy(held)
    n=next(i for i,h in enumerate(z) if h["orientation_pool_status"]=="exclude_ambiguous")
    z[n]["orientation_pool_status"]="candidate_primary"
    rejects(members,base,extra,z,"AMBIGUOUS_ORIENTATION_FORCED_TO_KNOWN_STATE")
    z=copy.deepcopy(held)
    n=next(i for i,h in enumerate(z) if h["discovery_overlap"]=="yes")
    z[n]["orientation_pool_status"]="candidate_primary"
    rejects(members,base,extra,z,"DISCOVERY_OVERLAPPING_TAXON_USED_TO_SCORE_PORTABILITY")
    z=copy.deepcopy(held);z.append(copy.deepcopy(z[0]))
    rejects(members,base,extra,z,"EAST_ASIAN_HELDOUT_DUPLICATE_NAME")
    z=copy.deepcopy(base);z[0]["paper_taxon_concept"]="Cirsium fabricated_species"
    rejects(members,z,extra,held,"JAPAN38_PAPER_AND_BOTANICAL_CONCEPT_DISAGREE")
    return {"status":"PASS_SYNTHETIC_SOURCE_LEAKAGE_AND_MEMORY_CLAIM_GATES",
            "n_source_corruption_negative_controls":7,
            "real_phylogenetic_memory_estimates_produced":0,
            "real_taxon_labels_verified":19+11,
            "taxon_label_forecasts_are_not_phylogenetic_memory":True}

def main():
    p=argparse.ArgumentParser()
    for key in ("japan38-members","nmns-base","nmns-extension","east-asia-heldout"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--out",type=Path)
    args=p.parse_args()
    inputs=[]
    for path,required in (
      (args.japan38_members,MEMBER_COLUMNS),
      (args.nmns_base,SRC_COLUMNS),
      (args.nmns_extension,SRC_COLUMNS),
      (args.east_asia_heldout,HOLD_COLUMNS)
    ):
        cols,rows=read(path)
        demand(required.issubset(cols),"ORIGINAL_EAZAMI_TRAIT_HEADER_CHANGED")
        inputs.append(rows)
    membership,seed,extension,external=inputs
    result={
        "version":"cirsium_three_question_state_portability_v1",
        "date":"2026-10-19",
        "authoritative_EAzami_commit":FROZEN_EAZAMI_SHA,
        "negative_controls":synthetic_negative_controls(
            membership,seed,extension,external),
        "existing_real_taxon_state_audit":audit(
            membership,seed,extension,external),
    }
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",
                            encoding="utf-8")
    x=result["existing_real_taxon_state_audit"]
    print(json.dumps({
        "Japanese_scored_orientation":x["discovery"]["orientation_scored_taxa"],
        "Eastern_candidate_orientation":x["external_candidate_pool"]["orientation_predeclared_eligible_taxa"],
        "Japan_marginal_heldout_Brier":x["taxon_label_only_sanity_check"]["train_Japan_upward_probability_heldout_Brier"],
        "Q1_memory":x["Q1_phylogenetic_memory"]["inferential_status"],
        "check":result["negative_controls"]["status"],
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
