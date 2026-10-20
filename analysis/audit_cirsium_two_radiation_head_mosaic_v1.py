#!/usr/bin/env python3
"""Real-source Cirsium HEAD geometry mosaics in two published taxon panels.

Question: do orientation and outward-phyllary posture deterministically form
one indivisible observed categorical head package? This is a taxon-coded
descriptive falsification of one-to-one state mapping, not experimental
causation, rate estimation, independent evolutionary gains, or genetics.

Source panels deliberately have DIFFERENT ONTOLOGIES and missingness:
 - North American 2023 study: erect(0)/pending(1), ascending(0) vs
   spreading(1)/reflexed(2) outer phyllaries, 63 original phylogeny tips.
 - EAzami Japanese NMNS authority: upward_or_erect / downward_or_nodding,
   appressed/ascending vs spreading/recurved (only directionally unambiguous).
 Mixed 'ascending_or_recurved' is UNKNOWN, never recoded as negative.
 Secondary-history JPN_06 dipsacolepis is kept SEPARATE from dominant sample.
 Exact authoritative botanical source text not new individual measurements.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import re
from pathlib import Path

SOURCE={
 "north_america":{
   "repo":"carol-siniscalchi/cirsium",
   "commit":"c925d4af86160064448f794a9715f70b917a2786",
   "paper_doi":"10.1086/724310",
   "traits_blob":"fa127f0587b676e6e004ea121424d26069586bcb",
   "tree_blob":"6417a38793427bc7c8562dc0e2cab59b74d8c316",
 },
 "japan":{
   "repo":"zuizui0223/EAzami",
   "commit":"c9870e91604ed1514344ddc8f765eae5237dbf27",
   "japan38_membership_blob":"426d7e5c015d0bd92a22d1cd2f4eba0d39fd49a8",
   "nmns_seed_blob":"c4da867491e913119d78e4b07df800965a0f021d",
   "nmns_extension_blob":"346381eeb1326c654cec22888467eb446ab058dd",
 }
}
CELL_CODES={
 (0,0):"erect_not_outward",
 (0,1):"erect_outward",
 (1,0):"nodding_not_outward",
 (1,1):"nodding_outward",
}
JP_ORIENTATION={
 "upward_or_erect":0,
 "downward_or_nodding":1,
 "unknown":None,
 "upward_or_ascending":None,
 "source_conflict_index_downward_detail_erect":None,
}
JP_OUTER={
 "appressed":0,
 "ascending":0,
 "appressed_or_ascending":0,
 "spreading":1,
 "recurved":1,
 "spreading_or_recurved":1,
 "unknown":None,
 "ascending_or_recurved":None,
}
SECONDARY_JAPAN={"JPN_06","JPN_15"}
MIN_JAPAN_COMPLETE=7


def check(ok,code):
    if not ok:raise ValueError(code)


def blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()


def read_verified(path, expected, record_kind):
    raw=Path(path).read_bytes()
    check(blob(raw)==expected,"ORIGINAL_SOURCE_GIT_BLOB_CHANGED_"+record_kind)
    return raw.decode("utf-8-sig")


def csv_rows(data, identity):
    rows=list(csv.DictReader(io.StringIO(data)))
    check(bool(rows),"EMPTY_"+identity)
    check(all(None not in z and z for z in rows),"INVALID_CSV_RECORD_"+identity)
    return rows


def unique_rows(rows,key):
    codes=[r.get(key) for r in rows]
    check(all(codes),"MISSING_SOURCE_TAXON")
    check(len(codes)==len(set(codes)),"DUPLICATE_SOURCE_TAXON")
    return {r[key]:r for r in rows}


def na_state(row):
    a=row.get("HeadPosition")
    b=row.get("OuterPhyllaryPosition")
    check(a in ("0","1"),"INVALID_NA_HEAD_ORIENTATION")
    check(b in ("0","1","2"),"INVALID_NA_OUTER_PHYLLARY")
    return int(a),int(b!="0")


def jp_state(row):
    orient=row.get("orientation_state")
    outer=row.get("phyllary_posture")
    check(orient in JP_ORIENTATION,"INVALID_JP_ORIENTATION_SOURCE_CODE")
    check(outer in JP_OUTER,"INVALID_JP_PHYLLARY_SOURCE_CODE")
    return JP_ORIENTATION[orient],JP_OUTER[outer]


def role(jpn_id):
    return "secondary" if jpn_id in SECONDARY_JAPAN else "dominant"


def fisher_exact_two_sided(a,b,c,d):
    """Hypergeometric probability-ordering two-sided Fisher, no SciPy."""
    check(all(isinstance(v,int) and v>=0 for v in (a,b,c,d))
          and a+b>0 and c+d>0,"INVALID_FISHER_TABLE")
    n=a+b+c+d
    row=a+b
    yes=a+c
    def hyper(k):
        return math.comb(yes,k)*math.comb(n-yes,row-k)/math.comb(n,row)
    reference=hyper(a)
    p=sum(hyper(k) for k in range(max(0,row-(n-yes)),min(row,yes)+1)
          if hyper(k)<=reference+1e-12)
    return min(1.,p)


def matrix_summary(records):
    """Measured categorical pair mapping, not independence or transition."""
    check(bool(records),"EMPTY_COHORT")
    check(all((a,b) in CELL_CODES for a,b in records),"INVALID_BINARY_HEAD_CELL")
    counts={v:sum((a,b)==k for a,b in records) for k,v in CELL_CODES.items()}
    all_four=all(counts.values())
    good_deletions=sum(
        all(sum(z==state for j,z in enumerate(records) if j!=i)>0
            for state in CELL_CODES)
        for i in range(len(records)))
    by_head={str(i):sorted({b for a,b in records if a==i}) for i in (0,1)}
    by_outer={str(i):sorted({a for a,b in records if b==i}) for i in (0,1)}
    # A deterministic mapping requires each observed input category to have
    # exactly one output category; missing input category is not a pass.
    fixed_outer=all(len(v)<=1 for v in by_head.values())
    fixed_head=all(len(v)<=1 for v in by_outer.values())
    erect_out=counts["erect_outward"]
    erect_not=counts["erect_not_outward"]
    nod_out=counts["nodding_outward"]
    nod_not=counts["nodding_not_outward"]
    return {
        "n_taxa":len(records),
        "four_cell_counts":counts,
        "n_observed_combinations":sum(bool(x) for x in counts.values()),
        "all_four_observed":all_four,
        "deterministic_outer_given_orientation":fixed_outer,
        "deterministic_orientation_given_outer":fixed_head,
        "orientation_only_can_uniquely_predict_outer":fixed_outer,
        "all_single_taxon_deletions_preserve_four":good_deletions==len(records),
        "leave_one_taxon_out_preserves_four_n":good_deletions,
        "descriptive_orientation_outer_fisher_p_two_sided":
            fisher_exact_two_sided(nod_out,nod_not,erect_out,erect_not)
            if (nod_out+nod_not)>0 and (erect_out+erect_not)>0 else None,
        "head_to_outer_observed_state_options":by_head,
        "outer_to_head_observed_state_options":by_outer,
        "warning":"nonrandom, source-coded tips, not independent lineages; Fisher p ignores phylogeny and sampling"
    }


def claim_ceiling():
    return {
      "observed_morphology_equals_adaptive_module":False,
      "four_tip_combinations_prove_independent_transitions":False,
      "jp_and_na_ontologies_are_identical_quantitative_angles":False,
      "same_head_visitor_access_fitness_verified":False,
      "genetic_independence_or_ancestral_trait_memory_demonstrated":False,
      "japan_secondary_and_dominant_origins_are_merged_as_one_independent_clade":False,
      "synonymous_expressions_are_same_measurement_protocol":False,
      "missing_japanese_source_trait_equals_negative":False,
      "whole_head_fully_nonintegrated":False,
    }


def analyze(na,jp_base,jp_extra,members,tree_text,source_check=True):
    north=unique_rows(na,"Species")
    japan=unique_rows(jp_base+jp_extra,"paper_japan_member_id")
    roster=unique_rows(members,"paper_japan_member_id")
    check(len(north)==63 and len(jp_base)==22 and len(jp_extra)==2
          and len(japan)==24 and len(roster)==38,
          "ORIGINAL_PANEL_SIZE_CHANGED")
    check(set(japan).issubset(roster),"NONJAPAN38_TAXON_LEAKED")
    check(all(v.get("authority_match_status")=="exact_authority_concept_match"
              for v in japan.values()),
          "SOURCE_AUTHORITY_CONCEPT_MISMATCH")
    names=re.findall(r"(Cirsium_[A-Za-z_]+)(?=:)",tree_text)
    check(len(names)==len(set(names))==63 and set(names)==set(north),
          "NA_63_NUCLEAR_TREE_TIPS_NOT_MATCHED")
    northern=[(taxon,na_state(row)) for taxon,row in sorted(north.items())]
    states=[(j,row,jp_state(row)) for j,row in sorted(japan.items())]
    scored=[(k,v[0],v[1],row) for k,row,v in states if None not in v]
    alljp=[(a,b) for _,a,b,_ in scored]
    dominant=[(a,b) for j,a,b,_ in scored if role(j)=="dominant"]
    secondary=[(a,b) for j,a,b,_ in scored if role(j)=="secondary"]
    check(len(scored)==MIN_JAPAN_COMPLETE and len(dominant)==6
          and len(secondary)==1,"SOURCE_JP_COMPLETE_HEAD_COHORT_CHANGED")
    nm=matrix_summary([v for _,v in northern])
    collapsed_names={"Cirsium_arizonicum_arizonicum",
                     "Cirsium_arizonicum_rothrockii"}
    check(collapsed_names.issubset(north) and
          na_state(north["Cirsium_arizonicum_arizonicum"])==
          na_state(north["Cirsium_arizonicum_rothrockii"]),
          "DUPLICATE_SUBSPECIES_OR_STATE_DRIFT")
    collapsed=matrix_summary([v for taxon,v in northern
        if taxon!="Cirsium_arizonicum_rothrockii"])
    check(collapsed["n_taxa"]==62 and
          collapsed["four_cell_counts"]=={
             "erect_not_outward":45,"erect_outward":12,
             "nodding_not_outward":2,"nodding_outward":3},
          "NA_SUBSPECIES_COLLAPSE_SENSITIVITY_CHANGED")
    jm=matrix_summary(alljp)
    dm=matrix_summary(dominant)
    sm=matrix_summary(secondary)
    check(nm["four_cell_counts"]=={
      "erect_not_outward":45,"erect_outward":13,
      "nodding_not_outward":2,"nodding_outward":3},
      "NA_HEAD_STATE_MATRIX_CHANGED")
    check(jm["four_cell_counts"]=={
      "erect_not_outward":3,"erect_outward":1,
      "nodding_not_outward":2,"nodding_outward":1},
      "JP_HEAD_STATE_MATRIX_CHANGED")
    check(dm["four_cell_counts"]=={
      "erect_not_outward":3,"erect_outward":0,
      "nodding_not_outward":2,"nodding_outward":1},
      "JP_DOMINANT_HEAD_STATE_MATRIX_CHANGED")
    results={
      "version":"cirsium_20261021_two_radiation_head_mosaic_real_sources_v1",
      "sources":SOURCE,
      "source_cohorts":{
        "north_america_trait_and_tree_tip_labels":len(north),
        "japan38_published_roster":len(roster),
        "japan_nmns_authority_coded":len(japan),
        "japan_orientation_unambiguous":sum(a is not None for _,_,(a,b) in states),
        "japan_outer_direction_unambiguous":sum(b is not None for _,_,(a,b) in states),
        "japan_both_axes_unambiguous":len(scored),
        "japan_secondary_history_in_cohort":[
             j for j,a,b,row in scored if role(j)=="secondary"],
        "japan_source_missing_ontology_not_biological_absence":True,
      },
      "morphological_comparison":{
        "north_america_63":nm,
        "north_america_second_arizonicum_subspecies_removed_62":collapsed,
        "japan_authority_all_7":jm,
        "japan_dominant_radiation_6":dm,
        "japan_secondary_history_1":sm,
      },
      "japan_actual_matched_taxon_codes":[{
         "paper_member_id":j,
         "taxon":row["paper_taxon_concept"],
         "source_orientation":row["orientation_state"],
         "source_outer_phyllary_posture":row["phyllary_posture"],
         "orientation_binary":a,"outward_phyllary_binary":b,
         "evolutionary_radiation_role":role(j),
      } for j,a,b,row in scored],
      "japan_excluded_unscorable_taxa":[{
         "paper_member_id":j,
         "original_orientation_state":row["orientation_state"],
         "original_outer_phyllary_state":row["phyllary_posture"],
         "missing_orientation":a is None,
         "missing_outer_state":b is None,
      } for j,row,(a,b) in states if a is None or b is None],
      "ontology":{
        "north_america_orientation":"HeadPosition: erect=0; pending=1",
        "japan_orientation":"NMNS strict upward_or_erect=0; downward_or_nodding=1",
        "north_america_outer":"OuterPhyllaryPosition ascending=0; spreading or reflexed=1",
        "japan_outer":"appressed,ascending,appressed_or_ascending=0; spreading,recurved,spreading_or_recurved=1",
        "held_unknown":"Japanese ascending_or_recurved; source_conflict_index_downward_detail_erect; upward_or_ascending; literal unknown remain NA",
        "comparability":"COARSENED_DIRECTIONAL_ANALOG_ONLY_NOT_IDENTICAL_PHYLLARY_ANGLE_OR_TRUE_SPINE",
      },
      "interpretation":{
        "observed_main_result":"The two morphological axes cannot be represented as one deterministic binary-state package in either the 63-tip NA or the 7-tip Japanese authority sample: each orientation occurs with both outer postures.",
        "na_four_combinations_robust_to_any_one_tip_deletion":nm["all_single_taxon_deletions_preserve_four"],
        "japan_four_combinations_require_secondary_history_JPN_06":True,
        "japan_dominant_only_has_three_combinations":dm["n_observed_combinations"]==3,
        "direct_phylogenetic_correlated_selection_test":False,
        "na_original_pollination_syndrome_unused_in_this_morphology_comparison":True,
        "fitness_interaction_or_pollen_deposition_observed":False,
        "japan_7_coverage_not_an_independent_random_sample":True,
      },
      "claim_ceiling":claim_ceiling(),
    }
    return results


def main():
    p=argparse.ArgumentParser()
    for arg in ("na-traits","na-tree","jp-roster","jp-base","jp-extension","out"):
        p.add_argument("--"+arg,type=Path,required=True)
    args=p.parse_args()
    N=SOURCE["north_america"]; J=SOURCE["japan"]
    na=csv_rows(read_verified(args.na_traits,N["traits_blob"],"NA_TRAITS"),"NA_TRAITS")
    tree=read_verified(args.na_tree,N["tree_blob"],"NA_TREE")
    roster=csv_rows(read_verified(args.jp_roster,J["japan38_membership_blob"],"JPN38"),"JPN38")
    jp_base=csv_rows(read_verified(args.jp_base,J["nmns_seed_blob"],"NMNS_BASE"),"NMNS_BASE")
    jp_ext=csv_rows(read_verified(args.jp_extension,J["nmns_extension_blob"],"NMNS_EXTENSION"),"NMNS_EXTENSION")
    result=analyze(na,jp_base,jp_ext,roster,tree)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
       "status":"SOURCE_PINNED_TWO_REAL_CIRSIUM_TAXON_PANELS",
       "north_america_4_states":result["morphological_comparison"]["north_america_63"]["four_cell_counts"],
       "japan_all7_states":result["morphological_comparison"]["japan_authority_all_7"]["four_cell_counts"],
       "japan_dominant6_states":result["morphological_comparison"]["japan_dominant_radiation_6"]["four_cell_counts"],
       "NA_63_leave1_retains_all4":result["interpretation"]["na_four_combinations_robust_to_any_one_tip_deletion"],
       "both_head_states_known_JP":result["source_cohorts"]["japan_both_axes_unambiguous"],
    },indent=2))


if __name__=="__main__":
    main()
