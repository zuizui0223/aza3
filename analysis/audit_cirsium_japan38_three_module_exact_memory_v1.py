#!/usr/bin/env python3
"""Exact tip-label permutation comparison for three Cirsium capitulum traits.

REPRODUCES: EAzami@c9870e91604ed1514344ddc8f765eae5237dbf27
    real Comp1061 Japan38 nuclear substitutions/site tree + NMNS authority
    trait states and concept-to-nuclear-tip map, source SHA256 checked.

Tests pairwise patristic distance against the 0/1 indicator that two
phylogenetic tips have DIFFERENT strict known trait states.
EXACT exhaustive reassignment of observed category frequencies across
the n scored tips, not the invalid 171-pair i.i.d. approximation:
    orientation: 19 tips, 13 U / 6 D -> 27132 assignments
    stickiness: 13 tips, 7 sticky / 6 nonsticky-or-nearly -> 1716
    phyllary posture: only 6 unambiguously SINGLETON-coded taxa
                      (4 ascending, 1 spreading, 1 appressed) -> 30

The phyllary n=6 screen is very low coverage and is always exploratory.
Taxa missing physical measurement or source-ambiguous words are excluded,
never imputed as a biological absence. The 2025 compatibility ML phylogram
has substitution branch lengths but no dated ancestral persistence.
No external clade-memory portability, adaptation, or genetic re-use claim.
"""
from __future__ import annotations
import argparse
import csv
import itertools
import json
import math
from collections import Counter
from pathlib import Path

from audit_cirsium_japan38_orientation_phylogenetic_memory_v1 import (
    prepare, rank_average, spearman_rank_with_binary, test_synthetic_geometry,
    REF, STATES as ORIENT, EXPECTED_SHA,
)

TRAITS={
    "orientation":{
        "source":"orientation_state",
        "allowed":ORIENT,
        "expected_counts":{"upward_or_erect":13,"downward_or_nodding":6},
        "minimum_taxa_for_inferential_test":10,
    },
    "stickiness":{
        "source":"stickiness_state",
        "allowed":{"sticky":1,"nonsticky_or_nearly_nonsticky":0},
        "expected_counts":{"sticky":7,"nonsticky_or_nearly_nonsticky":6},
        "minimum_taxa_for_inferential_test":10,
    },
    "phyllary_posture_STRICT_single_state_only":{
        "source":"phyllary_posture",
        "allowed":{"appressed":0,"ascending":1,"spreading":2,"recurved":3},
        "expected_counts":{"appressed":1,"ascending":4,"spreading":1},
        "minimum_taxa_for_inferential_test":10,
    },
}

def insist(test,code):
    if not test:raise ValueError(code)

def read(path):
    with open(path,encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def exact_reassignments(names,count):
    """Yield each distinct category allocation once, with fixed counts."""
    categories=sorted(count)
    free=tuple(range(len(names)))
    def rec(j,available,current):
        if j==len(categories)-1:
            if len(available)!=count[categories[j]]:
                raise AssertionError("Final exact state allocation invalid")
            final=current[:]
            for i in available:final[i]=categories[j]
            yield final
            return
        for taken in itertools.combinations(available,count[categories[j]]):
            x=current[:]
            for i in taken:x[i]=categories[j]
            selected=set(taken)
            keep=tuple(i for i in available if i not in selected)
            yield from rec(j+1,keep,x)
    yield from rec(0,free,[None]*len(names))

def statistic(tree,tips,states):
    n=len(tips)
    assert len(tips)==len(states)
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    distances=[tree.distance(tips[i],tips[j]) for i,j in pairs]
    insist(all(math.isfinite(d) and d>0 for d in distances),
           "REAL_TREE_NONPOSITIVE_OR_MISSING_DISTANCE")
    ranks=rank_average(distances)
    known=list(states)
    observed=[int(known[i]!=known[j]) for i,j in pairs]
    rho=spearman_rank_with_binary(ranks,observed)
    insist(rho is not None,"NO_TRAIT_CONTRAST_FOR_RHO")
    means={
      "same":sum(d for d,y in zip(distances,observed) if not y)/
             sum(1 for y in observed if not y),
      "different":sum(d for d,y in zip(distances,observed) if y)/
                  sum(observed),
    }
    counts=Counter(states)
    possible=math.factorial(n)
    for v in counts.values():possible//=math.factorial(v)
    total=0
    more=0
    same_tail=0
    null=[]
    for ys in exact_reassignments(tips,counts):
        y_diff=[int(ys[i]!=ys[j]) for i,j in pairs]
        r=spearman_rank_with_binary(ranks,y_diff)
        insist(r is not None,"EXACT_NULL_WITHOUT_CATEGORY_VARIATION")
        total+=1
        null.append(r)
        more+=abs(r)>=abs(rho)-1e-12
        same_tail+=r>=rho-1e-12
    insist(total==possible,"EXACT_FIXED_PREVALENCE_ENUMERATION_INCOMPLETE")
    null.sort()
    loo=[]
    for idx in range(n):
        sub=[k for k,(i,j) in enumerate(pairs)
             if i!=idx and j!=idx]
        r=spearman_rank_with_binary(
            rank_average([distances[k] for k in sub]),
            [observed[k] for k in sub])
        loo.append(r)
    return {
        "n_scored_taxa":n,
        "state_counts":dict(sorted(counts.items())),
        "n_nonindependent_pairs":len(pairs),
        "n_discordant_pairs":sum(observed),
        "n_fixed_prevalence_exact_assignments":total,
        "patristic_vs_discrete_mismatch_Spearman_rho":rho,
        "exact_two_sided_tip_label_permutation_p":more/total,
        "exact_positive_direction_tip_label_p":same_tail/total,
        "mean_patristic_same_state":means["same"],
        "mean_patristic_different_state":means["different"],
        "between_minus_within_state_patristic":means["different"]-means["same"],
        "exact_label_null_025":null[int(.025*(total-1))],
        "exact_label_null_975":null[int(.975*(total-1))],
        "leave_one_tip_out_rho_range":[min(loo),max(loo)],
        "leave_one_tip_out_direction_agreement":sum(r*rho>0 for r in loo)/len(loo),
    }

def test_synthetic():
    import io
    from Bio import Phylo
    t=Phylo.read(io.StringIO("((A:1,B:1):2,(C:1,D:1):2);"),"newick")
    q=statistic(t,["A","B","C","D"],[0,0,1,1])
    insist(abs(q["patristic_vs_discrete_mismatch_Spearman_rho"]-1)<1e-12,
           "SYNTHETIC_EXACT_TREE_NOT_CORRELATED")
    insist(q["n_fixed_prevalence_exact_assignments"]==6,
           "SYNTHETIC_4_TIP_ASSIGNMENTS_NOT_ENUMERATED")
    insist(abs(q["exact_two_sided_tip_label_permutation_p"]-1/3)<1e-12,
           "SYNTHETIC_EXACT_P_NOT_ENUMERATED")
    return {
        "status":"PASS_SYNTHETIC_EXACT_PERMUTATION_AND_TREE_DISTANCE_TESTS",
        "four_tip_orientation_exact_p":q["exact_two_sided_tip_label_permutation_p"],
        "biological_field_data_added":0,
    }

def main():
    p=argparse.ArgumentParser()
    for flag in ("tree","acceptance","concept-map","nmns-base","nmns-extension"):
        p.add_argument("--"+flag,type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    base,extension=read(a.nmns_base),read(a.nmns_extension)
    concepts=read(a.concept_map)
    acc=json.loads(a.acceptance.read_text(encoding="utf-8"))
    tree,orient,unresolved,sha=prepare(
        a.tree,acc,concepts,base,extension)
    by_id={x["paper_japan_member_id"]:x for x in concepts}
    src={x["paper_japan_member_id"]:x for x in (base+extension)}
    assert sha==EXPECTED_SHA
    modules={}
    for module,conf in TRAITS.items():
        subset=[]
        for cid,row in sorted(src.items()):
            state=row[conf["source"]]
            if state not in conf["allowed"]:
                continue
            cm=by_id[cid]
            if cm["trait_asr_primary_allowed"]!="true":
                continue
            code=cm["tip_ids"].split("|")
            insist(len(code)==1 and cm["n_biological_samples"]=="1",
                   "REPLICATED_NUCLEAR_TREE_CONCEPT_NOT_RESOLVED")
            subset.append((cid,code[0],state))
        counts=Counter(x[2] for x in subset)
        insist(dict(counts)==conf["expected_counts"],
               "FROZEN_CIRSIUM_SOURCE_STATE_COUNTS_CHANGED_"+module)
        out=statistic(tree,[x[1] for x in subset],[x[2] for x in subset])
        out.update({
           "status":("DESCRIPTIVE_EXACT_NULL_SINGLE_CANONICAL_TREE"
                     if len(subset)>=conf["minimum_taxa_for_inferential_test"]
                     else "EXPLORATORY_FEWER_THAN_TEN_UNAMBIGUOUS_SOURCE_TAXA"),
           "taxa":[{"concept":cid,"tip":tip,"state":state}
                   for cid,tip,state in subset],
           "trait_ontology":"Strict singleton NMNS categorical description; ambiguous sets not force coded",
        })
        modules[module]=out
    insist(modules["orientation"]["n_scored_taxa"]==19
           and modules["stickiness"]["n_scored_taxa"]==13
           and modules["phyllary_posture_STRICT_single_state_only"]["n_scored_taxa"]==6,
           "SOURCE_19_13_6_MODULE_COVERAGE_CHANGED")
    results={
      "version":"cirsium_Japan38_actual_three_head_module_memory_exact_null_v1",
      "as_of":"2026-10-19",
      "original_EAzami_commit":REF,
      "canonical_phylogram_SHA256":sha,
      "tree_source":"236-QC-locus Comp1061 compatibility reconstruction (not unavailable Moreyra final 350)",
      "branch_length_units":"substitutions/site, not dated evolutionary time",
      "synthetic_test":test_synthetic(),
      "actual_source_module_results":modules,
      "trait_module_primary_conclusion_ceiling":
       "Real within-Japan trait/distance signals only. Different source coverages and phylogeny uncertainty prohibit directly equating differences with evolutionary rates, trait-genetic resources or cross-family portability.",
      "external_head_state_memory_portability":
       "HOLD_11_EAST_ASIAN_CANDIDATES_WITH_INCOMPLETE_MATCHED_NUCLEAR_BRANCH_LENGTH_ENSEMBLE",
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(results,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
       module:{
         "n":v["n_scored_taxa"],
         "states":v["state_counts"],
         "rho":v["patristic_vs_discrete_mismatch_Spearman_rho"],
         "exact_two_sided_p":v["exact_two_sided_tip_label_permutation_p"],
         "assignments":v["n_fixed_prevalence_exact_assignments"],
         "status":v["status"],
       } for module,v in modules.items()
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
