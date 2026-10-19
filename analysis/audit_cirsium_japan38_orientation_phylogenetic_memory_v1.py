#!/usr/bin/env python3
"""Source-backed Japan38 capitulum orientation phylogenetic-memory diagnostic.

A real Comp1061 236-QC-locus *compatibility reconstruction* with branch
lengths in substitutions/site is paired with 19 strict NMNS authority-coded
orientations. EAzami source commit c9870e91604ed1514344ddc8f765eae5237dbf27.

Q1 within-Japan question only: do opposite upward/downward taxon labels
occur farther apart on this SPECIFIC nuclear phylogram than same-state tips?

Spearman(pairwise patristic distance, binary orientation dissimilarity),
compared with fixed-prevalence label permutations on the fixed tree.
Pairs are NONINDEPENDENT: never compute a naive n_pairs-based p-value.
Permutation p accounts for shared tips, conditional on this tree and 19
selected authority states. Does NOT test portability to other lineages,
OU half-life, dated evolutionary memory, hereditary trait causation, or
the evolution of partner-dependent head fitness.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter
from pathlib import Path

from Bio import Phylo

REF="c9870e91604ed1514344ddc8f765eae5237dbf27"
STATES={"upward_or_erect":1,"downward_or_nodding":0}
EXPECTED_SHA="ba592d3a9f52b5481314ab20e3202048505cdf72a3fc15134d5dca261b9d9b7e"
PERM_SEED=20261019

def insist(flag,code):
    if not flag:raise ValueError(code)

def csvrows(path,needed):
    with open(path,newline="",encoding="utf-8-sig") as f:
        r=csv.DictReader(f)
        insist(set(needed).issubset(r.fieldnames or []),
               "REAL_CIRSIUM_SOURCE_SCHEMA_DRIFT")
        return list(r)

def unique(rows,key,code):
    found={}
    for item in rows:
        v=(item.get(key) or "").strip()
        insist(v and v not in found,code)
        found[v]=item
    return found

def rank_average(values):
    values=[float(x) for x in values]
    order=sorted(range(len(values)),key=lambda i:values[i])
    res=[0.0]*len(values)
    n=len(values)
    i=0
    while i<n:
        j=i+1
        while j<n and values[order[j]]==values[order[i]]:
            j+=1
        rank=(i+j+1)/2  # 1-based average rank of i..j-1
        for k in range(i,j):res[order[k]]=rank
        i=j
    return res

def spearman_rank_with_binary(ranked_dist, comparison):
    n=len(comparison)
    if n<2:return None
    mu=sum(ranked_dist)/n
    p=sum(comparison)/n
    xx=sum((v-mu)**2 for v in ranked_dist)
    yy=sum((v-p)**2 for v in comparison)
    if xx<=0 or yy<=0:return None
    return sum((v-mu)*(y-p) for v,y in zip(ranked_dist,comparison))/math.sqrt(xx*yy)

def test_synthetic_geometry():
    import io
    simple=Phylo.read(io.StringIO("((A:1,B:1):2,(C:1,D:1):2);"),"newick")
    dd=[simple.distance(a,b) for a,b in
        [("A","B"),("A","C"),("A","D"),("B","C"),("B","D"),("C","D")]]
    labels=[0,1,1,1,1,0]
    rho=spearman_rank_with_binary(rank_average(dd),labels)
    insist(abs(rho-1)<1e-12,"PAIRWISE_RANK_MEMORY_FIXTURE_FAILED")
    insist(dd==[2.0,6.0,6.0,6.0,6.0,2.0],
           "PATRISTIC_DISTANCE_SYNTHETIC_FIXTURE_FAILED")
    insist(rank_average([1,1,2,2])==[1.5,1.5,3.5,3.5],
           "TIED_RANKS_NOT_AVERAGED")
    return {"status":"PASS_SYNTHETIC_TREE_RANK_AND_SHARED_TIP_PERMUTATION_FIXTURES",
            "exact_synthetic_tree_rho":rho,
            "real_biological_data_created":0}

def prepare(tree_path,tree_manifest,concept_rows,seed,extension):
    raw=Path(tree_path).read_bytes()
    sha=hashlib.sha256(raw).hexdigest()
    insist(sha==EXPECTED_SHA and
           tree_manifest["tree_sha256"]==EXPECTED_SHA,
           "CANONICAL_JAPAN38_TREE_HASH_MISMATCH")
    insist(tree_manifest["branch_length_semantics"]==
           "substitutions/site; not absolute time",
           "TREE_LENGTH_WRONGLY_ASSUMED_CALENDAR_TIME")
    insist(tree_manifest["focal_biological_samples"]==39
           and tree_manifest["tree_tip_count"]==40
           and tree_manifest["root_tip"]=="OUTGROUP_saff",
           "SOURCE_TREE_ACCEPTANCE_CHANGED")
    tree=Phylo.read(str(tree_path),"newick")
    terminals={t.name:t for t in tree.get_terminals()}
    insist(len(terminals)==40 and "OUTGROUP_saff" in terminals,
           "ORIGINAL_TREE_TIP_LIST_NOT_RECOVERED")
    concepts=unique(concept_rows,"paper_japan_member_id",
                    "CANONICAL_CONCEPT_ID_DUPLICATE_OR_EMPTY")
    insist(len(concepts)==38,"JAPAN38_CONCEPT_REGISTRY_NOT_38")
    scored=unique(seed+extension,"paper_japan_member_id",
                  "BOTANICAL_SOURCE_SCORED_DUPLICATE_ID")
    insist(len(seed)==22 and len(extension)==2 and len(scored)==24,
           "FROZEN_BOTANICAL_AUTHORITY_SCOPE_DRIFT")

    items=[]
    other_states=Counter()
    for cid,row in sorted(scored.items()):
        insist(cid in concepts,"BOTANICAL_TAXON_NOT_IN_SOURCE_TREE")
        insist(row["authority_match_status"]==
               "exact_authority_concept_match",
               "BOTANICAL_AUTHORITY_SOURCE_NOT_EXACT")
        st=row["orientation_state"]
        if st not in STATES:
            other_states[st]+=1
            continue
        cm=concepts[cid]
        insist(cm["trait_asr_primary_allowed"]=="true",
               "TAXON_HAS_BLOCKED_IDENTITY_OR_ANCESTRY")
        tips=cm["tip_ids"].split("|")
        insist(len(tips)==1 and cm["n_biological_samples"]=="1",
               "DUPLICATE_NUCLEAR_INDIVIDUALS_FOR_SINGLE_CONCEPT")
        insist(tips[0] in terminals,"AUTHORITATIVE_SOURCE_TIP_MISSING_FROM_TREE")
        items.append({
            "taxon_id":cid,"tip":tips[0],
            "state":st,"y":STATES[st],
            "source_scope":row["source_scope"],
            "source_url":row["source_url"],
        })
    counts=Counter(z["state"] for z in items)
    insist(len(items)==19
           and counts=={"upward_or_erect":13,"downward_or_nodding":6},
           "STRICT_19_SOURCE_ORIENTATION_BALANCE_CHANGED")
    insist("JPN_20" not in {r["taxon_id"] for r in items}
           and "JPN_31" not in {r["taxon_id"] for r in items},
           "SOURCE_NONMONOPHYLY_OR_TAXON_IDENTITY_CONFLICT_INCLUDED")
    insist(other_states["upward_or_ascending"]==1
           and other_states["source_conflict_index_downward_detail_erect"]==1,
           "ORIENTATION_AUTHORITY_CONFLICT_DISAPPEARED_OR_FORCED")
    # The exact number of source-unknown labels may change with
    # documentation additions; preserve rather than recode.
    return tree,items,dict(other_states),sha

def run_statistic(tree,items,permutations,seed):
    insist(1000<=permutations<=100000,"PERMUTATIONS_REPRODUCIBILITY_RANGE")
    tip=[x["tip"] for x in items]
    y=[x["y"] for x in items]
    n=len(y)
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    distances=[float(tree.distance(tip[i],tip[j])) for i,j in pairs]
    insist(all(math.isfinite(d) and d>0 for d in distances),
           "TREE_HAS_MISSING_OR_ZERO_PAIRWISE_BRANCH_LENGTH")
    rdist=rank_average(distances)
    mismatch=[int(y[i]!=y[j]) for i,j in pairs]
    insist(sum(mismatch)==13*6,"ORIENTATION_DISCORDANT_PAIR_BALANCE_DRIFT")
    rho=spearman_rank_with_binary(rdist,mismatch)
    insist(rho is not None,"UNCALCULABLE_EMPIRICAL_MEMORY_RHO")
    same=[d for d,m in zip(distances,mismatch) if not m]
    opposite=[d for d,m in zip(distances,mismatch) if m]
    mean_same=sum(same)/len(same)
    mean_opposite=sum(opposite)/len(opposite)

    rng=random.Random(seed)
    shuffled=y[:]
    null=[]
    ge_two=ge_pos=0
    for _ in range(permutations):
        rng.shuffle(shuffled)
        perm=[int(shuffled[i]!=shuffled[j]) for i,j in pairs]
        t=spearman_rank_with_binary(rdist,perm)
        insist(t is not None,"NULL_PERMUTATION_RHO_CONSTANT")
        null.append(t)
        ge_two+=abs(t)>=abs(rho)-1e-14
        ge_pos+=t>=rho-1e-14
    null.sort()
    loo=[]
    for leave in range(n):
        keep=[k for k,(i,j) in enumerate(pairs)
              if i!=leave and j!=leave]
        q=spearman_rank_with_binary(
            rank_average([distances[k] for k in keep]),
            [mismatch[k] for k in keep])
        loo.append({"excluded_taxon":items[leave]["taxon_id"],
                    "rho":q,"remaining_taxa":n-1})
    return {
        "status":"REAL_JAPAN38_ORIENTATION_PHYLOGRAM_DISTANCE_ASSOCIATION_NOT_LINEAGE_TRANSFER",
        "quantities":{
            "number_taxa":n,
            "n_up":sum(y),
            "n_down":n-sum(y),
            "n_nonindependent_pairs":len(pairs),
            "n_discordant_pairs":len(opposite),
            "n_same_state_pairs":len(same),
            "mean_patristic_same_orientation_substitutions_per_site":mean_same,
            "mean_patristic_opposite_orientation_substitutions_per_site":mean_opposite,
            "opposite_minus_same_patristic_distance":mean_opposite-mean_same,
            "spearman_patristic_vs_binary_orientation_difference":rho,
            "fixed_state_prevalence_permutation_two_sided_p":
                (ge_two+1)/(permutations+1),
            "fixed_state_prevalence_permutation_positive_tail_p":
                (ge_pos+1)/(permutations+1),
            "permutations":permutations,
            "seed":seed,
            "null_rho_025":null[int(.025*(permutations-1))],
            "null_rho_median":null[permutations//2],
            "null_rho_975":null[int(.975*(permutations-1))],
            "leave_one_taxon_out_rho_min":min(x["rho"] for x in loo),
            "leave_one_taxon_out_rho_max":max(x["rho"] for x in loo),
            "leave_one_taxon_out_same_sign_fraction":
                sum(x["rho"]*rho>0 for x in loo)/len(loo),
        },
        "leave_one_out":loo,
        "interpretation_rules":[
            "Positive rho = opposite orientation labels tend to be separated by larger nuclear phylogram distance, on this frozen tree only.",
            "Negative rho = opposite labels tend to occur at shorter nuclear phylogram distance, on this frozen tree only.",
            "No rho value estimates memory across another family or independently dated evolutionary persistence.",
            "Permutation shuffles only observed labels on fixed, unevenly sampled tips, not phylogenetic branch lengths.",
            "N=19 taxon labels, NOT 171 independently replicated evolutionary transitions.",
            "One reconstructed ML topology and unresolved head angle category measurement bias remain potential sensitivities.",
            "Multiple labelled tips can share past species histories, including hybridization and incomplete lineage sorting.",
        ],
        "not_a_test_of":["head orientation functional advantage",
                         "true spine or insect contact",
                         "ancestral pigment genetic retention",
                         "adaptation or causal evolutionary return",
                         "clade-to-clade transfer",
                         "OU phylogenetic memory half-life"],
    }

def main():
    p=argparse.ArgumentParser()
    for key in ("tree","acceptance","concept-map","nmns-base","nmns-extension"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--permutations",type=int,default=20000)
    p.add_argument("--seed",type=int,default=PERM_SEED)
    p.add_argument("--out",type=Path)
    args=p.parse_args()
    acc=json.loads(args.acceptance.read_text(encoding="utf-8"))
    cm=csvrows(args.concept_map,
               {"paper_japan_member_id","tip_ids","trait_asr_primary_allowed",
                "n_biological_samples"})
    base=csvrows(args.nmns_base,
                 {"paper_japan_member_id","orientation_state","source_scope",
                  "source_url","authority_match_status"})
    extension=csvrows(args.nmns_extension,
                      {"paper_japan_member_id","orientation_state","source_scope",
                       "source_url","authority_match_status"})
    tree,items,unknown,sha=prepare(
        args.tree,acc,cm,base,extension)
    report={
        "analysis_id":"cirsium_japan38_authentic_orientation_memory_v1",
        "as_of":"2026-10-19",
        "source_repository":"zuizui0223/EAzami",
        "source_ref":REF,
        "original_frozen_tree_sha256":sha,
        "tip_definition":"Strict NMNS binary anatomical capitulum orientation codes, NOT image vertical",
        "branch_length_unit":"substitutions/site, NOT calendar years",
        "source_discarded_other_orientation_labels":unknown,
        "source_scored_taxa":items,
        "synthetic_tests":test_synthetic_geometry(),
        "empirical_Japan38_statistic":run_statistic(
            tree,items,args.permutations,args.seed),
        "external_East_Asia_memory_transfer":"HOLD_MISSING_MATCHED_BRANCH_LENGTH_TOPOLOGY",
    }
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",
                            encoding="utf-8")
    q=report["empirical_Japan38_statistic"]["quantities"]
    print(json.dumps({
        "status":report["empirical_Japan38_statistic"]["status"],
        "n_taxa":q["number_taxa"],
        "nonindependent_pairs":q["n_nonindependent_pairs"],
        "rho":q["spearman_patristic_vs_binary_orientation_difference"],
        "two_sided_tip_label_permutation_p":
            q["fixed_state_prevalence_permutation_two_sided_p"],
        "leave_one_out_rho_range":[
            q["leave_one_taxon_out_rho_min"],
            q["leave_one_taxon_out_rho_max"]],
        "external_memory_test":"HOLD"
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
