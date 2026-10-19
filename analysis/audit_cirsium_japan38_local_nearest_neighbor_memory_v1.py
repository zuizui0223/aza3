#!/usr/bin/env python3
"""Out-of-tip Cirsium head-state predictability when GLOBAL phylogenetic rho≈0.

Check a distinct falsifiable within-Japan estimand: will a source-scored
taxon's closest OTHER taxon on the original nuclear phylogram forecast its
binary/multistate NMNS capitulum trait better than leave-one-tip-out (LOO)
marginal state frequencies? Use proper Brier scoring, not accuracy cherry-pick.

All fixed-prevalence label reallocations are enumerated EXACTLY on the
ORIGINAL Cirsium nuclear phylogram, 27,132 / 1,716 / 30 assignments. This
does not assess cross-lineage transfer or genetic retention/re-expression.
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
    prepare, EXPECTED_SHA, REF
)
from audit_cirsium_japan38_three_module_exact_memory_v1 import (
    TRAITS, exact_reassignments, test_synthetic
)

def insist(flag,code):
    if not flag: raise ValueError(code)

def read(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def design(tree,subset):
    tips=[x[1] for x in subset]
    n=len(tips)
    insist(n>=3,"TOO_FEW_TIPS_FOR_NEAREST_NEIGHBOR")
    distances=[[0. if i==j else float(tree.distance(tips[i],tips[j]))
                for j in range(n)] for i in range(n)]
    nearest=[]
    for i in range(n):
        eligible=[(d,j) for j,d in enumerate(distances[i]) if j!=i]
        smallest=min(d for d,j in eligible)
        insist(smallest>0 and math.isfinite(smallest),"MISSING_NUCLEAR_DISTANCE")
        tie=[j for d,j in eligible if abs(d-smallest)<=1e-11]
        nearest.append(tie)
    return nearest, distances

def half_multiclass_brier(y,probs,cats):
    """Normalized multiclass Brier, reduces to binary (y-p)^2."""
    insist(cats,"EMPTY_CATEGORY_UNIVERSE")
    return sum((probs.get(c,0.)-int(c==y))**2 for c in cats)/2

def two_LOO_briers(states,neighbors):
    categories=sorted(set(states))
    n=len(states)
    insist(len(categories)>1,"SOURCE_STATES_NO_VARIATION")
    details=[]
    nearest_scores=[]
    global_scores=[]
    nn_exact=0
    nn_ties=0
    for i,y in enumerate(states):
        others=[states[j] for j in range(n) if j!=i]
        neighbor_states=[states[j] for j in neighbors[i]]
        m=len(neighbor_states)
        if m>1:nn_ties+=1
        local=Counter(neighbor_states)
        baseline=Counter(others)
        p_nn={c:local[c]/m for c in categories}
        p_global={c:baseline[c]/(n-1) for c in categories}
        nnb=half_multiclass_brier(y,p_nn,categories)
        globb=half_multiclass_brier(y,p_global,categories)
        nearest_scores.append(nnb)
        global_scores.append(globb)
        nn_exact+=int(len(neighbor_states)==1 and neighbor_states[0]==y)
        details.append({
            "tip_index":i,"observed_category":y,
            "nn_prob_true_category":p_nn[y],
            "loo_global_prob_true_category":p_global[y],
            "nn_Brier":nnb,"global_Brier":globb,
            "number_equidistant_nearest_tips":len(neighbors[i]),
        })
    near=sum(nearest_scores)/n
    globalbase=sum(global_scores)/n
    return {
        "n_taxa":n,"n_NN_exact_correct_single_nearest":nn_exact,
        "n_taxa_with_equidistant_nearest_ties":nn_ties,
        "nearest_Brier":near,
        "leave_one_out_global_marginal_Brier":globalbase,
        "global_minus_NN_Brier_gain":globalbase-near,
        "correctness":"positive Brier gain means nearer phylogenetic neighbors predicted held-out tip state better",
        "per_tip_Brier":details,
    }

def exact_assessment(tree,subset):
    near,dist=design(tree,subset)
    states=[x[2] for x in subset]
    observed=two_LOO_briers(states,near)
    obs=observed["global_minus_NN_Brier_gain"]
    assignments=0
    null_ge=0
    null_abs=0
    null=[]
    for assignment in exact_reassignments(states,Counter(states)):
        gain=two_LOO_briers(assignment,near)["global_minus_NN_Brier_gain"]
        assignments+=1
        null_ge+=gain>=obs-1e-12
        null_abs+=abs(gain)>=abs(obs)-1e-12
        null.append(gain)
    null.sort()
    known=Counter(states)
    assert assignments==math.factorial(len(states))//math.prod(
        math.factorial(v) for v in known.values()
    )
    observed.update({
        "n_exact_fixed_label_allocations":assignments,
        "exact_one_sided_positive_predictive_gain_p":null_ge/assignments,
        "absolute_gain_extremeness_diagnostic_not_primary_p":
            null_abs/assignments,
        "null_gain_025":null[int(.025*(assignments-1))],
        "null_gain_median":null[assignments//2],
        "null_gain_975":null[int(.975*(assignments-1))],
        "taxa":[{"taxon":x[0],"tip":x[1],"state":x[2],
                 "nearest_taxa":[subset[j][0] for j in near[i]],
                 "nearest_distances_substitutions_per_site":[dist[i][j]
                                                            for j in near[i]]}
                for i,x in enumerate(subset)],
        "source_not_repeated_biological_individuals":True,
        "source_taxon_availability_is_nonrandom":True,
        "not_cross_clade_portability":True,
    })
    return observed

def synthetic_tests():
    import io
    from Bio import Phylo
    t=Phylo.read(io.StringIO("((A:1,B:1):2,(C:1,D:1):2);"),"newick")
    x=exact_assessment(t,[("A","A","U"),("B","B","U"),
                          ("C","C","D"),("D","D","D")])
    assert x["n_exact_fixed_label_allocations"]==6
    assert abs(x["nearest_Brier"])<1e-12
    assert x["leave_one_out_global_marginal_Brier"]>0
    assert abs(x["exact_one_sided_positive_predictive_gain_p"]-1/3)<1e-12
    # Deliberately reverse the one tip; genuine local predictability should
    # not be automatically manufactured by nearest-neighbor scoring.
    y=two_LOO_briers(["U","D","U","D"],design(
        t,[("A","A","U"),("B","B","D"),
           ("C","C","U"),("D","D","D")])[0])
    assert y["global_minus_NN_Brier_gain"]<0
    return {
        "status":"PASS_SYNTHETIC_NEAREST_NEIGHBOR_CV_AND_EXACT_FIXED_LABEL_NULL",
        "positive_four_tip_gain":x["global_minus_NN_Brier_gain"],
        "bad_four_tip_gain":y["global_minus_NN_Brier_gain"],
        "negative_or_missing_real_study_result_inferred":False
    }

def main():
    ap=argparse.ArgumentParser()
    for arg in ("tree","acceptance","concept-map","nmns-base","nmns-extension"):
        ap.add_argument("--"+arg,type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    base,extension=read(a.nmns_base),read(a.nmns_extension)
    maps=read(a.concept_map)
    acc=json.loads(a.acceptance.read_text(encoding="utf-8"))
    tree,_,unknown,sha=prepare(a.tree,acc,maps,base,extension)
    insist(sha==EXPECTED_SHA,"NUCLEAR_TREE_SOURCE_HASH_DRIFT")
    indexed={r["paper_japan_member_id"]:r for r in maps}
    source={r["paper_japan_member_id"]:r for r in base+extension}
    modules={}
    for name,conf in TRAITS.items():
        subset=[]
        for cid,row in sorted(source.items()):
            state=row[conf["source"]]
            if state not in conf["allowed"]:continue
            source_tax=indexed[cid]
            if source_tax["trait_asr_primary_allowed"]!="true":continue
            tips=source_tax["tip_ids"].split("|")
            insist(len(tips)==1 and source_tax["n_biological_samples"]=="1",
                   "INDEPENDENT_TIP_SET_UNKNOWN")
            subset.append((cid,tips[0],state))
        insist(dict(Counter(x[2] for x in subset))==
               conf["expected_counts"],"SOURCE_NUCLEAR_CIRSIUM_CATEGORIES_CHANGED")
        result=exact_assessment(tree,subset)
        result["status"]=(
            "EXPLORATORY_LOW_TAXON_COVERAGE"
            if len(subset)<10
            else "DESCRIPTIVE_SOURCE_TAXON_CV_SINGLE_TREE"
        )
        modules[name]=result
    receipt={
        "version":"cirsium_Japan38_same_trait_local_predictability_v1",
        "date":"2026-10-19",
        "source_repository":"zuizui0223/EAzami",
        "source_commit":REF,
        "canonical_source_sha256":sha,
        "type":"Real taxon-level nuclear nearest-neighbor prediction on original Japanese ML tree",
        "null":"all fixed-prevalence trait label reallocations, exact enumeration",
        "test_synthetic":synthetic_tests(),
        "real_Cirsium_results":modules,
        "interpretation_ceiling":[
            "All metrics are within-Japan single-tree leave-one-tip-out, not held-out geographic lineages or transfer across taxa families.",
            "Nearest-neighbor comparison is source-selected by categorical state availability.",
            "Shared phylogenetic branches/tips are correlated; exact label null addresses this within selected tree, not unmodeled gene tree heterogeneity.",
            "A negative nearest-neighbor gain cannot prove no latent ancestral genetic toolkit or re-expression.",
            "3 traits have different valid taxon coverage and phyllary posture n=6 exploratory.",
        ],
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",
                     encoding="utf-8")
    print(json.dumps({
       k:{"n":v["n_taxa"],
          "nn_Brier":v["nearest_Brier"],
          "LOO_baseline_Brier":v["leave_one_out_global_marginal_Brier"],
          "Brier_gain":v["global_minus_NN_Brier_gain"],
          "p_positive":v["exact_one_sided_positive_predictive_gain_p"],
          "n_assignments":v["n_exact_fixed_label_allocations"]}
       for k,v in modules.items()
    },ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
