#!/usr/bin/env python3
"""Source-authenticated North-American Cirsium capitulum joint-evolution audit.

Two *observed* head components: erect/pending head position, and
ascending vs spreading/reflexed OUTER phyllaries (NOT true spines).
Compare a nested independent binary-CTMC model to an eight-rate dependent
four-state CTMC on the *same pinned ultrametric tree*; also test the actual
joint overlap against conditional phylogenetic Mk simulations fixing the
other measured phenotype and each observed marginal frequency.

Tree/taxon coding: Siniscalchi et al. 2023 International Journal of Plant
Sciences, DOI 10.1086/724310; original public author GitHub commit
c925d4af86160064448f794a9715f70b917a2786. Not independent new
field observations; four-state rate differences do not establish selection.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import numpy as np
from Bio import Phylo
from scipy.linalg import expm
from scipy.optimize import minimize

AUTHOR="carol-siniscalchi/cirsium"
COMMIT="c925d4af86160064448f794a9715f70b917a2786"
BLOBS={
    "traits.csv":"fa127f0587b676e6e004ea121424d26069586bcb",
    "cirsium.phyluce95p.ultrametric.tre":"6417a38793427bc7c8562dc0e2cab59b74d8c316",
}
SEED=20261021
STATES=("erect_not_outward", "erect_outward", "pending_not_outward", "pending_outward")
EDGES=((0,1),(1,0),(2,3),(3,2),(0,2),(2,0),(1,3),(3,1))
# Parameters independent: outward gain, outward loss, head pending gain,
# head pending loss (branch length units = substitutions/site).
INDEPENDENT_EDGE_PARAMETER=(0,1,0,1,2,3,2,3)


def require(condition, label):
    if not condition: raise ValueError(label)


def source_blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()


def fetch_local_exact(path,name):
    raw=Path(path).read_bytes()
    require(source_blob(raw)==BLOBS[name],"SOURCE_GIT_BLOB_DRIFT_"+name)
    return raw.decode("utf-8-sig")


def parse_source(trait_csv,tree_newick):
    rows=list(csv.DictReader(io.StringIO(trait_csv)))
    require(len(rows)==63, "SOURCE_N63_CHANGED")
    names=[r.get("Species") for r in rows]
    require(len(set(names))==63 and all(names),"DUPLICATED_OR_MISSING_TIP")
    states={}
    for r in rows:
        require(r.get("HeadPosition") in {"0","1"},"SOURCE_HEAD_STATE_DRIFT")
        require(r.get("OuterPhyllaryPosition") in {"0","1","2"},"SOURCE_OUTER_STATE_DRIFT")
        states[r["Species"]]=(int(r["HeadPosition"]),int(r["OuterPhyllaryPosition"]!="0"))
    tree=Phylo.read(io.StringIO(tree_newick),"newick")
    tips=[t.name for t in tree.get_terminals()]
    require(len(tips)==63 and len(set(tips))==63 and set(tips)==set(names),"SOURCE_TREE_TIPS_DRIFT")
    depths=tree.depths()
    heights=[depths[t] for t in tree.get_terminals()]
    require(max(heights)-min(heights)<.001,"SOURCE_NOT_ULTRAMETRIC")
    for z in tree.find_clades():
        require(z.branch_length is None or
                (math.isfinite(z.branch_length) and z.branch_length>=0),"INVALID_PHYLOGENY_BRANCH")
    counted=np.bincount([2*h+o for h,o in states.values()],minlength=4).tolist()
    require(counted==[45,13,2,3],"SOURCE_HEAD_PAIR_COUNTS_DRIFT")
    require(states["Cirsium_vulgare"]==(0,0),"ORIGINAL_OUTGROUP_HEAD_STATES_DRIFT")
    return tree,states,{"tip_counts":dict(zip(STATES,counted)),
                        "n_pending":5,"n_outer_outward":16,
                        "observed_pending_and_outward":3,
                        "tree_depth_range":[min(heights),max(heights)]}


def binary_loglik(tree,states,gain,loss,root_state=0):
    require(root_state in (0,1) and gain>0 and loss>0,"INVALID_BINARY_PROCESS")
    r=gain+loss
    pi1=gain/r
    pi0=loss/r
    cache={}
    for n in tree.find_clades(order="postorder"):
        if n.is_terminal():
            require(n.name in states and states[n.name] in (0,1),"INVALID_MARGINAL_SOURCE_STATE")
            cache[n]=np.array([1.,0.] if states[n.name]==0 else [0.,1.])
        else:
            v=np.ones(2)
            for ch in n.clades:
                ee=math.exp(-r*(ch.branch_length or 0.))
                P=np.array([[pi0+pi1*ee,pi1*(1-ee)],[pi0*(1-ee),pi1+pi0*ee]])
                v *= P @ cache[ch]
            cache[n]=v
    p=cache[tree.root][root_state]
    return float(math.log(p)) if p>0 else -math.inf


def fit_binary(tree,states,root_state=0):
    def obj(v):
        z=binary_loglik(tree,states,math.exp(v[0]),math.exp(v[1]),root_state)
        return -z if math.isfinite(z) else 1e100
    starts=[(-2.,-1.),(-4.,0.),(-1.,-1.),(0.,0.)]
    fits=[minimize(obj,np.array(start),method="L-BFGS-B",bounds=[(-9,3),(-9,3)],
                   options={"maxiter":220,"ftol":5e-11}) for start in starts]
    success=[f for f in fits if f.success and math.isfinite(f.fun)]
    require(success,"BINARY_MK_FIT_FAILED")
    best=min(success,key=lambda z:z.fun)
    return {"root_state":root_state,"gain":math.exp(best.x[0]),
            "loss":math.exp(best.x[1]),"loglik":float(-best.fun),
            "converged":True,"rates_on_bounds_possible":bool(
              np.any(abs(best.x-np.array([-9.,-9.]))<.001) or
              np.any(abs(best.x-np.array([3.,3.]))<.001))}


def q_matrix(rates,dependent=False):
    rates=np.asarray(rates,dtype=float)
    require(rates.shape==((8,) if dependent else (4,)),"INVALID_RATE_DIMENSIONS")
    require(bool(np.all(np.isfinite(rates))) and bool(np.all(rates>0)),"INVALID_POSITIVE_RATES")
    Q=np.zeros((4,4),dtype=float)
    for i,(a,b) in enumerate(EDGES):
        Q[a,b]=rates[i] if dependent else rates[INDEPENDENT_EDGE_PARAMETER[i]]
    np.fill_diagonal(Q,-Q.sum(axis=1))
    require(bool(np.allclose(Q.sum(axis=1),0,atol=1e-12)),"CTMC_ROW_SUM_BROKEN")
    return Q


def four_loglik(tree,states,rates,dependent=False,root_state=0):
    require(root_state in range(4),"INVALID_FOUR_ROOT_STATE")
    Q=q_matrix(rates,dependent)
    # Numerical likelihood scaling reduces extreme-rate underflow.
    likelihood={}
    logscale=0.
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            h,o=states[node.name]
            v=np.zeros(4)
            v[2*h+o]=1.
            likelihood[node]=v
        else:
            v=np.ones(4)
            for child in node.clades:
                P=expm(Q*(child.branch_length or 0.))
                v *= P @ likelihood[child]
            s=float(v.sum())
            if s<=0 or not math.isfinite(s):return -math.inf
            likelihood[node]=v/s
            logscale+=math.log(s)
    p=likelihood[tree.root][root_state]
    return float(math.log(p)+logscale) if p>0 else -math.inf


def fit_pair_models(tree,states,orientation_fit,outer_fit):
    indrates=[outer_fit["gain"],outer_fit["loss"],
              orientation_fit["gain"],orientation_fit["loss"]]
    ind= four_loglik(tree,states,indrates,dependent=False)
    separate=orientation_fit["loglik"]+outer_fit["loglik"]
    require(math.isclose(ind,separate,abs_tol=1e-8),"INDEPENDENT_JOINT_LIKELIHOOD_IDENTITY_BROKEN")
    def expand(z):return np.array([z[0],z[1],z[0],z[1],z[2],z[3],z[2],z[3]])
    starts=[np.log(expand(indrates)),np.log(expand(indrates))+np.array([.2,-.2,.15,-.15,.2,-.2,.15,-.15])]
    def obj(x):
        ll=four_loglik(tree,states,np.exp(x),dependent=True)
        return -ll if math.isfinite(ll) else 1e100
    f=[minimize(obj,start,method="L-BFGS-B",bounds=[(-9,3)]*8,
                options={"maxiter":250,"ftol":3e-10}) for start in starts]
    success=[z for z in f if z.success and math.isfinite(z.fun)]
    require(success,"DEPENDENT_MK_FIT_FAILED")
    best=min(success,key=lambda z:z.fun)
    dependent_ll=max(ind, -float(best.fun)) # never assert improvement from a bad local optimum
    require(dependent_ll>=ind-1e-7,"NESTED_DEPENDENT_LOWER_LIKELIHOOD")
    lr=2*(dependent_ll-ind)
    aicind=-2*ind+2*4
    aicdep=-2*dependent_ll+2*8
    return {"root_state":"erect_not_outward_FIXED_00",
            "independent":{"n_rates":4,"rate_order":["outer_gain","outer_loss","head_pending_gain","head_pending_loss"],"rates":indrates,"loglik":ind,"AIC":aicind},
            "dependent":{"n_rates":8,"rate_order":["erect_close_to_outer","erect_outer_to_close","pending_close_to_outer","pending_outer_to_close","erect_close_to_pending","pending_close_to_erect","erect_outer_to_pending","pending_outer_to_erect"],"rates":list(map(float,np.exp(best.x))),"loglik":dependent_ll,"AIC":aicdep,"optimizer_converged":True,"some_rates_at_search_bound":bool(np.any(np.abs(best.x+9)<.015)|np.any(np.abs(best.x-3)<.015))},
            "lrt_descriptive":lr,"delta_AIC_dependent_minus_independent":aicdep-aicind,
            "chi_square_lrt_reference_p_omitted_due_to_sparse_cells_and_boundary":True,
            "parametric_bootstrap_lrt_not_performed":True,
            "na_63_code_tips_are_not_phylogenetically_independent_replication":True}


def conditional_tip_overlap(tree,var_states,fixed_states,gain,loss,n_sims=200000,
                            seed=SEED,root_state=0,min_accepted=300):
    """One trait evolves on original tree; other original trait vector fixed.

    Conditions simulated marginal prevalence to equal observed. This is a
    source-topology and root/rate-assumption sensitivity, not evidence of
    causal selection or covariance in real functional head components.
    """
    names=[x.name for x in tree.get_terminals()]
    require(set(names)==set(var_states)==set(fixed_states),"CTMC_TIP_NAMES_MISMATCH")
    require(root_state in (0,1),"INVALID_NULL_ROOT")
    k=int(sum(var_states.values()))
    require(0<k<len(names) and 1000<=n_sims<=1000000,"INVALID_SIMULATION_REQUEST")
    require(gain>0 and loss>0,"NONPOSITIVE_MARGINAL_RATES")
    rng=np.random.default_rng(seed)
    states={tree.root:np.full(n_sims,root_state,dtype=np.uint8)}
    parents={child:n for n in tree.find_clades(order="preorder") for child in n.clades}
    r=gain+loss;pi1=gain/r;pi0=loss/r
    for node in tree.find_clades(order="preorder"):
        if node is tree.root:continue
        parent=parents[node]
        e=math.exp(-r*(node.branch_length or 0))
        p01=pi1*(1-e);p11=pi1+pi0*e
        u=rng.random(n_sims)
        states[node]=np.where(states[parent]>0,u<p11,u<p01).astype(np.uint8)
    mat=np.vstack([states[x] for x in tree.get_terminals()])
    mask=(mat.sum(axis=0)==k)
    n=int(mask.sum())
    require(n>=min_accepted,"INSUFFICIENT_CONDITIONAL_MK_ACCEPTED_DRAWS")
    fixed=np.array([fixed_states[t] for t in names],dtype=np.uint8)
    overlap=(fixed[:,None]*mat[:,mask]).sum(axis=0)
    obs=sum(int(var_states[t])*int(fixed_states[t]) for t in names)
    h=np.bincount(overlap.astype(np.int32),minlength=min(k,int(fixed.sum()))+1)
    upper=int((overlap>=obs).sum());lower=int((overlap<=obs).sum())
    pup=(upper+1)/(n+1);plow=(lower+1)/(n+1)
    return {"observed_overlap":obs,"observed_var_positives":k,
            "observed_fixed_positives":int(fixed.sum()),
            "draws":n_sims,"accepted_exact_marginal":n,
            "root_state":root_state,"seed":seed,
            "null_overlap_mean":float(overlap.mean()),
            "null_overlap_sd":float(overlap.std()),
            "null_overlap_quantile_025":float(np.quantile(overlap,.025)),
            "null_overlap_quantile_975":float(np.quantile(overlap,.975)),
            "overlap_histogram":list(map(int,h)),
            "p_greater_equal":pup,"p_less_equal":plow,
            "two_inclusive_tails_doubled_min":min(1.,2*min(pup,plow)),
            "Monte_Carlo_SE_for_greater_equal":math.sqrt(pup*(1-pup)/(n+1)),
            "condition":"PRESERVE_TREE_AND_OTHER_TRAIT_TIP_PATTERN_AND_SIMULATED_TRAIT_POSITIVE_COUNT",
            "not_a_test_of_causal_head_function_or_genetic_independence":True}


def audit(tree,states,n_sims):
    orient={n:ho[0] for n,ho in states.items()}
    outer={n:ho[1] for n,ho in states.items()}
    source={"n_tips":len(states),"n_head_pending":sum(orient.values()),
            "n_outer_outward":sum(outer.values()),
            "n_joint_positive":sum(orient[n]*outer[n] for n in states)}
    independent={}
    for root in (0,1):
        mo=fit_binary(tree,orient,root)
        mb=fit_binary(tree,outer,root)
        independent[str(root)]={"head_orientation_Mk":mo,"outer_phyllary_Mk":mb,
             "head_evolved_outer_fixed":conditional_tip_overlap(tree,orient,outer,mo["gain"],mo["loss"],n_sims,seed=SEED+10+root,root_state=root),
             "outer_evolved_head_fixed":conditional_tip_overlap(tree,outer,orient,mb["gain"],mb["loss"],n_sims,seed=SEED+20+root,root_state=root)}
    combined=fit_pair_models(tree,states,independent["0"]["head_orientation_Mk"],independent["0"]["outer_phyllary_Mk"])
    return {"observed":source,
            "independence_reference":independent,
            "four_state_joint_model_comparison":combined,
            "interpretation":{
               "primary_test":"Tree CTMC parametric conditional null conditioned on observed count of one marginal tip state and fixing other marginal tip map",
               "root_model_sensitivity":"Root fixed to 0 and to 1 for each binary component separately. Joint four-state root 00 only.",
               "AIC_is_not_a_selection_causality_test":True,
               "observed_cell_counts_from_source":{"erect_not_outward":45,"erect_outward":13,"pending_not_outward":2,"pending_outward":3},
               "real_source_head_measurements_are_species_codes_not_same_individuals":True,
               "no_independent_genetic_module_mapping":True,
               "no_pollinator_or_seed_predator_benefits_measured":True,
               "original_phylogeny_uncertainty_not_propagated":True,
               "source_tip_taxon_sampling_nonrandom":True,
               "evolutionary_ancestral_reexpression_unidentifiable":True}}


def collapse_duplicated_arizonicum_subspecies(tree,states):
    """Taxonomically conservative tip sensitivity, not a new independent tree."""
    sibling=("Cirsium_arizonicum_arizonicum",
             "Cirsium_arizonicum_rothrockii")
    require(all(s in states for s in sibling),"ARIZONICUM_SOURCE_TIPS_MISSING")
    require(states[sibling[0]]==states[sibling[1]]==(0,1),
            "ARIZONICUM_SOURCE_MORPHOTYPE_DRIFT")
    t=copy.deepcopy(tree)
    t.prune(target=sibling[1])
    s={k:v for k,v in states.items() if k!=sibling[1]}
    require(len(t.get_terminals())==len(s)==62 and
            set(z.name for z in t.get_terminals())==set(s),
            "PRUNED_TIP_ALIGNMENT_FAILED")
    require(np.bincount([2*h+o for h,o in s.values()],minlength=4).tolist()==
            [45,12,2,3],"PRUNED_HEAD_COUNT_DRIFT")
    return t,s


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--traits",type=Path,required=True)
    ap.add_argument("--tree",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--simulations",type=int,default=200000)
    ap.add_argument("--subspecies-simulations",type=int,default=120000)
    args=ap.parse_args()
    traits=fetch_local_exact(args.traits,"traits.csv")
    newick=fetch_local_exact(args.tree,"cirsium.phyluce95p.ultrametric.tre")
    tree,tipstates,sinfo=parse_source(traits,newick)
    pruned_tree,pruned_states=collapse_duplicated_arizonicum_subspecies(
        tree,tipstates)
    reduced=audit(pruned_tree,pruned_states,args.subspecies_simulations)
    result={"version":"cirsium_na_capitulum_orientation_outer_mk_independence_v1",
            "as_of":"2026-10-21",
            "source":{"repo":AUTHOR,"commit":COMMIT,"paper_doi":"10.1086/724310",
                      "git_source_blobs":BLOBS,"branch_length_units":"source_ultrametric_units_not_independently_dated"},
            "source_checks":sinfo,
            "real_analysis":audit(tree,tipstates,args.simulations),
            "subspecies_62_tip_sensitivity":{
                "removed_taxon":"Cirsium_arizonicum_rothrockii",
                "remaining_taxa":62,
                "source_code_copy_is_not_independent_evolutionary_origin":True,
                "analysis":reduced,
            }}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    a=result["real_analysis"]
    print(json.dumps({"source":sinfo,"independent_vs_dependent":a["four_state_joint_model_comparison"],
         "conditioned_overlap_by_root":{root:{direction:{k:v[k] for k in ("accepted_exact_marginal","observed_overlap","null_overlap_mean","two_inclusive_tails_doubled_min")}
            for direction,v in item.items() if direction in ("head_evolved_outer_fixed","outer_evolved_head_fixed")}
            for root,item in a["independence_reference"].items()}},indent=2))

if __name__=="__main__": main()
