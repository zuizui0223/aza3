#!/usr/bin/env python3
"""Original-source North American Cirsium HEAD module vs pollination coding.

Taxa, morphotypes and an ultrametric tree from Siniscalchi et al. (2023),
International Journal of Plant Sciences, DOI 10.1086/724310, git-pinned
author repository carol-siniscalchi/cirsium.

Outcome is ORIGINAL "Hummingbird" vs "Insect" classification. It may
be an authored POLLINATION SYNDROME assembled partly from corolla/head
traits, not field-observed pollen deposition! Associations are descriptive.
Phylogenetic simulation holds this tree and measured head states fixed and
simulates INDEPENDENT Mk evolution of the original pollination classification.
It is a model-based sensitivity, NOT identification of correlated selection.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
from Bio import Phylo
from scipy.optimize import minimize
from scipy.stats import fisher_exact

SOURCE_REPO = "carol-siniscalchi/cirsium"
SOURCE_COMMIT = "c925d4af86160064448f794a9715f70b917a2786"
SOURCE_DOI = "10.1086/724310"
SOURCE_BLOBS = {
    "traits.csv": "fa127f0587b676e6e004ea121424d26069586bcb",
    "poll.csv": "641ad617b92857a4db0744f799c7ac66ba483ce7",
    "cirsium.phyluce95p.ultrametric.tre": "6417a38793427bc7c8562dc0e2cab59b74d8c316",
}
HEAD_TRAITS = {
    "outer_phyllaries_spreading_or_reflexed":
        ("OuterPhyllaryPosition", {"0","1","2"}, {"1","2"}),
    "inner_phyllaries_non_green":
        ("InnerPhyllaryColor", {"0","1","2","3"}, {"1","2","3"}),
    "head_pendulous":
        ("HeadPosition", {"0","1"}, {"1"}),
    "glutinous_ridge_present":
        ("GlutinousRidge", {"0","1"}, {"1"}),
}
ADDITIONAL_TRAITS = {
    "outer_phyllaries_reflexed":
        ("OuterPhyllaryPosition", {"0","1","2"}, {"2"}),
    "red_corolla_positive_control":
        ("CorollaColor", {"0","1","2","3"}, {"0"}),
}
PRIMARY_ORDER = tuple(HEAD_TRAITS)
N_ORIGINAL_TAXA = 63
N_BIRDS = 7
EXPECTED_SUBSPECIES = ("Cirsium_arizonicum_arizonicum",
                       "Cirsium_arizonicum_rothrockii")
SEED = 20261020


def assert_true(ok, message):
    if not ok:
        raise ValueError(message)


def git_blob(raw):
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def read_pinned_text(path, original_name):
    raw=Path(path).read_bytes()
    assert_true(git_blob(raw)==SOURCE_BLOBS[original_name],
                "SOURCE_GIT_BLOB_CHANGED_"+original_name)
    return raw.decode("utf-8-sig")


def read_csv_records(src, identity):
    reader=csv.DictReader(io.StringIO(src))
    rows=list(reader)
    assert_true(bool(rows),"EMPTY_SOURCE_"+identity)
    validate_unique_taxa(rows, identity)
    return rows


def validate_unique_taxa(rows, key):
    values=[row.get(key) for row in rows]
    assert_true(all(values),"MISSING_TAXON")
    assert_true(len(values)==len(set(values)),"DUPLICATE_TAXON")


def trait_positive(row, trait_name):
    if trait_name in HEAD_TRAITS:
        key,allowed,positive=HEAD_TRAITS[trait_name]
    elif trait_name in ADDITIONAL_TRAITS:
        key,allowed,positive=ADDITIONAL_TRAITS[trait_name]
    else:
        raise ValueError("UNKNOWN_HEAD_TRAIT_NAME")
    v=row.get(key)
    if v not in allowed:
        raise ValueError("UNKNOWN_HEAD_TRAIT_CODE_"+trait_name)
    return v in positive


def fisher_binary(bird_yes,bird_no,insect_yes,insect_no):
    vals=(bird_yes,bird_no,insect_yes,insect_no)
    assert_true(all(isinstance(z,(int,np.integer)) and z>=0 for z in vals)
                and bird_yes+bird_no>0 and insect_yes+insect_no>0,
                "INVALID_FISHER_TAXON_COUNTS")
    p=float(fisher_exact([[bird_yes,bird_no],
                          [insect_yes,insect_no]],alternative="two-sided").pvalue)
    odds=(float(bird_yes*insect_no/(bird_no*insect_yes))
          if bird_no*insect_yes else
          (0.0 if bird_yes*insect_no==0 and insect_yes>0 else None))
    return {
        "bird_positive":int(bird_yes),
        "bird_negative":int(bird_no),
        "insect_positive":int(insect_yes),
        "insect_negative":int(insect_no),
        "bird_fraction":bird_yes/(bird_yes+bird_no),
        "insect_fraction":insect_yes/(insect_yes+insect_no),
        "bird_minus_insect_prevalence":
            bird_yes/(bird_yes+bird_no)-insect_yes/(insect_yes+insect_no),
        "odds_ratio":odds,
        "p_two_sided":p,
        "p_is_iid_label_test_not_phylogenetically_independent":True,
    }


def bh_adjust(pvals):
    assert_true(bool(pvals) and all(math.isfinite(x) and 0<=x<=1
                                    for x in pvals),"INVALID_P_VALUES")
    n=len(pvals)
    z=sorted(range(n),key=lambda i:pvals[i])
    q=[0.0]*n
    running=1.0
    for rank in range(n,0,-1):
        k=z[rank-1]
        running=min(running, pvals[k]*n/rank)
        q[k]=running
    return q


def verify_tip_set(tree,taxon_set):
    names=[c.name for c in tree.get_terminals()]
    assert_true(bool(names) and len(names)==len(set(names))
                and set(names)==set(taxon_set),"TREE_TIP_SET_MISMATCH")
    for c in tree.find_clades(order="preorder"):
        d=0.0 if c.branch_length is None else float(c.branch_length)
        assert_true(math.isfinite(d) and d>=0,"INVALID_TREE_BRANCH_LENGTH")
    return len(names)


def parsimony_cost(tree, tip_states, root_state=None):
    """Min number binary state changes, *not* number of independent gains."""
    def visit(node):
        if node.is_terminal():
            if node.name not in tip_states or tip_states[node.name] not in (0,1):
                raise ValueError("MISSING_BINARY_TIP_STATE")
            k=tip_states[node.name]
            return [0 if k==0 else math.inf, 0 if k==1 else math.inf]
        child=[visit(z) for z in node.clades]
        return [sum(min(c[0],c[1]+1) for c in child),
                sum(min(c[0]+1,c[1]) for c in child)]
    v=visit(tree.root)
    if root_state not in (None,0,1):
        raise ValueError("INVALID_ROOT_STATE")
    return int(min(v) if root_state is None else v[root_state])


def binary_mk_loglik(tree,states,a,b):
    """Felsenstein binary CTMC: root fixed to 0 (insect-coded outgroup)."""
    if a<=0 or b<=0:raise ValueError("MK_RATES_MUST_BE_POSITIVE")
    cache={}
    rate=a+b
    pi1=a/rate
    pi0=b/rate
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            if node.name not in states or states[node.name] not in (0,1):
                raise ValueError("MISSING_SOURCE_TIP_POLLINATOR")
            v=states[node.name]
            cache[node]=(1.0,0.0) if v==0 else (0.0,1.0)
        else:
            l0,l1=1.0,1.0
            for child in node.clades:
                exp=math.exp(-rate*(child.branch_length or 0.0))
                p00=pi0+pi1*exp
                p01=pi1*(1-exp)
                p10=pi0*(1-exp)
                p11=pi1+pi0*exp
                c0,c1=cache[child]
                l0 *= p00*c0+p01*c1
                l1 *= p10*c0+p11*c1
            cache[node]=(l0,l1)
    p=cache[tree.root][0]
    return math.log(p) if p>0 else -math.inf


def fit_independent_binary_mk(tree,states):
    """Fit pollinator Mk only; head states are not used in fitting rates."""
    initial=[(-2,-1),(-3,0),(-1,-1),(0,0)]
    def objective(v):
        ll=binary_mk_loglik(tree,states,math.exp(v[0]),math.exp(v[1]))
        return -ll if math.isfinite(ll) else 1e250
    fits=[minimize(objective,v,method="L-BFGS-B",bounds=[(-9,3),(-9,3)])
          for v in initial]
    fits=[v for v in fits if v.success and math.isfinite(v.fun)]
    if not fits:
        raise ValueError("MK_PARAMETER_FIT_NOT_CONVERGED")
    best=min(fits,key=lambda f:f.fun)
    a,b=(float(math.exp(v)) for v in best.x)
    return {
        "insect_to_bird_rate":a,
        "bird_to_insect_rate":b,
        "root_condition":"INSECT_STATE_0_FIXED",
        "fit_source":"POLLINATOR_TIP_LABELS_ONLY",
        "maximum_loglik":float(-best.fun),
        "minimizer_converged":True,
    }


def simulate_mk_conditional(tree,tip_states,phenotypes,a,b,
                            n_sims=250000,seed=SEED,min_accepted=1000):
    """Independent-head Mk null, conditional on observed total bird tips.

    Under null pollination syndrome evolves by homogeneous binary Mk on the
    frozen tree with root insect and its independently fitted rates.
    Original morphology remains FIXED at tips; a tail tests concentration of
    coded head phenotypes among bird-coded tips under that conditional null.
    Model and classification correctness are separate assumptions.
    """
    taxa=[z.name for z in tree.get_terminals()]
    if set(taxa)!=set(tip_states) or set(taxa)!=set(phenotypes):
        raise ValueError("MK_SOURCE_TIP_ALIGNMENT_FAILURE")
    k=int(sum(tip_states.values()))
    if not (0<k<len(taxa)) or n_sims<100 or n_sims>1000000:
        raise ValueError("MK_SIMULATION_REQUEST_INVALID")
    size={len(v) for v in phenotypes.values()}
    if len(size)!=1 or not next(iter(size)):
        raise ValueError("MK_PHENOTYPE_MATRIX_INVALID")
    r=np.random.default_rng(seed)
    states={tree.root:np.zeros(n_sims,dtype=np.uint8)}
    rate=a+b
    pi1=a/rate
    pi0=b/rate
    for node in tree.find_clades(order="preorder"):
        if node is tree.root:
            continue
        # clade instances have identity-hash semantics in Biopython
        parent=None
        # a single precomputed parent map is much cheaper than nested searches
        if node not in parent_lookup:
            raise ValueError("TREE_PARENT_LOOKUP_DRIFT")
        parent=parent_lookup[node]
        e=math.exp(-rate*(node.branch_length or 0.0))
        p01=pi1*(1-e)
        p11=pi1+pi0*e
        u=r.random(n_sims)
        states[node]=np.where(states[parent]!=0,u<p11,u<p01).astype(np.uint8)
    actual=np.vstack([states[leaf] for leaf in tree.get_terminals()])
    mask=actual.sum(axis=0)==k
    n_accepted=int(mask.sum())
    if n_accepted<min_accepted:
        raise ValueError("INSUFFICIENT_CONDITIONED_MK_DRAWS")
    out={}
    for j in range(next(iter(size))):
        x=np.array([phenotypes[name][j] for name in taxa],dtype=np.uint8)
        if not np.isin(x,[0,1]).all():
            raise ValueError("INVALID_BINARY_HEAD_FEATURE")
        observed=sum(tip_states[taxa[i]]*int(x[i]) for i in range(len(taxa)))
        hits=np.matmul(x.astype(np.int16),actual[:,mask].astype(np.int16))
        count=int((hits>=observed).sum())
        probability=(count+1)/(n_accepted+1)
        out["trait_"+str(j)]={
            "bird_head_positive_observed":observed,
            "upper_tail_probability":probability,
            "n_upper_tail_draws":count,
            "conditional_monte_carlo_se":float(math.sqrt(
                probability*(1-probability)/(n_accepted+1))),
        }
    return {
        "n_simulated_trees":n_sims,
        "n_conditioned_draws":n_accepted,
        "observed_bird_tips":k,
        "seed":seed,
        "null_condition":"SIMULATED_BINARY_CTMC_CONDITION_ON_EXACT_OBSERVED_BIRD_TIP_COUNT",
        "phenotypes":out,
        "is_identification_of_natural_selection":False,
        "is_inference_of_mechanistic_spine_or_pollination_success":False,
    }


def claim_limits():
    return {
        "same_head_morphology_and_fitness_measured":False,
        "direct_bird_pollination_assayed":False,
        "adaptation_or_convergent_selection_identified":False,
        "independent_origins_of_correlated_trait_changes_identified":False,
        "hummingbird_label_independent_of_flower_colour_classification_verified":False,
        "true_spine_barrier_geometry_measured":False,
        "species_or_subspecies_are_independent_evolutionary_replicates":False,
        "generalize_to_east_asian_ryukyu_thistles":False,
        "phylogenetic_Mk_model_proves_discrete_evolution_history":False,
        "source_has_original_filled_achene_or_bird_pollen_success":False,
    }


def source_data(traits_path,poll_path,tree_path):
    ts=read_csv_records(read_pinned_text(traits_path,"traits.csv"),"Species")
    pp=read_csv_records(read_pinned_text(poll_path,"poll.csv"),"species")
    assert_true(len(ts)==len(pp)==N_ORIGINAL_TAXA,"ORIGINAL_63_TAXA_DRIFT")
    names={z["Species"] for z in ts}
    assert_true(names=={z["species"] for z in pp},"TRAIT_POLLINATOR_NAME_MISMATCH")
    zmap={z["species"]:z["pollinator"] for z in pp}
    assert_true(set(zmap.values())=={"Hummingbird","Insect"},
                "UNEXPECTED_POLLINATOR_CLASS")
    assert_true(sum(v=="Hummingbird" for v in zmap.values())==N_BIRDS,
                "ORIGINAL_7_HUMMINGBIRD_TAXA_DRIFT")
    tree=Phylo.read(io.StringIO(read_pinned_text(
        tree_path,"cirsium.phyluce95p.ultrametric.tre")),"newick")
    assert_true(verify_tip_set(tree,names)==N_ORIGINAL_TAXA,"ORIGINAL_TREE_DRIFT")
    depth=tree.depths()
    heights=[float(depth[tip]) for tip in tree.get_terminals()]
    assert_true(max(heights)-min(heights)<1e-3,"SOURCE_TREE_NOT_ULTRAMETRIC")
    return ts,zmap,tree,{"n_traits":len(ts),"n_pollinator":len(pp),
                        "n_tree_tips":len(tree.get_terminals()),
                        "n_hummingbird_syndrome":N_BIRDS,
                        "n_insect_syndrome":len(pp)-N_BIRDS,
                        "tree_depth_range":[min(heights),max(heights)]}


def run_analysis(ts,poll,tree,n_sims=250000):
    by={z["Species"]:z for z in ts}
    names=sorted(by)
    outcomes={}
    primary_p=[]
    for trait in list(HEAD_TRAITS)+list(ADDITIONAL_TRAITS):
        y={tax:trait_positive(by[tax],trait) for tax in names}
        b1=sum(y[s] for s in names if poll[s]=="Hummingbird")
        b0=sum(not y[s] for s in names if poll[s]=="Hummingbird")
        i1=sum(y[s] for s in names if poll[s]=="Insect")
        i0=sum(not y[s] for s in names if poll[s]=="Insect")
        out=fisher_binary(b1,b0,i1,i0)
        if trait in HEAD_TRAITS:
            primary_p.append(out["p_two_sided"])
            out["analysis_set"]="FOUR_DECLARED_HEAD_TRAITS"
        else:
            out["analysis_set"]="ANCILLARY_SYNDROME_OR_OUTER_SUBSTATE"
        outcomes[trait]=out
    for name,q in zip(HEAD_TRAITS,bh_adjust(primary_p)):
        outcomes[name]["BH_q_within_four_head_traits"]=q
    # Morphotype redundancy: one source taxon is a second subspecies entry.
    aa,rr=EXPECTED_SUBSPECIES
    assert_true(aa in by and rr in by and
                poll[aa]==poll[rr]=="Hummingbird" and
                all(by[aa][k]==by[rr][k] for k in
                    ("OuterPhyllaryPosition","InnerPhyllaryColor",
                     "HeadPosition","GlutinousRidge")),
                "ARIZONICUM_SUBSPECIES_SOURCE_STATE_CHANGED")
    taxa_one_subspecies=[s for s in names if s!=rr]
    collapse={}
    leave_one={}
    for trait in HEAD_TRAITS:
        def calc(selected):
            yy={s:trait_positive(by[s],trait) for s in selected}
            return fisher_binary(
                sum(yy[s] for s in selected if poll[s]=="Hummingbird"),
                sum(not yy[s] for s in selected if poll[s]=="Hummingbird"),
                sum(yy[s] for s in selected if poll[s]=="Insect"),
                sum(not yy[s] for s in selected if poll[s]=="Insect"))
        collapse[trait]=calc(taxa_one_subspecies)
        leave_one[trait]={
            s:calc([z for z in names if z!=s])
            for s in names if poll[s]=="Hummingbird"
        }
    binary_poll={s:int(poll[s]=="Hummingbird") for s in names}
    binary_phen={
        name:{s:int(trait_positive(by[s],name)) for s in names}
        for name in HEAD_TRAITS
    }
    mk=fit_independent_binary_mk(tree,binary_poll)
    # exact tip identity, orders from original HEAD_TRAITS dict insertion
    matrix={s:[binary_phen[t][s] for t in HEAD_TRAITS] for s in names}
    global parent_lookup
    parent_lookup={child:node for node in tree.find_clades(order="preorder")
                   for child in node.clades}
    sim=simulate_mk_conditional(
        tree,binary_poll,matrix,
        a=mk["insect_to_bird_rate"],b=mk["bird_to_insect_rate"],
        n_sims=n_sims,seed=SEED,
        min_accepted=1000 if n_sims>=100000 else 25)
    sim["trait_order"]=list(HEAD_TRAITS)
    for j,t in enumerate(HEAD_TRAITS):
        sim["phenotypes"][t]=sim["phenotypes"].pop("trait_"+str(j))
        sim["phenotypes"][t]["observed_trait_two_sided_fisher_p"] = outcomes[t]["p_two_sided"]
    return {
       "head_trait_binary_associations":outcomes,
       "subspecies_pair_sensitivity":{
           "removed_second_subspecies":rr,
           "taxon_count_after_collapse":62,
           "hummingbird_taxa_after_collapse":6,
           "associations":collapse,
       },
       "leave_one_bird_taxon_out":leave_one,
       "phylogenetic_sensitivity":{
           "minimum_change_count_for_pollinator_with_insect_root":
               parsimony_cost(tree,binary_poll,root_state=0),
           "minimum_head_state_changes_unconstrained_root":{
               n:parsimony_cost(tree,s) for n,s in binary_phen.items()
           },
           "parsimony_change_count_not_identified_convergent_gain_count":True,
           "independent_pollinator_Mk_fit":mk,
           "conditional_Mk_simulation":sim,
       },
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--traits",type=Path,required=True)
    p.add_argument("--poll",type=Path,required=True)
    p.add_argument("--tree",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--simulations",type=int,default=250000)
    a=p.parse_args()
    traits,poll,tree,metadata=source_data(a.traits,a.poll,a.tree)
    result={
        "version":"aza3_20261020_na_63_head_syndrome_phylogenetic_sensitivity_v1",
        "source":{
            "author_repo":SOURCE_REPO,"git_commit":SOURCE_COMMIT,
            "published_paper_doi":SOURCE_DOI,"published_year":2023,
            "git_blob_sha1":SOURCE_BLOBS,
            "source_classification_not_field_service":True,
            **metadata,
        },
        "analysis":run_analysis(traits,poll,tree,n_sims=a.simulations),
        "claim_limits":claim_limits(),
        "primary_hypothesis":"Do outer phyllary posture and other capitulum features show a specific taxon-state association with the original bird-vs-insect syndrome after phylogenetic and subspecies caution?",
        "interpretation_ceiling":"Comparative phylogenetic coding of head phenotype and pollinator syndromes. Not measured pollination efficiency, structural defence, correlated adaptive evolution, independent trait origin co-localization or genetic LEGO. Hummingbird classification might use the same flower characters.",
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":"SOURCE_AUTHENTICATED_63_TIP_MORPHOLOGY_COMPARISON",
        "pollination_changes_parsimony":result["analysis"]["phylogenetic_sensitivity"]["minimum_change_count_for_pollinator_with_insect_root"],
        "head_phenotypes":{name:{
            "bird_yes":x["bird_positive"],
            "insect_yes":x["insect_positive"],
            "fisher_p":x["p_two_sided"],
            "bh_q":x.get("BH_q_within_four_head_traits"),
        } for name,x in result["analysis"]["head_trait_binary_associations"].items()},
        "phylo_conditional_acceptance":result["analysis"]["phylogenetic_sensitivity"]["conditional_Mk_simulation"]["n_conditioned_draws"],
        "phylo_tail":{k:v["upper_tail_probability"] for k,v in result["analysis"]["phylogenetic_sensitivity"]["conditional_Mk_simulation"]["phenotypes"].items()},
    },indent=2))


if __name__=="__main__":
    main()
