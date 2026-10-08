#!/usr/bin/env python3
"""Prospectively frozen four-context check of Azami capitulum-only correlations.

Uses taxon medians preserved by frozen 36-pair replay and a separately
committed 12-test contract. No animal roles are observed in this dataset.

Null is exchangeable taxon memberships, NOT phylogenetic independent
evolution or randomized guild/fitness environment.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def bh(p):
    v=np.asarray(p,float)
    order=np.argsort(v)
    q=np.minimum.accumulate(
        (v[order]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]
    out=np.empty(len(v),dtype=float)
    out[order]=np.minimum(q,1.0)
    return out


def pearson(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    xx=x-x.mean(); yy=y-y.mean()
    den=np.sqrt((xx@xx)*(yy@yy))
    if den<=1e-12: return float("nan")
    return float(xx@yy/den)


def batch_r(mask, x, y):
    """Correlation within each boolean membership assignment, vectorized."""
    m=np.asarray(mask,float)
    x=np.asarray(x,float); y=np.asarray(y,float)
    n=m.sum(axis=1)
    sx=m@x; sy=m@y
    xc=m@(x*x)-sx*sx/n
    yc=m@(y*y)-sy*sy/n
    cc=m@(x*y)-sx*sy/n
    return cc/np.sqrt(np.maximum(1e-30,xc*yc))


def boot_corr_ci(x,y,rng,reps):
    n=len(x)
    idx=rng.integers(0,n,size=(reps,n))
    xa=x[idx]; ya=y[idx]
    xa-=xa.mean(axis=1,keepdims=True)
    ya-=ya.mean(axis=1,keepdims=True)
    r=np.sum(xa*ya,axis=1)/np.sqrt(np.maximum(1e-30,np.sum(xa*xa,axis=1)*np.sum(ya*ya,axis=1)))
    return [float(z) for z in np.quantile(r,[.025,.975])]


def analyze(replay, contract, nperm, nboot, seed):
    p=replay["taxon_level_context_profiles"]
    names=p["taxon_labels"]
    values=p["scalar_construct_medians"]
    context=p["matched_environment_and_image_quality"]
    n=len(names)
    assert n==42 and replay["n_complete_observations"]==1734
    assert all(len(v)==n for v in [*values.values(),*context.values()])
    tests=[]
    for split in contract["splits"]:
        field=split["field"]
        v=np.array(context[field],float)
        if split["id"]=="hemisphere_proxy":
            threshold=0.0
        else:
            threshold=float(np.median(v))
        high=v>=threshold
        nhigh=int(high.sum()); nlow=n-nhigh
        if min(nhigh,nlow)<8:
            tests.extend({"status":"NOT_IDENTIFIABLE_SMALL_GROUP","split":split["id"],"pair":pair}
                         for pair in contract["pairs"])
            continue
        rng=np.random.default_rng(seed)
        permutations=np.argsort(rng.random((nperm,n)),axis=1)
        perm_membership=np.zeros((nperm,n),dtype=bool)
        perm_membership[np.arange(nperm)[:,None],permutations[:,:nhigh]]=True
        for left,right in contract["pairs"]:
            x=np.array(values[left],float); y=np.array(values[right],float)
            lo=pearson(x[~high],y[~high]); hi=pearson(x[high],y[high])
            diff=float(np.arctanh(np.clip(hi,-.999999,.999999))-
                       np.arctanh(np.clip(lo,-.999999,.999999)))
            nullhigh=batch_r(perm_membership,x,y)
            nulllow=batch_r(~perm_membership,x,y)
            nulldiff=np.arctanh(np.clip(nullhigh,-.999999,.999999))-np.arctanh(np.clip(nulllow,-.999999,.999999))
            pv=float((1+np.count_nonzero(np.abs(nulldiff)>=abs(diff)))/(nperm+1))
            rngb=np.random.default_rng(seed+11)
            ci_low=boot_corr_ci(x[~high].copy(),y[~high].copy(),rngb,nboot)
            ci_high=boot_corr_ci(x[high].copy(),y[high].copy(),rngb,nboot)
            tests.append({
                "split":split["id"],"context_variable":field,
                "split_threshold":threshold,
                "low_n":nlow,"high_n":nhigh,
                "left":left,"right":right,"low_r":lo,"high_r":hi,
                "fisher_z_difference":diff,
                "conditional_taxon_shuffle_p_two_sided":pv,
                "low_boot_ci95":ci_low,"high_boot_ci95":ci_high,
                "both_groups_same_correlation_sign":bool(hi*lo>0),
                "role":"DESCRIPTIVE_CONTEXT_SENSITIVITY_NOT_ECOLOGICAL_REGIME",
                "scope_warning":"Taxon median sampling geography, not actual native region, effective guild pressure or random trait origin."
            })
    valid=[v for v in tests if "conditional_taxon_shuffle_p_two_sided" in v]
    corrected=bh([v["conditional_taxon_shuffle_p_two_sided"] for v in valid])
    for t,q in zip(valid,corrected):
        t["bh_q_12"]=float(q)
    return {
        "version":"capitulum_conditional_coupling_12test_result_v1",
        "status":"REPLAYED_CONTEXT_SENSITIVITY_WITHOUT_ADAPTIVE_MECHANISM",
        "n_taxa":n,"n_scored_pairs":len(contract["pairs"]),
        "n_declared_split_rules":len(contract["splits"]),
        "n_declared_tests":len(contract["pairs"])*len(contract["splits"]),
        "n_completed":len(valid),"n_bh_q_below_0_05":sum(t["bh_q_12"]<.05 for t in valid),
        "permutations":nperm,"taxon_bootstrap_draws_per_stratum":nboot,"seed":seed,
        "source_contract":contract,
        "tests":tests,
        "scope":{
            "negative_control":"photo_sharpness demonstrates technical heterogeneity may mimic ecological regimes",
            "this_is_not":"an interaction of experimentally tested floral traits on individual reproductive fitness",
            "western_eastern_taxon_split":"longitude median of sampled photos, not a native biogeographic region or independent evolutionary lineage",
            "genealogy":"No phylogenetic covariance is preserved by the taxon-shuffling null",
            "ecology":"No plant–insect interactions or fitness outcomes are present in these numbers",
            "null_result":"A null detected interaction does not prove equal effects or invariance across real herbivore/pollinator regimes."
        }
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--replay",type=Path,required=True)
    parser.add_argument("--contract",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--permutations",type=int,default=9999)
    parser.add_argument("--bootstrap",type=int,default=3000)
    parser.add_argument("--seed",type=int,default=20261008)
    a=parser.parse_args()
    if a.permutations<99 or a.bootstrap<100: parser.error("Resampling counts too small")
    replay=json.loads(a.replay.read_text())
    contract=json.loads(a.contract.read_text())
    r=analyze(replay,contract,a.permutations,a.bootstrap,a.seed)
    if r["n_completed"]!=12: raise ValueError("Context subgroup counts insufficient for declared contrasts")
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2)+"\n")
    print(json.dumps({"status":r["status"],"n_completed":r["n_completed"],
        "n_bh_significant":r["n_bh_q_below_0_05"],
        "contrasts": [{k:v[k] for k in ("split","left","right","low_n","high_n","low_r","high_r","conditional_taxon_shuffle_p_two_sided","bh_q_12")} for v in r["tests"]]},indent=2))

if __name__=="__main__":
    main()
