#!/usr/bin/env python3
"""Observed matched organ covariance: 3 Cirsium species x sex, fixed sampled sites.

Separates raw sex/taxon-centred phenotype correlation from a whole-head
size/capitulum-count/allometry sensitivity. Neither identifies developmental
mechanism, heritability, evolutionary change, or causal fitness selection.
"""
from __future__ import annotations

import argparse
import collections
import json
import math
import random
from pathlib import Path

import numpy as np
from estimate_cirsium_cross_organ_variability_v1 import load_rows, SEED

FIELDS=("Lf_L","Inv_L","Cor_L","Cap_W","123Lf_L","Cap_N")
PAIRS=(("leaf_involucre","Lf_L","Inv_L"),
       ("involucre_floret","Inv_L","Cor_L"),
       ("leaf_floret","Lf_L","Cor_L"))


def residual_correlations(groups, controls):
    arr={name:[] for name in FIELDS}
    for group in groups.values():
        for name in FIELDS:
            values=np.array([math.log1p(x[name]) if name=="Cap_N"
                             else math.log(x[name]) for x in group],dtype=float)
            arr[name].extend(values-values.mean())
    arr={k:np.asarray(v,dtype=float) for k,v in arr.items()}
    X=np.column_stack([arr[k] for k in controls]) if controls else None
    def resid(col):
        y=arr[col]
        return y if X is None else y-X@np.linalg.lstsq(X,y,rcond=None)[0]
    res={name:resid(name) for name in ("Lf_L","Inv_L","Cor_L")}
    result={label:float(np.corrcoef(res[a],res[b])[0,1])
            for label,a,b in PAIRS}
    result["difference"]=result["involucre_floret"]-result["leaf_involucre"]
    return result


def pct(xs,p):
    return float(np.quantile(np.asarray(xs,dtype=float),p))


def analyze(source, draws):
    rows, groups=load_rows(source)
    assert len(rows)==129 and sum(len(x) for x in groups.values())==114
    control=("Cap_W","123Lf_L","Cap_N")
    observed={
        "demeaned_only":residual_correlations(groups,()),
        "capitulum_width_adjusted":residual_correlations(groups,("Cap_W",)),
        "head_allometry_and_plant_proxies_adjusted":
            residual_correlations(groups,control),
    }
    rng=random.Random(SEED)
    boots=[]
    for _ in range(draws):
        rerows={k:[v[rng.randrange(len(v))] for i in range(len(v))]
                for k,v in groups.items()}
        boots.append(residual_correlations(rerows,control))
    ci={k:[pct([r[k] for r in boots],.025),pct([r[k] for r in boots],.975)]
        for k in boots[0]}
    return {"version":"aza3_cirsium_leaf_capitulum_covariance_v1",
       "source_doi":"10.1007/s00606-023-01854-2",
       "source_hash_verified_with_parent_script":True,
       "n_nonhybrid":114, "n_species_sex_strata":6,
       "scale":"log length residuals within each species x sex stratum, pooled over individuals",
       "controls":{
          "Cap_W":"log capitulum width",
          "123Lf_L":"log length of upper 3 leaves (vegetative size proxy)",
          "Cap_N":"log(1 + number of capitula)"},
       "models":observed,
       "paired_plant_bootstrap_draws":draws,"seed":SEED,
       "full_adjustment_ci95":ci,
       "bootstrap_fraction_difference_positive":
            sum(x["difference"]>0 for x in boots)/draws,
       "strong_claim_boundary":"Phenotypic integration difference is observed, but cannot distinguish common organ allometry, genetic pleiotropy, developmental constraint, adaptive selection or plant status. Conditioning on head width may overcontrol/collider bias; no causal mediation claim.",
       "population_boundary":"3 original parental taxa, sites fixed; no independence across phylogenetic branches"}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--xlsx",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--bootstrap",type=int,default=3000)
    args=parser.parse_args()
    r=analyze(args.xlsx,args.bootstrap)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":"PASS","models":r["models"],
          "allometry_CI":r["full_adjustment_ci95"]},indent=2))


if __name__=="__main__":
    main()
