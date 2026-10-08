#!/usr/bin/env python3
"""Paired Cirsium organ allometric buffering, conditional fixed-strata reanalysis.

Requires original published XLSX (DOI 10.1007/s00606-023-01854-2,
Online Resource 6), imported by artifact_tool in sibling load_rows.
Does not identify developmental canalization or selection.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
from estimate_cirsium_cross_organ_variability_v1 import load_rows, SEED

PROXIES={
    "upper_3_leaf_total":"123Lf_L",
    "upper_internode_distance":"3Lf_D",
    "capitulum_count":"Cap_N",
}
ENDPOINTS={
    "middle_leaf":"Lf_L",
    "involucre":"Inv_L",
    "corolla":"Cor_L",
}

def estimate(groups, predictor):
    xx=[]
    yy={k:[] for k in ENDPOINTS}
    for plants in groups.values():
        x=np.array([np.log1p(p[predictor]) if predictor=="Cap_N"
                    else np.log(p[predictor]) for p in plants],dtype=float)
        x-=x.mean()
        xx.extend(x)
        for biological_part, code in ENDPOINTS.items():
            y=np.log([p[code] for p in plants])
            yy[biological_part].extend(y-y.mean())
    x=np.asarray(xx,dtype=float)
    denom=x @ x
    assert denom > 0
    slopes={k:float(x @ np.asarray(y,dtype=float)/denom) for k,y in yy.items()}
    slopes.update({
        "leaf_minus_involucre":slopes["middle_leaf"]-slopes["involucre"],
        "leaf_minus_corolla":slopes["middle_leaf"]-slopes["corolla"],
        "involucre_minus_corolla":slopes["involucre"]-slopes["corolla"],
    })
    return slopes

def analyze(xlsx, n_boot=4000):
    original,groups=load_rows(xlsx)
    assert len(original)==129 and len(groups)==6 and sum(map(len,groups.values()))==114
    rng=np.random.default_rng(SEED)
    models={}
    for label,predictor in PROXIES.items():
        obs=estimate(groups,predictor)
        samples=[]
        for _ in range(n_boot):
            bootstrap={
               k:[v[idx] for idx in rng.integers(0,len(v),len(v))]
               for k,v in groups.items()
            }
            samples.append(estimate(bootstrap,predictor))
        ci={k:[float(v) for v in np.quantile([r[k] for r in samples],[.025,.975])]
            for k in obs}
        models[label]={
          "input_character":predictor,
          "observed_slopes":obs,
          "paired_plant_ci95":ci,
          "prob_leaf_slope_larger_than_corolla":
             float(np.mean([r["leaf_minus_corolla"]>0 for r in samples]))
        }
    return {
      "data":"Published Cirsium individual morphometrics: 3 parental taxa x sex",
      "source_doi":"10.1007/s00606-023-01854-2",
      "individuals":114,"strata":6,"bootstrap":n_boot,"seed":SEED,
      "sensitivity_models":models,
      "claim_ceiling":"Within-stratum log length slopes conditional on published individual morphology; not experimentally controlled reaction norms, canalization, genetic covariance, or selection",
      "critical_counterexample":"Plant size measured as organ-specific leaf totals or stem distances; allometry alone cannot identify causal ecological mechanism"
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--xlsx",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=4000)
    a=p.parse_args()
    if a.bootstrap<100:p.error("Need >= 100 resamples")
    r=analyze(a.xlsx,a.bootstrap)
    a.out.parent.mkdir(exist_ok=True,parents=True)
    a.out.write_text(json.dumps(r,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"n":r["individuals"],
                      "proxies":{k:v["observed_slopes"] for k,v in r["sensitivity_models"].items()}},indent=2))

if __name__=="__main__":main()
