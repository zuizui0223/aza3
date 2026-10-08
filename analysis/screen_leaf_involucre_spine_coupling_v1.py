#!/usr/bin/env python3
"""Exploratory leaf thorn × head projection screen; separate from Azami GEB.

This tests a *proxy* comparison, not the same anatomical spine trait:
  - leaf spine: terminal spine on longest middle-leaf lobe (Spi_L, mm);
  - involucre: HALF of the excess head width caused by bract tips
                (Ibr_W - Inv_W)/2, mm; NOT botanical spine length.
The alternative relative measures Lf_Sp and Inv_Sp-1 have different
denominators and are sensitivity only. Inference is conditional, not
evidence of evolutionary independence/equivalence or adaptive defence.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr
from estimate_cirsium_cross_organ_variability_v1 import load_rows, SEED

PAIRS = (
    ("leaf_spine_mm", "bract_width_extension_mm"),
    ("leaf_normalized_spine", "bract_projection_ratio"),
    ("involucre_length", "corolla_length"),
)

def get_value(r, name):
    if name == "leaf_spine_mm": return r["Spi_L"]
    if name == "bract_width_extension_mm":
        return (r["Ibr_W"] - r["Inv_W"]) / 2.0
    if name == "leaf_normalized_spine": return r["Lf_Sp"]
    if name == "bract_projection_ratio": return r["Inv_Sp"] - 1.0
    if name == "involucre_length": return r["Inv_L"]
    if name == "corolla_length": return r["Cor_L"]
    raise KeyError(name)

def group_rank_cor(groups, left, right):
    aa, bb = [], []
    for observations in groups.values():
        x = np.array([get_value(r,left) for r in observations], dtype=float)
        y = np.array([get_value(r,right) for r in observations], dtype=float)
        # Pairwise rank association, standardized for stratum n (n ranges 10–30).
        rx = rankdata(x,method="average") / len(observations)
        ry = rankdata(y,method="average") / len(observations)
        aa.extend(rx-rx.mean())
        bb.extend(ry-ry.mean())
    return float(np.corrcoef(aa,bb)[0,1])

def infer(source, boot_n=4000, permutations=4999):
    all_rows, groups = load_rows(source)
    assert len(all_rows)==129 and sum(len(x) for x in groups.values())==114
    assert len(groups)==6
    obs={f"{a}__{b}":group_rank_cor(groups,a,b) for a,b in PAIRS}
    boot={key:[] for key in obs}
    rng=np.random.default_rng(SEED)
    for _ in range(boot_n):
        resampled={k:[v[i] for i in rng.integers(0,len(v),len(v))]
                   for k,v in groups.items()}
        for a,b in PAIRS:
            boot[f"{a}__{b}"].append(group_rank_cor(resampled,a,b))
    ci={k:[float(np.percentile(v,2.5)),float(np.percentile(v,97.5))]
        for k,v in boot.items()}
    null=[]
    for _ in range(permutations):
        u,v=[],[]
        for rs in groups.values():
            a=np.array([get_value(r,"leaf_spine_mm") for r in rs])
            b=np.array([get_value(r,"bract_width_extension_mm") for r in rs])
            ra=rankdata(a)/len(rs)
            rb=rankdata(b)/len(rs)
            u.extend(ra-ra.mean())
            v.extend(rb[rng.permutation(len(rs))]-rb.mean())
        null.append(float(np.corrcoef(u,v)[0,1]))
    observed=obs["leaf_spine_mm__bract_width_extension_mm"]
    p=(1+sum(abs(v)>=abs(observed) for v in null))/(permutations+1)
    strata=[]
    for (taxon,sex),rs in sorted(groups.items()):
        x=np.array([get_value(r,"leaf_spine_mm") for r in rs])
        y=np.array([get_value(r,"bract_width_extension_mm") for r in rs])
        a=np.array([get_value(r,"leaf_normalized_spine") for r in rs])
        b=np.array([get_value(r,"bract_projection_ratio") for r in rs])
        strata.append({
          "taxon":taxon,"sex":sex,"n":len(rs),
          "raw_proxy_spearman":float(spearmanr(x,y).statistic),
          "relative_proxy_spearman":float(spearmanr(a,b).statistic),
          "zero_extension_count":int(sum(y==0))
        })
    # Tests: near-zero pooled association is NOT equivalence and does not
    # guarantee directionally homogeneous strata.
    return {
      "source_doi":"10.1007/s00606-023-01854-2",
      "status":"EXPLORATORY_POST_HOC_SENSITIVITY_ONLY",
      "n_parental_individuals":114,
      "strata_n":6,
      "rank_association":obs,
      "paired_plant_bootstrap_ci95":ci,
      "n_bootstrap":boot_n,
      "within_stratum_pairing_null":p,
      "n_permutation":permutations,
      "strata":strata,
      "observed_nonpositive_bract_extension_count":sum(z["zero_extension_count"] for z in strata),
      "claim_boundary":[
        "External tip-width extension is not individual phyllary spine length, hardness, or antiherbivore efficacy.",
        "Different relative normalizations and zero width differences may attenuate correlations.",
        "Near-zero observed r/p-value does not establish statistical independence or equivalence.",
        "Species and sex are only six fixed strata, not independent evolutionary transitions.",
        "Causal genetic/reproductive adaptation not measured; this screen cannot reclassify frozen Azami GEB results."
      ]
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--xlsx",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result=infer(args.xlsx)
    args.out.parent.mkdir(exist_ok=True,parents=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"n":result["n_parental_individuals"],
                      "rank_association":result["rank_association"],
                      "ci95":result["paired_plant_bootstrap_ci95"]},indent=2))

if __name__=="__main__":
    main()
