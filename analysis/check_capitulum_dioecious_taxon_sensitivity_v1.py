#!/usr/bin/env python3
"""Stress-test head-only Azami cross-taxon pair correlations removing dioecious Cirsium arvense.

This is NOT an image-level reproductive-sex classifier. C. arvense is
imperfectly dioecious; removal tests only reliance on that taxon's median.
Within-species unlabelled sex and phenology remain unresolved.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np

SELECTED=(
 ("presentation_angle","head_elongation"),
 ("floral_chroma","projection_prominence"),
 ("presentation_angle","projection_prominence"),
)
COV=("latitude","chelsa_bio12","chelsa_rsds_mean",
     "min_dimension","sharpness","mask_quality")

def cor(x,y): return float(np.corrcoef(x,y)[0,1])

def conditional(x,y,c):
    q=np.asarray(c,float)
    q=(q-q.mean(0))/np.where(q.std(0)>0,q.std(0),1)
    z=np.column_stack((np.ones(len(q)),q))
    a=x-z@np.linalg.lstsq(z,x,rcond=None)[0]
    b=y-z@np.linalg.lstsq(z,y,rcond=None)[0]
    return cor(a,b)

def compute(d):
    q=d["taxon_level_context_profiles"]
    names=q["taxon_labels"]
    X=q["scalar_construct_medians"]
    env=q["matched_environment_and_image_quality"]
    assert len(names)==42 and names.count("Cirsium arvense")==1
    Z=np.column_stack([env[c] for c in COV])
    result={}
    for left,right in SELECTED:
        x=np.asarray(X[left],float);y=np.asarray(X[right],float)
        r={}
        for name,exclude in (
            ("all_42",None),
            ("remove_Cirsium_arvense","Cirsium arvense"),
            ("remove_Cirsium_japonicum","Cirsium japonicum")
        ):
            mask=np.array([t!=exclude for t in names])
            r[name]={"n_taxa":int(mask.sum()),
                     "raw_r":cor(x[mask],y[mask]),
                     "six_covariate_partial_r":conditional(x[mask],y[mask],Z[mask])}
        result[left+"__"+right]=r
    assert result["presentation_angle__head_elongation"]["remove_Cirsium_arvense"]["raw_r"]<-.5
    assert result["floral_chroma__projection_prominence"]["remove_Cirsium_arvense"]["raw_r"]>.5
    return {
        "status":"PASS_UNLABELLED_SEX_NEGATIVE_CONTROL_LIMITED",
        "source":"Azami 42 taxa, 1734 visible heads; C. arvense imperfectly dioecious; source http://doi.org/10.1111/j.1365-2745.2010.01678.x",
        "paired_phenotype_sensitivity":result,
        "known_limitations":[
            "Within-C. arvense males/females may differ in head proportions; deleting a species does not resolve image sex or phenology.",
            "The six covariance controls are post-hoc, not a validated causal adjustment.",
            "Cirsium japonicum removal is additional species influence control, not reproductive sex control.",
            "No fitness outcome, true spine, directly verified access route or correlated evolution follows from these associations."
        ]
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--replay",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    d=json.loads(args.replay.read_text())
    output=compute(d)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps(output,indent=2))
if __name__=="__main__":
    main()
