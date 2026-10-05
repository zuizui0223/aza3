#!/usr/bin/env python3
"""Test whether capitulum constructs occupy module-specific ecological fingerprints.

This is a frozen-data synthesis of Azami v3 biological constructs. It does not
reinterpret the original endpoint-level multiplicity family and does not claim
selection or adaptation.

Primary estimand:
  mean cosine similarity among construct pairs within the predeclared biological
  modules minus mean cosine similarity among construct pairs between modules,
  using six normalized environmental-block RMS effect-magnitude signatures.

The null is evaluated exactly across all distinct assignments of module labels
with the observed module sizes (1,3,2,3): 1680 assignments.
"""
from __future__ import annotations
import argparse, csv, itertools, json, math
from pathlib import Path

CORE = [
    "presentation_angle",
    "floral_lightness",
    "floral_chroma",
    "floral_hue",
    "head_elongation",
    "head_compactness",
    "involucre_form",
    "projection_prominence",
    "projection_pattern",
]
MODULES = {
    "presentation_angle": "presentation",
    "floral_lightness": "colour",
    "floral_chroma": "colour",
    "floral_hue": "colour",
    "head_elongation": "head_form",
    "head_compactness": "head_form",
    "involucre_form": "involucre_armature",
    "projection_prominence": "involucre_armature",
    "projection_pattern": "involucre_armature",
}
ENV_BLOCKS = {
    "thermal": ["chelsa_bio01", "chelsa_bio04"],
    "hydric": ["chelsa_bio12", "chelsa_bio15"],
    "radiative_atmospheric": ["chelsa_rsds_mean", "chelsa_vpd_mean"],
    "mechanical": ["chelsa_sfcwind_mean"],
    "growing_season_water": ["chelsa_gsp"],
    "resource_productivity": ["chelsa_npp"],
}
LABELS = ["presentation"] + ["colour"]*3 + ["head_form"]*2 + ["involucre_armature"]*3

def load_axis(path: Path):
    rows=list(csv.DictReader(path.open(encoding="utf-8-sig")))
    d={}
    for r in rows:
        c=r["construct_id"]; p=r["predictor"]
        if c not in CORE:
            continue
        d[(c,p)] = float(r["effect_magnitude"])
    return d

def fingerprints(path: Path):
    d=load_axis(path)
    out={}
    raw={}
    for c in CORE:
        vals={}
        for b,ps in ENV_BLOCKS.items():
            xs=[d[(c,p)] for p in ps]
            vals[b]=math.sqrt(sum(x*x for x in xs)/len(xs))
        n=math.sqrt(sum(v*v for v in vals.values()))
        raw[c]=vals
        out[c]={b:(v/n if n else 0.0) for b,v in vals.items()}
    return raw,out

def cosine(a,b):
    keys=list(ENV_BLOCKS)
    num=sum(a[k]*b[k] for k in keys)
    da=math.sqrt(sum(a[k]**2 for k in keys))
    db=math.sqrt(sum(b[k]**2 for k in keys))
    return num/(da*db) if da and db else float("nan")

def pair_table(fp):
    rows=[]
    for i,a in enumerate(CORE):
        for b in CORE[i+1:]:
            rows.append({
                "construct_a":a,
                "construct_b":b,
                "module_a":MODULES[a],
                "module_b":MODULES[b],
                "same_module":MODULES[a]==MODULES[b],
                "cosine_similarity":cosine(fp[a],fp[b]),
            })
    return rows

def stat_for_labels(pair_rows, mapping):
    wi=[]; be=[]
    for r in pair_rows:
        v=r["cosine_similarity"]
        if mapping[r["construct_a"]]==mapping[r["construct_b"]]:
            wi.append(v)
        else:
            be.append(v)
    return sum(wi)/len(wi)-sum(be)/len(be), sum(wi)/len(wi), sum(be)/len(be)

def unique_label_assignments():
    # exact distinct partitions with fixed sizes: 1 presentation, 3 colour, 2 head, 3 involucre
    idx=range(len(CORE))
    for pres in itertools.combinations(idx,1):
        rem1=[i for i in idx if i not in pres]
        for colour in itertools.combinations(rem1,3):
            rem2=[i for i in rem1 if i not in colour]
            for head in itertools.combinations(rem2,2):
                inv=[i for i in rem2 if i not in head]
                m={}
                for i in pres:m[CORE[i]]="presentation"
                for i in colour:m[CORE[i]]="colour"
                for i in head:m[CORE[i]]="head_form"
                for i in inv:m[CORE[i]]="involucre_armature"
                yield m

def module_centroids(fp):
    out={}
    for m in sorted(set(MODULES.values())):
        cs=[c for c in CORE if MODULES[c]==m]
        vals={b:sum(fp[c][b] for c in cs)/len(cs) for b in ENV_BLOCKS}
        n=math.sqrt(sum(v*v for v in vals.values()))
        out[m]={b:(v/n if n else 0.0) for b,v in vals.items()}
    return out

def summarize(path: Path, scale: str):
    raw,fp=fingerprints(path)
    pairs=pair_table(fp)
    obs,wi,be=stat_for_labels(pairs,MODULES)
    null=[]
    for mapping in unique_label_assignments():
        null.append(stat_for_labels(pairs,mapping)[0])
    assert len(null)==1680
    p=(sum(x>=obs-1e-15 for x in null))/len(null)  # exact, observed included
    cent=module_centroids(fp)
    centroid_pairs=[]
    mods=sorted(cent)
    for i,a in enumerate(mods):
        for b in mods[i+1:]:
            centroid_pairs.append({
                "module_a":a,"module_b":b,
                "cosine_similarity":cosine(cent[a],cent[b]),
                "cosine_distance":1-cosine(cent[a],cent[b]),
            })
    rows=[]
    for c in CORE:
        row={"construct_id":c,"module":MODULES[c]}
        row.update({f"raw_{b}":raw[c][b] for b in ENV_BLOCKS})
        row.update({f"norm_{b}":fp[c][b] for b in ENV_BLOCKS})
        row["strongest_block"]=max(raw[c],key=raw[c].get)
        rows.append(row)
    return {
        "scale":scale,
        "constructs":9,
        "modules":4,
        "within_pairs":sum(r["same_module"] for r in pairs),
        "between_pairs":sum(not r["same_module"] for r in pairs),
        "mean_within_module_cosine":wi,
        "mean_between_module_cosine":be,
        "within_minus_between_cosine":obs,
        "exact_label_assignments":len(null),
        "exact_one_sided_p":p,
        "null_min":min(null),"null_max":max(null),
        "signature_rows":rows,
        "pairwise":pairs,
        "module_centroids":cent,
        "module_centroid_pairs":centroid_pairs,
    }

def write_csv(path,rows):
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--among",type=Path,required=True)
    ap.add_argument("--within",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    args=ap.parse_args()
    args.out_dir.mkdir(parents=True,exist_ok=True)
    a=summarize(args.among,"among_taxon")
    w=summarize(args.within,"within_taxon")
    report={
        "analysis_id":"aza3_fig1_ecological_modularity_v1",
        "primary_scale":"among_taxon",
        "primary_estimand":"mean within-module cosine similarity minus mean between-module cosine similarity of six-block normalized ecological fingerprints",
        "module_labels":"predeclared biological modules from Azami v3; independent of environmental outcomes",
        "environment_blocks":ENV_BLOCKS,
        "among_taxon":{k:v for k,v in a.items() if k not in {"signature_rows","pairwise"}},
        "within_taxon":{k:v for k,v in w.items() if k not in {"signature_rows","pairwise"}},
        "claim_boundary":"Environmental association fingerprint modularity only. Not selection, adaptation, plasticity, developmental modularity, genetic modularity, or causal ecological demand.",
    }
    (args.out_dir/"fig1_ecological_modularity_report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    write_csv(args.out_dir/"fig1_ecological_fingerprints_among.csv",a["signature_rows"])
    write_csv(args.out_dir/"fig1_ecological_fingerprints_within.csv",w["signature_rows"])
    write_csv(args.out_dir/"fig1_pairwise_cosine_among.csv",a["pairwise"])
    write_csv(args.out_dir/"fig1_pairwise_cosine_within.csv",w["pairwise"])
    print(json.dumps(report,indent=2))

if __name__=="__main__":
    main()
