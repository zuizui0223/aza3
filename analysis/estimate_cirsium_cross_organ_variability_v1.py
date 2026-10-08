#!/usr/bin/env python3
"""Independent Cirsium organ-variation analysis using published XLSX with artifact_tool.

No modification to the frozen azami (GEB) manuscript or original XLSX.
Published source: Michálková et al. 2023, 10.1007/s00606-023-01854-2,
Online Resource 6 (CC BY 4.0).  Conditional plant-level bootstrap only.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import random
import statistics
from pathlib import Path

from artifact_tool import Blob, SpreadsheetFile

EXPECTED_SHA = "c06423cc6f9800ba4775a4428afa072ba1b3dc6fb48ccd05fdb7e42a2f03ee79"
TAXA = ["Cirsium acaulon", "Cirsium bertolonii - Monte Sagro", "Cirsium erisithales"]
SEED = 20261008
PAIRS = [
    ("leaf_length_involucre_length", "Lf_L", "Inv_L"),
    ("leaf_width_involucre_width", "Lf_W", "Inv_W"),
    ("leaf_aspect_involucre_aspect", "Lf_Sh", "Inv_Sh"),
]


def cv(a):
    if len(a) < 3 or statistics.mean(a) <= 0:
        raise ValueError("CV undefined for insufficient or non-positive values")
    return statistics.stdev(a) / statistics.mean(a)


def log_sd(a):
    if any(x <= 0 for x in a):
        raise ValueError("Log-SD requires positive values")
    return statistics.stdev([math.log(x) for x in a])


def quantile(arr, p):
    data = sorted(arr)
    k = (len(data) - 1) * p
    lo, hi = math.floor(k), math.ceil(k)
    return data[lo] + (data[hi] - data[lo]) * (k - lo)


def geo_mean(vals):
    return math.exp(statistics.mean(math.log(x) for x in vals))


def load_rows(path):
    if hashlib.sha256(Path(path).read_bytes()).hexdigest() != EXPECTED_SHA:
        raise ValueError("Source XLSX SHA differs from frozen publisher file")
    wb = SpreadsheetFile.import_xlsx(Blob.load(str(path)))
    sheet = wb.worksheets.get_item_at(0)
    if sheet.name != "S5 Morphometry":
        raise ValueError("Unexpected sheet name")
    matrix = sheet.get_range("A9:AE181").values
    codes = matrix[1]
    rows = []
    for r in matrix[2:]:
        if not r[2] or str(r[3]) not in {"f", "h"}:
            continue
        row = {"taxon": str(r[1]), "plant_id": str(r[2]), "sex": str(r[3])}
        row.update({str(codes[i]): float(r[i]) if isinstance(r[i], (int,float)) else None
                    for i in range(4, 31)})
        rows.append(row)
    if len(rows) != 129 or len({r["plant_id"] for r in rows}) != 129:
        raise ValueError("Published cohort count or individual identities differ")
    groups = {(t, sex): [r for r in rows if r["taxon"] == t and r["sex"] == sex]
              for t in TAXA for sex in ("f","h")}
    if sorted(map(len, groups.values())) != [10,10,10,24,30,30]:
        raise ValueError("Parental source cohort unexpectedly changed")
    return rows, groups


def pair(d, x, y):
    a = [r[x] for r in d]
    b = [r[y] for r in d]
    return {"n":len(d),"leaf_cv":cv(a),"involucre_cv":cv(b),
            "cv_ratio":cv(a)/cv(b),"logsd_ratio":log_sd(a)/log_sd(b)}


def gradient(d):
    v = {k: cv([r[k] for r in d]) for k in ("Lf_L","Inv_L","Cor_L")}
    return {"leaf_cv":v["Lf_L"],"involucre_cv":v["Inv_L"],
            "corolla_cv":v["Cor_L"],"ordered":v["Lf_L"]>v["Inv_L"]>v["Cor_L"],
            "leaf_involucre":v["Lf_L"]/v["Inv_L"],
            "involucre_corolla":v["Inv_L"]/v["Cor_L"],
            "leaf_corolla":v["Lf_L"]/v["Cor_L"]}


def analyze(source, draws):
    all_rows, groups = load_rows(source)
    rng = random.Random(SEED)
    outputs = {}
    for label, a, b in PAIRS:
        observed = {f"{k[0]}|{k[1]}": pair(v,a,b) for k,v in groups.items()}
        boot = []
        for _ in range(draws):
            ss = [pair([v[rng.randrange(len(v))] for _ in range(len(v))],a,b)
                  for k,v in groups.items()]
            boot.append((geo_mean([z["cv_ratio"] for z in ss]),
                         geo_mean([z["logsd_ratio"] for z in ss])))
        outputs[label] = {
            "n_nonhybrid_plants":114,
            "all_six_strata_cv_ratio_gt_one":all(z["cv_ratio"]>1 for z in observed.values()),
            "mean_cv_ratio":geo_mean([z["cv_ratio"] for z in observed.values()]),
            "cv_boot_ci95":[quantile([z[0] for z in boot],.025),
                            quantile([z[0] for z in boot],.975)],
            "mean_logsd_ratio":geo_mean([z["logsd_ratio"] for z in observed.values()]),
            "logsd_boot_ci95":[quantile([z[1] for z in boot],.025),
                               quantile([z[1] for z in boot],.975)],
            "strata":observed,
        }
    # Reseed, because this hypothesis was a separately evaluated 3-tier check.
    rng = random.Random(SEED)
    original = {f"{t}|{sex}": gradient(d) for (t,sex),d in groups.items()}
    vals = collections.defaultdict(list)
    all6_count = 0
    for _ in range(draws):
        m = [gradient([d[rng.randrange(len(d))] for _ in range(len(d))])
             for k,d in groups.items()]
        all6_count += int(all(z["ordered"] for z in m))
        for key in ("leaf_involucre","involucre_corolla","leaf_corolla"):
            vals[key].append(geo_mean([z[key] for z in m]))
    gradients = {
        key: {"cv_ratio":geo_mean([z[key] for z in original.values()]),
              "bootstrap_ci95":[quantile(vals[key],.025), quantile(vals[key],.975)]}
        for key in ("leaf_involucre","involucre_corolla","leaf_corolla")
    }
    return {
        "source_doi":"10.1007/s00606-023-01854-2",
        "source_xlsx_sha256":EXPECTED_SHA,
        "n_total_published":len(all_rows), "n_nonhybrid_plants":114,
        "species":TAXA, "n_strata":6, "bootstrap_draws":draws, "seed":SEED,
        "cohort_notice":"Parental plants only; sex stratified; hybrid plants n=15 excluded",
        "organ_pair_contrasts":outputs,
        "leaf_involucre_corolla_gradient":{
            "all_six_strata_ordered":all(z["ordered"] for z in original.values()),
            "bootstrap_probability_all_six_ordered":all6_count/draws,
            "cv_ratios":gradients, "strata":original
        },
        "claim_ceiling":"Measured within-lineage phenotypic variability, NOT selection, genetic constraints, evolutionary rates or causal ecological adaptation",
        "uncertainty_boundary":"Plant resampling conditional on 3 taxa and observed localities, not phylogenetic/site resampling."
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--bootstrap",type=int,default=3000)
    args = ap.parse_args()
    if args.bootstrap < 100: ap.error("Bootstrap draws must be at least 100")
    result = analyze(args.xlsx,args.bootstrap)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":"PASS","n":result["n_nonhybrid_plants"],
          "ratios":result["leaf_involucre_corolla_gradient"]["cv_ratios"]},
          ensure_ascii=False,indent=2))

if __name__ == "__main__":
    main()
