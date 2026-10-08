#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math, random
from pathlib import Path

def sign_flip_p(ds, rng, sign_draws):
    obs=sum(ds)
    ge=1
    total=sign_draws+1
    for _ in range(sign_draws):
        s=sum(d if rng.random()<0.5 else -d for d in ds)
        if s>=obs: ge+=1
    return ge/total

def power(delta,n,reps,seed,sign_draws):
    rng=random.Random(seed + int(delta*10000) + n*100000)
    hit=0
    for _ in range(reps):
        ds=[rng.gauss(delta,1.0) for _ in range(n)]
        if sum(ds)>0 and sign_flip_p(ds,rng,sign_draws)<0.05:
            hit+=1
    return hit/reps

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--reps",type=int,default=10000)
    ap.add_argument("--seed",type=int,default=20261006)
    ap.add_argument("--usable-fraction",type=float,default=None)
    ap.add_argument("--sign-draws",type=int,default=1999)
    ap.add_argument("--output-json",required=True,type=Path)
    ap.add_argument("--output-csv",required=True,type=Path)
    args=ap.parse_args()

    deltas=[0.20,0.30,0.40,0.50,0.70]
    ns=[20,24,30,36,40,50,60]
    rows=[]
    for d in deltas:
        for n in ns:
            p=power(d,n,args.reps,args.seed,args.sign_draws)
            collected=None
            if args.usable_fraction is not None:
                if not (0 < args.usable_fraction <= 1):
                    raise SystemExit("usable fraction must be in (0,1]")
                collected=math.ceil(n/args.usable_fraction)
            rows.append({
                "delta":d,"validation_usable_n":n,"power":p,
                "usable_fraction":args.usable_fraction,
                "validation_collected_n_if_fraction_given":collected
            })

    primary=[r for r in rows if abs(r["delta"]-0.40)<1e-12 and r["power"]>=0.80]
    secondary=[r for r in rows if abs(r["delta"]-0.50)<1e-12 and r["power"]>=0.90]
    min_primary=min((r["validation_usable_n"] for r in primary),default=None)
    min_secondary=min((r["validation_usable_n"] for r in secondary),default=None)

    out={
        "result_version":"aza3_fig3_estimand_level_power_v1",
        "simulation_reps":args.reps,
        "seed":args.seed,\n        "monte_carlo_sign_flip_draws_per_power_replicate":args.sign_draws,
        "standardized_delta_grid":deltas,
        "validation_n_grid":ns,
        "primary_requirement":{"delta":0.40,"power_gte":0.80,"minimum_validation_usable_n":min_primary},
        "secondary_benchmark":{"delta":0.50,"power_gte":0.90,"minimum_validation_usable_n":min_secondary},
        "usable_fraction":args.usable_fraction,
        "boundary":"Estimand-level validation-discrimination power only. Power uses Monte Carlo random-sign nulls for computational feasibility; the final held-out scorer retains its frozen exact/Monte-Carlo rule. TRAIN genomic-discovery adequacy must be assessed separately from genotype-only pilot quantities."
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    with args.output_csv.open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
