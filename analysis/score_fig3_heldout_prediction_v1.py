#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math, random
from pathlib import Path

TRUE={"1","true","t","yes","y"}
FALSE={"0","false","f","no","n"}

def bval(x):
    s=(x or "").strip().casefold()
    if s in TRUE: return True
    if s in FALSE: return False
    return None

def fnum(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None

def sign_flip_p(ds, exact_max_n=20, mc_draws=1000000, seed=20261006):
    n=len(ds)
    obs=sum(ds)
    if n==0:
        return None, "none"
    if n<=exact_max_n:
        ge=0
        total=1<<n
        for mask in range(total):
            s=0.0
            for i,d in enumerate(ds):
                s += d if (mask>>i)&1 else -d
            if s >= obs - 1e-15:
                ge += 1
        return ge/total, f"exact_{total}_sign_flips"
    rng=random.Random(seed)
    ge=1
    total=mc_draws+1
    for _ in range(mc_draws):
        s=sum(d if rng.random()<0.5 else -d for d in ds)
        if s >= obs:
            ge += 1
    return ge/total, f"monte_carlo_{mc_draws}_plus_observed_seed_{seed}"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()

    rows=list(csv.DictReader(args.input.open(encoding="utf-8-sig")))
    if not rows:
        raise SystemExit("empty registry")

    eligible=[]
    excluded=[]
    for r in rows:
        reason=[]
        if (r.get("split") or "").strip().upper()!="VALIDATION":
            reason.append("not_validation")
        if bval(r.get("phenotype_labels_opened")) is not True:
            reason.append("phenotype_labels_not_opened")
        if (r.get("taxonomic_confidence") or "").strip().casefold()!="high":
            reason.append("taxonomic_confidence_not_high")
        if bval(r.get("reference_stable")) is not True:
            reason.append("reference_not_stable")
        if (r.get("technical_exclusion_reason") or "").strip():
            reason.append("technical_exclusion")
        nums=[fnum(r.get(k)) for k in (
            "model_F_joint_log_score","model_R_joint_log_score",
            "colour_F_log_score","colour_R_log_score",
            "orientation_F_log_score","orientation_R_log_score"
        )]
        if any(v is None for v in nums):
            reason.append("missing_log_score")
        if reason:
            excluded.append({"individual_id":r.get("individual_id",""),"reason":";".join(reason)})
            continue
        f_joint,r_joint,c_f,c_r,o_f,o_r=nums
        eligible.append({
            "individual_id":r.get("individual_id",""),
            "d_joint":r_joint-f_joint,
            "d_colour":c_r-c_f,
            "d_orientation":o_r-o_f,
            "mosaic_registered":bval(r.get("prospective_mosaic_registered")) is True,
            "colour_correct":bval(r.get("colour_prediction_correct")) is True,
            "orientation_correct":bval(r.get("orientation_prediction_correct")) is True,
        })

    ds=[r["d_joint"] for r in eligible]
    n=len(ds)
    delta=sum(ds)
    colour_delta=sum(r["d_colour"] for r in eligible)
    orientation_delta=sum(r["d_orientation"] for r in eligible)
    p,method=sign_flip_p(ds)

    registered=[r for r in eligible if r["mosaic_registered"]]
    correct_mosaics=[
        r for r in registered
        if r["colour_correct"] and r["orientation_correct"]
    ]

    if n < 20:
        decision="NOT_IDENTIFIABLE"
        reason="untouched_validation_n_below_20"
    else:
        r_predictive=(delta>0 and p is not None and p<0.05 and colour_delta>=0 and orientation_delta>=0)
        if r_predictive and len(correct_mosaics)>=3:
            decision="FIG3_R_STRONG"
            reason="heldout_joint_prediction_plus_at_least_3_correct_prospective_mosaics"
        elif r_predictive:
            decision="FIG3_R_PARTIAL"
            reason="heldout_prediction_supports_R_but_mosaic_kill_shot_not_met"
        elif delta<=0:
            decision="FIG3_F_SUPPORTIVE"
            reason="shared_state_model_equal_or_better_on_joint_heldout_log_score"
        else:
            decision="NOT_IDENTIFIABLE"
            reason="R_direction_positive_but_frozen_inference_criterion_not_met"

    out={
        "result_version":"aza3_fig3_heldout_prediction_score_v1",
        "eligible_validation_n":n,
        "excluded_n":len(excluded),
        "delta_joint":delta,
        "mean_d_joint":(delta/n if n else None),
        "median_d_joint":(sorted(ds)[n//2] if n and n%2==1 else ((sorted(ds)[n//2-1]+sorted(ds)[n//2])/2 if n else None)),
        "colour_delta_total":colour_delta,
        "orientation_delta_total":orientation_delta,
        "sign_flip_one_sided_p":p,
        "sign_flip_method":method,
        "prospective_mosaic_registered_n":len(registered),
        "prospective_mosaic_correct_both_n":len(correct_mosaics),
        "decision":decision,
        "decision_reason":reason,
        "excluded":excluded,
        "claim_boundary":"This scorer evaluates the frozen held-out prediction contract only. Adaptation, causality, haplotype age and structural-variant reuse require separate evidence."
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
