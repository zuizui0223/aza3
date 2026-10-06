#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math, random, statistics
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

def intval(x):
    try:
        return int(x)
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

def med(xs):
    return statistics.median(xs) if xs else None

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
        if bval(r.get("construct_admission_frozen")) is not True:
            reason.append("construct_admission_not_frozen")
        if (r.get("taxonomic_confidence") or "").strip().casefold()!="high":
            reason.append("taxonomic_confidence_not_high")
        if bval(r.get("reference_stable")) is not True:
            reason.append("reference_not_stable")
        if (r.get("technical_exclusion_reason") or "").strip():
            reason.append("technical_exclusion")

        nums=[fnum(r.get(k)) for k in (
            "model_F_joint_log_score","model_M_joint_log_score","model_R_joint_log_score"
        )]
        ac=intval(r.get("admitted_construct_count"))
        am=intval(r.get("admitted_multi_trait_module_count"))
        if any(v is None for v in nums):
            reason.append("missing_joint_log_score")
        if ac is None or am is None:
            reason.append("missing_construct_admission_count")

        if reason:
            excluded.append({"individual_id":r.get("individual_id",""),"reason":";".join(reason)})
            continue

        f_joint,m_joint,r_joint=nums
        module_deltas={}
        for module in ("colour","head_form","involucre"):
            mv=fnum(r.get(f"module_{module}_M_log_score"))
            rv=fnum(r.get(f"module_{module}_R_log_score"))
            if mv is not None and rv is not None:
                module_deltas[module]=rv-mv

        eligible.append({
            "individual_id":r.get("individual_id",""),
            "d_RM":r_joint-m_joint,
            "d_MF":m_joint-f_joint,
            "d_RF":r_joint-f_joint,
            "module_deltas":module_deltas,
            "admitted_construct_count":ac,
            "admitted_multi_trait_module_count":am,
            "mosaic_registered":bval(r.get("prospective_mosaic_registered")) is True,
            "colour_correct":bval(r.get("colour_prediction_correct")) is True,
            "orientation_correct":bval(r.get("orientation_prediction_correct")) is True,
        })

    n=len(eligible)
    d_rm=[r["d_RM"] for r in eligible]
    d_mf=[r["d_MF"] for r in eligible]
    d_rf=[r["d_RF"] for r in eligible]

    delta_rm=sum(d_rm)
    delta_mf=sum(d_mf)
    delta_rf=sum(d_rf)

    p_rm,method_rm=sign_flip_p(d_rm)
    p_mr,method_mr=sign_flip_p([-x for x in d_rm])
    p_mf,method_mf=sign_flip_p(d_mf)
    p_fm,method_fm=sign_flip_p([-x for x in d_mf])
    p_fr,method_fr=sign_flip_p([-x for x in d_rf])

    construct_counts={r["admitted_construct_count"] for r in eligible}
    module_counts={r["admitted_multi_trait_module_count"] for r in eligible}
    construct_count=(next(iter(construct_counts)) if len(construct_counts)==1 else None)
    multi_module_count=(next(iter(module_counts)) if len(module_counts)==1 else None)
    full_module_admission=(
        construct_count is not None and multi_module_count is not None
        and construct_count>=5 and multi_module_count>=2
    )

    module_delta_totals={}
    for module in ("colour","head_form","involucre"):
        vals=[r["module_deltas"][module] for r in eligible if module in r["module_deltas"]]
        if len(vals)==n and n>0:
            module_delta_totals[module]=sum(vals)
    stable_nonnegative_modules=sum(v>=0 for v in module_delta_totals.values())

    registered=[r for r in eligible if r["mosaic_registered"]]
    correct_mosaics=[r for r in registered if r["colour_correct"] and r["orientation_correct"]]

    if n<20:
        decision="NOT_IDENTIFIABLE"
        reason="untouched_validation_n_below_20"
    else:
        r_over_m=(delta_rm>0 and p_rm is not None and p_rm<0.05)
        m_over_r=(delta_rm<0 and p_mr is not None and p_mr<0.05)
        m_over_f=(delta_mf>0 and p_mf is not None and p_mf<0.05)
        f_over_m=(delta_mf<0 and p_fm is not None and p_fm<0.05)
        f_over_r=(delta_rf<0 and p_fr is not None and p_fr<0.05)

        if r_over_m:
            if full_module_admission and stable_nonnegative_modules>=2 and len(correct_mosaics)>=3:
                decision="FIG3_R_STRONG"
                reason="R_beats_module_baseline_with_multi_module_support_and_mosaic_kill_shot"
            else:
                decision="FIG3_R_PARTIAL"
                reason="R_beats_module_baseline_but_full_module_or_mosaic_requirement_not_met"
        elif m_over_r and m_over_f:
            decision="FIG3_MODULE_INHERITANCE_SUPPORTIVE"
            reason="module_model_beats_both_reconfigurable_and_whole_organ_models"
        elif f_over_m and f_over_r:
            decision="FIG3_WHOLE_ORGAN_SUPPORTIVE"
            reason="whole_organ_model_beats_module_and_reconfigurable_models"
        else:
            decision="NOT_IDENTIFIABLE"
            reason="no_model_meets_frozen_predictive_superiority_rule"

    out={
        "result_version":"aza3_fig3_FMR_heldout_prediction_score_v2",
        "eligible_validation_n":n,
        "excluded_n":len(excluded),
        "admitted_construct_count":construct_count,
        "admitted_multi_trait_module_count":multi_module_count,
        "full_module_admission":full_module_admission,
        "delta_RM":delta_rm,
        "delta_MF":delta_mf,
        "delta_RF":delta_rf,
        "mean_d_RM":(delta_rm/n if n else None),
        "median_d_RM":med(d_rm),
        "p_R_over_M_one_sided":p_rm,
        "p_M_over_R_one_sided":p_mr,
        "p_M_over_F_one_sided":p_mf,
        "p_F_over_M_one_sided":p_fm,
        "p_F_over_R_one_sided":p_fr,
        "sign_flip_methods":{
            "R_over_M":method_rm,"M_over_R":method_mr,"M_over_F":method_mf,
            "F_over_M":method_fm,"F_over_R":method_fr
        },
        "module_R_minus_M_log_score_totals":module_delta_totals,
        "nonnegative_complete_module_count":stable_nonnegative_modules,
        "prospective_mosaic_registered_n":len(registered),
        "prospective_mosaic_correct_both_n":len(correct_mosaics),
        "decision":decision,
        "decision_reason":reason,
        "excluded":excluded,
        "claim_boundary":"Held-out F-M-R prediction only. Adaptation, causal variants, haplotype age and SV reuse require separate evidence."
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
