#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, json, math, random, statistics
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
    except Exception: return None

def intval(x):
    try: return int(x)
    except Exception: return None

def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024*1024), b""): h.update(chunk)
    return h.hexdigest()

def sign_flip_p(ds, exact_max_n=20, mc_draws=1000000, seed=20261006):
    n=len(ds); obs=sum(ds)
    if n==0: return None,"none"
    if n<=exact_max_n:
        ge=0; total=1<<n
        for mask in range(total):
            s=sum(d if (mask>>i)&1 else -d for i,d in enumerate(ds))
            if s>=obs-1e-15: ge+=1
        return ge/total,f"exact_{total}_sign_flips"
    rng=random.Random(seed); ge=1; total=mc_draws+1
    for _ in range(mc_draws):
        s=sum(d if rng.random()<0.5 else -d for d in ds)
        if s>=obs: ge+=1
    return ge/total,f"monte_carlo_{mc_draws}_plus_observed_seed_{seed}"

def med(xs): return statistics.median(xs) if xs else None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--scores",required=True,type=Path)
    ap.add_argument("--prediction-registry",required=True,type=Path)
    ap.add_argument("--unlock-receipt",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()

    pred_hash=sha256(args.prediction_registry)
    receipt=json.loads(args.unlock_receipt.read_text(encoding="utf-8"))
    if receipt.get("status")!="UNLOCKED_FOR_SCORING":
        raise SystemExit("unlock receipt status is not UNLOCKED_FOR_SCORING")
    if (receipt.get("prediction_registry_sha256") or "").strip().lower()!=pred_hash:
        raise SystemExit("prediction registry SHA256 does not match unlock receipt")
    frozen_at=receipt.get("prediction_registry_frozen_at")
    opened_at=receipt.get("validation_phenotypes_opened_at")
    if not frozen_at or not opened_at or str(opened_at)<=str(frozen_at):
        raise SystemExit("invalid prediction-freeze / phenotype-unlock chronology")

    preds=list(csv.DictReader(args.prediction_registry.open(encoding="utf-8-sig")))
    scores=list(csv.DictReader(args.scores.open(encoding="utf-8-sig")))
    if not preds or not scores: raise SystemExit("empty prediction or scoring registry")
    pred_by_id={r["individual_id"]:r for r in preds}
    if len(pred_by_id)!=len(preds): raise SystemExit("duplicate individual_id in prediction registry")

    eligible=[]; excluded=[]
    for r in scores:
        iid=(r.get("individual_id") or "").strip(); reason=[]; p=pred_by_id.get(iid)
        if p is None: reason.append("missing_preunblinding_prediction")
        if (r.get("prediction_registry_sha256") or "").strip().lower()!=pred_hash:
            reason.append("score_row_prediction_hash_mismatch")
        if (r.get("taxonomic_confidence") or "").strip().casefold()!="high": reason.append("taxonomic_confidence_not_high")
        if bval(r.get("reference_stable")) is not True: reason.append("reference_not_stable")
        if (r.get("technical_exclusion_reason") or "").strip(): reason.append("technical_exclusion")
        nums=[fnum(r.get(k)) for k in ("model_F_joint_log_score","model_M_joint_log_score","model_R_joint_log_score")]
        if any(v is None for v in nums): reason.append("missing_joint_log_score")
        if reason:
            excluded.append({"individual_id":iid,"reason":";".join(reason)}); continue
        ac=intval(p.get("admitted_construct_count")); am=intval(p.get("admitted_multi_trait_module_count"))
        if ac is None or am is None:
            excluded.append({"individual_id":iid,"reason":"missing_preunblinding_construct_admission_count"}); continue
        f_joint,m_joint,r_joint=nums
        module_deltas={}
        for module in ("colour","head_form","involucre"):
            mv=fnum(r.get(f"module_{module}_M_log_score")); rv=fnum(r.get(f"module_{module}_R_log_score"))
            if mv is not None and rv is not None: module_deltas[module]=rv-mv
        eligible.append({
            "individual_id":iid,"d_RM":r_joint-m_joint,"d_MF":m_joint-f_joint,"d_RF":r_joint-f_joint,
            "module_deltas":module_deltas,"admitted_construct_count":ac,"admitted_multi_trait_module_count":am,
            "mosaic_registered":bval(p.get("prospective_mosaic_registered")) is True,
            "colour_correct":bval(r.get("colour_prediction_correct")) is True,
            "orientation_correct":bval(r.get("orientation_prediction_correct")) is True
        })

    n=len(eligible); d_rm=[x["d_RM"] for x in eligible]; d_mf=[x["d_MF"] for x in eligible]; d_rf=[x["d_RF"] for x in eligible]
    delta_rm=sum(d_rm); delta_mf=sum(d_mf); delta_rf=sum(d_rf)
    p_rm,m_rm=sign_flip_p(d_rm); p_mr,m_mr=sign_flip_p([-x for x in d_rm])
    p_mf,m_mf=sign_flip_p(d_mf); p_fm,m_fm=sign_flip_p([-x for x in d_mf]); p_fr,m_fr=sign_flip_p([-x for x in d_rf])
    cc={x["admitted_construct_count"] for x in eligible}; mm={x["admitted_multi_trait_module_count"] for x in eligible}
    construct_count=next(iter(cc)) if len(cc)==1 else None; multi_module_count=next(iter(mm)) if len(mm)==1 else None
    full_module_admission=(construct_count is not None and multi_module_count is not None and construct_count>=5 and multi_module_count>=2)
    module_totals={}
    for module in ("colour","head_form","involucre"):
        vals=[x["module_deltas"][module] for x in eligible if module in x["module_deltas"]]
        if len(vals)==n and n>0: module_totals[module]=sum(vals)
    stable_nonnegative_modules=sum(v>=0 for v in module_totals.values())
    registered=[x for x in eligible if x["mosaic_registered"]]
    correct=[x for x in registered if x["colour_correct"] and x["orientation_correct"]]

    if n<20:
        decision="NOT_IDENTIFIABLE"; reason="untouched_validation_n_below_20"
    else:
        r_over_m=delta_rm>0 and p_rm is not None and p_rm<0.05
        m_over_r=delta_rm<0 and p_mr is not None and p_mr<0.05
        m_over_f=delta_mf>0 and p_mf is not None and p_mf<0.05
        f_over_m=delta_mf<0 and p_fm is not None and p_fm<0.05
        f_over_r=delta_rf<0 and p_fr is not None and p_fr<0.05
        if r_over_m:
            if full_module_admission and stable_nonnegative_modules>=2 and len(correct)>=3:
                decision="FIG3_R_STRONG"; reason="R_beats_M_with_multi_module_support_and_prospective_mosaics"
            else:
                decision="FIG3_R_PARTIAL"; reason="R_beats_M_but_full_module_or_mosaic_requirement_not_met"
        elif m_over_r and m_over_f:
            decision="FIG3_MODULE_INHERITANCE_SUPPORTIVE"; reason="M_beats_R_and_F"
        elif f_over_m and f_over_r:
            decision="FIG3_WHOLE_ORGAN_SUPPORTIVE"; reason="F_beats_M_and_R"
        else:
            decision="NOT_IDENTIFIABLE"; reason="no_model_meets_frozen_predictive_superiority_rule"

    out={
      "result_version":"aza3_fig3_FMR_heldout_prediction_score_v3",
      "prediction_registry_sha256":pred_hash,"prediction_registry_frozen_at":frozen_at,"validation_phenotypes_opened_at":opened_at,
      "eligible_validation_n":n,"excluded_n":len(excluded),"admitted_construct_count":construct_count,
      "admitted_multi_trait_module_count":multi_module_count,"full_module_admission":full_module_admission,
      "delta_RM":delta_rm,"delta_MF":delta_mf,"delta_RF":delta_rf,"mean_d_RM":delta_rm/n if n else None,"median_d_RM":med(d_rm),
      "p_R_over_M_one_sided":p_rm,"p_M_over_R_one_sided":p_mr,"p_M_over_F_one_sided":p_mf,
      "p_F_over_M_one_sided":p_fm,"p_F_over_R_one_sided":p_fr,
      "sign_flip_methods":{"R_over_M":m_rm,"M_over_R":m_mr,"M_over_F":m_mf,"F_over_M":m_fm,"F_over_R":m_fr},
      "module_R_minus_M_log_score_totals":module_totals,"nonnegative_complete_module_count":stable_nonnegative_modules,
      "prospective_mosaic_registered_n":len(registered),"prospective_mosaic_correct_both_n":len(correct),
      "decision":decision,"decision_reason":reason,"excluded":excluded,
      "claim_boundary":"Held-out F-M-R prediction only. Adaptation, causal variants, haplotype age and SV reuse require separate evidence."
    }
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__": main()
