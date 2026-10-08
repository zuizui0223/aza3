#!/usr/bin/env python3
"""Stage-separated parasitoid timing and local next-generation Cirsium head damage.

This is a deliberately minimal conditional model, not a fitted ecosystem
or genetic-selection model. Late parasitoids change local recruitment AFTER
the current head's seed loss is already irreversible.

A0 = I / (1 - K*(1-q0)) is the equilibrium **potential seed-feeder
arrivals** in a baseline patch with baseline habitat/architecture.
When a barrier is introduced at t=0, newly established larvae = r * A1[t].
Potential arrivals next generation:
 A1[t+1] = I + K * r * (1-q1) * A1[t]
with K a local next-generation replacement/retention coefficient. Baseline
arrivals A0[t] remain at the original equilibrium. Current seed damage is
D0=d*A0 and D1[t]=d*r*A1[t], with *zero effective parasitoid killing
before irreversible seed damage in both conditions*.

The stable-equilibrium reversal D1*/D0*>1 requires:
   r*K*(q0-q1) > 1-r,
provided 0<r<=1, 0<=q<=1, I>0 and both reproduction multipliers <1.
The sign of this PATCH-LEVEL long-run ecological feedback does not
identify same-plant fitness, genetic selection or adaptive trait evolution.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path


def finite_between(x, lower, upper, name):
    if not math.isfinite(x) or not (lower <= x <= upper):
        raise ValueError(f"{name} must be finite and in [{lower},{upper}]")


def evaluate(*, establishment_ratio: float, baseline_late_kill: float,
             barrier_late_kill: float, local_replacement: float,
             immigration: float = 1.0, damage_per_established: float = 1.0,
             generations: int = 12):
    r = establishment_ratio
    q0 = baseline_late_kill
    q1 = barrier_late_kill
    K = local_replacement
    I = immigration
    d = damage_per_established
    finite_between(r, 0, 1, "establishment_ratio")
    finite_between(q0, 0, 1, "baseline_late_kill")
    finite_between(q1, 0, 1, "barrier_late_kill")
    for n,v in (("local_replacement",K),("immigration",I),
                ("damage_per_established",d)):
        if not math.isfinite(v) or (v<=0 if n in {"immigration", "damage_per_established"} else v<0):
            raise ValueError(f"{n} must be finite and nonnegative; immigration and damage positive")
    if not isinstance(generations,int) or generations<1 or generations>200:
        raise ValueError("generations must be 1..200")

    base_multiplier=K*(1-q0)
    defended_multiplier=K*r*(1-q1)
    stable=(base_multiplier<1 and defended_multiplier<1)
    if not stable:
        return {
          "status": "NO_STATIONARY_COMPARISON_SUPERCRITICAL_LINEAR_MODEL",
          "base_multiplier": base_multiplier,
          "defended_multiplier": defended_multiplier,
          "model_limit": "Density dependence and movement are omitted; no numerical fitness extrapolation.",
          "biological_data": "NONE",
        }

    control_arrivals=I/(1-base_multiplier)
    defended_arrivals=control_arrivals
    records=[]
    for t in range(generations+1):
        d0=d*control_arrivals
        d1=d*r*defended_arrivals
        records.append({
            "generation": t,
            "control_potential_arrivals": control_arrivals,
            "defended_potential_arrivals": defended_arrivals,
            "control_current_head_damage": d0,
            "defended_current_head_damage": d1,
            "defended_to_control_damage_ratio": d1/d0 if d0>0 else None,
        })
        defended_arrivals=I+defended_multiplier*defended_arrivals

    steady_ratio=r*(1-base_multiplier)/(1-defended_multiplier)
    reversal=steady_ratio>1+1e-12
    threshold=(1-r)/(r*K) if r>0 and K>0 else None
    sufficient_condition=(r*K*(q0-q1) > (1-r)+1e-12)
    assert reversal == sufficient_condition
    assert abs(records[0]["defended_to_control_damage_ratio"]-r)<1e-12
    first_cross=next((v["generation"] for v in records
                      if v["defended_to_control_damage_ratio"]>1+1e-12),None)
    return {
        "status": "CONDITIONAL_STABLE_ECOLOGICAL_PATCH_MODEL_NOT_FIT",
        "assumptions": [
            "late killing prevents feeder recruitment but not already irreversible current seed damage",
            "same external immigration each generation for control and defended patch",
            "fixed local replacement and density-independent survival; no plant depletion or saturation",
            "same potential seed-feeder arrival pool just before barrier introduced",
            "resident next generation returns to the same kind of defended patch",
            "constant per-larva seed damage, no pollen/flower production or barrier cost",
        ],
        "base_multiplier":base_multiplier,
        "defended_multiplier":defended_multiplier,
        "first_generation_relative_damage":r,
        "stable_equilibrium_relative_damage":steady_ratio,
        "stable_equilibrium_damage_reversal":reversal,
        "necessary_and_sufficient_late_kill_difference_threshold":
            threshold,
        "observed_late_kill_difference":q0-q1,
        "first_cross_over_generation_within_simulation":first_cross,
        "trajectory":records,
        "biological_data":"NONE",
        "inferential_boundary":"Patch-level enemy carryover, NOT individual thistle fitness or adaptation.",
    }


def test_synthetic():
    a=evaluate(establishment_ratio=.7,baseline_late_kill=.6,
               barrier_late_kill=.1,local_replacement=1.2,
               generations=20)
    assert a["status"].startswith("CONDITIONAL_STABLE")
    assert abs(a["base_multiplier"]-.48)<1e-12
    assert abs(a["defended_multiplier"]-.756)<1e-12
    assert abs(a["stable_equilibrium_relative_damage"]-1.4918032786885245)<1e-10
    assert a["stable_equilibrium_damage_reversal"]
    assert a["first_cross_over_generation_within_simulation"]==2
    assert abs(a["first_generation_relative_damage"]-.7)<1e-12
    no_return=evaluate(establishment_ratio=.7,baseline_late_kill=.6,
                       barrier_late_kill=.1,local_replacement=0.,
                       generations=20)
    assert no_return["first_cross_over_generation_within_simulation"] is None
    assert abs(no_return["stable_equilibrium_relative_damage"]-.7)<1e-12
    equal_late=evaluate(establishment_ratio=.7,baseline_late_kill=.6,
                        barrier_late_kill=.6,local_replacement=1.2)
    assert not equal_late["stable_equilibrium_damage_reversal"]
    no_eggs=evaluate(establishment_ratio=0,baseline_late_kill=.6,
                     barrier_late_kill=0,local_replacement=1.2)
    assert no_eggs["stable_equilibrium_relative_damage"]==0
    unstable=evaluate(establishment_ratio=1,baseline_late_kill=0,
                      barrier_late_kill=0,local_replacement=2)
    assert unstable["status"].startswith("NO_STATIONARY")
    for x in [(-.1,.4,.4,1),(.7,1.1,.4,1),(.7,.4,.4,-1)]:
        try:
            evaluate(establishment_ratio=x[0], baseline_late_kill=x[1],
                     barrier_late_kill=x[2],local_replacement=x[3])
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid ecological parameter accepted")

    # Check both boundary cases: exact equality and just across it.
    # Choose r=0.7, K=1.2, q0=.6, then q1 = .6 - .3/.84.
    boundary_q1=.6-(.3/.84)
    boundary=evaluate(establishment_ratio=.7,baseline_late_kill=.6,
                      barrier_late_kill=boundary_q1,local_replacement=1.2)
    assert abs(boundary["stable_equilibrium_relative_damage"]-1)<1e-12
    return {
        "status":"PASS_SYNTHETIC_DELAYED_ENEMY_FEEDBACK_TESTS",
        "illustrative_first_reversal_generation":a["first_cross_over_generation_within_simulation"],
        "immediate_ratio":a["first_generation_relative_damage"],
        "steady_ratio":a["stable_equilibrium_relative_damage"],
        "no_local_retention_prevents_reversal":not no_return["stable_equilibrium_damage_reversal"],
        "boundary_ratio":boundary["stable_equilibrium_relative_damage"],
        "empirical_observations":0,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    verified=test_synthetic()
    example=evaluate(establishment_ratio=.7,baseline_late_kill=.6,
                     barrier_late_kill=.1,local_replacement=1.2,
                     generations=12)
    response={
        "version":"capitulum_late_enemy_multigeneration_feedback_v1",
        "asof":"2026-10-08",
        "test_result":verified,
        "illustrative_only":example,
        "direct_evidence": {
            "source":"Vanbergen et al. 2006, Journal of Animal Ecology, 10.1111/j.1365-2656.2006.01099.x",
            "support":"T. conura larvae feed June–July, P. elevatus oviposition peaks early/mid-August; species/grazing interaction differences.",
            "does_not_support":"Exact individual head seed-loss timing; post-loss q measured; local host return; morphology effect; q0=.6 or q1=.1 numerical inputs.",
        },
        "never_conclude":[
            "observed insect generation feedback in the Cirsium field",
            "long-run fitness for a biennial plant that dies after flowering",
            "parasitoid adult emergence equals pre-seed-damage seed rescue",
            "head shape/spines adapted via apparent inter-generational delayed effect",
            "species with bivoltine insects follow exactly annual same-stage timing",
        ],
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(response,ensure_ascii=False,indent=2)+"\n",
                        encoding="utf-8")
    print(json.dumps({
        "test":verified["status"],
        "late_parasitoid_reversal":example["stable_equilibrium_damage_reversal"],
        "first_cross_generation":example["first_cross_over_generation_within_simulation"],
        "empirical_rows":0,
    },indent=2))

if __name__=="__main__":
    main()
