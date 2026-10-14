#!/usr/bin/env python3
"""Conditional Cirsium capitulum parasitoid-sign gate (NOT biological data).

Separate (i) established egg-laying herbivores, (ii) parasitoids that kill a
host BEFORE irreversible seed damage, and (iii) koinobiont parasitoids whose
hosts CONTINUE TO FEED and may consume more than unparasitized larvae.

D_i = E0 * r_i * d * (1-q_i) * (1 + p_i*(g-1))
r_0=1, r_1=r.
q_i: verified predamage host killing fraction, not parasitoid emergence.
p_i: fraction of NOT-predamage-killed hosts under late parasitism.
g: relative per-host SEED damage by later-parasitized vs non-parasitized host.
g is common to both arms solely in this narrow hypothetical model.

No currently accessible focal Cirsium field records measure these parameters.
Xi et al. 2015 DOI 10.1111/1365-2656.12361 contains C. setosum
observational cases, but its randomized field experiment was on Saussurea
nigrescens. Do NOT reuse a difference in mean head damage as per-host g.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

def probability(x, name):
    if not isinstance(x, (int,float)) or not math.isfinite(x) or not 0 <= x <= 1:
        raise ValueError(f"{name} must be finite in [0,1]")

def evaluate(*, r, q0, q1, p0, p1, g, E0=10., d=1.):
    for name,x in (("r",r),("q0",q0),("q1",q1),("p0",p0),("p1",p1)):
        probability(x,name)
    for name,x in (("g",g),("E0",E0),("d",d)):
        if not isinstance(x,(int,float)) or not math.isfinite(x) or x<=0:
            raise ValueError(f"{name} must be finite positive")
    # Keep observed damage totals and exactly named denominators independent.
    D0=E0*d*(1-q0)*(1+p0*(g-1))
    D1=E0*r*d*(1-q1)*(1+p1*(g-1))
    out={
        "status":"CONDITIONAL_TWO_CHANNEL_ACCOUNTING_ONLY",
        "baseline_damage":D0,"barrier_damage":D1,
        "barrier_minus_control_seed_damage":D1-D0,
        "barrier_reduces_current_seed_damage":D1 < D0-1e-12,
        "barrier_increases_current_seed_damage":D1 > D0+1e-12,
        "ratio":D1/D0 if D0>0 else None,
        "baseline_zero_damage_no_ratio":D0==0,
        "theoretical_model_only":True,
    }
    if D0>0:
        # Independent symbolic ratio check for numerical identity.
        A=r*(1-q1)/(1-q0)
        late=(1+p1*(g-1))/(1+p0*(g-1))
        assert abs(out["ratio"]-A*late) <= 1e-10*max(1.,A*late)
        out["early_kill_and_establishment_factor"]=A
        out["feeding_response_factor"]=late
    else:
        out["early_kill_and_establishment_factor"]=None
        out["feeding_response_factor"]=None
    return out

def synthetic_tests():
    # Restrictive case from earlier program: fewer eggs, but much less early
    # parasitoid killing. With no late feeding effect, defence backfires.
    neutral=evaluate(r=.7,q0=.6,q1=.2,p0=.6,p1=.1,g=1)
    feeding=evaluate(r=.7,q0=.6,q1=.2,p0=.6,p1=.1,g=2)
    assert abs(neutral["ratio"]-1.4)<1e-12
    assert abs(feeding["ratio"]-.9625)<1e-12
    assert neutral["barrier_increases_current_seed_damage"]
    assert feeding["barrier_reduces_current_seed_damage"]
    # Closed-form threshold at which parasitism-stimulated host feeding
    # exactly cancels lost EARLY natural enemy control in the illustration.
    A=neutral["early_kill_and_establishment_factor"]
    t=(A-1)/(.6-A*.1)
    crossover=evaluate(r=.7,q0=.6,q1=.2,p0=.6,p1=.1,g=1+t)
    assert abs(crossover["ratio"]-1)<1e-12
    assert abs(1+t-1.8695652173913042)<1e-12
    # Late parasitism may increase head damage even when direct attack
    # completely fails to rescue seed production. Host feeding matters.
    late_more=evaluate(r=.7,q0=0,q1=0,p0=.6,p1=.1,g=2)
    assert abs(late_more["ratio"]-.48125)<1e-12
    # Conversely, an enemy that lowers per-host consumption produces the
    # conventional enemy-shielding penalty even without predamage killing.
    protects=evaluate(r=.9,q0=0,q1=0,p0=.9,p1=.1,g=.5)
    assert protects["barrier_increases_current_seed_damage"]
    # Baseline early removal of all herbivores cannot yield a relative ratio.
    zero=evaluate(r=.7,q0=1,q1=.2,p0=.6,p1=.1,g=2)
    assert zero["ratio"] is None and zero["baseline_zero_damage_no_ratio"]
    # Confirm repeated nested observations NEVER become a data estimate.
    n_grid=0
    for r in (0.,.5,.7,1.):
        for q0,q1 in ((0.,0.),(.6,.2),(.9,.9)):
            for p0,p1 in ((0.,0.),(.6,.1),(.1,.9)):
                for g in (.5,1.,2.,5.):
                    v=evaluate(r=r,q0=q0,q1=q1,p0=p0,p1=p1,g=g)
                    assert v["baseline_damage"]>=0 and v["barrier_damage"]>=0
                    n_grid+=1
    for bad in ({"r":-1}, {"p1":1.1}, {"q0":float("nan")}, {"g":0}):
        args={"r":.7,"q0":.6,"q1":.2,"p0":.6,"p1":.1,"g":2}
        args.update(bad)
        try:evaluate(**args)
        except ValueError:pass
        else:raise AssertionError("INVALID_DYNAMICAL_PARAMETER_ACCEPTED")
    return {
        "status":"PASS_SYNTHETIC_EARLY_KILL_VS_LATE_HOST_FEEDING",
        "grid_combinations":n_grid,
        "neutral_host_feeding_reversal_ratio":neutral["ratio"],
        "greater_host_feeding_cancels_reversal_ratio":feeding["ratio"],
        "analytical_host_feeding_multiplier_crossover":1+t,
        "late_feeding_only_damage_ratio":late_more["ratio"],
        "empirical_field_coefficients_fitted":0,
        "not_a_Cirsium_experiment":True,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    result={
        "version":"capitulum_parasitoid_host_feeding_sign_v1",
        "date":"2026-10-08",
        "test_result":synthetic_tests(),
        "literature":{
            "study":"Xi, Eisenhauer & Sun (2015), Journal of Animal Ecology",
            "doi":"10.1111/1365-2656.12361",
            "Cirsium_species_role":"C. setosum INCLUDED ONLY IN OBSERVATIONAL FIVE-SPECIES COMPARISON",
            "field_experiment_species":"Saussurea nigrescens",
            "field_experimental_head_damage_means":{
               "tephritid_only":8.04,"tephritid_plus_parasitoid":14.83
            },
            "do_not_substitute_for_g":"Means per capitulum are NOT known per-host causal g for Cirsium, and the exposure/selection scheme differs."
        },
        "outcome_gates":[
            "Identify parasitoid strategy/host per independent head; no insect order-only assignment.",
            "Link parasitoid timing and host feeding before/after irreversible achene damage.",
            "Measure establishment E, verified predamage killing q, continuing-host parasitism p and damage d separately.",
            "Measure genuine spine/phyllary/orientation access before claiming architectural gating.",
            "Verify filled/viable achenes, pollen limitation, host count and head stage.",
            "No adaptation/correlated evolution from this model without heritable variation and replicated selection.",
        ],
        "current_focal_Cirsium_fitness_effect":"NOT_IDENTIFIABLE_NO_REAL_GUILD_STAGE_OR_VIABLE_SEED_DATA",
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["test_result"]["status"],
                     "grid":result["test_result"]["grid_combinations"],
                     "field_coefficients":0},indent=2))
if __name__=="__main__":
    main()
