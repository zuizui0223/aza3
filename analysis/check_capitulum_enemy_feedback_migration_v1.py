#!/usr/bin/env python3
"""Diffusion/migration robustness of the Cirsium delayed enemy-shielding thought experiment.

SYNTHETIC, NOT FITTED. Two equally sized habitat patches with baseline (0)
or more armoured (1) flower heads. q0/q1 = effective parasitoid killing
AFTER irreversible current achene damage, lowering next-generation recruitment
only. r = reduced initial seed-feeder establishment on armoured heads.

A0[t+1] = I + (1-m) k0 A0[t] + m k1 A1[t]
A1[t+1] = I + m k0 A0[t] + (1-m) k1 A1[t]
k0=K(1-q0), k1=K*r(1-q1), 0<=m<=0.5, I>0

For 0<=k0,k1<1 (sufficient for stable equilibria):
D1*/D0* = r*[1-(1-2m)k0]/[1-(1-2m)k1].
A local patch-level reversal requires
    (1-2m) r K (q0-q1) > 1-r.
This is NOT a claim about plant genotype invasion fitness or selection.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

def _prob(x, name, lower=0.0, upper=1.0):
    if not isinstance(x, (float, int)) or not math.isfinite(x) or not lower <= x <= upper:
        raise ValueError(f"{name} outside [{lower}, {upper}]")

def evaluate(*, r, q0, q1, K, mixing, I=1.0, d=1.0, generations=12):
    for key, value, lo, hi in (
        ("r", r, 0, 1), ("q0", q0, 0, 1), ("q1", q1, 0, 1),
        ("mixing", mixing, 0, .5),
    ):
        _prob(value, key, lo, hi)
    if not isinstance(K, (int, float)) or not math.isfinite(K) or K<0:
        raise ValueError("K must be finite nonnegative")
    if not all(isinstance(x, (int,float)) and math.isfinite(x) and x>0 for x in (I,d)):
        raise ValueError("I and d must be finite positive")
    if not isinstance(generations, int) or not 1<=generations<=200:
        raise ValueError("generations 1..200")
    k0=K*(1-q0)
    k1=K*r*(1-q1)
    # Clearly delimit the sufficient domain rather than silently extrapolating
    # an unstable linear feedback or unknown density dependence.
    if k0>=1 or k1>=1:
        return {
            "status":"HOLD_OUTSIDE_STABLE_COMPONENT_DOMAIN",
            "k0":k0,"k1":k1,
            "no_equilibrium_fitness_inference":True,
        }
    u=1-2*mixing
    ratio_closed=r*(1-u*k0)/(1-u*k1)
    # Independently invert I-B, then check the two equilibrium equations.
    a=1-(1-mixing)*k0
    b=-(mixing*k1)
    c=-(mixing*k0)
    e=1-(1-mixing)*k1
    det=a*e-b*c
    if det<=0:
        raise AssertionError("Unexpected determinant in stable component domain")
    A0star=I*(e-b)/det
    A1star=I*(a-c)/det
    assert abs(A1star*r/A0star-ratio_closed)<1e-10
    assert abs(A0star-(I+(1-mixing)*k0*A0star+mixing*k1*A1star))<1e-10
    assert abs(A1star-(I+mixing*k0*A0star+(1-mixing)*k1*A1star))<1e-10
    if K>0 and r>0 and q0>q1:
        m_threshold=(1-(1-r)/(r*K*(q0-q1)))/2
        feasible_threshold=(max(0.,min(.5,m_threshold))
                            if m_threshold>0 else None)
    else:
        feasible_threshold=None
    reversal=ratio_closed>1+1e-12
    inequality=u*r*K*(q0-q1)>(1-r)+1e-12
    assert inequality==reversal
    # Starting from uniform baseline before one patch obtains the barrier.
    A0=A1=I/(1-k0)
    series=[]
    for t in range(generations+1):
        series.append({
            "generation":t, "A_control":A0, "A_armoured":A1,
            "damage_control":A0*d, "damage_armoured":r*A1*d,
            "relative_damage_armoured_to_control":r*A1/A0,
        })
        A0,A1=(I+(1-mixing)*k0*A0+mixing*k1*A1,
               I+mixing*k0*A0+(1-mixing)*k1*A1)
    assert abs(series[0]["relative_damage_armoured_to_control"]-r)<1e-12
    return {
        "status":"SYNTHETIC_STABLE_TWO_PATCH_MODEL",
        "migration_fraction_per_offspring_generation":mixing,
        "k0":k0,"k1":k1,
        "relative_patch_damage_at_equilibrium":ratio_closed,
        "equilibrium_patch_reversal":reversal,
        "critical_mixing_if_positive":feasible_threshold,
        "control_arrivals_equilibrium":A0star,
        "armoured_arrivals_equilibrium":A1star,
        "first_reversal_within_simulated_generations":next(
            (x["generation"] for x in series if
             x["relative_damage_armoured_to_control"]>1+1e-12),None),
        "trajectory":series,
        "empirical_observations":0,
        "not_identified":["genetic selection","phylogenetic evolution",
                          "herbivore dispersal kernel","actual parasite kill timing",
                          "viable achenes","pollination","host density dependence"],
    }

def synthetic_tests():
    args={"r":.7,"q0":.6,"q1":.1,"K":1.2}
    zero=evaluate(**args,mixing=0)
    mixed=evaluate(**args,mixing=.5)
    assert abs(zero["relative_patch_damage_at_equilibrium"]-1.4918032786885247)<1e-12
    assert abs(mixed["relative_patch_damage_at_equilibrium"]-.7)<1e-12
    assert zero["equilibrium_patch_reversal"]
    assert not mixed["equilibrium_patch_reversal"]
    threshold=zero["critical_mixing_if_positive"]
    assert abs(threshold-1/7)<1e-12
    boundary=evaluate(**args,mixing=threshold)
    assert abs(boundary["relative_patch_damage_at_equilibrium"]-1)<1e-12
    below=evaluate(**args,mixing=.1)
    above=evaluate(**args,mixing=.15)
    assert below["equilibrium_patch_reversal"] and not above["equilibrium_patch_reversal"]
    assert all(abs(evaluate(**args,mixing=m)["relative_patch_damage_at_equilibrium"]-.7)<1e-12
               for m in (.5,))
    for x in (
        {"r":.7,"q0":.6,"q1":.6,"K":1.2,"mixing":0},
        {"r":.7,"q0":.6,"q1":.1,"K":0,"mixing":0},
    ):
        result=evaluate(**x)
        assert not result["equilibrium_patch_reversal"]
        assert result["critical_mixing_if_positive"] is None
    for q in (-.1,1.01):
        try:
            evaluate(r=.7,q0=q,q1=.1,K=1.2,mixing=.2)
        except ValueError:
            pass
        else:
            raise AssertionError("Bad mortality admitted")
    for m in (-.1,.51):
        try:
            evaluate(**args,mixing=m)
        except ValueError:
            pass
        else:
            raise AssertionError("Bad mixing admitted")
    unstable=evaluate(r=1,q0=0,q1=0,K=1.2,mixing=.1)
    assert unstable["status"].startswith("HOLD_")
    # Random-looking deterministic grid: direct and matrix solutions agree
    # inside the stable component domain, even when migration is intermediate.
    grid=0
    for r in (.2,.5,.7,1.):
        for K in (0.,.5,1.2):
            for q0,q1 in ((.6,.1),(.3,.3),(.1,.6)):
                for m in (0.,.1,.25,.5):
                    out=evaluate(r=r,q0=q0,q1=q1,K=K,mixing=m)
                    if out["status"]=="SYNTHETIC_STABLE_TWO_PATCH_MODEL":
                        grid+=1
    assert grid>=80
    return {
        "status":"PASS_SYNTHETIC_MIGRATION_DESTROYS_PATCH_REVERSAL",
        "grid_stable_parameter_sets_checked":grid,
        "threshold_exchange_fraction":threshold,
        "unmixed_equilibrium_damage_ratio":zero["relative_patch_damage_at_equilibrium"],
        "mixing_10pct_ratio":below["relative_patch_damage_at_equilibrium"],
        "mixing_15pct_ratio":above["relative_patch_damage_at_equilibrium"],
        "fully_mixed_ratio":mixed["relative_patch_damage_at_equilibrium"],
        "empirical_estimates":0,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    receipt={
        "name":"capitulum_enemy_escape_migration_robustness_v1",
        "date":"2026-10-08",
        "version":"minimal_two_equally_sized_patches",
        "synthetic_tests":synthetic_tests(),
        "illustrative_mixing_sensitivity":[
            {
              "m":m,
              "ratio":evaluate(r=.7,q0=.6,q1=.1,K=1.2,mixing=m)["relative_patch_damage_at_equilibrium"],
            }
            for m in (0.,.05,.1,.14,.15,.2,.5)
        ],
        "literature_context":{
            "known":"Enemy-free-space and plant-mediated parasitoid accessibility are well-studied, e.g. Peterson et al. 2016 DOI 10.3389/fpls.2016.01794",
            "Cirsium_primary":"Vanbergen et al. 2006 DOI 10.1111/j.1365-2656.2006.01099.x",
            "unobserved":"No study cited here directly measures thistle-phyllary-induced late parasitoid suppression together with returning herbivore dispersal and genotype-specific fitness",
        },
        "ecological_boundary":"patch damage ratio != trait invasion fitness; model assumes two habitat types with equal resources and immigration.",
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":receipt["synthetic_tests"]["status"],
                      "threshold":receipt["synthetic_tests"]["threshold_exchange_fraction"],
                      "real_data":0},indent=2))

if __name__=="__main__":
    main()
