#!/usr/bin/env python3
"""Why Cirsium head egg/visit totals cannot identify attraction versus defence.

Poisson thinning: arrivals ~ Poisson(T*lambda); reach conditional on
arrivals Bernoulli(p); egg-laying conditional on reach Bernoulli(s).
Then final egg events ~ Poisson(T*lambda*p*s). Without intermediate
observations there is a continuum of observationally equivalent
models with different biological mechanisms.

This is a standard probabilistic identity, NOT a discovered biological
law, fitted Cirsium result, selection proof, or evidence of exact
Poisson processes in natural arthropod communities.
"""
from __future__ import annotations

import json
import argparse
import math
from pathlib import Path


def rate(T, arrival_lambda, reach_p, oviposit_p):
    if not (T>0 and arrival_lambda>=0 and
            0<=reach_p<=1 and 0<=oviposit_p<=1 and
            all(math.isfinite(x) for x in
                (T,arrival_lambda,reach_p,oviposit_p))):
        raise ValueError("invalid effort/exposure/access conditional probability")
    return T*arrival_lambda*reach_p*oviposit_p


def poisson_pmf(k, mu):
    if not (isinstance(k,int) and k>=0 and mu>=0):
        raise ValueError("PMF input invalid")
    if mu==0:
        return 1.0 if k==0 else 0.0
    return math.exp(k*math.log(mu)-mu-math.lgamma(k+1))


def evidence_counterfactual(T, A, E, Y):
    """A is detected approaches, E verified head-zone reaches,
       Y confirmed eggs; cannot use if video not full effort."""
    if not(all(isinstance(k,int) for k in (A,E,Y)) and
            A>=E>=Y>=0 and T>0):
        raise ValueError("Invalid nested counts or denominator")
    return {
      "counted_approaches":A,
      "counted_zone_reaches":E,
      "counted_egg_events":Y,
      "approach_rate_per_minute":A/T,
      "reach_given_approach":E/A if A>0 else None,
      "egg_given_reach":Y/E if E>0 else None,
      "scientific_limit":"Conditional *descriptive* estimates, not selection; only with full audited approach-zone detection and independent guild IDs."
    }


def synthetic():
    T=1
    scenarioA={"arrival_lambda":4.0,"reach_p":0.25,"oviposit_p":0.5}
    scenarioB={"arrival_lambda":2.0,"reach_p":0.5,"oviposit_p":0.5}
    means=[rate(T,**z) for z in (scenarioA,scenarioB)]
    assert means==[.5,.5]
    for y in range(12):
        assert abs(poisson_pmf(y,means[0])-poisson_pmf(y,means[1]))<1e-14
    # The same 'observed egg count' can be caused by 2x arrivals or 2x entry;
    # observing head approaches discriminates these scenarios.
    assert scenarioA["arrival_lambda"]!=scenarioB["arrival_lambda"]
    accepted=evidence_counterfactual(T=2,A=8,E=2,Y=1)
    assert accepted["approach_rate_per_minute"]==4
    assert accepted["reach_given_approach"]==.25
    assert accepted["egg_given_reach"]==.5
    assert evidence_counterfactual(T=2,A=0,E=0,Y=0)["reach_given_approach"] is None
    invalid=[(2,5,6,2),(2,4,3,5),(0,4,3,1)]
    for T,A,E,Y in invalid:
        try:evidence_counterfactual(T,A,E,Y)
        except ValueError:pass
        else: raise AssertionError("allowed invalid count chain")
    return {
      "status":"PASS_POISSON_THINNING_NONIDENTIFIABILITY",
      "synthetic_equal_expected_egg_means":means,
      "same_final_egg_count_distribution_under_distinct_access_mechanisms":True,
      "different_expected_approach_counts":True,
      "positive_and_zero_verified_effort_examples":2,
      "invalid_nested_count_cases_rejected":3,
      "not_biological_data":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result={
        "status":"STRUCTURAL_IDENTIFIABILITY_PROOF_NOT_REAL_FIELD_RESULT",
        "published_identity":"Poisson thinning and binomial conditional decomposition (standard statistics)",
        "equations":{
           "approaches":"A ~ Poisson(T * lambda), if access-zone effort fully audited",
           "reaches_given_approaches":"R | A ~ Binomial(A, p_entry)",
           "eggs_given_reaches":"Y | R ~ Binomial(R, p_egg)",
           "marginal_eggs":"Y ~ Poisson(T * lambda * p_entry * p_egg)"
        },
        "mechanistic_ambiguity":[
          "Less oviposition can result from reduced arrival (floral cue or alternate-host availability), reduced access (head orientation × phyllary/spine), or reduced acceptance after access (host-race phenology/physiology).",
          "Even a true head-level phenotype association plus fewer eggs cannot isolate 'defence' if intermediate arrival and head-zone reach are missing.",
          "Neighboring other-Cirsium open buds/heads change the arrival process lambda independently of observed focal head defence traits.",
          "Adult seed-feeder versus parasitoid late-stage access are separate host-stage transitions; final parasitism is not a surrogate for early seed rescue.",
          "The factorization is observational unless morphologic perturbations are randomized and other covariates/observation errors are independently addressed."
        ],
        "observation_requirements":{
           "T":"Verified, coverage-complete and guild-calibrated video minutes; negative and zero bouts retained.",
           "lambda":"Observed pre-entry independently identified approaches per same head and stage.",
           "p_entry":"Approach-to-disc/involucre and true morphological penetration outcomes, counting blocked attempted entries.",
           "p_egg":"Oviposition confirmed conditional on physical head reach, including zero success.",
           "alternate_hosts":"Measured buds/open heads of alternative local Cirsium, with spatial radius and same-stage timing.",
           "fitness":"Separate confirmed pollen transfer, viable seeds, predamage parasitoid death, sex function and herbivore damage."
        },
        "synthetic_tests":synthetic()
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"receipt":result["synthetic_tests"]},indent=2))


if __name__=="__main__":main()
