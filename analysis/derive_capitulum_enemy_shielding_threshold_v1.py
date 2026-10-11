#!/usr/bin/env python3
"""An exact *conditional* boundary for thistle seed-feeder enemy shielding.

Let established damaging larvae E; q is effective mortality BEFORE any
irreversible seed loss (not final parasitoid prevalence); d is loss per
damaging larva after that mortality. Then D=E*(1-q)*d.

For equal d and E1/E0=r<1, damage can INCREASE under added defense only
if q0 > 1-r and q1 < 1 - (1-q0)/r. This is a deterministic identity
within the stated restricted model, not a fitted ecological estimate.

More general d1/d0 adds a multiplier; background pollination potential,
resource costs, density dependence, after-attack damage, compensatory
reproduction, parasitoid-induced oviposition avoidance and movement are
outside this identity. Do not interpret this as plant total fitness.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def damage(E, q, d=1.0):
    if not (math.isfinite(E) and E>=0 and math.isfinite(q) and
            0<=q<=1 and math.isfinite(d) and d>=0):
        raise ValueError("E>=0, q in [0,1], d>=0, all finite")
    return E*(1-q)*d


def exact(r, q0, q1, d_ratio=1.0):
    if not (0<r<=1 and 0<=q0<=1 and 0<=q1<=1 and d_ratio>0
            and all(map(math.isfinite,(r,q0,q1,d_ratio)))):
        raise ValueError("Requires 0<r<=1, 0<=q0,q1<=1, d_ratio>0")
    if q0==1:
        return {
          "damage0":0.0,
          "relative_damage_ratio":None,
          "absolute_new_damage_positive":q1<1,
          "interpretation":"Baseline damage is zero; ratio undefined; new damage positive if some damaging larvae remain.",
          "no_ratio_claim":True
        }
    ratio=r*((1-q1)/(1-q0))*d_ratio
    critical=1-(1-q0)/(r*d_ratio)
    return {
      "egg_establishment_ratio":r,
      "baseline_predamage_kill":q0,
      "barrier_predamage_kill":q1,
      "damage_per_surviving_larva_ratio":d_ratio,
      "relative_damage_ratio":ratio,
      "increased_damage":ratio>1+1e-12,
      "unchanged_damage":abs(ratio-1)<=1e-12,
      "critical_max_q1_exclusive":critical,
      "necessary_baseline_q0_for_any_reversal":1-r*d_ratio,
      "could_ever_reverse_for_some_q1_in_0_to_1":critical>0,
      "interpretation":"Conditional damage, NOT viable seed total or pollen/seed-fitness; stage-specific q is unobserved in original Cirsium GloBI."
    }


def synthetic_tests():
    # Added barrier 30% reduction in establishment, but later less successful
    # predamage parasitoid killing can increase expected damage.
    a=exact(.7,.6,.2)
    assert abs(a["relative_damage_ratio"]-1.4)<1e-12
    assert a["increased_damage"]
    assert abs(a["critical_max_q1_exclusive"]-3/7)<1e-12
    # Important empirical feasibility threshold: when pre-damage killing is
    # at most the fraction by which barriers reduce eggs, paradox impossible.
    b=exact(.7,.2,0)
    assert not b["increased_damage"]
    assert not b["could_ever_reverse_for_some_q1_in_0_to_1"]
    c=exact(.7,.3,0)
    assert c["unchanged_damage"]
    # A high percentage of parasitoid emergence after damage is NOT an
    # early q. This sanity check intentionally supplies ONLY the early q.
    late_parasitism_prevalence=.8
    predamage_effective_q0=.04
    d=exact(.7,predamage_effective_q0,0)
    assert d["relative_damage_ratio"] < 1
    # Invariance of dimensionless ratio to arbitrary absolute head count.
    for E in (1,7,13,100):
        assert abs(damage(.7*E,.2)/damage(E,.6)-1.4)<1e-12
    # Reversal of seed damage != reversal of total plant fitness if
    # pollination potential changes (not estimated by this model).
    example_seed_potential_0=80.
    example_seed_potential_1=120.
    d0=damage(10,.6,1)
    d1=damage(7,.2,1)
    assert d1>d0 and (example_seed_potential_1-d1)>(example_seed_potential_0-d0)
    # Another boundary: a cost paid for a barrier could reverse total
    # plant fitness without any herbivore or parasitoid response.
    assert (70-damage(7,.6)) < (80-damage(10,.6))
    try:
        exact(.7,1.2,.1)
    except ValueError:pass
    else: raise AssertionError("Invalid mortality admitted")
    return {
      "status":"PASS_CONDITIONAL_IDENTITY",
      "number_of_test_families":7,
      "hypothetical_damage_reversal_ratio":a["relative_damage_ratio"],
      "necessary_predamage_killing_rate_if_eggs_drop_30pct":"baseline q0 must exceed 0.30 (assuming equal d)",
      "late_parasitism_not_interchangeable_with_early_host_killing":late_parasitism_prevalence != predamage_effective_q0,
      "not_a_biological_fitness_estimate":True
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    results={
       "title":"Capitulum seed-damage sign reversal is conditional on *early* effective enemy killing",
       "date":"2026-10-08",
       "empirical_status":"THEORETICAL_TEST_ONLY_NOT_OBSERVED_TRAIT_ADAPTATION",
       "necessary_under_equal_per_larva_damage":
         "If barrier reduces established eggs from E0 to E1=r E0, q0 > 1-r is necessary for a later parasitoid-exclusion penalty to outweigh egg reduction.",
       "sufficient_full_condition":
         "r * (1-q1)/(1-q0) * (d1/d0) > 1, with 0<=q<1 representing only parasitoid killing BEFORE irreversible seed damage.",
       "illustrative_scenarios":[
         {"scenario":"strong_effective_early_parasitism",**exact(.7,.6,.2)},
         {"scenario":"early_parasitism_too_weak_for_reversal",**exact(.7,.2,0)},
         {"scenario":"boundary_equal",**exact(.7,.3,0)},
         {"scenario":"late_parasitism_is_biologically_irrelevant_for_seed_rescue",**exact(.7,.04,0)},
         {"scenario":"damage_per_larva_also_changes",**exact(.7,.6,.2,.8)},
       ],
       "critical_distinction":[
         "Egg-laying attempt, confirmed egg, established larva and final surviving seed-feeder adult are distinct denominators.",
         "Original published parasitoid incidence in a head can be high even if all attacks occur after seed damage: it then cannot justify high q.",
         "Actual mean viable achene W involves pollen/ovule potential, seed production and resource/defense costs in addition to herbivore damage D; D1>D0 does not imply W1<W0.",
         "Even if D1>D0 and W1<W0 are measured under intervention, community resource availability, alternate host flowering and actual head/sex stages must be controlled.",
         "E, q, d are not estimated by GloBI relationships, head image covariation, or 2023 leaf/flower morphometrics."
       ],
       "unobserved_inputs_required_for_focal_Cirsium":{
         "E":"Pre-existing eggs or established damaging larvae, measured before later enemy-gating manipulation",
         "q0_q1":"Verified mortality of individually linked larval hosts BEFORE irreversible achene damage in the same treatment contrast",
         "d0_d1":"Damage per surviving larva given head size, floret count and stage",
         "W":"Ultimate filled or viable achenes in authorized samples; pollen delivery/potential and male function separately"
       },
       "synthetic_test_receipt":synthetic_tests(),
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(results,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":results["synthetic_test_receipt"]["status"],
          "critical":results["necessary_under_equal_per_larva_damage"],
          "n_scenarios":len(results["illustrative_scenarios"])},indent=2))


if __name__=="__main__":
    main()
