#!/usr/bin/env python3
"""Prospective phenotype-to-function equivalence gate for Cirsium.

No user field, RAD, pollen, or seed rows are ingested. All scenario values printed
by this script are explicit hypothetical demonstrations, never empirical effects.
"""
from __future__ import annotations
import json
import math
from statistics import NormalDist

ENDPOINTS = ("conspecific_pollen", "direct_seed_damage", "filled_achenes")
REQUIRED_GATES = (
    "identified_independent_mothers",
    "joint_shape_function_outcome_ids",
    "matched_flowering_stage_and_sex",
    "comparable_site_and_realized_guild",
    "known_video_effort",
    "predeclared_equivalence_margins",
)

def precision_halfwidth(n_mothers: int, maternal_difference_sd: float,
                        confidence: float = 0.90) -> float:
    """Normal-approx CI half-width; NOT power, does not account for site clustering."""
    if (not isinstance(n_mothers, int) or n_mothers < 2 or
        not math.isfinite(maternal_difference_sd) or maternal_difference_sd <= 0 or
        not 0 < confidence < 1):
        raise ValueError("Invalid prospective independent-mother precision parameters")
    z = NormalDist().inv_cdf((1 + confidence) / 2)
    return z * maternal_difference_sd / math.sqrt(n_mothers)

def classify_evidence(*, gates: dict[str, bool], distinct_morphologies: bool,
                      intervals: dict[str, tuple[float, float]],
                      margins: dict[str, float]) -> dict:
    """Classify already valid mother-cluster intervals, not infer those intervals.

    Each interval must be a 90% CI for the predeclared group contrast to use
    the usual alpha=0.05 TOST decision. Intervals are NOT computed from raw heads.
    Positive/negative mechanistic responses cannot be added as if comparable.
    """
    missing = [key for key in REQUIRED_GATES if gates.get(key) is not True]
    if missing:
        return {"status": "NOT_IDENTIFIABLE", "missing_gates": missing}
    if not distinct_morphologies:
        return {"status": "NO_DISTINCT_MORPHOLOGY_CONTRAST"}
    if set(intervals) != set(ENDPOINTS) or set(margins) != set(ENDPOINTS):
        return {"status": "NOT_IDENTIFIABLE", "missing_gates": ["three_joint_endpoints"]}
    equivalent = []
    clearly_different = []
    for outcome in ENDPOINTS:
        pair = intervals[outcome]
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            raise ValueError("Malformed endpoint interval")
        low, high = pair
        margin = margins[outcome]
        if not all(isinstance(z, (int,float)) and math.isfinite(z)
                   for z in (low, high, margin)):
            raise ValueError("Nonfinite interval or equivalence margin")
        if low > high or margin <= 0:
            raise ValueError("Reversed interval or nonpositive equivalence margin")
        equivalent.append(low > -margin and high < margin)
        clearly_different.append(low > margin or high < -margin)

    if all(equivalent):
        status = "CONDITIONAL_MULTI_ENDPOINT_EQUIVALENCE_COMPATIBLE"
    elif equivalent[2] and any(clearly_different[:2]):
        status = "SAME_NET_FITNESS_COMPONENT_DIFFERENT_PATHWAYS"
    elif clearly_different[2]:
        status = "MATERnal_REPRODUCTIVE_OUTCOME_DIFFERENT"
    else:
        status = "INDETERMINATE_NOT_EQUIVALENT_BY_DEFAULT"

    return {
        "status": status,
        "equivalent_endpoint_names": [k for k, ok in zip(ENDPOINTS, equivalent) if ok],
        "clearly_different_endpoint_names": [k for k, ok in zip(ENDPOINTS, clearly_different) if ok],
        "inference_limit": ("A morphological natural contrast is observational; "
            "identical maternal seed output does not prove same mechanism, "
            "trait selection, heritability, or evolutionary convergence."),
    }

def demo() -> dict:
    ready = dict.fromkeys(REQUIRED_GATES, True)
    margins = dict.fromkeys(ENDPOINTS, 0.1)  # HYPOTHETICAL common standardized scale
    case_equiv = classify_evidence(
        gates=ready, distinct_morphologies=True, margins=margins,
        intervals={"conspecific_pollen": (-.05,.05),
                   "direct_seed_damage": (-.04,.06),
                   "filled_achenes": (-.06,.04)})
    case_routes = classify_evidence(
        gates=ready, distinct_morphologies=True, margins=margins,
        intervals={"conspecific_pollen": (.20,.32),
                   "direct_seed_damage": (.18,.30),
                   "filled_achenes": (-.04,.04)})
    case_broad = classify_evidence(
        gates=ready, distinct_morphologies=True, margins=margins,
        intervals={name: (-.30,.30) for name in ENDPOINTS})
    blocked = classify_evidence(
        gates=dict.fromkeys(REQUIRED_GATES, False), distinct_morphologies=True,
        margins=margins, intervals={name: (-.04,.04) for name in ENDPOINTS})
    return {
        "status": "SYNTHETIC_PROSPECTIVE_DESIGN_ONLY",
        "not_empirical": True,
        "scenario_assumptions": {
            "precision_scale": "mother-level paired difference, hypothetical standardized units",
            "ci_type": "normal approximate 90 percent interval",
            "hypothetical_equivalence_margin": 0.1,
            "mothers": [24,60,120],
            "maternal_sd": [0.4,0.7],
            "not_power_or_sample_size_recommendation": True,
            "excluded_variance": "year, site, guild, unequal group counts, phenology, exposure and head clustering",
        },
        "ci_halfwidths": [
            {"n_independent_mothers": n, "hypothetical_maternal_sd": sd,
             "ci90_halfwidth": round(precision_halfwidth(n, sd), 6),
             "zero_effect_equivalence_possible_if_margin_point_one":
                 precision_halfwidth(n, sd) < 0.1}
            for sd in (0.4,0.7) for n in (24,60,120)
        ],
        "synthetic_evidence_examples": {
            "three_components_narrow": case_equiv,
            "same_final_but_different_routes": case_routes,
            "all_broad": case_broad,
            "no_joint_head_data": blocked,
        },
        "empirical_fitness_equivalence_authorized": False,
        "RAD_ancestral_reconstruction_completed": False,
    }

if __name__ == "__main__":
    print(json.dumps(demo(), indent=2, ensure_ascii=False))
