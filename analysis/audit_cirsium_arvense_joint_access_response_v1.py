#!/usr/bin/env python3
"""Original Cirsium arvense four-arm *joint* seed/insect response, 2022.

Source-authenticated Karpate et al. 2025 PREPRINT DOI 10.21203/rs.3.rs-6634172/v1,
original author's immutable seeds.csv from 2022 field experiment
Yogita25/Multitrophic_Interactions commit 8a151c62ed069f775a16b9bc5fbdce209cbb8691.

Source grain: 615 head rows => 612 eligible => 164 inferred plant proxies
(site x block x arm) => 45 site blocks, 36 with all FOUR intervention arms.
Within-block *two-outcome* response is determined at inferred-plant grain,
not assumed insect-to-seed lineage, and not an identified ecological mediation.

The organza bag excludes legitimate pollinators AND can exclude antagonists,
while predator exclusion is physical stem barrier and manual removals.
The plastic predator model is no reintroduction of real predators.
Results are descriptive of an existing source study; not new same-head
capitulum orientation/spine, genetic LEGO or adaptive landscape fitness.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path

from audit_cirsium_2022_four_arm_head_outcomes_v1 import (
    SOURCE_COMMIT, SOURCE_GIT_BLOB, SOURCE_REPO, ARCHIVE_DOI, PAPER_DOI,
    SITES, ARMS, parse_rows, prep, quantile, mean, check,
)

METRICS = (
    "viable_fraction", "emerged_insects_per_head",
    "observed_damaged_fraction", "viable_per_head",
)
COMPARISONS = (
    ("bag_minus_control", "Pollinator Exclusion", "Control"),
    ("pred_exclusion_minus_control", "Predator Exclusion", "Control"),
    ("pred_model_minus_control", "Predator Model", "Control"),
    ("pred_model_minus_pred_exclusion", "Predator Model", "Predator Exclusion"),
    ("bag_minus_pred_exclusion", "Pollinator Exclusion", "Predator Exclusion"),
)
SEED = 20261023
DRAWS = 10000
MASK = 0xFFFFFFFF


def xorshift32(seed=SEED):
    check(isinstance(seed, int) and 0 < seed < 2**32,
          "INVALID_XORSHIFT_SEED")
    state = seed
    while True:
        state ^= (state << 13) & MASK
        state ^= state >> 17
        state ^= (state << 5) & MASK
        state &= MASK
        yield state


def sign_counts(arr):
    return {
        "positive": sum(z > 1e-12 for z in arr),
        "negative": sum(z < -1e-12 for z in arr),
        "zero": sum(abs(z) <= 1e-12 for z in arr),
    }


def direct_response(block, comp, metric):
    return block["arms"][comp[1]][metric] - block["arms"][comp[2]][metric]


def joint_counts(complete, comp):
    c = Counter()
    for b in complete:
        v = direct_response(b, comp, "viable_fraction")
        i = direct_response(b, comp, "emerged_insects_per_head")
        va = "viable_up" if v > 1e-12 else (
            "viable_down" if v < -1e-12 else "viable_equal")
        ia = "emerged_up" if i > 1e-12 else (
            "emerged_down" if i < -1e-12 else "emerged_equal")
        c[va + "__" + ia] += 1
    check(sum(c.values()) == len(complete), "LOST_MULTIVARIATE_BLOCK")
    return dict(sorted(c.items()))


def summarize(complete):
    check(len(complete) >= 2 and len(set(b["site"] for b in complete)) >= 2,
          "REQUIRE_MULTIPLE_ORIGINAL_SITES_AND_BLOCKS")
    ordered = sorted(complete, key=lambda x: (x["site"], x["block"]))
    by_site = {s: [z for z in ordered if z["site"] == s] for s in SITES}
    check(all(by_site[s] for s in SITES),
          "SOURCE_ALL_FOUR_SITES_REQUIRED")
    results = {}
    for comp in COMPARISONS:
        measures = {}
        for metric in METRICS:
            actual = [direct_response(z, comp, metric) for z in ordered]
            site_means = {
                s: mean([direct_response(z, comp, metric) for z in by_site[s]])
                for s in SITES
            }
            measures[metric] = {
                "mean": mean(actual),
                "sign": sign_counts(actual),
                "site_means": site_means,
                "site_equal_weight_mean": mean(list(site_means.values())),
                "leave_one_site_out": {
                    s: mean([direct_response(z, comp, metric)
                             for z in ordered if z["site"] != s])
                    for s in SITES
                },
            }
        results[comp[0]] = {
            "arm_a": comp[1],
            "arm_b": comp[2],
            "n_original_matched_blocks": len(ordered),
            "measures": measures,
            "joint_response_counts": joint_counts(ordered, comp),
        }
    return results, by_site


def bootstrap(results, by_site, draws=DRAWS, seed=SEED):
    check(isinstance(draws, int) and 100 <= draws <= 20000,
          "INVALID_BOOTSTRAP_DRAWS")
    rng = xorshift32(seed)
    distributions = {
        c[0]: {m: [] for m in METRICS} for c in COMPARISONS
    }
    for _ in range(draws):
        sampled = []
        for site in SITES:
            block_list = by_site[site]
            for _ in range(len(block_list)):
                sampled.append(block_list[next(rng) % len(block_list)])
        for comp in COMPARISONS:
            for metric in METRICS:
                distributions[comp[0]][metric].append(
                    mean([direct_response(b, comp, metric) for b in sampled])
                )
    for comp in COMPARISONS:
        for metric in METRICS:
            arr = distributions[comp[0]][metric]
            results[comp[0]]["measures"][metric][
                "site_stratified_bootstrap_ci95"
            ] = [quantile(arr, .025), quantile(arr, .975)]
    return results


def synthetic_controls():
    def arm(viable, emerged, damage=.1, viables=20):
        return {
            "viable_fraction": viable,
            "emerged_insects_per_head": emerged,
            "observed_damaged_fraction": damage,
            "viable_per_head": viables,
        }
    complete = []
    for i, site in enumerate(SITES):
        complete.append({
            "site": site, "block": str(i), "arms": {
                "Control": arm(.8, 3),
                "Pollinator Exclusion": arm(.2, 0),
                "Predator Exclusion": arm(.8, 4),
                "Predator Model": arm(.8, 2),
            }
        })
    got = joint_counts(
        complete, ("bag_minus_control", "Pollinator Exclusion", "Control"))
    check(got == {"viable_down__emerged_down": 4},
          "INVALID_JOINT_QUADRANT_COUNT")
    got = joint_counts(
        complete, ("pred_exclusion_minus_control",
                   "Predator Exclusion", "Control"))
    check(got == {"viable_equal__emerged_up": 4},
          "INVALID_PREDATOR_COMPARATOR_TEST")
    check(sign_counts([-1.0, .0, .3]) == {
        "positive": 1, "negative": 1, "zero": 1
    }, "SIGN_ACCOUNTING_TEST")
    g = xorshift32(SEED)
    check(next(g) == next(xorshift32(SEED)), "RNG_REPLAY_BROKEN")
    return "PASS_MUTUALLY_EXCLUSIVE_JOINT_SIGNS_AND_PREDATOR_NEGATIVE_CONTROL"


def analyze(raw, draws=DRAWS):
    rows, source_sha = parse_rows(raw, strict=True)
    eligible, plants, blocks, complete = prep(rows, strict=True)
    ordered = sorted(complete, key=lambda x: (x["site"], x["block"]))
    actual_site_n = {s: sum(z["site"] == s for z in ordered) for s in SITES}
    check(actual_site_n == {
        "campus": 12, "construction": 9, "ditch": 8, "rapseed": 7
    }, "SOURCE_SITE_BLOCK_DISTRIBUTION_CHANGED")
    results, by_site = summarize(ordered)
    results = bootstrap(results, by_site, draws, SEED)
    check(results["bag_minus_control"]["measures"]["viable_fraction"]
          ["sign"] == {"positive": 0, "negative": 36, "zero": 0},
          "UNEXPECTED_BAGGED_SEED_RESPONSE")
    check(results["bag_minus_control"]["joint_response_counts"] ==
          {"viable_down__emerged_down": 25,
           "viable_down__emerged_equal": 8,
           "viable_down__emerged_up": 3},
          "OBSERVED_BAGGED_JOINT_RESPONSE_CHANGED")
    check(results["pred_exclusion_minus_control"]["measures"]
          ["viable_fraction"]["sign"] ==
          {"positive": 17, "negative": 19, "zero": 0},
          "PREDATOR_EXCLUSION_SEED_RESPONSE_CHANGED")
    return {
        "version": "cirsium_arvense_four_arm_joint_response_python_v1",
        "analyzed_date": "2026-10-23",
        "source": {
            "study_year": 2022,
            "preprint_year": 2025,
            "doi": PAPER_DOI,
            "dataset_doi": ARCHIVE_DOI,
            "author_repo": SOURCE_REPO,
            "author_commit": SOURCE_COMMIT,
            "seeds_git_blob_sha": source_sha,
        },
        "source_counts": {
            "heads_all": len(rows),
            "eligible_head_rows": len(eligible),
            "ineligible_rows": len(rows) - len(eligible),
            "plant_proxy_groups": len(plants),
            "source_site_blocks": len(blocks),
            "complete_four_arm_blocks": len(complete),
            "site_complete_counts": actual_site_n,
        },
        "bootstrap": {
            "seed": SEED,
            "draws": draws,
            "sites_fixed": True,
            "site_stratified_original_blocks": True,
            "unit": ("site x source block x arm (plant proxy); resampled "
                     "complete site-blocks within each of the four original sites"),
            "rng": "xorshift32",
            "site_order": list(SITES),
            "sorting": "lexicographic site then block",
            "not_causal_confidence_interval": True,
        },
        "comparisons": results,
        "synthetic_controls": synthetic_controls(),
        "limits": {
            "author_original_plant_id_not_verified": True,
            "pollinator_bag_compound_access_and_microclimate_intervention": True,
            "predator_exclusion_stem_adhesive_plus_manual_removal": True,
            "plastic_spider_not_actual_predator_addition": True,
            "insect_emergence_not_verified_flower_benefit_or_seed_feeder_only": True,
            "no_true_cirsium_phyllary_spine_or_head_orientation_measurement": True,
            "not_demonstrated_partner_specific_adaptive_head_LEGO": True,
            "new_biological_field_observations": 0,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=DRAWS)
    args = parser.parse_args()
    result = analyze(args.seeds.read_bytes(), args.bootstrap)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "source_counts": result["source_counts"],
        "bag_minus_control": {
            m: result["comparisons"]["bag_minus_control"]["measures"][m]
            for m in ("viable_fraction", "emerged_insects_per_head")
        },
        "pred_exclusion_minus_control": {
            m: result["comparisons"]["pred_exclusion_minus_control"]["measures"][m]
            for m in ("viable_fraction", "emerged_insects_per_head")
        },
    }, indent=2))


if __name__ == "__main__":
    main()
