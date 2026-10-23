#!/usr/bin/env python3
"""Original-source same-BRANCH capitulum seed-production hurdle audit.

West & Louda (2021), AJB 10.1002/ajb2.1744; USDA CC0
10.15482/USDA.ADC/1522649. Head codes 200/210, 300/310 etc.
identify head position on a branch, NOT orientation angle, age or
phyllary/spine biomechanics. This study is descriptive and post hoc.
It does not identify adaptive branch architecture, pollination or
a causal pathway from insect feeding to mature seed production.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
from statistics import mean

from audit_cirsium_undulatum_within_branch_position_v1 import (
    BLOBS, YEARS, check, exact_sign_p, identity, position,
    quantile, raw_audit, read_locked, value, xorshift,
)

RNG_SEED = 20261023


def build_pairs(mothers, heads, dictionary):
    """Use ONLY a documented terminal and intermediate head on the same branch."""
    check('Apical head is ""100""' in dictionary
          and "between the terminals and stem" in dictionary,
          "SOURCE_HEAD_CODE_ONTOLOGY_NOT_VERIFIED")
    parents, clean, source_audit = raw_audit(mothers, heads)
    grouped = defaultdict(lambda: {"branch_terminal": [], "intermediate": []})
    for h in clean:
        kind = position(h)
        if kind == "apical":
            continue
        number = int(h["hd"])
        branch = number // 100
        check(number in (branch * 100, branch * 100 + 10),
              "HEAD_CODE_NOT_SINGLE_BRANCH_LOCATION")
        grouped[(h["year"], h["plt"], branch)][kind].append(h)

    matched = []
    unmatched_intermediate = []
    for (year, plant, branch), kinds in sorted(grouped.items()):
        t = kinds["branch_terminal"]
        i = kinds["intermediate"]
        check(len(t) <= 1 and len(i) <= 1,
              "DUPLICATED_LOCATION_ON_SINGLE_BRANCH")
        if i and not t:
            unmatched_intermediate.append(f"{year}|{plant}|{branch}")
        if not t or not i:
            continue
        terminal, inner = t[0], i[0]
        ts, iseed = value(terminal, "hd.gsd"), value(inner, "hd.gsd")
        td, injury = value(terminal, "av.dam"), value(inner, "av.dam")
        check(all(v is not None for v in (ts, iseed, td, injury)),
              "MISSING_PAIRED_OUTCOME_CANNOT_BE_ZERO_IMPUTED")
        check(all(v >= 0 and float(v).is_integer() for v in (ts, iseed)),
              "GOOD_SEED_COUNT_NOT_NONNEGATIVE_INTEGER")
        check(all(0 <= v <= 1 for v in (td, injury)),
              "DAMAGE_SCORE_OUTSIDE_ORIGINAL_SCALE")
        mother_id = f"{year}|{plant}"
        check(mother_id in parents, "UNKNOWN_ORIGINAL_PARENT")
        if ts > 0 and iseed > 0:
            route = "both"
        elif ts > 0:
            route = "terminal_only"
        elif iseed > 0:
            route = "intermediate_only"
        else:
            route = "neither"
        matched.append({
            "id": f"{mother_id}|{branch}",
            "mother_id": mother_id,
            "year": year,
            "arm": parents[mother_id]["trt"],
            "terminal_head_code": str(branch * 100),
            "intermediate_head_code": str(branch * 100 + 10),
            "terminal_good_seed": int(ts),
            "intermediate_good_seed": int(iseed),
            "terminal_damage": td,
            "intermediate_damage": injury,
            "good_seed_delta": iseed - ts,
            "damage_delta": injury - td,
            "route": route,
        })

    check(len(matched) == 20
          and Counter(x["year"] for x in matched) == {"2007": 11, "2008": 9},
          "ORIGINAL_20_SAME_BRANCH_MATCHES_CHANGED")
    check(unmatched_intermediate == ["2007|355|3"],
          "UNMATCHED_INTERMEDIATE_HEAD_CHANGED")
    check(len({x["id"] for x in matched}) == 20,
          "DUPLICATE_PAIRED_BRANCH_ID")
    return matched, source_audit, unmatched_intermediate


def sign_counts(values):
    return {
        "negative": sum(v < 0 for v in values),
        "positive": sum(v > 0 for v in values),
        "equal": sum(v == 0 for v in values),
    }


def yearly_results(pairs, n_bootstrap=5000):
    check(100 <= n_bootstrap <= 15000, "INVALID_RESAMPLING_BUDGET")
    # A branch is a measured contrast, NOT an independent mother.
    out = {}
    mother_deltas = {}
    for year in YEARS:
        sample = [b for b in pairs if b["year"] == year]
        by_mother = defaultdict(list)
        for b in sample:
            by_mother[b["mother_id"]].append(b["good_seed_delta"])
        m = {k: mean(v) for k, v in sorted(by_mother.items())}
        check(len(m) == 7, "EXACTLY_SEVEN_PAIRED_MOTHERS_PER_YEAR")
        mother_deltas[year] = list(m.values())
        seed = [b["good_seed_delta"] for b in sample]
        damage = [b["damage_delta"] for b in sample]
        both = [b["good_seed_delta"] for b in sample if b["route"] == "both"]
        sc = sign_counts(m.values())
        out[year] = {
            "n_paired_branches": len(sample),
            "n_paired_mothers": len(m),
            "paired_branch_seed_directions": sign_counts(seed),
            "paired_branch_damage_directions": sign_counts(damage),
            "paired_branch_seed_mean_delta": mean(seed),
            "paired_branch_damage_mean_delta": mean(damage),
            "branch_seed_routes": {
                cat: sum(b["route"] == cat for b in sample)
                for cat in ("neither", "terminal_only", "intermediate_only", "both")
            },
            "n_terminal_good_seed_positive": sum(
                b["terminal_good_seed"] > 0 for b in sample
            ),
            "n_intermediate_good_seed_positive": sum(
                b["intermediate_good_seed"] > 0 for b in sample
            ),
            "mean_seed_delta_conditioned_on_both_positive":
                mean(both) if both else None,
            "paired_mother_seed_deltas": m,
            "paired_mother_seed_directions": sc,
            "paired_mother_mean_seed_delta": mean(m.values()),
            "posthoc_mother_exact_two_sided_sign_p": exact_sign_p(
                sc["positive"], sc["negative"]
            ),
        }
    rng = xorshift(RNG_SEED)
    draws = {y: [] for y in YEARS}
    for _ in range(n_bootstrap):
        for y in YEARS:
            arr = mother_deltas[y]
            draws[y].append(mean(arr[next(rng) % len(arr)] for _i in arr))
    for y in YEARS:
        out[y]["paired_mother_bootstrap_ci95"] = [
            quantile(draws[y], 0.025), quantile(draws[y], 0.975)
        ]
    out["bootstrap"] = {
        "draws": n_bootstrap,
        "seed": RNG_SEED,
        "unit": "original mother with 1+ exact same-branch pair",
        "interpretation": "descriptive conditional, not causal confidence interval",
    }
    return out


def analyze(mothers, heads, dictionary, n_bootstrap=5000):
    pairs, original, unmatched = build_pairs(mothers, heads, dictionary)
    years = yearly_results(pairs, n_bootstrap)
    return {
        "version": "aza3_undulatum_nested_branch_seed_hurdle_v2",
        "source": {
            "paper_doi": "10.1002/ajb2.1744",
            "dataset_doi": "10.15482/USDA.ADC/1522649",
            "git_blob_sources": BLOBS,
            "head_number_semantics": "same 100s block is a documented branch location",
        },
        "source_audit": original,
        "n_exact_branch_pairs": len(pairs),
        "n_original_mothers_with_exact_branch_pairs": len({
            x["mother_id"] for x in pairs
        }),
        "unmatched_intermediate_branch_ids": unmatched,
        "paired_branches": pairs,
        "by_year": years,
        "limits": {
            "branch_not_independent_of_original_mother": True,
            "only_original_mothers_with_both_codes_eligible": True,
            "unrecorded_head_not_assumed_biological_zero": True,
            "posthoc_sign_test_not_confirmatory_and_original_six_test_BH_q_2008_0p09375": True,
            "head_code_not_flowering_age_or_gravity_orientation": True,
            "cannot_separate_pollen_resource_phenology_insect_causes": True,
            "no_genetic_lego_or_evolutionary_regain_test": True,
            "new_field_observations": 0,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mothers", required=True, type=Path)
    parser.add_argument("--heads", required=True, type=Path)
    parser.add_argument("--dictionary", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--bootstrap", type=int, default=5000)
    a = parser.parse_args()
    mothers = read_locked(
        a.mothers, "mother_plants.csv",
        ("plt", "year", "dam", "ins", "trt")
    )
    heads = read_locked(
        a.heads, "flower_heads.csv",
        ("plt", "year", "dam", "ins", "hd", "av.dam", "hd.gsd")
    )
    dictionary = read_locked(a.dictionary, "source_dictionary.txt")
    z = analyze(mothers, heads, dictionary, a.bootstrap)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(z, ensure_ascii=False, indent=2,
                                allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "n_exact_branch_pairs": z["n_exact_branch_pairs"],
        "n_original_mothers_with_exact_branch_pairs":
            z["n_original_mothers_with_exact_branch_pairs"],
        "by_year": {
            y: {k: v for k, v in z["by_year"][y].items()
                if k != "paired_mother_seed_deltas"}
            for y in YEARS
        },
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
