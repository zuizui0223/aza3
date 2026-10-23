#!/usr/bin/env python3
"""Real C. undulatum source: nested branch sample selection and year-method confounding.

West & Louda 2021, AJB 10.1002/ajb2.1744; USDA original data DOI
10.15482/USDA.ADC/1522649. The 2012 West dissertation, chapter 2
(University of Nebraska bioscidiss/41), independently details changed
insecticide identities/schedules and apical damage methods across years.
The 2007 and 2008 year labels therefore do NOT identify environment-only
effects. Only 14 of 106 original maternal plants have both branch-code types.
No causal selection effect, adaptive head value, visitor mechanism or
genetic modularity is identified by this descriptive reanalysis.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from math import comb
import json
from pathlib import Path
from statistics import mean, median

from audit_cirsium_undulatum_nested_branch_hurdle_v2 import (
    BLOBS, YEARS, build_pairs, check, identity, read_locked, value,
)

ARMS = ("ndi", "ndni", "di", "dni")
ROUTES = ("neither", "terminal_only", "intermediate_only", "both")

# Human-verified EXTERNAL documentary provenance (not computable from USDA CSV).
# West 2012 ch. 2, pages printed 61-64; PDF zero-based pages 68-72 and
# figure 2.2 legend PDF page 87. These are not newly inferred measurements.
METHODS_BY_YEAR = {
    "2007": {
        "insecticide_active_ingredient": "acephate",
        "application_interval_days": "7-10",
        "application_label_solution_percent": 1.0,
        "apical_damage": "insect_cage_plus_drill_with_girdling_if_needed",
        "flowerhead_census_interval_weeks": 2,
    },
    "2008": {
        "insecticide_active_ingredient": "bifenthrin",
        "application_interval_days": "14",
        "application_label_solution_percent": 0.06,
        "apical_damage": "insect_assisted_plus_mechanical_only_replicates",
        "flowerhead_census_interval_weeks": 4,
    },
}

def summarize(rows):
    vals = list(rows)
    return {
        "n_observed": len(vals),
        "mean": mean(vals) if vals else None,
        "median": median(vals) if vals else None,
        "min": min(vals) if vals else None,
        "max": max(vals) if vals else None,
    }


def exact_discordant_symmetry(terminal_only, intermediate_only):
    """Two-sided binomial sign reference; meaningful only for independent pairs."""
    n = terminal_only + intermediate_only
    if n == 0:
        return 1.0
    return min(1.0, 2 * sum(comb(n, i) for i in
                         range(min(terminal_only, intermediate_only) + 1))
               / (2 ** n))


def select_one_branch_per_mother(pairs, year, rule):
    """Fixed smallest/largest ORIGINAL branch code; not temporal head age."""
    check(rule in ("smallest_code", "largest_code"), "INVALID_CODE_CHOICE_RULE")
    per_mother = defaultdict(list)
    for p in pairs:
        if p["year"] == year:
            per_mother[p["mother_id"]].append(p)
    selected = []
    for key, rows in sorted(per_mother.items()):
        ordered = sorted(rows, key=lambda x: int(x["terminal_head_code"]))
        selected.append(ordered[0 if rule == "smallest_code" else -1])
    check(len(selected) == 7, "SEVEN_MOTHERS_EXPECTED_PER_YEAR")
    routes = Counter(z["route"] for z in selected)
    discordant = {
        "terminal_only": routes["terminal_only"],
        "intermediate_only": routes["intermediate_only"],
    }
    return {
        "rule": rule,
        "one_original_branch_code_pair_per_mother": True,
        "no_head_code_to_flowering_age_inference": True,
        "n_distinct_mothers": len(selected),
        "selected_mother_branch_ids": [x["id"] for x in selected],
        "route_counts": {r: routes[r] for r in ROUTES},
        "discordant": discordant,
        "two_sided_exact_discordant_symmetry_reference":
            exact_discordant_symmetry(routes["terminal_only"],
                                      routes["intermediate_only"]),
        "not_randomized_branch_placement": True,
        "not_confirmatory_or_unselected_population_inference": True,
    }


def mother_stage(m):
    seeds = value(m, "plt.gsd")
    check(seeds is not None and seeds >= 0 and seeds.is_integer(),
          "MOTHER_GOOD_SEED_COUNT_NOT_VERIFIED")
    d = value(m, "r.diam")
    check(d is None or d >= 0, "INVALID_ROSETTE_DIAMETER")
    return seeds, d


def group_mothers(mothers, matched_ids):
    byyear = {}
    for year in YEARS:
        plants = [m for m in mothers if m["year"] == year]
        matched = [m for m in plants if identity(m) in matched_ids]
        remainder = [m for m in plants if identity(m) not in matched_ids]
        def grp(subset):
            dat = [mother_stage(x) for x in subset]
            return {
                "n_original_mothers": len(subset),
                "author_good_seed_total": summarize(s for s, _ in dat),
                "source_rosette_diameter_nonmissing": summarize(
                    d for _, d in dat if d is not None),
                "n_good_seed_positive": sum(s > 0 for s, _ in dat),
            }
        arms = {}
        for arm in ARMS:
            sample = [m for m in plants if m["trt"] == arm]
            selected = [m for m in sample if identity(m) in matched_ids]
            not_selected = [m for m in sample if identity(m) not in matched_ids]
            check(len(sample) == len(selected) + len(not_selected),
                  "ARM_SELECTION_DENOMINATOR_CHANGED")
            arms[arm] = {
                "n_original_mothers": len(sample),
                "n_with_matching_branch_codes": len(selected),
                "share_selected": len(selected) / len(sample),
                "selected": grp(selected),
                "not_selected": grp(not_selected),
            }
        check(sum(x["n_original_mothers"] for x in arms.values())
              == len(plants), "ARM_COUNTS_DO_NOT_SUM_TO_ORIGINAL_MOTHERS")
        selected_stats = grp(matched)
        excluded_stats = grp(remainder)
        a = selected_stats["author_good_seed_total"]["mean"]
        b = excluded_stats["author_good_seed_total"]["mean"]
        byyear[year] = {
            "all_original_mothers": len(plants),
            "mothers_with_both_branch_code_positions": len(matched),
            "mothers_without_both_branch_code_positions": len(remainder),
            "fraction_with_both": len(matched) / len(plants),
            "selected": selected_stats,
            "not_selected": excluded_stats,
            "selected_divided_by_unselected_mean_good_seed": (
                a / b if b else None),
            "by_year_treatment_arm": arms,
        }
    check(byyear["2007"]["all_original_mothers"] == 39
          and byyear["2008"]["all_original_mothers"] == 67,
          "EXACT_YEAR_MOTHER_DENOMINATORS_CHANGED")
    check(all(byyear[y]["mothers_with_both_branch_code_positions"] == 7
              for y in YEARS), "FOURTEEN_PAIRED_SOURCE_MOTHERS_CHANGED")
    return byyear


def analyze(mothers, heads, dictionary):
    check(len(mothers) == 106 and len(heads) == 323,
          "SOURCE_DENOMINATORS_UNEXPECTED")
    check(len({identity(m) for m in mothers}) == 106,
          "DUPLICATE_ORIGINAL_MOTHER_ID")
    pairs, source_audit, unpaired = build_pairs(mothers, heads, dictionary)
    paired_ids = {p["mother_id"] for p in pairs}
    check(len(paired_ids) == 14 and len(pairs) == 20,
          "EXACT_PAIR_COUNTS_CHANGED")
    by_year = group_mothers(mothers, paired_ids)
    route_year = {}
    for y in YEARS:
        yy = [p for p in pairs if p["year"] == y]
        route_year[y] = {
            "all_actual_recorded_branch_code_pairs": len(yy),
            "n_distinct_original_mothers": len({
                p["mother_id"] for p in yy}),
            "all_pair_route_counts_descriptive_not_independent":
                {r: sum(p["route"] == r for p in yy) for r in ROUTES},
            "lowest_code_one_per_mother": select_one_branch_per_mother(
                pairs, y, "smallest_code"),
            "highest_code_one_per_mother": select_one_branch_per_mother(
                pairs, y, "largest_code"),
        }
    return {
        "version": "aza3_undulatum_year_protocol_selection_gate_v3",
        "original_sources": {
            "dataset_doi": "10.15482/USDA.ADC/1522649",
            "2021_paper_doi": "10.1002/ajb2.1744",
            "source_git_blob_shas": BLOBS,
            "dissertation_2012_chapter2":
                "https://digitalcommons.unl.edu/bioscidiss/41/",
            "independently_reviewed_pdf_pages_zero_based": [68, 69, 70, 71, 87],
            "figure_2_2_subordinate_heads_below_terminal_on_same_branch":
                True,
        },
        "manual_external_source_method_provenance": METHODS_BY_YEAR,
        "original_source_audit": source_audit,
        "n_unmatched_intermediate_branch_codes": len(unpaired),
        "original_denominator_selection": by_year,
        "nested_branch_binary_route_sensitivity": route_year,
        "limits": {
            "year_environment_protocol_is_perfectly_conflated": True,
            "insecticide_chemicals_not_equivalent_due_to_dose_label": True,
            "rosette_timing_in_current_export_not_ascertained": True,
            "conditioning_on_observed_both_head_types_is_post_outcome_selection": True,
            "same_plant_head_pairs_are_correlated": True,
            "smallest_code_is_not_flowering_order": True,
            "2012_method_metadata_is_externally_manually_verified_not_CSV_encoded": True,
            "dissertation_figure_code_correspondence_not_per_head_photo_join": True,
            "separate_mechanisms_pollen_injury_resources_not_identified": True,
            "observed_good_seed_not_germination_or_lifetime_fitness": True,
            "no_supported_evolutionary_adaptation_genetic_lego_or_trait_regain": True,
            "new_field_observations": 0,
        },
    }


def main():
    q = argparse.ArgumentParser()
    q.add_argument("--mothers", required=True, type=Path)
    q.add_argument("--heads", required=True, type=Path)
    q.add_argument("--dictionary", required=True, type=Path)
    q.add_argument("--out", required=True, type=Path)
    args = q.parse_args()
    mothers = read_locked(args.mothers, "mother_plants.csv",
                          ("year", "plt", "trt", "plt.gsd", "r.diam"))
    heads = read_locked(args.heads, "flower_heads.csv",
                        ("year", "plt", "dam", "ins", "hd", "hd.gsd", "av.dam"))
    dictionary = read_locked(args.dictionary, "source_dictionary.txt")
    z = analyze(mothers, heads, dictionary)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(z, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    print(json.dumps({
        "original_denominator_selection": {
            y: {k: v for k, v in z["original_denominator_selection"][y].items()
                if k != "by_year_treatment_arm"} for y in YEARS},
        "nested_route_sensitivity": z["nested_branch_binary_route_sensitivity"],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
