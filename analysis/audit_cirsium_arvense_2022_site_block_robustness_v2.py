#!/usr/bin/env python3
"""Site/block robustness of actual 2022 Cirsium arvense treatment outcomes.

Source: Karpate et al. 2025 PREPRINT 10.21203/rs.3.rs-6634172/v1;
author open archive 10.5281/zenodo.15364152, CC BY 4.0.
seeds.csv is pinned to Yogita25/Multitrophic_Interactions @
8a151c62ed069f775a16b9bc5fbdce209cbb8691, exact Git blob SHA checked
through the existing original-source reader.

Compare organza BAG pollinator exclusion vs open control on the original
site x block pairs. This is an additional robustness audit on observed seed
bookkeeping, NOT a causal decomposition of pollination, seed feeder arrival,
real plant spine/angle, heritable traits, or adaptive fitness.
"""
from __future__ import annotations

import argparse
import collections
import json
import math
from pathlib import Path

from reanalyse_cirsium_arvense_2022_real_head_seed_bookkeeping_v1 import (
    source_csv, parse_head_rows, plant_stats, SITES,
    GIT_BLOBS, SOURCE_COMMIT, SOURCE_REPO, DOI_DATA, DOI_PAPER,
)

SITES_ORDER = ("campus", "construction", "ditch", "rapseed")
ARMS = {"Control", "Pollinator Exclusion"}
METRICS = {
    "viable_fraction": "viable_fraction_all",
    "direct_damaged_fraction": "directly_recorded_damaged_fraction_all",
    "missing_signed_fraction": "missing_residual_signed_fraction_all",
    "inferred_eaten_signed_fraction": "inferred_eaten_signed_fraction_all",
    "emerged_insects_per_head": "emerged_insects_per_head_all",
}
RAW_COUNT_METRICS = {
    "viable_per_sampled_head": "Viable_Seeds",
    "damaged_per_sampled_head": "Damaged_Seeds",
    "florets_per_sampled_head": "Flowers",
}


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def mean(v):
    require(len(v) > 0, "EMPTY_SITE_OR_BLOCK_COMPARISON")
    return sum(v) / len(v)


def sign_counts(v):
    eps = 1e-12
    return {
        "positive": sum(x > eps for x in v),
        "negative": sum(x < -eps for x in v),
        "zero": sum(abs(x) <= eps for x in v),
    }


def per_group(rows):
    q = plant_stats(rows)
    for name, column in RAW_COUNT_METRICS.items():
        q[name] = sum(x[column] for x in rows) / len(rows)
    return q


def derive_pairs(rows, enforce_source_n=True):
    if enforce_source_n:
        require(len(rows) == 615, "ORIGINAL_SOURCE_ROW_COUNT_DRIFT")
    keep = [r for r in rows if r["site"] in SITES
            and r["treatment"] in ARMS]
    groups = collections.defaultdict(list)
    for r in keep:
        groups[(r["site"], r["block"], r["treatment"])].append(r)
    lookup = {k: per_group(v) for k, v in groups.items()}
    blocks = sorted({(s, b) for s, b, _ in groups})
    paired = []
    for site, block in blocks:
        control = lookup.get((site, block, "Control"))
        bag = lookup.get((site, block, "Pollinator Exclusion"))
        if control is None or bag is None:
            continue
        d = {"site": site, "block": block,
             "n_control_heads": control["n_heads"],
             "n_bag_heads": bag["n_heads"]}
        for name, field in METRICS.items():
            d[name] = bag[field] - control[field]
        for name in RAW_COUNT_METRICS:
            d[name] = bag[name] - control[name]
        require(math.isclose(
            d["inferred_eaten_signed_fraction"],
            d["direct_damaged_fraction"] + d["missing_signed_fraction"],
            rel_tol=0, abs_tol=1e-10), "SEED_RESIDUAL_COMPONENT_IDENTITY_BROKEN")
        paired.append(d)
    if enforce_source_n:
        require(len(paired) == 39, "PINNED_STUDY_39_MATCHED_BLOCKS_DRIFT")
        require(collections.Counter(x["site"] for x in paired) ==
                {"campus": 12, "construction": 10, "ditch": 9, "rapseed": 8},
                "PINNED_STUDY_SITE_BLOCK_COUNTS_DRIFT")
    return paired


def analyze(paired):
    require(len(paired) >= 2 and len(set(x["site"] for x in paired)) > 1,
            "AT_LEAST_TWO_BLOCKS_AND_SITES_REQUIRED")
    names = tuple(METRICS) + tuple(RAW_COUNT_METRICS)
    sites = sorted({x["site"] for x in paired})
    summary = {
        "n_matched_site_blocks": len(paired),
        "site_n": {site: sum(x["site"] == site for x in paired)
                   for site in sites},
        "site": {},
        "overall": {},
        "site_equal_weight": {},
        "leave_one_site_out": {},
        "leave_one_block_out_range": {},
        "response_definitions": {
            "viable_fraction": "sum author 'Viable_Seeds' / sum 'Flowers' per inferred plant",
            "direct_damaged_fraction": "sum visually Damaged_Seeds / sum Flowers per inferred plant",
            "missing_signed_fraction": "sum author signed accounting Missing_Seeds / sum Flowers",
            "inferred_eaten_signed_fraction": "direct_damaged_fraction + missing_signed_fraction (identity, NOT observed feeding)",
            "emerged_insects_per_head": "author Insects count / sampled heads, not verified seed-feeders",
            "viable_per_sampled_head": "sum author Viable_Seeds / number of sampled head rows",
            "damaged_per_sampled_head": "sum author Damaged_Seeds / number of sampled head rows",
            "florets_per_sampled_head": "sum author Flowers / number of sampled head rows",
        },
    }
    for site in sites:
        sub = [x for x in paired if x["site"] == site]
        summary["site"][site] = {
            name: {"mean_bag_minus_control": mean([x[name] for x in sub]),
                   "sign_counts": sign_counts([x[name] for x in sub])}
            for name in names
        }
    for name in names:
        summary["overall"][name] = {
            "mean_bag_minus_control": mean([x[name] for x in paired]),
            "sign_counts": sign_counts([x[name] for x in paired]),
        }
        summary["site_equal_weight"][name] = mean([
            summary["site"][s][name]["mean_bag_minus_control"] for s in sites
        ])
        summary["leave_one_site_out"][name] = {
            site: mean([x[name] for x in paired if x["site"] != site])
            for site in sites
        }
        loo = [mean([p[name] for j, p in enumerate(paired) if j != i])
               for i in range(len(paired))]
        summary["leave_one_block_out_range"][name] = [min(loo), max(loo)]
    return summary


def synthetic_tests():
    tiny = [
        {"site": "A", "block": "1", "viable_fraction": -0.2,
         "direct_damaged_fraction": -0.1, "missing_signed_fraction": 0.3,
         "inferred_eaten_signed_fraction": 0.2, "emerged_insects_per_head": -1,
         "viable_per_sampled_head": -2, "damaged_per_sampled_head": -1,
         "florets_per_sampled_head": 0},
        {"site": "B", "block": "2", "viable_fraction": -0.6,
         "direct_damaged_fraction": 0.2, "missing_signed_fraction": -0.1,
         "inferred_eaten_signed_fraction": 0.1, "emerged_insects_per_head": 0,
         "viable_per_sampled_head": -6, "damaged_per_sampled_head": 2,
         "florets_per_sampled_head": 4},
    ]
    s = analyze(tiny)
    assert math.isclose(s["overall"]["viable_fraction"]["mean_bag_minus_control"], -0.4)
    assert s["overall"]["viable_fraction"]["sign_counts"]["negative"] == 2
    assert math.isclose(s["site_equal_weight"]["direct_damaged_fraction"], .05)
    assert math.isclose(s["overall"]["inferred_eaten_signed_fraction"]["mean_bag_minus_control"], .15)
    assert s["leave_one_block_out_range"]["viable_fraction"] == [-.6, -.2]
    assert math.isclose(s["leave_one_site_out"]["direct_damaged_fraction"]["A"], .2)
    return "PASS_SYNTHETIC_SIGN_AND_BLOCK_SITE_HOLDOUT_TESTS"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    hdr, raw = source_csv(args.seeds, "seeds.csv")
    require({"Flowers", "Site_Code", "Block_Code", "Treatment_Code",
             "Missing_Seeds", "Viable_Seeds", "Damaged_Seeds"}.issubset(hdr),
            "PINNED_ORIGINAL_SEED_HEADER_DRIFT")
    rows = parse_head_rows(raw)
    out = {
        "version": "cirsium_arvense_2022_site_block_robustness_v2",
        "as_of": "2026-10-20",
        "source": {
            "repo": SOURCE_REPO, "commit": SOURCE_COMMIT,
            "seeds_git_blob_sha": GIT_BLOBS["seeds.csv"],
            "zenodo_doi": DOI_DATA, "paper_preprint_doi": DOI_PAPER,
            "study_year": 2022, "source_raw_head_rows": len(rows),
            "treatment_is_organza_BAG_not_pollination_only": True,
        },
        "synthetic_tests": synthetic_tests(),
        "actual_Cirsium_observed_site_block_sensitivity": analyze(
            derive_pairs(rows)),
        "biological_claim_ceiling":
            "A within-study descriptive robustness audit of two treatments on inferred plant/site-block keys. Signed Missing is an accounting residual; Eaten includes this residual. No isolated pollination path, verified seed-feeder larvae, controlled capitulum morphology, or causal guild-dependent adaptive fitness.",
        "own_Cirsium_field_plant_assays_added": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, ensure_ascii=False)
                        + "\n", encoding="utf-8")
    print(json.dumps({
        "source_raw_heads": len(rows), "matched_blocks": out[
          "actual_Cirsium_observed_site_block_sensitivity"]["n_matched_site_blocks"],
        "site_n": out["actual_Cirsium_observed_site_block_sensitivity"]["site_n"],
        "overall": out["actual_Cirsium_observed_site_block_sensitivity"]["overall"],
        "synthetic_tests": out["synthetic_tests"],
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
