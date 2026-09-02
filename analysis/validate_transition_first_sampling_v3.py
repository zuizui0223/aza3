#!/usr/bin/env python3
"""Validate the authoritative transition-first Chapter 3 sampling plan v3."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "data" / "planning" / "chapter3_transition_first_sampling_v3.json"
DOC = ROOT / "docs" / "CHAPTER3_TRANSITION_FIRST_SAMPLING_PLAN_V3.md"
BUILDER = ROOT / "analysis" / "build_nmns_transition_candidate_universe_v3.py"
README = ROOT / "README.md"


def load() -> dict:
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def validate_contract() -> dict:
    d = load()
    if d.get("contract_version") != "chapter3_transition_first_sampling_v3":
        raise AssertionError("v3 contract version drift")
    if d.get("authoritative_for_new_sampling") is not True:
        raise AssertionError("v3 is not authoritative for new sampling")

    u = d["candidate_universe"]
    if u["source_record_count_observed_2026_09_02"] != 161:
        raise AssertionError("NMNS source-record snapshot drift")
    if "not 161 guaranteed independent species" not in u["record_count_interpretation"]:
        raise AssertionError("authority records were promoted to species count")
    if u["moreyra_japan_subset"] != 38:
        raise AssertionError("Moreyra subset size drift")
    if u["moreyra_subset_role"] != "HYPOTHESIS_ORIGIN_SUBSET_NOT_SAMPLING_UNIVERSE":
        raise AssertionError("Japan38 was promoted back to the Chapter 3 sampling universe")

    overlap = d["overlap_policy"]
    if overlap["same_individual_multitrait_reuse"] != "REQUIRED_AND_EFFICIENT":
        raise AssertionError("efficient same-individual reuse was removed")
    if overlap["same_neighbourhood_as_independent_transition_replicate"] != "PROHIBITED":
        raise AssertionError("one neighbourhood can be double-counted as evolutionary replication")
    if overlap["p03_minimum_nonoverlapping_neighbourhoods"] != 4:
        raise AssertionError("P03 non-overlapping neighbourhood minimum drift")

    a1 = d["A1_fixed_focal_population_layer"]
    bank = a1["population_bank"]
    if (bank["primary_rad_per_population"], bank["reserve_per_population"], bank["physical_per_population"]) != (12, 3, 15):
        raise AssertionError("focal population bank drift")
    if a1["minimum_launch_total"] != {"populations": 11, "physical": 165, "primary_rad": 132}:
        raise AssertionError("minimum focal launch totals drift")
    if a1["target_total"] != {"populations": 14, "physical": 210, "primary_rad": 168}:
        raise AssertionError("target focal replication totals drift")
    target = {row["system"]: row["populations"] for row in a1["target_replication"]}
    expected = {
        "JPN36 Cirsium sieboldii": 4,
        "JPN06 Cirsium dipsacolepis": 3,
        "JPN15 Cirsium lineare": 3,
        "Cirsium brevicaule": 2,
        "Cirsium irumtiense": 2,
    }
    if target != expected:
        raise AssertionError("fixed focal population replication drift")

    pilot = d["A1_rad_assay_pilot"]
    if len(pilot["systems"]) != 5:
        raise AssertionError("focal RAD pilot must remain nested in five mandatory systems")
    if pilot["stage_A"]["libraries"] != 15 or pilot["stage_B"]["total_libraries"] != 25:
        raise AssertionError("focal RAD pilot library counts drift")
    if pilot["minimum_overlap_genotype_concordance"] != 0.95 or pilot["minimum_core_locus_recovery"] != 0.90:
        raise AssertionError("RAD technical thresholds drift")
    if "second stratum-specific assay" not in pilot["polyploid_rule"]:
        raise AssertionError("later polyploid entry no longer triggers a separate assay check")

    a2 = d["A2_backbone_augmentation"]
    if a2["initial_nonbackbone_taxon_cap"] != 12 or a2["placement_representatives_per_taxon"] != 2:
        raise AssertionError("initial target-capture placement batch drift")
    if a2["rad_before_stable_placement"] is not False:
        raise AssertionError("non-backbone taxa can enter RAD before nuclear placement")
    if not any("no close-relative claim based on geography alone" == x for x in a2["nomination_requirements"]):
        raise AssertionError("geography-only close-relative shortcut reopened")

    p03 = d["A3_transition_neighbourhoods"]["P03"]
    if p03["minimum_neighbourhoods"] != 4 or p03["must_be_nonoverlapping"] is not True:
        raise AssertionError("P03 replicate-neighbourhood rule drift")
    if "not posterior probability mass" not in p03["selection_basis"]:
        raise AssertionError("unweighted parsimony histories were incorrectly converted to probability mass")
    for needed in (
        "two populations of one species do not count as two evolutionary-transition slots",
        "JPN36 can count at most once",
    ):
        if needed not in p03["requirements"]:
            raise AssertionError(f"P03 anti-pseudoreplication rule missing: {needed}")

    added = d["added_neighbourhood_population_default"]
    if (added["populations_per_taxon"], added["primary_rad_per_population"], added["reserve_per_population"]) != (2, 10, 2):
        raise AssertionError("added-neighbourhood population default drift")
    if "morph is not identical to locality" not in added["within_species_polymorphism_rule"]:
        raise AssertionError("morph-locality confounding rule missing")

    if any(bool(v) for k, v in d["current_state"].items() if k != "physical_samples") or d["current_state"]["physical_samples"] != 0:
        raise AssertionError("v3 fail-closed state was opened")
    return d


def validate_narrative() -> None:
    doc = DOC.read_text(encoding="utf-8")
    required = [
        "Japan38 is not the sampling universe",
        "161 authority records",
        "four non-overlapping candidate transition neighbourhoods",
        "authority-wide trait screen → phylogenetic placement where needed → transition-neighbourhood selection → population RAD",
        "Geographic proximity alone cannot promote a taxon to `close relative`",
        "165",
        "210",
        "15 shallow libraries",
        "25 libraries",
        "up to 12 non-backbone authority concepts",
        "not posterior probabilities",
    ]
    missing = [x for x in required if x not in doc]
    if missing:
        raise AssertionError(f"v3 narrative missing: {missing}")

    builder = BUILDER.read_text(encoding="utf-8")
    for needed in ("DEFAULT_URL", "source_catchphrase_sha256", "record_count_is_species_count", "represented_in_moreyra_binomial_screen"):
        if needed not in builder:
            raise AssertionError(f"authority-universe builder missing: {needed}")

    readme = README.read_text(encoding="utf-8")
    if "CHAPTER3_TRANSITION_FIRST_SAMPLING_PLAN_V3.md" not in readme:
        raise AssertionError("README does not route to authoritative v3")
    if "Japan38 = hypothesis-origin subset" not in readme:
        raise AssertionError("README still treats Japan38 as sampling universe")


def main() -> int:
    for path in (CONTRACT, DOC, BUILDER, README):
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing v3 artifact: {path.relative_to(ROOT)}")
    d = validate_contract()
    validate_narrative()
    print("chapter3_transition_first_sampling_v3_valid=true")
    print(f"nmns_authority_records_snapshot={d['candidate_universe']['source_record_count_observed_2026_09_02']}")
    print(f"moreyra_hypothesis_origin_subset={d['candidate_universe']['moreyra_japan_subset']}")
    print(f"minimum_focal_primary_rad={d['A1_fixed_focal_population_layer']['minimum_launch_total']['primary_rad']}")
    print(f"target_focal_primary_rad={d['A1_fixed_focal_population_layer']['target_total']['primary_rad']}")
    print("new_sampling_authorized=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
