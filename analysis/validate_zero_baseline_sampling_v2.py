#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "data" / "planning" / "chapter3_zero_baseline_sampling_v2.json"
TARGETS = ROOT / "data" / "planning" / "chapter3_zero_baseline_population_targets_v2.csv"
DOC = ROOT / "docs" / "CHAPTER3_ZERO_BASELINE_SAMPLING_PLAN_V2.md"
README = ROOT / "README.md"


def load_contract() -> dict:
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def rows() -> list[dict[str, str]]:
    with TARGETS.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_contract() -> dict:
    d = load_contract()
    if d.get("contract_version") != "chapter3_zero_baseline_sampling_v2":
        raise AssertionError("zero-baseline contract version drift")
    if d.get("starting_material") != "ZERO_PHYSICAL_SAMPLES":
        raise AssertionError("zero physical sample premise drift")
    if d["species_scaffold"]["primary"] != "Moreyra et al. 2025 target-capture nuclear framework":
        raise AssertionError("species scaffold was reassigned to RAD")

    s0 = d["S0_trait_reconnaissance"]
    if (s0["target_concepts"], s0["target_populations_per_concept"], s0["target_flowering_individuals_per_population"]) != (38, 2, 5):
        raise AssertionError("S0 trait reconnaissance target drift")
    if s0["destructive_collection_required"] is not False:
        raise AssertionError("S0 was made destructively mandatory")

    s1 = d["S1_focal_launch"]
    if s1["default_population_bank"] != {"physical_individuals": 15, "primary_rad": 12, "reserves": 3}:
        raise AssertionError("S1 12+3 population bank drift")
    if s1["derived_minimum"] != {"physical_plants": 165, "initial_rad_individuals": 132}:
        raise AssertionError("S1 minimum total drift")
    if s1["derived_target_mainland_replication"]["physical_plants"] != 210 or s1["derived_target_mainland_replication"]["initial_rad_individuals"] != 168:
        raise AssertionError("S1 target total drift")

    pilot = d["S1_rad_pilot_nested"]
    if len(pilot["systems"]) != 5:
        raise AssertionError("focal RAD pilot must use five mandatory S1 systems")
    if pilot["stage_A"]["libraries"] != 15:
        raise AssertionError("S1 Stage A pilot must remain 15 libraries")
    if pilot["stage_B"]["primary_libraries"] != 20 or pilot["stage_B"]["technical_repeats"] != 5 or pilot["stage_B"]["libraries_total"] != 25:
        raise AssertionError("S1 Stage B pilot must remain 20+5=25")
    if pilot["technical_gates"]["minimum_overlap_genotype_concordance"] != 0.95 or pilot["technical_gates"]["minimum_core_locus_recovery"] != 0.90:
        raise AssertionError("focal pilot technical gate drift")
    if "does not qualify high-ploidy/polyploid" not in pilot["polyploid_boundary"]:
        raise AssertionError("focal pilot was silently generalized to polyploids")

    s2 = d["S2_transition_neighbourhood"]
    if not any(">=80%" in x for x in s2["nomination_algorithm"]):
        raise AssertionError("S2 uncertainty-cover rule drift")
    if s2["per_concept"]["primary_rad_per_concept"] != 12 or s2["per_concept"]["reserves_per_concept"] != 4:
        raise AssertionError("S2 per-concept localization panel drift")

    s3 = d["S3_m01_expansion"]
    if s3["added_physical"] != 60 or s3["added_initial_rad"] != 48:
        raise AssertionError("S3 added panel drift")
    if s3["final_m01_physical"] != 120 or s3["final_m01_initial_rad"] != 96:
        raise AssertionError("S3 final M01 bank drift")
    if "non-RAD confirmation assay" not in s3["selection_boundary"]:
        raise AssertionError("RAD-only selection claim reopened")

    if d["S4_optional_all_japan_rad"]["priority"] != "LAST_OPTIONAL":
        raise AssertionError("38-concept all-Japan RAD was moved back to first-line sampling")
    if any(d["current_state"].values()):
        raise AssertionError("zero-baseline future work was silently opened")
    return d


def validate_targets() -> list[dict[str, str]]:
    data = rows()
    minimum = [r for r in data if r["population_role"].startswith("MINIMUM") or r["population_role"].startswith("DISCOVERY")]
    if sum(int(r["physical_target"]) for r in minimum) != 165:
        raise AssertionError("minimum S1 physical target rows do not total 165")
    if sum(int(r["primary_rad"]) for r in minimum) != 132:
        raise AssertionError("minimum S1 RAD target rows do not total 132")
    conditional = [r for r in data if r["stage"] == "S3"]
    if len(conditional) != 4 or sum(int(r["physical_target"]) for r in conditional) != 60:
        raise AssertionError("S3 conditional expansion must remain four x15")
    if any(r["status"] != "CONDITIONAL_CLOSED" for r in conditional):
        raise AssertionError("S3 was opened before promotion")
    return data


def validate_narrative() -> None:
    doc = DOC.read_text(encoding="utf-8")
    for phrase in (
        "zero physical biological material",
        "165 physical plants, 132 initial RAD individuals",
        "15 shallow libraries",
        "25 libraries",
        ">=80% of cumulative transition-placement uncertainty",
        "S4 — optional all-Japan RAD sensitivity, last not first",
    ):
        if phrase not in doc:
            raise AssertionError(f"zero-baseline narrative missing: {phrase}")
    readme = README.read_text(encoding="utf-8")
    if "CHAPTER3_ZERO_BASELINE_SAMPLING_PLAN_V2.md" not in readme:
        raise AssertionError("README does not route to zero-baseline v2")
    if "38 concepts × 3 wild individuals = **114** all-Japan floor" in readme:
        raise AssertionError("README still presents 38x3 as the first sampling step")


def main() -> int:
    for p in (CONTRACT, TARGETS, DOC, README):
        if not p.exists() or p.stat().st_size == 0:
            raise AssertionError(f"missing zero-baseline file: {p.relative_to(ROOT)}")
    d = validate_contract()
    data = validate_targets()
    validate_narrative()
    print("chapter3_zero_baseline_sampling_v2_valid=true")
    print("starting_material=ZERO_PHYSICAL_SAMPLES")
    print("S1_minimum_physical=165")
    print("S1_minimum_initial_rad=132")
    print("focal_stage_A_libraries=15")
    print("focal_stage_B_libraries=25")
    print("S4_all_japan_rad=LAST_OPTIONAL")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
