from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COHORT = ROOT / "data" / "planning" / "chapter3_acsp_prospective_validation_cohort_v1.csv"
PREFILL = ROOT / "data" / "planning" / "chapter3_acsp_field_handoff_prefill_v1.csv"
CONTRACT = ROOT / "data" / "planning" / "chapter3_acsp_prospective_validation_contract_v1.json"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _anchor_class(count: int) -> str:
    if count == 0:
        return "ZERO_PRIMARY_ANCHOR"
    if count == 1:
        return "SINGLE_PRIMARY_ANCHOR"
    return "MULTIPLE_PRIMARY_ANCHORS"


def validate() -> dict[str, object]:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    cohort = _read_csv(COHORT)
    prefill = _read_csv(PREFILL)

    assert contract["status"] == "FROZEN_PRE_PATCH_PRE_FIELD"
    assert contract["field_outcomes_opened"] is False
    assert contract["candidate_patches_built"] is False
    assert contract["public_exact_coordinates_written"] is False
    assert contract["cohort_size"] == 13

    assert len(cohort) == 13
    assert len(prefill) == 13
    assert len({row["cohort_unit_id"] for row in cohort}) == 13
    assert len({row["species_binomial"] for row in cohort}) == 13
    assert len({row["aza3_slot_id"] for row in cohort}) == 13

    forbidden = ("latitude", "longitude", "decimaldegrees", "exact_coordinate")
    for field in cohort[0]:
        lowered = field.lower()
        assert not any(token in lowered for token in forbidden)
    for field in prefill[0]:
        lowered = field.lower()
        assert not any(token in lowered for token in ("latitude", "longitude", "decimaldegrees"))

    regimes = Counter(row["occurrence_problem_class"] for row in cohort)
    anchor_classes = Counter(row["anchor_replication_class"] for row in cohort)
    method_arms = Counter(row["method_arm"] for row in cohort)
    families = Counter(row["structural_feature_family"] for row in cohort)

    assert regimes == Counter(contract["regime_counts"])
    assert anchor_classes == Counter(contract["anchor_replication_counts"])
    assert method_arms == Counter(contract["method_arm_counts"])
    assert families == Counter(
        {
            "ALPINE_TOPOGRAPHIC_STRUCTURE": 4,
            "OPEN_GRASSLAND_STRUCTURE": 3,
            "WETLAND_MOISTURE_STRUCTURE": 2,
            "COASTAL_ISLAND_STRUCTURE": 2,
            "FOREST_EDGE_STRUCTURE": 1,
            "GENERAL_SPATIAL_BASELINE_ONLY": 1,
        }
    )

    cohort_by_unit: dict[str, dict[str, str]] = {}
    for row in cohort:
        unit = row["cohort_unit_id"]
        cohort_by_unit[unit] = row
        anchor_count = int(row["primary_unique_coordinate_count"])
        assert row["anchor_replication_class"] == _anchor_class(anchor_count)
        assert row["outcome_opened"] == "false"
        assert row["candidate_patch_status"] == "NOT_BUILT"
        assert row["field_performance_denominator"] == "true"

        expected_slot_id = f"{row['aza3_priority']}_{row['species_binomial'].replace(' ', '_')}_{row['selected_slot']}"
        assert row["aza3_slot_id"] == expected_slot_id

        comparators = set(row["comparators"].split("|"))
        arm = row["method_arm"]
        regime = row["occurrence_problem_class"]
        family = row["structural_feature_family"]
        if arm == "STRUCTURAL_LOCAL":
            assert regime == "LOCAL_CONTINUATION"
            assert anchor_count >= 1
            assert {"ANNULAR_NEAREST_KNOWN", "DETERMINISTIC_SPATIAL_BALANCE"} <= comparators
            assert family != "GENERAL_SPATIAL_BASELINE_ONLY"
        elif arm == "STRUCTURAL_SENTINEL":
            assert regime == "SENTINEL"
            assert anchor_count == 0
            assert {"VALIDATED_BROAD_ROBUST_SUPPORT", "DETERMINISTIC_SPATIAL_BALANCE"} <= comparators
            assert family != "GENERAL_SPATIAL_BASELINE_ONLY"
        elif arm == "SPATIAL_BASELINE_ONLY":
            assert regime == "LOCAL_CONTINUATION"
            assert family == "GENERAL_SPATIAL_BASELINE_ONLY"
            assert {"ANNULAR_NEAREST_KNOWN", "DETERMINISTIC_SPATIAL_BALANCE"} <= comparators
        else:
            raise AssertionError(f"unknown method arm: {arm}")

    assert {row["acsp_validation_unit_id"] for row in prefill} == set(cohort_by_unit)
    for row in prefill:
        source = cohort_by_unit[row["acsp_validation_unit_id"]]
        assert row["aza3_slot_id"] == source["aza3_slot_id"]
        assert row["priority"] == source["aza3_priority"]
        assert row["species_binomial"] == source["species_binomial"]
        assert row["sample_slot"] == source["selected_slot"]
        assert row["range_sector"] == source["range_sector"]
        assert row["discovery_regime"] == source["occurrence_problem_class"]
        assert row["structural_feature_family"] == source["structural_feature_family"]
        assert row["anchor_evidence_class"] == source["anchor_replication_class"]
        assert int(row["anchor_count"]) == int(source["primary_unique_coordinate_count"])
        assert row["prospective_method_arm"] == source["method_arm"]
        assert row["comparator_assignment"] == source["comparators"]

        assert row["acsp_patch_id"] == ""
        assert row["acsp_abstention_or_block_reason"] == ""
        assert row["current_occurrence_field_state"] == "NOT_OPENED"
        assert row["identity_verification_status"] == "NOT_EVALUATED"
        assert row["search_minutes"] == ""
        assert row["observer_count"] == ""
        assert row["traversed_length_m"] == ""
        assert row["searched_area_m2"] == ""
        assert row["tissue_acquired_secondary"] == "false"
        assert row["deidentified_locality_id"] == ""
        assert row["private_exact_site_record_status"] == "NOT_CREATED"
        assert row["aza3_site_freeze_status"] == "NOT_FROZEN"
        assert row["collection_permission_status"] == "NOT_AUTHORIZED"
        assert row["tissue_permission_status"] == "NOT_AUTHORIZED"

    assert contract["freeze_rules"]["outcome_dependent_taxon_replacement_allowed"] is False
    assert contract["freeze_rules"]["outcome_dependent_feature_family_switch_allowed"] is False
    assert contract["freeze_rules"]["outcome_dependent_method_arm_switch_allowed"] is False
    assert contract["freeze_rules"]["access_failure_is_biological_nondetection"] is False
    assert contract["freeze_rules"]["permission_block_is_biological_nondetection"] is False

    return {
        "status": "OK",
        "cohort_size": len(cohort),
        "regime_counts": dict(regimes),
        "anchor_replication_counts": dict(anchor_classes),
        "method_arm_counts": dict(method_arms),
        "structural_family_counts": dict(families),
        "candidate_patches_built": False,
        "field_outcomes_opened": False,
    }


if __name__ == "__main__":
    print(json.dumps(validate(), ensure_ascii=False, indent=2))
