#!/usr/bin/env python3
"""Export the coordinate-free ACSP range-sector provenance receipt from aza3.

The exporter is intentionally fail-closed. The public readiness registry may
exist while field/locality gates are incomplete; in that state --check-only
reports the blockers and normal export refuses to create a provenance receipt.

No exact coordinates, sensitive access instructions, permit documents, or
personal details are read or written by this adapter.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = ROOT / "data/planning/chapter3_acsp_range_sector_readiness_v1.csv"
SCHEMA_VERSION = "cirsium-fresh-sentinel-range-sector-provenance-v1"
STATUS = "PRE_GEOMETRY_RANGE_SECTOR_PROVENANCE_FROZEN"
CONTRACT_VERSION = "chapter3_exact_site_freeze_v8"
EXPECTED_UNITS = ("CIR02", "CIR06", "CIR12", "CIR13")
EXPECTED = {
    "CIR02": {
        "species_binomial": "Cirsium inundatum",
        "aza3_slot_id": "P1_Cirsium_inundatum_A",
        "range_sector_label": "青森県側の確認済み野生集団",
        "p02_required": False,
    },
    "CIR06": {
        "species_binomial": "Cirsium yezoalpinum",
        "aza3_slot_id": "P1_Cirsium_yezoalpinum_A",
        "range_sector_label": "知床山系の確認済み野生集団",
        "p02_required": False,
    },
    "CIR12": {
        "species_binomial": "Cirsium dipsacolepis",
        "aza3_slot_id": "P6_Cirsium_dipsacolepis_OWN",
        "range_sector_label": "P02で最初に検証される野生集団セクター",
        "p02_required": True,
    },
    "CIR13": {
        "species_binomial": "Cirsium lineare",
        "aza3_slot_id": "P6_Cirsium_lineare_OWN",
        "range_sector_label": "P02で最初に検証される野生集団セクター",
        "p02_required": True,
    },
}
REQUIRED_COLUMNS = (
    "cohort_unit_id",
    "species_binomial",
    "aza3_slot_id",
    "range_sector_label",
    "freeze_status",
    "current_occurrence_status",
    "conservation_review_status",
    "land_manager_gate_satisfied",
    "collection_permission_status",
    "tissue_collection_permission_status",
    "field_window_checked",
    "target_locality_id",
    "private_exact_site_record_exists",
    "p02_first_validated_wild_population_required",
    "p02_first_validated_wild_population_link_status",
    "range_sector_geometry_may_now_be_materialized",
)
TRUE = {"true", "1", "yes"}
FALSE = {"false", "0", "no"}
PERMISSION_READY = {"APPROVED", "NOT_REQUIRED"}
CONSERVATION_COMPLETE = {"CLEAR", "RESTRICTION_PRESENT", "NOT_APPLICABLE"}


def _bool(value: object, label: str) -> bool:
    text = str(value or "").strip().lower()
    if text in TRUE:
        return True
    if text in FALSE:
        return False
    raise ValueError(f"{label} must be an explicit boolean")


def load_registry(path: Path = DEFAULT_REGISTRY) -> dict[str, dict[str, str]]:
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"missing ACSP range-sector readiness registry: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != REQUIRED_COLUMNS:
            raise ValueError("ACSP range-sector readiness registry columns changed")
        rows = list(reader)
    if len(rows) != len(EXPECTED_UNITS):
        raise ValueError("ACSP range-sector readiness registry must contain exactly four rows")
    by_unit: dict[str, dict[str, str]] = {}
    for row in rows:
        unit = str(row.get("cohort_unit_id") or "").strip()
        if unit not in EXPECTED:
            raise ValueError(f"unexpected cohort_unit_id: {unit!r}")
        if unit in by_unit:
            raise ValueError(f"duplicate cohort_unit_id: {unit}")
        expected = EXPECTED[unit]
        for key in ("species_binomial", "aza3_slot_id", "range_sector_label"):
            if str(row.get(key) or "") != expected[key]:
                raise ValueError(f"{unit} {key} drifted from the ACSP handoff identity")
        p02_required = _bool(
            row.get("p02_first_validated_wild_population_required"),
            f"{unit}.p02_first_validated_wild_population_required",
        )
        if p02_required is not bool(expected["p02_required"]):
            raise ValueError(f"{unit} P02 dependency identity changed")
        by_unit[unit] = row
    if tuple(by_unit) != EXPECTED_UNITS:
        raise ValueError("registry row order must be CIR02,CIR06,CIR12,CIR13")
    return by_unit


def assess_registry(rows: dict[str, dict[str, str]]) -> dict[str, Any]:
    blockers: dict[str, list[str]] = {}
    for unit in EXPECTED_UNITS:
        row = rows[unit]
        reasons: list[str] = []
        if row["freeze_status"] != "FROZEN_FOR_FIELD_COLLECTION":
            reasons.append("freeze_status")
        if row["current_occurrence_status"] != "SUPPORTED_CURRENT":
            reasons.append("current_occurrence_status")
        if row["conservation_review_status"] not in CONSERVATION_COMPLETE:
            reasons.append("conservation_review_status")
        if not _bool(row["land_manager_gate_satisfied"], f"{unit}.land_manager_gate_satisfied"):
            reasons.append("land_manager_gate_satisfied")
        if row["collection_permission_status"] not in PERMISSION_READY:
            reasons.append("collection_permission_status")
        if row["tissue_collection_permission_status"] not in PERMISSION_READY:
            reasons.append("tissue_collection_permission_status")
        if not _bool(row["field_window_checked"], f"{unit}.field_window_checked"):
            reasons.append("field_window_checked")
        if not str(row["target_locality_id"]).strip():
            reasons.append("target_locality_id")
        if not _bool(row["private_exact_site_record_exists"], f"{unit}.private_exact_site_record_exists"):
            reasons.append("private_exact_site_record_exists")
        if not _bool(
            row["range_sector_geometry_may_now_be_materialized"],
            f"{unit}.range_sector_geometry_may_now_be_materialized",
        ):
            reasons.append("range_sector_geometry_may_now_be_materialized")

        p02_required = bool(EXPECTED[unit]["p02_required"])
        link_status = str(row["p02_first_validated_wild_population_link_status"]).strip()
        if p02_required:
            if link_status != "SATISFIED":
                reasons.append("p02_first_validated_wild_population_link_status")
        elif link_status != "NOT_APPLICABLE":
            reasons.append("p02_first_validated_wild_population_link_status")
        if reasons:
            blockers[unit] = reasons

    return {
        "schema_version": "chapter3-acsp-range-sector-readiness-assessment-v1",
        "status": "READY_FOR_ACSP_PROVENANCE_EXPORT" if not blockers else "BLOCKED_UPSTREAM_EXACT_SITE_FREEZE",
        "cohort_unit_ids": list(EXPECTED_UNITS),
        "ready_unit_count": len(EXPECTED_UNITS) - len(blockers),
        "blocked_unit_count": len(blockers),
        "blockers_by_unit": blockers,
        "p02_first_validated_wild_population_dependency_units": ["CIR12", "CIR13"],
        "exact_coordinates_read": False,
        "sensitive_access_instructions_read": False,
    }


def _git_snapshot_commit(registry: Path, *, repo_root: Path = ROOT) -> str:
    repo = Path(repo_root).resolve()
    path = Path(registry).resolve()
    try:
        relative = path.relative_to(repo).as_posix()
    except ValueError as exc:
        raise ValueError("registry must be inside the aza3 repository") from exc
    subprocess.run(["git", "ls-files", "--error-unmatch", "--", relative], cwd=repo, check=True, capture_output=True)
    status = subprocess.run(
        ["git", "status", "--porcelain", "--", relative],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if status:
        raise ValueError("registry has uncommitted changes; provenance must bind a committed snapshot")
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def build_provenance(
    rows: dict[str, dict[str, str]],
    *,
    upstream_snapshot_commit: str,
) -> dict[str, Any]:
    assessment = assess_registry(rows)
    if assessment["blocked_unit_count"]:
        raise ValueError("range-sector readiness registry is not fully frozen for ACSP export")
    if len(upstream_snapshot_commit) != 40 or any(ch not in "0123456789abcdef" for ch in upstream_snapshot_commit):
        raise ValueError("upstream snapshot commit must be a lowercase 40-hex SHA")

    units: dict[str, Any] = {}
    for unit in EXPECTED_UNITS:
        row = rows[unit]
        p02_required = bool(EXPECTED[unit]["p02_required"])
        units[unit] = {
            "species_binomial": row["species_binomial"],
            "aza3_slot_id": row["aza3_slot_id"],
            "range_sector_label": row["range_sector_label"],
            "freeze_status": "FROZEN_FOR_FIELD_COLLECTION",
            "current_occurrence_supported": True,
            "permission_gate_satisfied": True,
            "target_locality_id": str(row["target_locality_id"]).strip(),
            "private_exact_site_record_exists": True,
            "range_sector_geometry_may_now_be_materialized": True,
            "p02_first_validated_wild_population_required": p02_required,
            "p02_first_validated_wild_population_link_satisfied": p02_required,
        }
    return {
        "schema_version": SCHEMA_VERSION,
        "status": STATUS,
        "cohort_unit_ids": list(EXPECTED_UNITS),
        "aza3_exact_site_contract_version": CONTRACT_VERSION,
        "upstream_snapshot_commit": upstream_snapshot_commit,
        "unit_provenance": units,
        "exact_coordinates_included": False,
        "sensitive_access_instructions_included": False,
        "prospective_acsp_field_outcomes_opened": False,
        "acsp_field_outcomes_used_to_define_sector": False,
        "post_geometry_edits_allowed": False,
        "post_outcome_edits_allowed": False,
        "public_safe_to_commit": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--out-json", type=Path)
    args = parser.parse_args()

    rows = load_registry(args.registry)
    assessment = assess_registry(rows)
    if args.check_only:
        print(json.dumps(assessment, ensure_ascii=False, indent=2))
        return 0

    if args.out_json is None:
        raise SystemExit("--out-json is required unless --check-only is used")
    if assessment["blocked_unit_count"]:
        print(json.dumps(assessment, ensure_ascii=False, indent=2))
        return 2

    commit = _git_snapshot_commit(args.registry)
    value = build_provenance(rows, upstream_snapshot_commit=commit)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    if args.out_json.exists():
        raise SystemExit(f"refusing to overwrite existing provenance receipt: {args.out_json}")
    args.out_json.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": value["status"], "out_json": str(args.out_json)}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
