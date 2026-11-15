#!/usr/bin/env python3
"""Fail-closed read-only audit of potential DNA × morphology × ecological functions.

Counts documented linkages only. Does not validate real tissue, genotypes, authorization
or trait causation. Never prints sample IDs, exact sites or record-level metadata.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = {
    "dna": "data/intake/aza3_nature_individual_intake_v1.csv",
    "shape": "data/intake/orientation_aspect_projection_calibration_v1.csv",
    "events": "data/intake/capitulum_guild_access_event_ledger_v1.csv",
    "effort": "data/intake/capitulum_video_effort_denominator_v1.csv",
    "pollen": "data/intake/aza3_focal_head_pollen_assay_v1.csv",
    "seed": "data/intake/aza3_focal_head_seed_outcome_v1.csv",
}
SCHEMA = {
    "dna": {"individual_id","voucher_id","leaf_dna_tissue_id","population_id",
            "admission_status","tier1_library_id"},
    "shape": {"individual_id","population_id","capitulum_id",
              "head_gravity_angle_deg","minimum_phyllary_entry_gap_mm",
              "gravity_reference_visible","scale_reference_visible",
              "field_access_permission_id"},
    "events": {"individual_id","population_id","capitulum_id","record_id",
               "oviposition_confirmed","evidence_uri"},
    "effort": {"individual_id","population_id","capitulum_id",
               "evaluable_duration_s","valid_rate_denominator",
               "head_zone_fully_covered"},
    "pollen": {"individual_id","population_id","capitulum_id","floret_id",
               "conspecific_pollen_grains","pollen_assay_method",
               "head_link_verified","admission_status","authorization_key","evidence_uri"},
    "seed": {"individual_id","population_id","capitulum_id","achenes_assessed",
             "achenes_filled","achenes_damaged","head_link_verified",
             "admission_status","authorization_key","evidence_uri"},
}
YES = {"yes","true","1","confirmed","verified"}
ADMITTED = {"admitted","approved"}


def textval(row: dict[str,str], field: str) -> str:
    return (row.get(field) or "").strip()


def truthy(row: dict[str,str], field: str) -> bool:
    return textval(row,field).lower() in YES


def admitted(row: dict[str,str]) -> bool:
    return textval(row,"admission_status").lower() in ADMITTED


def read_rows(root: Path, kind: str) -> list[dict[str,str]]:
    path = root / FILES[kind]
    if not path.is_file():
        raise ValueError(f"Missing mandatory intake template: {FILES[kind]}")
    with path.open(encoding="utf-8-sig",newline="") as stream:
        reader = csv.DictReader(stream)
        missing = SCHEMA[kind] - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Invalid {kind} template; missing columns: {sorted(missing)}")
        return list(reader)


def key(row: dict[str,str]) -> tuple[str,str,str] | None:
    k = tuple(textval(row,col) for col in ("population_id","individual_id","capitulum_id"))
    return k if all(k) else None


def number(row: dict[str,str], col: str, *, integer: bool = False) -> float | int | None:
    raw = textval(row,col)
    if not raw:
        return None
    try:
        n = int(raw) if integer else float(raw)
    except ValueError as exc:
        raise ValueError(f"Invalid {col} value") from exc
    if not math.isfinite(float(n)) or n < 0:
        raise ValueError(f"Nonfinite or negative {col}")
    return n


def audit(root: Path) -> dict:
    tables={kind:read_rows(root,kind) for kind in FILES}
    # These are *documentary labels*, not QC-verified genotype calls.
    tissue={
        (textval(r,"population_id"),textval(r,"individual_id"))
        for r in tables["dna"]
        if textval(r,"population_id") and textval(r,"individual_id")
        and textval(r,"voucher_id") and textval(r,"leaf_dna_tissue_id") and admitted(r)
    }
    head_shapes=set()
    for r in tables["shape"]:
        k=key(r)
        if not k or not textval(r,"field_access_permission_id"):
            continue
        orientation=number(r,"head_gravity_angle_deg")
        gap=number(r,"minimum_phyllary_entry_gap_mm")
        if (orientation is not None and orientation <= 180 and gap is not None
            and truthy(r,"gravity_reference_visible")
            and truthy(r,"scale_reference_visible")):
            head_shapes.add(k)

    pollen_heads=set()
    pollen_keys=set()
    for r in tables["pollen"]:
        k=key(r)
        if not k or not admitted(r) or not truthy(r,"head_link_verified"):
            continue
        if not (textval(r,"authorization_key") and textval(r,"evidence_uri")
                and textval(r,"pollen_assay_method") and textval(r,"floret_id")):
            continue
        n=number(r,"conspecific_pollen_grains",integer=True)
        if n is None:
            continue
        floret_key=(*k,textval(r,"floret_id"),textval(r,"observation_bout_id"))
        if floret_key in pollen_keys:
            raise ValueError("Duplicate pollen assay floret/bout identity")
        pollen_keys.add(floret_key)
        pollen_heads.add(k)

    seed_heads=set()
    for r in tables["seed"]:
        k=key(r)
        if not k or not admitted(r) or not truthy(r,"head_link_verified"):
            continue
        if not (textval(r,"authorization_key") and textval(r,"evidence_uri")):
            continue
        n=number(r,"achenes_assessed",integer=True)
        filled=number(r,"achenes_filled",integer=True)
        damaged=number(r,"achenes_damaged",integer=True)
        if not n or filled is None or damaged is None:
            continue
        if filled>n or damaged>n:
            raise ValueError("Seed outcomes exceed examined achenes")
        if k in seed_heads:
            raise ValueError("Duplicate head-level final seed outcome")
        seed_heads.add(k)

    effort_heads=set()
    for r in tables["effort"]:
        k=key(r)
        if not k:
            continue
        dur=number(r,"evaluable_duration_s")
        if (dur is not None and dur>0 and truthy(r,"valid_rate_denominator")
            and truthy(r,"head_zone_fully_covered")):
            effort_heads.add(k)

    oviposition_heads=set()
    for r in tables["events"]:
        k=key(r)
        if k and truthy(r,"oviposition_confirmed") and textval(r,"evidence_uri"):
            oviposition_heads.add(k)

    dna_linked_shape={k for k in head_shapes if (k[0],k[1]) in tissue}
    linked=dna_linked_shape & pollen_heads & seed_heads & effort_heads
    results={
        "contract_version":"cirsium_ecological_solution_space_readiness_v1",
        "source_record_counts":{kind:len(rows) for kind,rows in tables.items()},
        "documented_dna_tissue_individuals":len(tissue),
        "calibrated_geometry_heads":len(head_shapes),
        "documented_true_pollen_assay_heads":len(pollen_heads),
        "documented_head_seed_outcomes":len(seed_heads),
        "documented_valid_video_effort_heads":len(effort_heads),
        "documented_positive_oviposition_heads":len(oviposition_heads),
        "documented_dna_shape_pollen_seed_effort_joint_head_candidates":len(linked),
        "genotype_sequence_qc_verified_by_this_audit":False,
        "field_permits_verified_by_this_audit":False,
        "functional_equivalence_or_pareto_inference_authorized":False,
        "status":("POTENTIAL_JOINT_ROWS_REQUIRE_QC_AND_REPLICATION"
                  if linked else "HOLD_NO_JOINT_HEAD_LEVEL_GENOTYPE_FUNCTION_FITNESS"),
        "claim_boundary":"Counts of public documentary links, not independent samples or fitness effects; no genomic or causal inference.",
    }
    return results


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=ROOT)
    parser.add_argument("--require-joined-head",action="store_true")
    args=parser.parse_args()
    r=audit(args.root)
    print(json.dumps(r,ensure_ascii=False,indent=2))
    return 2 if args.require_joined_head and not r["documented_dna_shape_pollen_seed_effort_joint_head_candidates"] else 0


if __name__=="__main__":
    raise SystemExit(main())
