#!/usr/bin/env python3
"""Validate the executable RAD-seq pilot design and library allocation."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "data" / "planning" / "chapter3_radseq_pilot_execution_v1.json"
CANDIDATE_PATH = ROOT / "data" / "planning" / "chapter3_radseq_pilot_protocol_candidates_v1.csv"
ALLOCATION_PATH = ROOT / "data" / "planning" / "chapter3_radseq_pilot_library_allocation_v1.csv"
ANCHOR_PATH = ROOT / "data" / "planning" / "chapter3_radseq_pilot_anchor_ledger_v1.csv"
DOC_PATH = ROOT / "docs" / "CHAPTER3_RADSEQ_PILOT_EXECUTION_V1.md"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def load_contract() -> dict:
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def validate_contract() -> dict:
    d = load_contract()
    if d.get("contract_version") != "chapter3_radseq_pilot_execution_v1":
        raise AssertionError("RAD pilot execution contract version drift")
    candidates = d["candidate_protocols"]
    expected_ids = ["RAD_C01_AVAII_MSPI", "RAD_C02_ECORI_MSPI", "RAD_C03_PSTI_MSPI"]
    if [row["protocol_id"] for row in candidates] != expected_ids:
        raise AssertionError("RAD pilot enzyme candidates drift")
    if any(row["status"] != "SCREEN_CANDIDATE_NOT_WINNER" for row in candidates):
        raise AssertionError("an enzyme candidate was silently promoted before Stage A")
    if any(row["size_window_status"] != "SCREEN_ONLY_NOT_PRODUCTION_FROZEN" for row in candidates):
        raise AssertionError("400-700 pilot window was silently frozen as production")

    a = d["stage_A_screen"]
    if (a["anchor_systems"], a["templates_per_anchor"], a["candidate_protocols_per_template"], a["planned_libraries"]) != (8, 1, 3, 24):
        raise AssertionError("Stage A 8 x 1 x 3 = 24 design drift")
    if "equal-read downsampling" not in a["sequencing_rule"]:
        raise AssertionError("Stage A protocol comparison is not depth-normalized")
    if "technical Pareto" not in a["advancement_rule"]:
        raise AssertionError("Stage A winner can be chosen by biological outcomes")

    b = d["stage_B_reproducibility"]
    if (b["anchor_systems"], b["biological_templates_per_anchor"], b["advanced_protocols"], b["primary_libraries"], b["technical_repeat_libraries"], b["planned_libraries_total"]) != (8, 2, 2, 32, 8, 40):
        raise AssertionError("Stage B 32 primary + 8 repeat = 40 design drift")
    gate = b["technical_concordance_gate"]
    if gate["minimum_overlap_genotype_concordance"] != 0.95:
        raise AssertionError("technical genotype concordance gate drift")
    if gate["minimum_core_locus_recovery"] != 0.90:
        raise AssertionError("technical core-locus recovery gate drift")
    if "four repeats assigned to TOP1 and four to TOP2" not in b["repeat_balance"]:
        raise AssertionError("technical repeat protocol balance drift")
    if "cancel a single all-Japan SNP matrix" not in b["decision_rule"]:
        raise AssertionError("single-protocol failure no longer cancels pooled matrix")

    sat = d["pilot_depth_saturation"]
    if sat["downsampling_fractions"] != [0.2, 0.4, 0.6, 0.8, 1.0]:
        raise AssertionError("pilot downsampling grid drift")
    if "less than 10 percent additional usable loci" not in sat["plateau_rule"]:
        raise AssertionError("production read target lost plateau rule")

    stacks = d["stacks_pilot_rule"]
    if stacks["m_grid"] != [2,3,4,5,6] or stacks["M_grid"] != [1,2,3,4,5,6]:
        raise AssertionError("pilot Stacks grid drift")
    if "within focal or population-like subsets" not in stacks["within_population_parameter_optimization"]:
        raise AssertionError("r80 was incorrectly generalized across Japan38")
    if "6/8, 7/8, 8/8" not in stacks["cross_concept_gate"]:
        raise AssertionError("cross-concept anchor occupancy gate drift")

    if any(d["current_state"].values()):
        raise AssertionError("RAD pilot current state was opened without new evidence")
    return d


def validate_candidate_registry() -> list[dict[str, str]]:
    rows = read_csv(CANDIDATE_PATH)
    if [r["protocol_id"] for r in rows] != ["RAD_C01_AVAII_MSPI", "RAD_C02_ECORI_MSPI", "RAD_C03_PSTI_MSPI"]:
        raise AssertionError("protocol candidate registry drift")
    if any(r["stage_a_status"] != "SCREEN_CANDIDATE_NOT_WINNER" for r in rows):
        raise AssertionError("candidate registry contains a preselected winner")
    return rows


def validate_allocation() -> list[dict[str, str]]:
    rows = read_csv(ALLOCATION_PATH)
    if len(rows) != 64:
        raise AssertionError(f"pilot allocation must contain 64 library slots, observed {len(rows)}")
    stage_a = [r for r in rows if r["pilot_stage"] == "STAGE_A"]
    stage_b = [r for r in rows if r["pilot_stage"] == "STAGE_B"]
    repeats = [r for r in rows if r["pilot_stage"] == "STAGE_B_REPEAT"]
    if (len(stage_a), len(stage_b), len(repeats)) != (24, 32, 8):
        raise AssertionError("pilot library allocation counts drift")
    anchors = {f"RAD_A{i:02d}" for i in range(1,9)}
    if {r["anchor_id"] for r in stage_a} != anchors or {r["anchor_id"] for r in repeats} != anchors:
        raise AssertionError("Stage A/repeat anchor coverage drift")
    for anchor in anchors:
        protocols = {r["protocol_slot"] for r in stage_a if r["anchor_id"] == anchor}
        if protocols != {"RAD_C01_AVAII_MSPI", "RAD_C02_ECORI_MSPI", "RAD_C03_PSTI_MSPI"}:
            raise AssertionError(f"Stage A candidate coverage drift: {anchor}")
    if sum(r["protocol_slot"] == "TOP1" for r in repeats) != 4 or sum(r["protocol_slot"] == "TOP2" for r in repeats) != 4:
        raise AssertionError("Stage B repeats are not balanced 4/4 across TOP1/TOP2")
    if any(not r["technical_repeat_of"] for r in repeats):
        raise AssertionError("technical repeat lacks declared original library")
    return rows


def validate_anchor_alignment() -> None:
    anchors = read_csv(ANCHOR_PATH)
    if {r["anchor_id"] for r in anchors} != {f"RAD_A{i:02d}" for i in range(1,9)}:
        raise AssertionError("pilot execution lost alignment to anchor ledger")


def validate_narrative() -> None:
    text = DOC_PATH.read_text(encoding="utf-8")
    required = [
        "24 shallow pilot libraries",
        "40",
        ">= 0.95",
        ">= 0.90",
        "6/8, 7/8 and 8/8",
        "400–700 bp",
        "different restriction protocols are not concatenated into one SNP matrix",
        "Stage A authorized: **false**",
    ]
    missing = [x for x in required if x not in text]
    if missing:
        raise AssertionError(f"RAD pilot narrative missing: {missing}")


def main() -> int:
    for path in (CONTRACT_PATH, CANDIDATE_PATH, ALLOCATION_PATH, ANCHOR_PATH, DOC_PATH):
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing or empty RAD pilot file: {path.relative_to(ROOT)}")
    d = validate_contract()
    candidates = validate_candidate_registry()
    allocation = validate_allocation()
    validate_anchor_alignment()
    validate_narrative()
    print("chapter3_radseq_pilot_execution_valid=true")
    print(f"candidate_protocols={len(candidates)}")
    print(f"stage_a_libraries={sum(r['pilot_stage']=='STAGE_A' for r in allocation)}")
    print(f"stage_b_primary={sum(r['pilot_stage']=='STAGE_B' for r in allocation)}")
    print(f"stage_b_repeats={sum(r['pilot_stage']=='STAGE_B_REPEAT' for r in allocation)}")
    print(f"genotype_concordance_gate={d['stage_B_reproducibility']['technical_concordance_gate']['minimum_overlap_genotype_concordance']}")
    print("production_open=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
