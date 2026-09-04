#!/usr/bin/env python3
"""Fail-closed validation of the Chapter 3 RAD-seq end-to-end contract."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "data" / "planning" / "chapter3_radseq_end_to_end_v1.json"
ANCHOR_PATH = ROOT / "data" / "planning" / "chapter3_radseq_pilot_anchor_ledger_v1.csv"
LIBRARY_INTAKE_PATH = ROOT / "data" / "intake" / "chapter3_radseq_library_intake_v1.csv"
DOC_PATH = ROOT / "docs" / "CHAPTER3_RADSEQ_END_TO_END_V1.md"
README_PATH = ROOT / "README.md"

LIBRARY_INTAKE_FIELDS = [
    "individual_id", "rad_tissue_id", "extraction_id", "extraction_batch",
    "dna_concentration_ng_ul", "dna_input_ng", "dna_purity_state", "dna_integrity_state",
    "rad_protocol_id", "enzyme_pair_id", "size_window_id", "sample_index", "rad_library_id",
    "library_batch", "pcr_cycles", "technical_duplicate_group", "sequencing_run_id",
    "lane_or_partition", "read_configuration", "raw_read_pairs", "retained_read_pairs",
    "library_complexity_state", "median_locus_depth", "usable_locus_count", "missingness_rate",
    "technical_replicate_concordance", "library_qc_status", "exclusion_reason",
]


def load_contract() -> dict:
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_anchor_ledger() -> list[dict[str, str]]:
    rows = read_rows(ANCHOR_PATH)
    if [row["anchor_id"] for row in rows] != [f"RAD_A{i:02d}" for i in range(1, 9)]:
        raise AssertionError("RAD pilot anchors must remain ordered RAD_A01-RAD_A08")
    if sum(int(row["biological_templates"]) for row in rows) != 16:
        raise AssertionError("RAD pilot must retain 16 biological templates")
    if any(row["counts_toward_core_biological_n"] != "false" for row in rows):
        raise AssertionError("assay pilot material was promoted into biological sample n")
    required = {
        "JPN_03 Cirsium alpicola",
        "JPN_05 Cirsium aomorense",
        "JPN_06 Cirsium dipsacolepis",
        "JPN_15 Cirsium lineare",
        "JPN_35 Cirsium nipponicum",
        "JPN_36 Cirsium sieboldii",
        "Cirsium brevicaule",
        "Cirsium irumtiense",
    }
    if {row["taxon_concept"] for row in rows} != required:
        raise AssertionError("RAD pilot anchor concepts drift")
    return rows


def validate_empty_library_intake() -> None:
    with LIBRARY_INTAKE_PATH.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != LIBRARY_INTAKE_FIELDS:
            raise AssertionError("RAD library intake schema drift")
        rows = list(reader)
    if rows:
        raise AssertionError("RAD library intake must remain empty before production RAD authorization")
    forbidden = {"latitude", "longitude", "exact_locality", "permit_document"}
    if forbidden.intersection(LIBRARY_INTAKE_FIELDS):
        raise AssertionError("sensitive locality fields entered RAD library intake")


def validate_contract() -> dict:
    d = load_contract()
    if d.get("contract_version") != "chapter3_radseq_end_to_end_v1":
        raise AssertionError("RAD contract version drift")
    if "population-ancestry" not in d.get("role", ""):
        raise AssertionError("RAD role drifted away from population ancestry")

    products = d["primary_products"]
    if "secondary strict-shared-locus" not in products["R2_all_japan_sensitivity"]:
        raise AssertionError("all-Japan RAD was promoted beyond a secondary sensitivity")
    if "no pooled diploid matrix" not in products["R3_mixed_ploidy_overlay"]:
        raise AssertionError("mixed ploidy was forced into one diploid matrix")
    if "RAD-only selection outliers are not final E3 evidence" not in products["R4_m01_selection_support"]:
        raise AssertionError("RAD-only M01 selection was promoted")

    pilot = d["pilot"]
    if pilot["required"] is not True or pilot["anchor_biological_templates"] != 16:
        raise AssertionError("RAD pilot gate drift")
    if "at least three candidate" not in pilot["stage_A_enzyme_screen"]["rule"]:
        raise AssertionError("enzyme-pair screen was removed")
    if pilot["stage_B_reproducibility_screen"]["technical_repeat_minimum"] != 8:
        raise AssertionError("pilot technical-repeat minimum drift")
    if "prohibit concatenating different RAD protocols" not in pilot["stratification_fallback"]:
        raise AssertionError("protocol-stratification boundary drift")

    lab = d["dna_and_library"]
    if lab["whole_genome_amplification_allowed"] is not False:
        raise AssertionError("whole-genome amplification was opened in the primary RAD protocol")
    if lab["paired_end_required"] is not True:
        raise AssertionError("paired-end production requirement drift")
    if "5-10 percent" not in lab["production_technical_duplicates"]:
        raise AssertionError("production technical-duplicate rule drift")
    if "no focal population may be confounded" not in lab["batch_randomization"]:
        raise AssertionError("batch/population anti-confounding rule drift")

    bio = d["bioinformatics"]
    if bio["primary_assembly"] != "Stacks 2 de novo assembly because no single lineage-matched nuclear reference is assumed across Japan38":
        raise AssertionError("primary RAD assembly drift")
    stacks = bio["stacks_parameter_freeze"]
    if stacks["m_grid"] != [2, 3, 4, 5, 6] or stacks["M_grid"] != [1, 2, 3, 4, 5, 6]:
        raise AssertionError("Stacks parameter grid drift")
    if "r80-style" not in stacks["optimization_basis"]:
        raise AssertionError("outcome-blind r80 optimization was removed")
    if "before inspecting P01-P05" not in stacks["freeze_timing"]:
        raise AssertionError("Stacks parameters can be retuned after trait outcomes")

    matrices = d["analysis_matrices"]
    if "secondary sensitivity/network" not in matrices["all_japan_strict_matrix"]["claim_ceiling"]:
        raise AssertionError("all-Japan RAD claim ceiling drift")
    mixed = matrices["mixed_ploidy_matrices"]
    if "polyRAD" not in mixed["polyploid"]:
        raise AssertionError("ploidy-aware probabilistic route missing")
    if "no single diploid genotype caller" not in mixed["cross_ploidy"]:
        raise AssertionError("cross-ploidy diploid-caller prohibition missing")

    m01 = d["m01_selection_boundary"]
    if "RAD-only cross-species FST" not in " ".join(m01["rad_not_sufficient"]):
        if "two-species FST outlier" not in m01["rad_not_sufficient"]:
            raise AssertionError("RAD-only cross-species selection shortcut reopened")
    if "non-RAD confirmation route" not in m01["promotion_rule"]:
        raise AssertionError("M01 E3 no longer requires independent candidate-region confirmation")
    if "restriction-site dropout" not in m01["dropout_audit"]:
        raise AssertionError("RAD dropout audit missing from M01 E3")

    state = d["current_state"]
    for key in (
        "rad_pilot_completed",
        "enzyme_pair_frozen",
        "size_window_frozen",
        "production_read_target_frozen",
        "single_protocol_viability_passed",
        "production_rad_authorized",
        "m01_rad_only_selection_claim_authorized",
    ):
        if state[key] is not False:
            raise AssertionError(f"RAD fail-closed state was opened: {key}")
    return d


def validate_narrative() -> None:
    doc = DOC_PATH.read_text(encoding="utf-8")
    required = [
        "RAD-seq is not one universal answer",
        "16 anchor templates",
        "Single-protocol viability gate",
        "Stacks 2 de novo",
        "m = 2–6",
        "M = 1–6",
        "polyRAD",
        "RAD alone is **not sufficient**",
        "targeted capture, amplicon resequencing or qualified low-pass/whole-genome resequencing",
        "no RAD-only selection claim",
    ]
    missing = [x for x in required if x not in doc]
    if missing:
        raise AssertionError(f"RAD narrative missing: {missing}")
    readme = README_PATH.read_text(encoding="utf-8")
    if "CHAPTER3_RADSEQ_END_TO_END_V1.md" not in readme:
        raise AssertionError("README does not route to RAD end-to-end design")
    if "chapter3_radseq_library_intake_v1.csv" not in readme:
        raise AssertionError("README does not route to RAD library intake schema")


def main() -> int:
    for path in (CONTRACT_PATH, ANCHOR_PATH, LIBRARY_INTAKE_PATH, DOC_PATH, README_PATH):
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing or empty RAD contract file: {path.relative_to(ROOT)}")
    anchors = validate_anchor_ledger()
    validate_empty_library_intake()
    d = validate_contract()
    validate_narrative()
    print("chapter3_radseq_end_to_end_valid=true")
    print(f"pilot_anchors={len(anchors)}")
    print(f"pilot_biological_templates={sum(int(r['biological_templates']) for r in anchors)}")
    print(f"primary_assembly={d['bioinformatics']['primary_assembly']}")
    print("rad_library_records_admitted=0")
    print("production_rad_authorized=false")
    print("m01_rad_only_selection_claim_authorized=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
