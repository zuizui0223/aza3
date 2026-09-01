#!/usr/bin/env python3
"""Fail-closed validation of the independent M01 mechanistic demonstration lane."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
M01_PATH = ROOT / "data" / "planning" / "chapter3_mechanistic_demonstration_v1.json"
PRIORITY_PATH = ROOT / "data" / "planning" / "chapter3_sampling_priorities_v1.csv"
M01_DOC_PATH = ROOT / "docs" / "M01_FLORAL_PIGMENTATION_MECHANISTIC_DEMONSTRATION_V1.md"
SCOPE_PATH = ROOT / "docs" / "CHAPTER3_SCOPE_AND_HANDOFF_V1.md"
README_PATH = ROOT / "README.md"

M01_LITERATURE_DOI = "10.1186/s12870-026-08097-6"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_core_priority_independence() -> list[dict[str, str]]:
    with PRIORITY_PATH.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if [row["priority_id"] for row in rows] != ["P01", "P02", "P03", "P04", "P05"]:
        raise AssertionError("M01 must not alter or extend the frozen P01-P05 core priority list")
    return rows


def validate_m01_contract() -> dict:
    m01 = load_json(M01_PATH)
    if m01.get("contract_version") != "chapter3_mechanistic_demonstration_v1":
        raise AssertionError("M01 contract version drift")
    if m01.get("demonstration_id") != "M01":
        raise AssertionError("mechanistic demonstration must remain M01")

    role = m01.get("role", {})
    if role.get("core_priority_status") != "INDEPENDENT_NOT_P06":
        raise AssertionError("M01 was promoted into the frozen P01-P05 priority list")
    if "does not reorder, replace, validate or invalidate" not in role.get("core_independence", ""):
        raise AssertionError("M01 core-independence boundary drift")
    if "worked example" not in role.get("not_primary_subject", ""):
        raise AssertionError("M01 focal taxa were promoted into the dissertation's primary taxonomic subject")

    anchor = m01.get("literature_anchor", {})
    if anchor.get("doi") != M01_LITERATURE_DOI:
        raise AssertionError("M01 literature anchor drift")
    if anchor.get("admission_status") != "REFERENCE_PREMISE_ONLY":
        raise AssertionError("M01 literature anchor was promoted to focal evidence")
    if not any("species delimitation" in item for item in anchor.get("not_estimands", [])):
        raise AssertionError("M01 began re-testing basic species delimitation as its primary aim")
    if not any("published phylogenomic backbone" in item for item in anchor.get("not_estimands", [])):
        raise AssertionError("few-locus Sanger data were promoted over the published phylogenomic backbone")

    worked = m01.get("worked_example", {})
    if worked.get("focal_concepts") != ["Cirsium brevicaule", "Cirsium irumtiense"]:
        raise AssertionError("M01 worked-example taxa drift")
    if "without presupposing loss or regain" not in worked.get("primary_question", ""):
        raise AssertionError("M01 began presupposing regain or loss")

    histories = m01.get("competing_histories", [])
    if [row.get("history_id") for row in histories] != ["H1", "H2", "H3"]:
        raise AssertionError("M01 must retain competing histories H1-H3")
    if any(row.get("prior_status") != "COMPETING_NOT_PREFERRED" for row in histories):
        raise AssertionError("an M01 history was made preferred before focal evidence")
    if "lost or strongly reduced" not in histories[0].get("statement", ""):
        raise AssertionError("H1 loss history drift")
    if "gained or regained" not in histories[1].get("statement", ""):
        raise AssertionError("H2 secondary-gain/regain history drift")
    if "neither a simple H1 nor H2" not in histories[2].get("statement", ""):
        raise AssertionError("H3 unresolved/complex-history lane drift")

    ladder = m01.get("evidence_ladder", [])
    if [row.get("level") for row in ladder] != ["E0", "E1", "E2", "E3", "E4"]:
        raise AssertionError("M01 evidence ladder must remain E0-E4")
    if "no loss, regain" not in ladder[0].get("maximum_claim", ""):
        raise AssertionError("E0 claim ceiling drift")
    if "no molecular mechanism or selection claim" not in ladder[1].get("maximum_claim", ""):
        raise AssertionError("E1 claim ceiling drift")
    if "not a causal regulatory variant" not in ladder[2].get("maximum_claim", ""):
        raise AssertionError("E2 expression ceiling drift")
    if "does not identify pollinators" not in ladder[3].get("maximum_claim", ""):
        raise AssertionError("E3 selective-agent ceiling drift")
    if "direct agent-specific ecological evidence" not in ladder[4].get("maximum_claim", ""):
        raise AssertionError("E4 selective-agent boundary drift")

    sampling = m01.get("sampling_contract", {})
    if sampling.get("field_manipulation_required") is not False:
        raise AssertionError("M01 incorrectly requires field manipulation")
    if sampling.get("pollinator_observation_required") is not False:
        raise AssertionError("M01 incorrectly requires pollinator observation")
    if "two-species FST contrast alone is not admitted" not in sampling.get(
        "selection_lane_population_requirement", ""
    ):
        raise AssertionError("M01 population-replication boundary drift")

    selection = m01.get("selection_contract", {})
    if "neutral population structure and demographic history" not in selection.get("background_null", ""):
        raise AssertionError("M01 demographic null drift")
    if "genome-size difference" not in selection.get("mapping_and_orthology_gate", ""):
        raise AssertionError("M01 mapping/orthology gate drift")
    if "genomics alone cannot identify" not in selection.get("selective_agent_boundary", ""):
        raise AssertionError("M01 selective-agent boundary drift")
    shortcuts = selection.get("prohibited_shortcuts", [])
    for required in (
        "two-species FST outlier alone",
        "branch-specific dN/dS alone",
        "post-hoc threshold tuning after seeing candidate genes",
    ):
        if required not in shortcuts:
            raise AssertionError(f"M01 prohibited shortcut missing: {required}")

    state = m01.get("current_state", {})
    if state.get("own_m01_biological_records_admitted") != 0:
        raise AssertionError("M01 own data were admitted without a new intake contract")
    for key in (
        "tissue_collection_authorized",
        "field_manipulation_authorized",
        "pollinator_agent_claim_authorized",
        "regain_claim_authorized",
        "selection_claim_authorized",
    ):
        if state.get(key) is not False:
            raise AssertionError(f"M01 fail-closed state was opened: {key}")
    return m01


def validate_narrative() -> None:
    m01_doc = M01_DOC_PATH.read_text(encoding="utf-8")
    required_m01 = [
        "Cirsium brevicaule",
        "Cirsium irumtiense",
        "without presupposing regain",
        "H1",
        "H2",
        "H3",
        "E0",
        "E1",
        "E2",
        "E3",
        "E4",
        "genomics alone cannot identify the selective agent",
        "not P06",
    ]
    missing = [needle for needle in required_m01 if needle not in m01_doc]
    if missing:
        raise AssertionError(f"M01 narrative missing: {missing}")

    scope = SCOPE_PATH.read_text(encoding="utf-8")
    required_scope = [
        "M01",
        "not P06",
        "does not invalidate P01-P05",
        "two-species FST",
        "genomics alone",
    ]
    missing = [needle for needle in required_scope if needle not in scope]
    if missing:
        raise AssertionError(f"Chapter 3 scope missing M01 boundary language: {missing}")

    readme = README_PATH.read_text(encoding="utf-8")
    required_readme = ["M01", "not P06", "mechanistic demonstration", "P01-P05"]
    missing = [needle for needle in required_readme if needle not in readme]
    if missing:
        raise AssertionError(f"README missing M01 architecture: {missing}")


def main() -> int:
    for path in (M01_PATH, PRIORITY_PATH, M01_DOC_PATH, SCOPE_PATH, README_PATH):
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing or empty M01 file: {path.relative_to(ROOT)}")
    priorities = validate_core_priority_independence()
    m01 = validate_m01_contract()
    validate_narrative()
    print("chapter3_m01_mechanistic_demonstration_valid=true")
    print(f"core_priorities={len(priorities)}")
    print(f"mechanistic_demonstration={m01['demonstration_id']}")
    print("own_m01_biological_records_admitted=0")
    print("regain_claim_authorized=false")
    print("selection_claim_authorized=false")
    print("pollinator_agent_claim_authorized=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
