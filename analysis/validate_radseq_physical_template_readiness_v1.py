#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "data" / "planning" / "chapter3_radseq_physical_template_inventory_v1.csv"
ANCHORS = ROOT / "data" / "planning" / "chapter3_radseq_pilot_anchor_ledger_v1.csv"
DOC = ROOT / "docs" / "CHAPTER3_RADSEQ_PHYSICAL_TEMPLATE_READINESS_V1.md"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_inventory() -> list[dict[str, str]]:
    rows = read_csv(INVENTORY)
    if len(rows) != 16:
        raise AssertionError(f"physical inventory must contain 16 template rows, observed {len(rows)}")
    expected = {(f"RAD_A{i:02d}", role) for i in range(1, 9) for role in ("TEMPLATE_1", "TEMPLATE_2")}
    observed = {(r["anchor_id"], r["template_role"]) for r in rows}
    if observed != expected:
        raise AssertionError("physical template inventory lost the 8 anchors x 2 templates structure")

    anchor_rows = read_csv(ANCHORS)
    taxon_by_anchor = {r["anchor_id"]: r["taxon_concept"] for r in anchor_rows}
    for row in rows:
        if row["taxon_concept"] != taxon_by_anchor[row["anchor_id"]]:
            raise AssertionError(f"taxon mismatch for {row['anchor_id']} {row['template_role']}")

    stage_a = [r for r in rows if r["template_role"] == "TEMPLATE_1"]
    stage_b = [r for r in rows if r["template_role"] == "TEMPLATE_2"]
    if len(stage_a) != 8 or len(stage_b) != 8:
        raise AssertionError("Stage A/Stage B template balance drift")
    if any(r["stage_requirement"] != "STAGE_A_REQUIRED" for r in stage_a):
        raise AssertionError("TEMPLATE_1 ceased to be mandatory for Stage A")
    if any(r["stage_requirement"] != "STAGE_B_REQUIRED" for r in stage_b):
        raise AssertionError("TEMPLATE_2 ceased to be mandatory for Stage B")

    # The repository has not admitted physical holdings yet. Public accessions must never be treated as a tube in hand.
    if any(r["physical_material_status"] != "UNKNOWN_NOT_VERIFIED" for r in rows):
        raise AssertionError("physical possession was asserted without an inventory update supported by new evidence")
    if any(r["physical_sample_id"] for r in rows):
        raise AssertionError("physical sample IDs were invented before actual inventory admission")
    if any(r["stage_a_eligible"] != "false" or r["stage_b_eligible"] != "false" for r in rows):
        raise AssertionError("RAD pilot eligibility opened before physical template verification")
    if any(r["primary_blocker"] != "PHYSICAL_SAMPLE_NOT_CONFIRMED" for r in rows):
        raise AssertionError("initial physical-readiness blocker drift")

    joined = " ".join(r["candidate_acquisition_route"] for r in rows)
    if "public sequence references do not count as physical pilot DNA" not in joined and "public sequence records do not establish physical possession" not in joined:
        raise AssertionError("public-reference versus physical-material boundary missing")
    return rows


def validate_narrative() -> None:
    text = DOC.read_text(encoding="utf-8")
    for required in (
        "public sequence/reference exists",
        "physical tissue or DNA template is actually available",
        "all 16 planned biological templates start as `UNKNOWN_NOT_VERIFIED`",
        "Stage A remains unauthorized until all eight `TEMPLATE_1` rows",
        "TEMPLATE_2 must be biologically independent of TEMPLATE_1",
        "does not authorize collection",
    ):
        if required not in text:
            raise AssertionError(f"physical readiness narrative missing: {required}")


def main() -> int:
    for path in (INVENTORY, ANCHORS, DOC):
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing or empty physical readiness file: {path.relative_to(ROOT)}")
    rows = validate_inventory()
    validate_narrative()
    print("chapter3_radseq_physical_template_readiness_valid=true")
    print(f"template_rows={len(rows)}")
    print("verified_stage_a_templates=0")
    print("stage_A_authorized=false")
    print("stage_B_authorized=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
