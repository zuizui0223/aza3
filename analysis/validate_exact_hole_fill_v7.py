#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "data/planning/chapter3_exact_hole_fill_summary_v7.json"
DOC = ROOT / "docs/CHAPTER3_EXACT_BACKBONE_HOLE_FILL_LEDGER_V7.md"
FILES = {
    1: [ROOT / "data/planning/chapter3_holefill_priority1_empty_blocks_v7.csv"],
    2: [ROOT / "data/planning/chapter3_holefill_priority2_broken_tips_v7.csv"],
    3: [
        ROOT / "data/planning/chapter3_holefill_priority3_kaganoazami_block_v7.csv",
        ROOT / "data/planning/chapter3_holefill_priority3_sawaazami_block_v7.csv",
        ROOT / "data/planning/chapter3_holefill_priority3_norikura_series_v7.csv",
        ROOT / "data/planning/chapter3_holefill_priority3_hamaazami_block_v7.csv",
        ROOT / "data/planning/chapter3_holefill_priority3_yamaazami_block_v7.csv",
    ],
    4: [ROOT / "data/planning/chapter3_holefill_priority4_partial_blocks_v7.csv"],
    5: [ROOT / "data/planning/chapter3_holefill_priority5_local_gaps_v7.csv"],
    6: [ROOT / "data/planning/chapter3_holefill_priority6_add_own_trait_link_v7.csv"],
    7: [ROOT / "data/planning/chapter3_holefill_priority7_tree_filled_v7.csv"],
}
SECTOR_FILES = {
    1: ROOT / "data/planning/chapter3_holefill_priority1_sample_sectors_v7.csv",
    2: ROOT / "data/planning/chapter3_holefill_priority2_sample_sectors_v7.csv",
    3: ROOT / "data/planning/chapter3_holefill_priority3_sample_sectors_v7.csv",
    4: ROOT / "data/planning/chapter3_holefill_priority4_sample_sectors_v7.csv",
    5: ROOT / "data/planning/chapter3_holefill_priority5_sample_sectors_v7.csv",
    6: ROOT / "data/planning/chapter3_holefill_priority6_sample_sectors_v7.csv",
    7: ROOT / "data/planning/chapter3_holefill_priority7_sample_sectors_v7.csv",
}
EXTENSION = ROOT / "data/planning/chapter3_holefill_taxonomy_extension_v7.csv"
CONFLICTS = ROOT / "data/planning/chapter3_holefill_taxonomy_conflicts_v7.csv"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def assert_no_exact_coordinate(text: str, species: str) -> None:
    if re.search(r"[-+]?\d{1,3}\.\d{3,}", text or ""):
        raise AssertionError(f"exact coordinate prematurely frozen for {species}")


def validate_two_sample_sector_file(path: Path, priority: int, expected_species: set[str], expected_n: int, expected_samples: int) -> None:
    rows = read_rows(path)
    if len(rows) != expected_n:
        raise AssertionError(f"P{priority} sector ledger count drift: {len(rows)}")
    if {r["species_binomial"] for r in rows} != expected_species:
        raise AssertionError(f"P{priority} sector species do not match exact hole set")
    if sum(int(r["new_tree_samples_needed"]) for r in rows) != expected_samples:
        raise AssertionError(f"P{priority} sector sample total drift")
    for row in rows:
        for field in ("nmns_distribution_source", "sample_A_sector", "sample_B_sector", "sector_design"):
            if not row.get(field, "").strip():
                raise AssertionError(f"P{priority} sector field {field} blank for {row['species_binomial']}")
        if row.get("locality_freeze_gate") != "CURRENT_OCCURRENCE_PLUS_PERMISSION_REQUIRED":
            raise AssertionError(f"P{priority} exact-locality gate drift for {row['species_binomial']}")
        assert_no_exact_coordinate(row["sample_A_sector"], row["species_binomial"])
        assert_no_exact_coordinate(row["sample_B_sector"], row["species_binomial"])


def validate_singleton_complements(path: Path, expected_species: set[str]) -> None:
    rows = read_rows(path)
    if len(rows) != 26 or {r["species_binomial"] for r in rows} != expected_species:
        raise AssertionError("P6 singleton-complement ledger drift")
    if sum(int(r["new_tree_samples_needed"]) for r in rows) != 26:
        raise AssertionError("P6 singleton-complement sample total drift")
    for row in rows:
        for field in ("nmns_distribution_source", "own_sample_sector_rule", "sector_design"):
            if not row.get(field, "").strip():
                raise AssertionError(f"P6 field {field} blank for {row['species_binomial']}")
        if row.get("public_locality_audit_gate") != "REQUIRED_BEFORE_FINAL_SECTOR_FREEZE":
            raise AssertionError(f"P6 public-locality audit gate drift for {row['species_binomial']}")
        if row.get("locality_freeze_gate") != "CURRENT_OCCURRENCE_PLUS_PERMISSION_REQUIRED":
            raise AssertionError(f"P6 exact-locality gate drift for {row['species_binomial']}")
        assert_no_exact_coordinate(row["own_sample_sector_rule"], row["species_binomial"])


def validate_p7_optional(path: Path, expected_species: set[str]) -> None:
    rows = read_rows(path)
    if len(rows) != 1 or {r["species_binomial"] for r in rows} != expected_species:
        raise AssertionError("P7 optional-trait-link ledger drift")
    row = rows[0]
    if int(row["tree_new_samples_needed"]) != 0 or int(row["optional_trait_link_samples"]) != 1:
        raise AssertionError("P7 tree/trait-link sample counts drift")
    if row.get("public_locality_audit_gate") != "REQUIRED_BEFORE_OPTIONAL_SECTOR_FREEZE":
        raise AssertionError("P7 public-locality audit gate drift")
    if row.get("locality_freeze_gate") != "CURRENT_OCCURRENCE_PLUS_PERMISSION_REQUIRED":
        raise AssertionError("P7 exact-locality gate drift")
    assert_no_exact_coordinate(row.get("optional_own_trait_link_sector_rule", ""), row["species_binomial"])


def main() -> int:
    paths = [SUMMARY, DOC, EXTENSION, CONFLICTS] + list(SECTOR_FILES.values()) + [p for ps in FILES.values() for p in ps]
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing v7 hole-fill artifact: {path.relative_to(ROOT)}")

    s = json.loads(SUMMARY.read_text(encoding="utf-8"))
    if s.get("contract_version") != "chapter3_exact_hole_fill_v7":
        raise AssertionError("v7 summary version drift")
    src = s["source_snapshot"]
    if (src["nmns_authority_records"], src["nmns_unique_species_binomials"], src["published_or_no_new_label_core_species"]) != (161, 154, 128):
        raise AssertionError("NMNS species-screen snapshot drift")
    if (src["unpublished_new_name_species"], src["provisional_name_species"]) != (25, 1):
        raise AssertionError("taxonomy extension counts drift")

    expected_species = {1: 33, 2: 6, 3: 44, 4: 16, 5: 2, 6: 26, 7: 1}
    expected_samples = {1: 66, 2: 12, 3: 88, 4: 32, 5: 4, 6: 26, 7: 0}
    seen: set[str] = set()
    by_priority: dict[int, list[dict[str, str]]] = {}
    for priority, paths_for_priority in FILES.items():
        rows = [row for path in paths_for_priority for row in read_rows(path)]
        by_priority[priority] = rows
        if len(rows) != expected_species[priority]:
            raise AssertionError(f"priority {priority} species count drift: {len(rows)}")
        if sum(int(row["new_tree_samples_needed"]) for row in rows) != expected_samples[priority]:
            raise AssertionError(f"priority {priority} sample count drift")
        for row in rows:
            sp = row["species_binomial"]
            if sp in seen:
                raise AssertionError(f"species appears in multiple core priorities: {sp}")
            seen.add(sp)

    if len(seen) != 128:
        raise AssertionError(f"published-core ledger must contain 128 unique species, got {len(seen)}")
    if sum(expected_samples.values()) != 228:
        raise AssertionError("internal sample arithmetic error")

    for priority in (1, 2, 3, 4, 5):
        validate_two_sample_sector_file(
            SECTOR_FILES[priority],
            priority,
            {r["species_binomial"] for r in by_priority[priority]},
            expected_species[priority],
            expected_samples[priority],
        )
    validate_singleton_complements(SECTOR_FILES[6], {r["species_binomial"] for r in by_priority[6]})
    validate_p7_optional(SECTOR_FILES[7], {r["species_binomial"] for r in by_priority[7]})

    p2 = {row["species_binomial"] for row in by_priority[2]}
    p3 = {row["species_binomial"] for row in by_priority[3]}
    p6 = {row["species_binomial"] for row in by_priority[6]}
    if "Cirsium sieboldii" not in p2:
        raise AssertionError("C. sieboldii must remain a broken-tip replacement drawn from P01 bank")
    if not {"Cirsium brevicaule", "Cirsium irimtiense"}.issubset(p3):
        raise AssertionError("M01 species must remain Hamaazami-block missing tips")
    if not {"Cirsium dipsacolepis", "Cirsium lineare"}.issubset(p6):
        raise AssertionError("P02 focal species must remain public-single + own-trait-link complements")

    ext = read_rows(EXTENSION)
    if len(ext) != 26 or sum(int(r["new_samples_if_extension_admitted"]) for r in ext) != 51:
        raise AssertionError("taxonomy extension lane drift")
    if sum(r["slot_priority"] == "P1_PROVISIONAL_NAME" for r in ext) != 1:
        raise AssertionError("expected exactly one provisional-name unit")
    if sum(r["slot_priority"] == "P2_UNPUBLISHED_NEW_NAME" for r in ext) != 25:
        raise AssertionError("expected exactly 25 unpublished-new-name units")

    conflicts = read_rows(CONFLICTS)
    if {r["member_id"] for r in conflicts} != {"JPN_29", "JPN_31", "JPN_33"}:
        raise AssertionError("taxonomy-gate concept set drift")

    doc = DOC.read_text(encoding="utf-8")
    for term in ("128", "228", "Priority 1", "33", "Priority 2", "six broken/weak", "25 `新称`", "JPN29"):
        if term not in doc:
            raise AssertionError(f"v7 narrative missing {term!r}")

    print("chapter3_exact_hole_fill_v7_valid=true")
    print("published_core_species=128")
    print("new_tree_samples=228")
    print("two_sample_sector_species=101")
    print("two_sample_sector_samples=202")
    print("priority6_singleton_complements=26")
    print("priority7_optional_trait_link=1")
    print("extension_species=26")
    print("extension_new_samples_if_admitted=51")
    print("taxonomy_gate_concepts=3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
