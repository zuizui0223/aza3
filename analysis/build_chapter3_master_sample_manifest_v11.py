#!/usr/bin/env python3
"""Build the authoritative public-safe Chapter 3 one-row-per-physical-source manifest.

The manifest de-duplicates the nationwide tree representation slots against nested
P02 RAD and M01 discovery plants. It deliberately stores no exact coordinates or
sensitive micro-locality text.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "data" / "planning"
TREE = P / "chapter3_exact_site_freeze_queue_v8.csv"
P02 = P / "chapter3_p02_recommended_sample_manifest_v1.csv"
M01 = P / "m01_operational_population_ledger_v1.csv"
BLOCKS = P / "chapter3_master_sample_blocks_v11.csv"
OUT = P / "chapter3_master_sample_manifest_v11.csv"
SUMMARY = P / "chapter3_master_sample_summary_v11.json"

FIELDS = [
    "master_sample_id",
    "taxon",
    "geographic_unit",
    "primary_program",
    "program_roles",
    "stage_status",
    "tree_slot_id",
    "tree_priority",
    "required_for_tree",
    "source_slot",
    "population_slot",
    "acquisition_mode",
    "authorization_gate",
    "planned_material",
    "minimum_core",
    "deidentified_locality_id",
    "public_exact_coordinates_allowed",
    "notes",
]

TAXON_EQUIV = {
    # Current NMNS operational spelling in the tree ledger versus focal M01 spelling.
    "Cirsium irumtiense": "Cirsium irimtiense",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h))


def norm_taxon(name: str) -> str:
    return TAXON_EQUIV.get(name, name)


def public_safe(text: str) -> str:
    # Inputs are already range-sector/deidentified planning data. Guard against
    # accidentally introducing coordinate-bearing fields into the master output.
    low = text.lower()
    if "latitude" in low or "longitude" in low or "gps" in low:
        raise AssertionError("coordinate-like field leaked into public master manifest")
    return text


def blank_row() -> dict[str, str]:
    return {k: "" for k in FIELDS}


def tree_rows() -> list[dict[str, str]]:
    if not TREE.exists():
        raise AssertionError("run build_core_conservation_screen_v10.py and build_exact_site_freeze_queue_v8.py first")
    rows = read_rows(TREE)
    if len(rows) != 229:
        raise AssertionError(f"tree queue must contain 229 rows, got {len(rows)}")
    out: list[dict[str, str]] = []
    for r in rows:
        o = blank_row()
        required = r["required_for_tree"] == "true"
        o.update({
            "master_sample_id": f"TREE::{r['slot_id']}",
            "taxon": r["species_binomial"],
            "geographic_unit": public_safe(r["range_sector"]),
            "primary_program": "NATIONWIDE_TREE",
            "program_roles": "NATIONWIDE_TREE",
            "stage_status": "ACTIVE_TREE_REQUIRED" if required else "OPTIONAL_TREE_TRAIT",
            "tree_slot_id": r["slot_id"],
            "tree_priority": r["priority"],
            "required_for_tree": "true" if required else "false",
            "source_slot": r["sample_slot"],
            "population_slot": "",
            "acquisition_mode": r["acquisition_mode"],
            "authorization_gate": r["freeze_status"],
            "planned_material": "Comp1061-compatible target capture or conservation-safe equivalent; live trait linkage when fresh material is used",
            "minimum_core": "yes" if required else "no",
            "deidentified_locality_id": r["target_locality_id"],
            "public_exact_coordinates_allowed": "false",
            "notes": r["exclusion_or_block_reason"],
        })
        out.append(o)
    return out


def choose_tree_row(rows: list[dict[str, str]], taxon: str, sector_token: str | None = None) -> dict[str, str]:
    nt = norm_taxon(taxon)
    candidates = [r for r in rows if norm_taxon(r["taxon"]) == nt and r["required_for_tree"] == "true"]
    if sector_token:
        hit = [r for r in candidates if sector_token.lower() in r["geographic_unit"].lower()]
        if hit:
            candidates = hit
    if len(candidates) != 1:
        raise AssertionError(f"expected one tree row for {taxon} token={sector_token!r}, got {len(candidates)}")
    return candidates[0]


def add_role(row: dict[str, str], role: str, note: str) -> None:
    roles = row["program_roles"].split("|") if row["program_roles"] else []
    if role not in roles:
        roles.append(role)
    row["program_roles"] = "|".join(roles)
    if note:
        row["notes"] = (row["notes"] + " | " + note).strip(" |").strip()


def append_p02(rows: list[dict[str, str]]) -> None:
    p02 = read_rows(P02)
    if len(p02) != 48:
        raise AssertionError(f"P02 manifest must contain 48 rows, got {len(p02)}")
    overlap_slots = {"JPN06": "J06-P1-01", "JPN15": "J15-P1-01"}
    for taxon_concept, planned in overlap_slots.items():
        source = next(r for r in p02 if r["planned_slot"] == planned)
        tree = choose_tree_row(rows, source["source_taxon_name"])
        add_role(tree, f"P02_RAD::{taxon_concept}", f"nested P02 sample {planned}; also fills own P6 tree complement")
        tree["population_slot"] = source["population_slot"]
        tree["planned_material"] += "; RAD tissue; gland/exudate record; cytotype/genome-size status"

    for r in p02:
        if r["planned_slot"] in overlap_slots.values():
            continue
        o = blank_row()
        core = r["minimum_core"] == "yes"
        o.update({
            "master_sample_id": f"P02::{r['planned_slot']}",
            "taxon": r["source_taxon_name"],
            "geographic_unit": "private authorized population slot; exact locality not public",
            "primary_program": "P02_RAD",
            "program_roles": f"P02_RAD::{r['taxon_concept']}",
            "stage_status": "ACTIVE_P02_MINIMUM" if core else "RECOMMENDED_P02_EXTENSION",
            "required_for_tree": "false",
            "source_slot": r["planned_slot"],
            "population_slot": r["population_slot"],
            "acquisition_mode": "AUTHORIZED_FRESH_POPULATION_SAMPLE",
            "authorization_gate": r["authorization_status"],
            "planned_material": "RAD tissue + voucher + same-individual stickiness/gland/orientation/phyllary + cytotype/genome-size",
            "minimum_core": "yes" if core else "no",
            "deidentified_locality_id": r["deidentified_locality_id"],
            "public_exact_coordinates_allowed": "false",
            "notes": "P02 own-data population ancestry sample; historical localities cannot be substituted for verified current populations",
        })
        rows.append(o)


def m01_stage(role: str) -> str:
    if role in {"DISCOVERY_AND_E3"}:
        return "ACTIVE_M01_DISCOVERY"
    if role == "E3_EXPANSION":
        return "CONDITIONAL_M01_E3"
    if role in {"E1_SISTER_CLADE_BRACKET", "E1_STATE_BALANCED_DEEP_BRACKET"}:
        return "MATERIAL_GAP_DEPENDENT_E1"
    if role == "OPTIONAL_MOLECULAR_POSITIVE_CONTROL":
        return "OPTIONAL_M01_CONTROL"
    raise AssertionError(f"unknown M01 role: {role}")


def append_m01(rows: list[dict[str, str]]) -> None:
    m01 = read_rows(M01)
    discovery_tree_token = {
        "M01_BREV_OKI": ("Cirsium brevicaule", "Okinawa"),
        "M01_BREV_AMAMI": ("Cirsium brevicaule", "Amami"),
        "M01_IRUM_MIYAKO": ("Cirsium irumtiense", "Miyako"),
        "M01_IRUM_ISHIGAKI": ("Cirsium irumtiense", "Ishigaki"),
    }
    for r in m01:
        n_primary = int(r["primary_individuals"])
        n_reserve = int(r["reserve_individuals"])
        if r["sample_unit_id"] in discovery_tree_token:
            taxon, token = discovery_tree_token[r["sample_unit_id"]]
            tree = choose_tree_row(rows, taxon, token)
            add_role(tree, f"M01_DISCOVERY::{r['sample_unit_id']}", "one M01 discovery individual also fills this P3 tree slot")
            tree["planned_material"] += "; M01 leaf DNA/photo/trait; nested pigment/RNA subset where assigned"
            start = 2
        else:
            start = 1

        stage = m01_stage(r["population_role"])
        for i in range(start, n_primary + 1):
            o = blank_row()
            o.update({
                "master_sample_id": f"M01::{r['sample_unit_id']}::I{i:02d}",
                "taxon": r["taxon_or_role"],
                "geographic_unit": public_safe(r["geographic_unit"]),
                "primary_program": "M01",
                "program_roles": f"M01::{r['sample_unit_id']}",
                "stage_status": stage,
                "required_for_tree": "false",
                "source_slot": f"I{i:02d}",
                "population_slot": r["sample_unit_id"],
                "acquisition_mode": r["acquisition_route"],
                "authorization_gate": r["site_lock_rule"],
                "planned_material": "photo/phenotype + leaf DNA" + (" + nested pigment/RNA subset" if stage == "ACTIVE_M01_DISCOVERY" else ""),
                "minimum_core": "yes" if stage == "ACTIVE_M01_DISCOVERY" else "no",
                "public_exact_coordinates_allowed": "false",
                "notes": r["notes"],
            })
            rows.append(o)

        for j in range(1, n_reserve + 1):
            o = blank_row()
            o.update({
                "master_sample_id": f"M01::{r['sample_unit_id']}::R{j:02d}",
                "taxon": r["taxon_or_role"],
                "geographic_unit": public_safe(r["geographic_unit"]),
                "primary_program": "M01",
                "program_roles": f"M01_RESERVE::{r['sample_unit_id']}",
                "stage_status": "RECOMMENDED_M01_RESERVE",
                "required_for_tree": "false",
                "source_slot": f"R{j:02d}",
                "population_slot": r["sample_unit_id"],
                "acquisition_mode": r["acquisition_route"],
                "authorization_gate": r["site_lock_rule"],
                "planned_material": "reserve leaf DNA/photo; reserve floral material only under matched protocol",
                "minimum_core": "no",
                "public_exact_coordinates_allowed": "false",
                "notes": "reserve is insurance, not an extra independent primary replicate by default",
            })
            rows.append(o)


def append_p01b(rows: list[dict[str, str]]) -> None:
    # Conditional function experiment. Do not presume tree/P01 history samples are eligible.
    for i in range(1, 25):
        o = blank_row()
        o.update({
            "master_sample_id": f"P01B::JPN36::I{i:02d}",
            "taxon": "Cirsium sieboldii",
            "geographic_unit": "one verified authorized population; exact locality private",
            "primary_program": "P01B_FUNCTION",
            "program_roles": "P01B_PHYLLARY_ACCESS_FUNCTION",
            "stage_status": "CONDITIONAL_P01B_FUNCTION",
            "required_for_tree": "false",
            "source_slot": f"I{i:02d}",
            "population_slot": "P01B_AUTHORIZED_POPULATION",
            "acquisition_mode": "LIVE_NONDESTRUCTIVE_EXPERIMENT",
            "authorization_gate": "P01 history retained + site-specific manipulation/conservation/terminal-collection/viability/device gates",
            "planned_material": "live experimental plant; one focal head; visitor/enemy/reproductive endpoints",
            "minimum_core": "no",
            "public_exact_coordinates_allowed": "false",
            "notes": "conditional follow-up only; do not assume overlap with tree representatives",
        })
        rows.append(o)


def validate(rows: list[dict[str, str]]) -> dict:
    if len(rows) != 451:
        raise AssertionError(f"full prospective master must contain 451 unique rows, got {len(rows)}")
    ids = [r["master_sample_id"] for r in rows]
    if len(ids) != len(set(ids)):
        raise AssertionError("master_sample_id must be unique")
    if any(r["public_exact_coordinates_allowed"] != "false" for r in rows):
        raise AssertionError("exact coordinates must remain prohibited in public master")

    counts: dict[str, int] = {}
    for r in rows:
        counts[r["stage_status"]] = counts.get(r["stage_status"], 0) + 1
    expected = {
        "ACTIVE_TREE_REQUIRED": 228,
        "OPTIONAL_TREE_TRAIT": 1,
        "ACTIVE_P02_MINIMUM": 30,
        "RECOMMENDED_P02_EXTENSION": 16,
        "ACTIVE_M01_DISCOVERY": 56,
        "RECOMMENDED_M01_RESERVE": 12,
        "CONDITIONAL_M01_E3": 60,
        "MATERIAL_GAP_DEPENDENT_E1": 18,
        "OPTIONAL_M01_CONTROL": 6,
        "CONDITIONAL_P01B_FUNCTION": 24,
    }
    if counts != expected:
        raise AssertionError(f"stage counts differ: {counts}")

    active_minimum = counts["ACTIVE_TREE_REQUIRED"] + counts["ACTIVE_P02_MINIMUM"] + counts["ACTIVE_M01_DISCOVERY"]
    p02_recommended = active_minimum + counts["RECOMMENDED_P02_EXTENSION"]
    with_m01_reserve = p02_recommended + counts["RECOMMENDED_M01_RESERVE"]
    with_optional_tree = with_m01_reserve + counts["OPTIONAL_TREE_TRAIT"]
    assert active_minimum == 314
    assert p02_recommended == 330
    assert with_m01_reserve == 342
    assert with_optional_tree == 343

    tree_nested = [r for r in rows if r["stage_status"] == "ACTIVE_TREE_REQUIRED" and "|" in r["program_roles"]]
    if len(tree_nested) != 6:
        raise AssertionError(f"expected 6 nested tree/focal rows, got {len(tree_nested)}")

    return {
        "contract_version": "chapter3_master_sample_manifest_v11",
        "full_prospective_rows": len(rows),
        "stage_counts": counts,
        "active_minimum_unique_physical_sources": active_minimum,
        "with_p02_recommended_depth": p02_recommended,
        "with_m01_discovery_reserves": with_m01_reserve,
        "with_optional_p7_trait_link": with_optional_tree,
        "conditional_after_active_core": {
            "m01_e3_expansion": counts["CONDITIONAL_M01_E3"],
            "m01_e1_material_gap_dependent": counts["MATERIAL_GAP_DEPENDENT_E1"],
            "m01_optional_controls": counts["OPTIONAL_M01_CONTROL"],
            "p01b_function": counts["CONDITIONAL_P01B_FUNCTION"],
        },
        "tree_required_slots": 228,
        "tree_optional_trait_link_slots": 1,
        "tree_focal_overlap_rows": 6,
        "exact_coordinates_allowed_in_public_manifest": False,
        "authoritative_block_source": str(BLOCKS.relative_to(ROOT)),
        "notes": [
            "314 is the active minimum unique physical-source count, not a target for 314 destructive wild collections.",
            "P02 and M01 tree representatives are nested and therefore not double counted.",
            "Historical/herbarium acquisitions are representation sources but cannot supply the full live-trait bundle.",
            "Conditional M01 and P01b rows are not current collection requirements.",
        ],
    }


def main() -> int:
    rows = tree_rows()
    append_p02(rows)
    append_m01(rows)
    append_p01b(rows)
    summary = validate(rows)
    with OUT.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("chapter3_master_sample_manifest_v11_built=true")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
