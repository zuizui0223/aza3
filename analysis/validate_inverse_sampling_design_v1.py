#!/usr/bin/env python3
"""Fail-closed validation of the claim-backward Chapter 3 sampling design."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DESIGN_PATH = ROOT / "data" / "planning" / "chapter3_inverse_sampling_design_v1.json"
DOC_PATH = ROOT / "docs" / "CHAPTER3_INVERSE_SAMPLING_DESIGN_V1.md"
README_PATH = ROOT / "README.md"
SCOPE_PATH = ROOT / "docs" / "CHAPTER3_SCOPE_AND_HANDOFF_V1.md"


def load_design() -> dict:
    return json.loads(DESIGN_PATH.read_text(encoding="utf-8"))


def validate_design() -> dict:
    d = load_design()
    if d.get("contract_version") != "chapter3_inverse_sampling_design_v1":
        raise AssertionError("inverse sampling contract version drift")
    if "claim-backward" not in d.get("design_principle", ""):
        raise AssertionError("inverse sampling lost claim-backward design")

    core = d["core_breadth_panel"]
    floor = core["all_japan_floor"]
    if (floor["concepts"], floor["biological_individuals_per_concept"], floor["total_individuals"]) != (38, 3, 114):
        raise AssertionError("all-Japan 38 x 3 floor drift")

    focal = {row["concept"]: row for row in core["focal_oversampling"]}
    expected = {
        "JPN_36 Cirsium sieboldii": (30, 40),
        "JPN_06 Cirsium dipsacolepis": (16, 24),
        "JPN_15 Cirsium lineare": (16, 24),
    }
    for concept, totals in expected.items():
        row = focal.get(concept)
        if row is None or (row["minimum_total_individuals"], row["recommended_total_individuals"]) != totals:
            raise AssertionError(f"core focal sample target drift: {concept}")
    if core["derived_totals"] != {
        "minimum_primary_individuals": 167,
        "recommended_primary_individuals": 193,
        "calculation": "114 all-Japan floor plus focal top-ups for JPN36, JPN06 and JPN15; P04 reuses these individuals and P05 initially reuses their calibrated images rather than adding fresh field samples",
    }:
        raise AssertionError("core derived totals drift")
    if "no additional fresh-collection quota" not in core["p05_rule"]:
        raise AssertionError("P05 silently opened new field sampling")

    m01 = d["m01_depth_panel"]
    m1 = m01["stage_m1_history_and_expression_discovery"]
    if len(m1["focal_population_roles"]) != 2:
        raise AssertionError("M01 discovery must retain two focal lineages")
    bank = m1["per_population_field_bank"]
    if (bank["primary_individuals"], bank["predeclared_technical_reserves"]) != (15, 3):
        raise AssertionError("M01 discovery field-bank target drift")
    rna = m1["rna_subset"]
    if (rna["selected_populations"], rna["collect_individuals_per_population"], rna["sequence_primary_individuals_per_population"], rna["planned_libraries_at_one_frozen_developmental_stage"]) != (4, 6, 5, 20):
        raise AssertionError("M01 RNA design drift")
    pigment = m1["pigment_subset"]
    if (pigment["primary_individuals_per_population"], pigment["total_focal_individuals"]) != (3, 12):
        raise AssertionError("M01 pigment subset drift")
    out = m1["ancestral_state_outgroup_panel"]
    if (out["minimum_nonfocal_concepts"], out["minimum_individuals_per_concept"], out["minimum_nonfocal_individuals"]) != (6, 3, 18):
        raise AssertionError("M01 ancestral-state outgroup panel drift")
    if "per observed colour state" not in out["polymorphism_rule"]:
        raise AssertionError("M01 reintroduced taxon-level colour compression")

    e3 = m01["stage_m2_population_genomic_selection"]
    final = e3["final_population_design"]
    if len(final["Cirsium brevicaule"]) != 4 or len(final["Cirsium irumtiense"]) != 4:
        raise AssertionError("M01 E3 must retain four populations per lineage")
    if (e3["primary_individuals_per_population"], e3["populations_per_lineage"], e3["primary_individuals_per_lineage"], e3["primary_individuals_total"]) != (15, 4, 60, 120):
        raise AssertionError("M01 E3 population-genomic target drift")
    if "RAD-only cross-species FST is not sufficient" not in e3["sequencing_platform_gate"]:
        raise AssertionError("M01 E3 allowed RAD-only two-species selection shortcut")

    gates = d["assay_staging"]
    if [row["gate"] for row in gates] != [
        "G0 identity/authorization",
        "G1 core breadth",
        "G2 M01 E1",
        "G3 M01 E2",
        "G4 M01 E3",
    ]:
        raise AssertionError("sampling gate sequence drift")
    if "only if E2 produces a replicable candidate mechanism" not in gates[-1]["action"]:
        raise AssertionError("E3 sequencing was opened before E2 promotion gate")

    state = d["current_state"]
    if state["sampling_authorized"] is not False or state["biological_records_admitted"] != 0:
        raise AssertionError("inverse sampling design opened biological sampling")
    if state["e3_population_genomic_sequencing_open"] is not False:
        raise AssertionError("E3 sequencing opened before evidence gate")
    return d


def validate_narrative() -> None:
    doc = DOC_PATH.read_text(encoding="utf-8")
    required = [
        "38 concepts × 3 individuals = **114",
        "minimum core primary n = 167",
        "20 primary RNA-seq libraries",
        "120 primary E3 individuals",
        "not the first assay",
        "two populations per lineage",
        "no bagging",
    ]
    missing = [x for x in required if x not in doc]
    if missing:
        raise AssertionError(f"inverse sampling narrative missing: {missing}")

    readme = README_PATH.read_text(encoding="utf-8")
    if "CHAPTER3_INVERSE_SAMPLING_DESIGN_V1.md" not in readme:
        raise AssertionError("README does not route to inverse sampling design")
    scope = SCOPE_PATH.read_text(encoding="utf-8")
    if "claim-backward sampling" not in scope:
        raise AssertionError("scope does not state claim-backward sampling architecture")


def main() -> int:
    for path in (DESIGN_PATH, DOC_PATH, README_PATH, SCOPE_PATH):
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing or empty inverse sampling file: {path.relative_to(ROOT)}")
    d = validate_design()
    validate_narrative()
    print("chapter3_inverse_sampling_design_valid=true")
    print(f"core_minimum_primary={d['core_breadth_panel']['derived_totals']['minimum_primary_individuals']}")
    print(f"core_recommended_primary={d['core_breadth_panel']['derived_totals']['recommended_primary_individuals']}")
    print(f"m01_rna_primary={d['m01_depth_panel']['stage_m1_history_and_expression_discovery']['rna_subset']['planned_libraries_at_one_frozen_developmental_stage']}")
    print(f"m01_e3_primary={d['m01_depth_panel']['stage_m2_population_genomic_selection']['primary_individuals_total']}")
    print("sampling_authorized=false")
    print("e3_population_genomic_sequencing_open=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
