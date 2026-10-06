#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "C_SIEBOLDII_COMBINATORIAL_SUBSTRATE_GATE_V1.md"
CON = ROOT / "data" / "contracts" / "c_sieboldii_combinatorial_substrate_gate_v1.json"
AUDIT = ROOT / "data" / "evidence" / "c_sieboldii_combinatorial_substrate_public_audit_v1.json"
LOCALITY = ROOT / "data" / "evidence" / "c_sieboldii_public_locality_prescreen_v1.json"
TEMPLATE = ROOT / "data" / "templates" / "c_sieboldii_combinatorial_census_v1.csv"
LOCALITY_LEDGER = ROOT / "data" / "planning" / "c_sieboldii_cs0_locality_priority_v1.csv"
FIELD_CARD = ROOT / "docs" / "C_SIEBOLDII_CS0_FIELD_CARD_V1.md"
TAXONOMIC = ROOT / "docs" / "C_SIEBOLDII_CS0_TAXONOMIC_DIAGNOSTIC_V1.md"
WGS = ROOT / "data" / "contracts" / "aza3_r1b_own_wgs_transferability_pilot_v1.json"
NATURE = ROOT / "data" / "contracts" / "aza3_nature_single_question_v2.json"
PANEL = ROOT / "data" / "planning" / "aza3_nature_genomic_sampling_panel_v1.csv"
README = ROOT / "README.md"

def need(text, token):
    if token not in text:
        raise AssertionError(f"missing: {token}")

def main():
    doc = DOC.read_text(encoding="utf-8")
    con = json.loads(CON.read_text(encoding="utf-8"))
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    locality = json.loads(LOCALITY.read_text(encoding="utf-8"))
    header = TEMPLATE.read_text(encoding="utf-8").splitlines()[0].split(",")
    locality_rows = list(csv.DictReader(LOCALITY_LEDGER.open(encoding="utf-8")))
    field_card = FIELD_CARD.read_text(encoding="utf-8")
    taxdoc = TAXONOMIC.read_text(encoding="utf-8")
    wgs = json.loads(WGS.read_text(encoding="utf-8"))
    nature = json.loads(NATURE.read_text(encoding="utf-8"))
    panel = list(csv.DictReader(PANEL.open(encoding="utf-8")))
    readme = README.read_text(encoding="utf-8")

    assert con["status"] == "REQUIRED_BEFORE_BIOLOGICAL_WGS_SAMPLE_SELECTION"
    assert con["focal_taxon"] == "Cirsium sieboldii"
    assert con["primary_pair"] == ["floral_colour", "anthesis_orientation"]
    assert con["orientation_geometry"]["descriptive_upright_max_deg"] == 60
    assert con["orientation_geometry"]["descriptive_nodding_min_deg"] == 120
    assert con["census"]["minimum_complete_primary_pair_A1_A2"] == 20
    assert con["variability_thresholds"]["primary_pair_required_occupied_cells_of_4"] == 3
    assert con["variability_thresholds"]["continuous_primary_pair_abs_spearman_max"] == 0.8

    # Taxonomic fail-closed rule.
    assert con["taxonomic_diagnostic"] == "docs/C_SIEBOLDII_CS0_TAXONOMIC_DIAGNOSTIC_V1.md"
    taxrule = con["primary_taxonomic_rule"]
    assert taxrule["required_confidence"] == "high"
    assert taxrule["medium_role"] == "sensitivity_only"
    assert taxrule["unresolved_role"] == "descriptive_only"
    assert "high-confidence focal C. sieboldii" in taxrule["occupancy_rule"]
    assert "taxonomically unresolved individuals cannot fill primary combinatorial cells" in con["hard_stops"]

    # Public evidence remains a prescreen, not the substrate result.
    assert audit["audit_decision"] == "PUBLIC_EVIDENCE_SUPPORTS_SPECIES_WIDE_MULTIPLICITY_BUT_NOT_A_VERIFIED_COMBINATORIAL_POPULATION_SUBSTRATE"
    assert "orientation" in audit["key_confound"]
    assert con["public_locality_prescreen"] == "data/evidence/c_sieboldii_public_locality_prescreen_v1.json"
    assert con["public_prescreen_status"] == "TARGETED_CS0_LOCALITIES_IDENTIFIED__SUBSTRATE_NOT_YET_VERIFIED"
    assert locality["status"] == "PUBLIC_PRESCREEN_SUPPORTS_TARGETED_CS0__NO_COMBINATORIAL_SUBSTRATE_YET"
    assert locality["white_form"]["type_locality"] == "Nagano Prefecture, Kirigamine, Kurumayama Moor"
    assert locality["decision"]["public_gate_class"] == "NOT_IDENTIFIABLE_FOR_COMBINATORIAL_SUBSTRATE"

    # Locality routing is frozen prospectively.
    assert con["locality_priority_ledger"] == "data/planning/c_sieboldii_cs0_locality_priority_v1.csv"
    assert con["field_card"] == "docs/C_SIEBOLDII_CS0_FIELD_CARD_V1.md"
    assert con["public_flowering_window"]["authoritative_species"] == "August-October"
    assert [x["rank"] for x in con["cs0_priority_localities"]] == [1, 2, 3]
    assert con["cs0_priority_localities"][0]["locality"].startswith("Kurumayama Moor")
    assert con["cs0_priority_localities"][1]["locality"].startswith("Kiyooka-Mukaiyama Wetland")
    assert con["cs0_priority_localities"][2]["locality"].startswith("Kashibaru Wetland")
    assert [r["site_id"] for r in locality_rows] == ["CS0_KURU", "CS0_KIYO", "CS0_KASHI"]
    assert [int(r["rank"]) for r in locality_rows] == [1, 2, 3]
    assert locality_rows[0]["primary_role"] == "white-form verification"
    assert locality_rows[1]["primary_role"] == "standardized anthesis-orientation census"
    assert locality_rows[2]["primary_role"] == "upright-flowering validation"

    # Biological WGS sample selection remains blocked.
    assert "combinatorial_substrate_gate" in wgs
    assert wgs["sample_selection_gate"]["required_before_biological_sample_selection"] is True
    assert "CS_GREEN" in wgs["sample_selection_gate"]["biological_use_open_if"]
    assert nature["existing_evidence"]["focal_same_population_combinatorial_substrate"] == "not_yet_verified"
    assert nature["focal_system_gate"] == "data/contracts/c_sieboldii_combinatorial_substrate_gate_v1.json"

    sie = [r for r in panel if r["panel_id"] == "SIE"]
    assert len(sie) == 1
    assert "CONDITIONAL" in sie[0]["geographic_design"]
    assert "Gate CS" in sie[0]["geographic_design"]
    assert "post-anthesis upright heads" in sie[0]["stop_rule"]

    for col in (
        "observation_id", "population_id", "local_set_id", "flower_stage",
        "colour_chroma", "colour_white_descriptive", "orientation_deg_gravity",
        "phyllary_rows", "stickiness_score", "taxonomic_confidence"
    ):
        assert col in header, col

    for token in (
        "0° = vertically upward",
        "180° = vertically downward",
        "A1–A2 only",
        "3 of the 4",
        "CS_GREEN",
        "CS_GREEN_LOCAL_SET",
        "CS_RED_PHENOLOGY",
        "Taxonomic identity gate",
        "taxonomic_confidence = high",
        "Do not sequence 40",
    ):
        need(doc, token)

    for token in (
        "Target **30 flowering individuals**",
        ">=20 complete A1/A2 individuals",
        "high-confidence focal",
        "Do not combine distant populations",
        "Do not choose Nature WGS individuals",
        "C_SIEBOLDII_CS0_TAXONOMIC_DIAGNOSTIC_V1.md",
    ):
        need(field_card, token)

    for token in (
        "Never use head orientation to confirm the taxon",
        "Cirsium austrokiusianum",
        "Cirsium hidapaludosum",
        "Primary CS_GREEN and CS_GREEN_LOCAL_SET calculations use high-confidence",
        "unresolved individuals do not contribute",
    ):
        need(taxdoc, token)

    need(readme, "C_SIEBOLDII_COMBINATORIAL_SUBSTRATE_GATE_V1.md")
    need(readme, "C_SIEBOLDII_CS0_FIELD_CARD_V1.md")

    print(json.dumps({
        "status": "ok",
        "gate": "CS",
        "public_substrate_status": "unverified",
        "primary_pair": "colour_x_anthesis_orientation",
        "primary_taxonomic_confidence": "high_only",
        "cs0_sites": [r["site_id"] for r in locality_rows],
        "biological_WGS_selection": "blocked_until_CS_classified",
        "full_C_sieboldii_panel": "blocked"
    }, indent=2))

if __name__ == "__main__":
    main()
