#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"/"C_SIEBOLDII_COMBINATORIAL_SUBSTRATE_GATE_V1.md"
CON=ROOT/"data"/"contracts"/"c_sieboldii_combinatorial_substrate_gate_v1.json"
AUDIT=ROOT/"data"/"evidence"/"c_sieboldii_combinatorial_substrate_public_audit_v1.json"
TEMPLATE=ROOT/"data"/"templates"/"c_sieboldii_combinatorial_census_v1.csv"
LOCALITY=ROOT/"data"/"evidence"/"c_sieboldii_public_locality_prescreen_v1.json"\nWGS=ROOT/"data"/"contracts"/"aza3_r1b_own_wgs_transferability_pilot_v1.json"
NATURE=ROOT/"data"/"contracts"/"aza3_nature_single_question_v2.json"
PANEL=ROOT/"data"/"planning"/"aza3_nature_genomic_sampling_panel_v1.csv"
LOCALITY_LEDGER=ROOT/"data"/"planning"/"c_sieboldii_cs0_locality_priority_v1.csv"\nFIELD_CARD=ROOT/"docs"/"C_SIEBOLDII_CS0_FIELD_CARD_V1.md"\nREADME=ROOT/"README.md"

def need(t,x):
    if x not in t:
        raise AssertionError(f"missing: {x}")

def main():
    d=DOC.read_text(encoding="utf-8")
    c=json.loads(CON.read_text(encoding="utf-8"))
    a=json.loads(AUDIT.read_text(encoding="utf-8"))
    loc=json.loads(LOCALITY.read_text(encoding="utf-8"))\n    w=json.loads(WGS.read_text(encoding="utf-8"))
    n=json.loads(NATURE.read_text(encoding="utf-8"))
    p=list(csv.DictReader(PANEL.open(encoding="utf-8")))
    header=TEMPLATE.read_text(encoding="utf-8").splitlines()[0].split(",")
    locality_rows=list(csv.DictReader(LOCALITY_LEDGER.open(encoding="utf-8")))\n    field_card=FIELD_CARD.read_text(encoding="utf-8")\n    readme=README.read_text(encoding="utf-8")

    assert c["status"]=="REQUIRED_BEFORE_BIOLOGICAL_WGS_SAMPLE_SELECTION"
    assert c["focal_taxon"]=="Cirsium sieboldii"
    assert c["primary_pair"]==["floral_colour","anthesis_orientation"]
    assert c["orientation_geometry"]["descriptive_upright_max_deg"]==60
    assert c["orientation_geometry"]["descriptive_nodding_min_deg"]==120
    assert c["census"]["minimum_complete_primary_pair_A1_A2"]==20
    assert c["variability_thresholds"]["primary_pair_required_occupied_cells_of_4"]==3
    assert c["variability_thresholds"]["continuous_primary_pair_abs_spearman_max"]==0.8

    assert a["audit_decision"]=="PUBLIC_EVIDENCE_SUPPORTS_SPECIES_WIDE_MULTIPLICITY_BUT_NOT_A_VERIFIED_COMBINATORIAL_POPULATION_SUBSTRATE"
    assert "orientation" in a["key_confound"]\n    assert c["public_locality_prescreen"]=="data/evidence/c_sieboldii_public_locality_prescreen_v1.json"\n    assert c["public_prescreen_status"]=="TARGETED_CS0_LOCALITIES_IDENTIFIED__SUBSTRATE_NOT_YET_VERIFIED"\n    assert [x["rank"] for x in c["cs0_priority_localities"]]==[1,2]\n    assert c["cs0_priority_localities"][0]["locality"].startswith("Kurumayama Moor")\n    assert c["cs0_priority_localities"][1]["locality"].startswith("Kiyooka-Mukaiyama Wetland")\n    assert loc["status"]=="PUBLIC_PRESCREEN_SUPPORTS_TARGETED_CS0__NO_COMBINATORIAL_SUBSTRATE_YET"\n    assert loc["white_form"]["type_locality"]=="Nagano Prefecture, Kirigamine, Kurumayama Moor"\n    assert loc["decision"]["public_gate_class"]=="NOT_IDENTIFIABLE_FOR_COMBINATORIAL_SUBSTRATE"\n    assert c["locality_priority_ledger"]=="data/planning/c_sieboldii_cs0_locality_priority_v1.csv"\n    assert c["field_card"]=="docs/C_SIEBOLDII_CS0_FIELD_CARD_V1.md"\n    assert c["public_flowering_window"]["authoritative_species"]=="August-October"\n    assert [r["site_id"] for r in locality_rows]==["CS0_KURU","CS0_KIYO"]\n    assert [int(r["rank"]) for r in locality_rows]==[1,2]\n    assert locality_rows[0]["primary_role"]=="white-form verification"\n    assert locality_rows[1]["primary_role"]=="standardized anthesis-orientation census"\n    for x in ("Target **30 flowering individuals**",">=20 complete A1/A2 individuals","Do not combine distant populations","Do not choose Nature WGS individuals"):\n        need(field_card,x)

    assert "combinatorial_substrate_gate" in w
    assert w["sample_selection_gate"]["required_before_biological_sample_selection"] is True
    assert "CS_GREEN" in w["sample_selection_gate"]["biological_use_open_if"]
    assert n["existing_evidence"]["focal_same_population_combinatorial_substrate"]=="not_yet_verified"
    assert n["focal_system_gate"]=="data/contracts/c_sieboldii_combinatorial_substrate_gate_v1.json"

    sie=[r for r in p if r["panel_id"]=="SIE"]
    assert len(sie)==1
    assert "CONDITIONAL" in sie[0]["geographic_design"]
    assert "Gate CS" in sie[0]["geographic_design"]
    assert "post-anthesis upright heads" in sie[0]["stop_rule"]

    for x in ("observation_id","population_id","local_set_id","flower_stage","colour_chroma",
              "colour_white_descriptive","orientation_deg_gravity","phyllary_rows","stickiness_score"):
        assert x in header, x

    for x in (
        "0° = vertically upward",
        "180° = vertically downward",
        "A1–A2 only",
        "3 of the 4",
        "CS_GREEN",
        "CS_GREEN_LOCAL_SET",
        "CS_RED_PHENOLOGY",
        "post-flowering upright heads",\n        "Priority 1 — Kurumayama Moor / Kirigamine, Nagano",\n        "Priority 2 — Kiyooka-Mukaiyama Wetland, Tsukude, Aichi",
        "do not sequence 40",
    ):
        need(d,x)

    need(readme,"C_SIEBOLDII_COMBINATORIAL_SUBSTRATE_GATE_V1.md")

    print(json.dumps({
        "status":"ok",
        "gate":"CS",
        "public_substrate_status":"unverified",
        "primary_pair":"colour_x_anthesis_orientation",
        "biological_WGS_selection":"blocked_until_CS_classified",
        "full_C_sieboldii_panel":"blocked"
    },indent=2))

if __name__=="__main__":
    main()
