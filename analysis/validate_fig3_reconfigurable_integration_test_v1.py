#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"/"FIG3_RECONFIGURABLE_INTEGRATION_TEST_V1.md"
CON=ROOT/"data"/"contracts"/"aza3_fig3_reconfigurable_integration_test_v1.json"
TPL=ROOT/"data"/"templates"/"fig3_heldout_prediction_registry_v1.csv"
NATURE=ROOT/"data"/"contracts"/"aza3_nature_single_question_v2.json"
CS=ROOT/"data"/"contracts"/"c_sieboldii_combinatorial_substrate_gate_v1.json"
WGS=ROOT/"data"/"contracts"/"aza3_r1b_own_wgs_transferability_pilot_v1.json"

def need(t,x):
    if x not in t:
        raise AssertionError(f"missing: {x}")

def main():
    d=DOC.read_text(encoding="utf-8")
    c=json.loads(CON.read_text(encoding="utf-8"))
    nature=json.loads(NATURE.read_text(encoding="utf-8"))
    cs=json.loads(CS.read_text(encoding="utf-8"))
    wgs=json.loads(WGS.read_text(encoding="utf-8"))
    header=next(csv.reader(TPL.open(encoding="utf-8")))

    assert c["status"]=="PROSPECTIVE__LOCKED_BEFORE_FOCAL_GENOMIC_DISCOVERY"
    assert c["question"]==nature["one_question"]
    assert c["prerequisites"]["gate_cs_allowed"]==["CS_GREEN","CS_GREEN_LOCAL_SET"]
    assert c["prerequisites"]["reference_gate"]=="OWN_WGS_GREEN"
    assert c["prerequisites"]["primary_taxonomic_confidence"]=="high"
    assert c["sample_design"]["technical_pilot_n"]==8
    assert c["sample_design"]["training_target_min"]==40
    assert c["sample_design"]["untouched_validation_min"]==20
    assert c["sample_design"]["strongest_route_min_total_usable"]==60
    assert c["sample_design"]["legacy_40_status"]=="DISCOVERY_ONLY_UNLESS_EXPANDED_FOR_UNTOUCHED_VALIDATION"
    assert c["models"]["F"]["label"]=="fixed_evolutionary_integration"
    assert c["models"]["R"]["label"]=="reconfigurable_integration"
    assert c["primary_estimand"]["id"]=="E1_UNTOUCHED_JOINT_PREDICTION"
    assert c["kill_shot_estimand"]["id"]=="E3_PROSPECTIVE_MOSAIC_PREDICTION"
    assert c["kill_shot_estimand"]["validation_mosaic_min"]==3
    assert set(c["decisions"])=={"FIG3_R_STRONG","FIG3_R_PARTIAL","FIG3_F_SUPPORTIVE","NOT_IDENTIFIABLE"}

    assert cs["primary_pair"]==["floral_colour","anthesis_orientation"]
    assert wgs["sample_selection_gate"]["required_before_biological_sample_selection"] is True

    required={
        "individual_id","validation_unit","split","phenotype_labels_opened",
        "model_F_joint_log_score","model_R_joint_log_score","d_i",
        "colour_F_log_score","colour_R_log_score",
        "orientation_F_log_score","orientation_R_log_score",
        "colour_genomic_state","orientation_genomic_state",
        "prospective_mosaic_registered","colour_prediction_correct",
        "orientation_prediction_correct","reference_stable",
        "taxonomic_confidence","technical_exclusion_reason"
    }
    assert required.issubset(set(header))

    for x in (
        "The decisive test is prediction, not merely non-overlapping GWAS peaks.",
        "at least 60 usable individuals",
        "at least 20 individuals",
        "rank-1",
        "Delta_joint",
        "sign-flip",
        "at least 3 high-confidence prospective mosaic",
        "FIG3_R_STRONG",
        "FIG3_F_SUPPORTIVE",
        "Do not reinterpret NOT_IDENTIFIABLE as support for Model F.",
    ):
        need(d,x)

    print(json.dumps({
        "status":"ok",
        "question":c["question"],
        "primary_estimand":c["primary_estimand"]["id"],
        "kill_shot":c["kill_shot_estimand"]["id"],
        "minimum_training":c["sample_design"]["training_target_min"],
        "minimum_untouched_validation":c["sample_design"]["untouched_validation_min"],
        "legacy_40":"discovery_only_without_expansion"
    },indent=2))

if __name__=="__main__":
    main()
