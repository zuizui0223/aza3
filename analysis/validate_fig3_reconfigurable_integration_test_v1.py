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

    assert c["contract_version"]=="aza3_fig3_reconfigurable_integration_test_v2"
    assert c["status"]=="PROSPECTIVE__LOCKED_BEFORE_FOCAL_GENOMIC_DISCOVERY"
    assert c["question"]==nature["one_question"]

    mods=c["frozen_phenotypic_modules"]
    assert mods["presentation"]==["presentation_angle"]
    assert mods["colour"]==["floral_lightness","floral_chroma","floral_hue"]
    assert mods["head_form"]==["head_elongation","head_compactness"]
    assert mods["involucre_armature"]==["involucre_form","projection_prominence","projection_pattern"]

    adm=c["construct_admission"]
    assert adm["strongest_route_multi_trait_modules_min"]==2
    assert adm["constructs_per_admitted_multi_trait_module_min"]==2
    assert adm["admitted_constructs_total_min"]==5

    assert c["prerequisites"]["gate_cs_allowed"]==["CS_GREEN","CS_GREEN_LOCAL_SET"]
    assert c["prerequisites"]["reference_gate"]=="OWN_WGS_GREEN"
    assert c["sample_design"]["technical_pilot_n"]==8
    assert c["sample_design"]["training_target_min"]==40
    assert c["sample_design"]["untouched_validation_min"]==20
    assert c["sample_design"]["strongest_route_min_total_usable"]==60
    assert c["sample_design"]["legacy_40_status"]=="DISCOVERY_ONLY_UNLESS_EXPANDED_FOR_UNTOUCHED_VALIDATION"

    assert c["models"]["F"]["label"]=="whole_organ_fixed_integration"
    assert c["models"]["M"]["label"]=="frozen_phenotypic_module_inheritance"
    assert c["models"]["R"]["label"]=="reconfigurable_integration"
    assert c["primary_estimand"]["id"]=="E1_R_VS_MODULE_HELDOUT_PREDICTION"
    assert c["baseline_estimand"]["id"]=="E1B_MODULE_VS_WHOLE_ORGAN"
    assert c["kill_shot_estimand"]["id"]=="E3_PROSPECTIVE_COLOUR_ORIENTATION_MOSAICS"
    assert c["kill_shot_estimand"]["correct_validation_mosaic_min"]==3
    assert c["kill_shot_estimand"]["preregistration_provenance"]==[
        "mosaic_registry_commit_sha","mosaic_registry_frozen_at_utc","mosaic_prediction_hash"
    ]
    assert set(c["decisions"])=={
        "FIG3_R_STRONG","FIG3_R_PARTIAL","FIG3_MODULE_INHERITANCE_SUPPORTIVE",
        "FIG3_WHOLE_ORGAN_SUPPORTIVE","NOT_IDENTIFIABLE"
    }

    assert cs["primary_pair"]==["floral_colour","anthesis_orientation"]
    assert wgs["sample_selection_gate"]["required_before_biological_sample_selection"] is True

    required={
        "individual_id","validation_unit","split","phenotype_labels_opened",
        "construct_admission_frozen","admitted_construct_count",
        "admitted_multi_trait_module_count",
        "model_F_joint_log_score","model_M_joint_log_score","model_R_joint_log_score",
        "d_MF","d_RM",
        "module_colour_M_log_score","module_colour_R_log_score",
        "module_head_form_M_log_score","module_head_form_R_log_score",
        "module_involucre_M_log_score","module_involucre_R_log_score",
        "colour_genomic_state","orientation_genomic_state",
        "prospective_mosaic_registered","mosaic_registry_commit_sha",
        "mosaic_registry_frozen_at_utc","mosaic_prediction_hash","colour_prediction_correct",
        "orientation_prediction_correct","reference_stable",
        "taxonomic_confidence","technical_exclusion_reason"
    }
    assert required.issubset(set(header))

    for x in (
        "Nature-scale support requires R to beat M, not merely F.",
        "at least two predeclared multi-trait modules",
        "at least five admitted constructs total",
        "colour × orientation mosaicism",
        "cannot by itself establish FIG3_R_STRONG",
        "Minimum strong-route total after QC: 60.",
        "Delta_RM",
        "FIG3_MODULE_INHERITANCE_SUPPORTIVE",
        "FIG3_WHOLE_ORGAN_SUPPORTIVE",
        "Do not reinterpret NOT_IDENTIFIABLE as support for F or M.",
        "same TRAIN-derived candidate genomic-block pool",
        "same effective parameter/regularization budget",
        "design floor, not a power-derived target",
    ):
        need(d,x)

    print(json.dumps({
        "status":"ok",
        "question":c["question"],
        "primary_comparison":"R_vs_frozen_phenotypic_module_inheritance",
        "baseline_comparison":"module_vs_whole_organ",
        "min_constructs":adm["admitted_constructs_total_min"],
        "min_multi_trait_modules":adm["strongest_route_multi_trait_modules_min"],
        "minimum_training":c["sample_design"]["training_target_min"],
        "minimum_untouched_validation":c["sample_design"]["untouched_validation_min"],
        "legacy_40":"discovery_only_without_expansion"
    },indent=2))

if __name__=="__main__":
    main()
