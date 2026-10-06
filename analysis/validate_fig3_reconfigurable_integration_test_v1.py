#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"/"FIG3_RECONFIGURABLE_INTEGRATION_TEST_V1.md"
CON=ROOT/"data"/"contracts"/"aza3_fig3_reconfigurable_integration_test_v1.json"
PRED=ROOT/"data"/"templates"/"fig3_preunblinding_prediction_registry_v1.csv"
UNLOCK=ROOT/"data"/"templates"/"fig3_validation_unlock_receipt_v1.json"
SCORE=ROOT/"data"/"templates"/"fig3_postunblinding_score_registry_v1.csv"
NATURE=ROOT/"data"/"contracts"/"aza3_nature_single_question_v2.json"
CS=ROOT/"data"/"contracts"/"c_sieboldii_combinatorial_substrate_gate_v1.json"
WGS=ROOT/"data"/"contracts"/"aza3_r1b_own_wgs_transferability_pilot_v1.json"

def need(t,x):
    if x not in t: raise AssertionError(f"missing: {x}")

def main():
    d=DOC.read_text(encoding="utf-8")
    c=json.loads(CON.read_text(encoding="utf-8"))
    nature=json.loads(NATURE.read_text(encoding="utf-8"))
    cs=json.loads(CS.read_text(encoding="utf-8"))
    wgs=json.loads(WGS.read_text(encoding="utf-8"))
    ph=next(csv.reader(PRED.open(encoding="utf-8")))
    sh=next(csv.reader(SCORE.open(encoding="utf-8")))
    unlock=json.loads(UNLOCK.read_text(encoding="utf-8"))

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
    assert "multi-trait modules" in c["primary_estimand"]["analysis_scope"]
    assert c["kill_shot_estimand"]["correct_validation_mosaic_min"]==3
    assert c["validation_blinding"]["scoring_rule"].startswith("Observed validation phenotypes are joined only after")
    assert c["model_fairness"]["candidate_block_pool"].startswith("identical TRAIN-derived")\n    assert c["module_baselines"]["M_local"]["role"].startswith("sensitivity baseline")\n    assert "non-negative held-out advantage over M_local" in c["module_baselines"]["strong_R_rule"]
    assert "same complexity budget" in c["model_fairness"]["R_effect_structure"]
    assert cs["primary_pair"]==["floral_colour","anthesis_orientation"]
    assert wgs["sample_selection_gate"]["required_before_biological_sample_selection"] is True

    pred_required={
      "individual_id","validation_unit","split","prediction_registry_frozen_at","prediction_model_version",
      "model_F_prediction_object","model_M_prediction_object","model_Mlocal_prediction_object","model_R_prediction_object",
      "colour_genomic_state","orientation_genomic_state","prospective_mosaic_registered",
      "reference_stable","taxonomic_confidence","technical_exclusion_reason",
      "admitted_construct_count","admitted_multi_trait_module_count"
    }
    score_required={
      "individual_id","validation_unit","prediction_registry_sha256","observed_phenotype_source",
      "model_F_joint_log_score","model_M_joint_log_score","model_Mlocal_joint_log_score","model_R_joint_log_score",
      "d_RM_i","d_MF_i","prospective_mosaic_registered","colour_prediction_correct",
      "orientation_prediction_correct","reference_stable","taxonomic_confidence","technical_exclusion_reason"
    }
    assert pred_required.issubset(set(ph))
    assert score_required.issubset(set(sh))
    assert unlock["status"]=="TEMPLATE__UNFILLED"
    assert unlock["prediction_registry_sha256"] is None

    for x in (
      "Nature-scale support requires R to beat M, not merely F.",
      "at least two predeclared multi-trait modules",
      "at least five admitted constructs total",
      "cannot by itself establish FIG3_R_STRONG",
      "Minimum strong-route total after QC: 60.",
      "same TRAIN-derived candidate genomic-block pool",
      "same effective parameter/regularization budget",
      "E1 scoring scope",
      "single-trait presentation module is excluded",
      "FIG3_MODULE_INHERITANCE_SUPPORTIVE",
      "FIG3_WHOLE_ORGAN_SUPPORTIVE",
      "Do not reinterpret NOT_IDENTIFIABLE as support for F or M."
    ): need(d,x)

    print(json.dumps({
      "status":"ok","question":c["question"],
      "primary_comparison":"R_vs_frozen_phenotypic_module_inheritance",
      "preunblinding_registry":"required_and_hash_locked",
      "minimum_training":c["sample_design"]["training_target_min"],
      "minimum_untouched_validation":c["sample_design"]["untouched_validation_min"]
    },indent=2))

if __name__=="__main__": main()
