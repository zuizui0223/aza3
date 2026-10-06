#!/usr/bin/env python3
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"/"AZA3_NATURE_SINGLE_QUESTION_V2.md"
CON=ROOT/"data"/"contracts"/"aza3_nature_single_question_v2.json"
FIG1=ROOT/"data"/"evidence"/"aza3_fig1_ecological_module_alignment_result_v1.json"
R1=ROOT/"data"/"contracts"/"aza3_reference_public_data_readiness_v1.json"
WGS=ROOT/"data"/"contracts"/"aza3_r1b_own_wgs_transferability_pilot_v1.json"
README=ROOT/"README.md"

def need(t,x):
    if x not in t:
        raise AssertionError(f"missing: {x}")

def main():
    d=DOC.read_text(encoding="utf-8")
    c=json.loads(CON.read_text(encoding="utf-8"))
    f=json.loads(FIG1.read_text(encoding="utf-8"))
    r=json.loads(R1.read_text(encoding="utf-8"))
    w=json.loads(WGS.read_text(encoding="utf-8"))
    readme=README.read_text(encoding="utf-8")

    assert c["status"]=="ACTIVE_CONCEPTUAL_TARGET__GENOMIC_ANSWER_PENDING"
    assert c["one_question"]=="Must an integrated complex organ evolve as an integrated unit?"
    assert c["working_label"]=="reconfigurable_integration"
    assert set(c["competing_models"])=={
        "whole_organ_fixed_integration",
        "phenotypic_module_inheritance",
        "reconfigurable_integration",
    }
    assert c["fig3_test_contract"]=="data/contracts/aza3_fig3_reconfigurable_integration_test_v1.json"
    assert "R_must_outpredict_frozen_phenotypic_module_model_for_FIG3_R_STRONG" in c["claim_boundaries"]
    assert "colour_orientation_separability_alone_is_not_reconfigurable_integration" in c["claim_boundaries"]

    assert f["decision"]=="FIXED_ECOLOGICAL_MODULE_ALIGNMENT_NOT_SUPPORTED"
    assert f["test"]["among_taxon"]["exact_one_sided_p"] > 0.05
    assert f["test"]["within_taxon"]["exact_one_sided_p"] > 0.05
    assert f["contrast_with_phenotypic_modularity"]["within_taxon_module_cohesion"]["permutation_p_one_sided"] < 0.05
    assert f["contrast_with_phenotypic_modularity"]["among_taxon_module_cohesion"]["permutation_p_one_sided"] < 0.05

    assert r["current_decision"]=="PROCEED_TO_OWN_WGS_TRANSFERABILITY_PILOT"
    assert r["reference_build_decision"]=="DO_NOT_BUILD_FOCAL_C_SIEBOLDII_REFERENCE_YET"
    assert w["hard_stop"]=="Do not sequence the full 298-individual Phase-A panel, and do not choose C. sieboldii Nature-test individuals, before Gate CS and this transferability pilot are classified."

    for x in (
        "Must an integrated complex organ evolve as an integrated unit?",
        "reconfigurable integration",
        "FIXED_ECOLOGICAL_MODULE_ALIGNMENT_NOT_SUPPORTED",
        "integration itself need not be the unit of inheritance",
        "do not search for another grouping that makes it significant",
        "cross-level mismatch itself is unprecedented",
        "C_SIEBOLDII_COMBINATORIAL_SUBSTRATE_GATE_V1.md",
        "The Nature-scale comparison is **R versus M**, not merely R versus F.",
        "at least five admitted constructs",
        "at least 20 untouched validation individuals",
        "Colour × anthesis orientation",
        "cannot establish FIG3_R_STRONG",
    ):
        need(d,x)

    # The doc must explicitly deny the strongest novelty overclaim.
    need(d,"must not claim that cross-level mismatch itself is unprecedented")

    need(readme,"AZA3_NATURE_SINGLE_QUESTION_V2.md")
    need(readme,"aza3_fig1_ecological_module_alignment_result_v1.json")

    print(json.dumps({
        "status":"ok",
        "nature_question":c["one_question"],
        "fig1_ecological_module_alignment":"not_supported",
        "phenotypic_module_cohesion":"supported",
        "genomic_answer":"pending",
        "next_gate":"Gate_CS_then_8_individual_own_WGS_transferability_pilot_then_F_M_R_prediction",
        "fig3_primary_comparison":"R_vs_frozen_phenotypic_module_M"
    },indent=2))

if __name__=="__main__":
    main()
