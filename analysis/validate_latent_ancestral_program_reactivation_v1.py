#!/usr/bin/env python3
from pathlib import Path
import csv,json

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"/"LATENT_ANCESTRAL_PROGRAM_REACTIVATION_V1.md"
HIST=ROOT/"docs"/"REACTIVATION_HISTORY_GATE_V1.md"
CON=ROOT/"data"/"contracts"/"aza3_latent_ancestral_program_reactivation_v1.json"
HCON=ROOT/"data"/"contracts"/"aza3_reactivation_history_gate_v1.json"
HAND=ROOT/"data"/"evidence"/"eazami_latent_pathway_handoff_v1.json"
LIT=ROOT/"data"/"evidence"/"latent_program_reactivation_literature_map_v1.csv"

def need(t,x):
    if x not in t: raise AssertionError(f"missing: {x}")

def main():
    d=DOC.read_text(encoding="utf-8")
    h=HIST.read_text(encoding="utf-8")
    c=json.loads(CON.read_text(encoding="utf-8"))
    hc=json.loads(HCON.read_text(encoding="utf-8"))
    hand=json.loads(HAND.read_text(encoding="utf-8"))
    lit=list(csv.DictReader(LIT.open(encoding="utf-8")))

    assert c["status"]=="MECHANISM_HYPOTHESIS__DOWNSTREAM_OF_FMR"
    assert "LATENT_ANCESTRAL_PROGRAM_REACTIVATION" in c["mechanism_classes"]
    assert c["reactivation_required_layers"]==[
      "ancestral_state_supported",
      "machinery_persists_through_loss",
      "loss_state_has_regulatory_or_developmental_suppression",
      "recurrent_state_restores_homologous_pathway_in_same_tissue_or_stage"
    ]
    assert c["nature_hierarchy"]["headline"].startswith("module_boundary_evolvability")
    assert hc["status"]=="REQUIRED_BEFORE_REGAIN_LANGUAGE"
    assert hc["current_system_status"]["takaoense"]=="H_REGAIN_COMPATIBLE"
    assert hc["current_system_status"]["arenicola"]=="H_NOT_IDENTIFIABLE"
    assert hand["source_status"]=="EAzami_hypothesis_not_demonstrated"
    assert hand["historical_identifiability"]["takaoense_six_morph_linked_samples"]["rooted_binary_resolutions"]==945
    assert hand["historical_identifiability"]["takaoense_six_morph_linked_samples"]["topologies_regain_required"]==270
    assert hand["historical_identifiability"]["takaoense_six_morph_linked_samples"]["topologies_no_regain_optimum_allowed"]==675
    systems={r["system"] for r in lit}
    assert {"Iochrominae_Solanaceae","Cypriniform_fishes","Trait_reversal_review_2026"}.issubset(systems)

    for x in (
      "Only route 3 is the strict suppressed-pathway re-expression hypothesis.",
      "It is that a pathway can remain evolutionarily available while being phenotypically silent.",
      "REACTIVATION_STRONG",
      "standing variation",
      "introgress",
      "non-homologous",
      "does not replace the F/M/R inheritance-unit test"
    ): need(d,x)

    for x in (
      "H_REGAIN_SUPPORTED",
      "H_REGAIN_COMPATIBLE",
      "H_RETICULATE_REACQUISITION",
      "Historical classification is frozen before opening focal molecular reactivation results.",
      "The Nature paper does not require reversal."
    ): need(h,x)

    print(json.dumps({
      "status":"ok",
      "mechanism":"latent_ancestral_program_reactivation",
      "takaoense_history":"regain_compatible_not_demonstrated",
      "arenicola_history":"not_identifiable",
      "nature_role":"downstream_mechanism_not_headline"
    },indent=2))

if __name__=="__main__": main()
