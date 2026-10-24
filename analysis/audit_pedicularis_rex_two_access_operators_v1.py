#!/usr/bin/env python3
"""Pedicularis rex: two ORIGINAL published ecological interfaces within one species.

Audits the 2015 *functional intervention* and the 2016 *corolla-exsertion
observational selection* at their correct and distinct experimental grains.

Sun & Huang 2015 DOI 10.1093/aobpla/plv019: water draining done by a scissors
hole (water plus cut), visit study Shama one population with 10 randomized
plots/year; seed outcome study six geographic populations. Note pollen was NOT
experimentally quantified for the 2015 visitation comparison.
Sun, Armbruster & Huang 2016 DOI 10.1093/aob/mcw097: naturally varying corolla
exsertion; stigma pollen and seed predation both increase with exsertion in
published separate fitted models; morphological joins retained for only SEVEN
populations (14 morphology/pollen, 12 seed outcomes). No causal exsertion
manipulation or water x exsertion full factorial experiment.

Source data are PUBLISHED TABLE COEFFICIENTS, not original individual records.
The xls supplement's URL is known but original raw individual records have
NOT been authenticated / read. This validator never pretends to fit new GLMs.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

EXPECTED_SITES={
    "Baishuitai","Sanba","Zhongdian","Deqin","Daxueshan","Shama"
}
EXPECTED_EFFECTS={
    "legitimate_bumblebee_visitation_Shama_only":.958,
    "nectar_robber_visitation_Shama_only":.951,
    "initial_seed_set_six_sites":.906,
}
MAX_AUTHENTIC_REPORT_P=0.0001

def need(ok,why):
    if not ok:
        raise ValueError(why)

def conservative_bh(p_bounds):
    """BH q using UPPER BOUNDS of reported 'p<0.0001', not fake p=0.
    Larger p-values cannot improve a q value; the result is a conservative
    upper bound on the real BH-adjusted p for a fixed six-test family.
    """
    need(len(p_bounds)==6 and all(0<=p<=1 for p in p_bounds),
         "INVALID_OR_INCOMPLETE_ORIGINAL_SIX_LOCATION_P_VALUES")
    sorted_index=sorted(range(len(p_bounds)),key=lambda i:(p_bounds[i],i))
    qs=[None]*6
    bound=1.0
    for n in range(5,-1,-1):
        i=sorted_index[n]
        bound=min(bound,p_bounds[i]*6/(n+1))
        qs[i]=bound
    return qs

def audit(evidence):
    need(evidence["version"]=="pedicularis_rex_two_operators_source_table_v1",
         "WRONG_REVIEWED_PRIMARY_EVIDENCE")
    s=evidence["sources"]
    a=s["drainage_2015"]
    b=s["exsertion_2016"]
    need(a["doi"]=="10.1093/aobpla/plv019" and
         a["pmcid"]=="PMC4392828" and
         b["doi"]=="10.1093/aob/mcw097" and
         b["pmcid"]=="PMC4970362",
         "TWO_PRIMARY_PAPERS_MISSING")
    need(a["seed_populations"]==6 and
         a["pollinator_video_populations"]==1 and
         a["pollinator_video_population"]=="Shama" and
         a["pollinator_video_population"] in EXPECTED_SITES,
         "FALSE_SIX_POPULATION_POLLINATOR_EFFECT")
    need(a["water_only_sham_performed"] if "water_only_sham_performed" in a
         else s.get("nonsense_test") is None,
         "UNEXPECTED_SHAM_FIELD")
    need("not sham-equivalent" in a["treatment"] or
         "not sham-equivalent" in a.get("experimental_units","") or
         evidence["synthesis"]["study_2015_water_cut_not_sham_water_only"],
         "SCISSOR_CUT_CONFOUNDED_WATER_ISOLATION_NOT_DISCLOSED")
    endpoint={x["label"]:x for x in a["tables_1_2"]}
    need(len(endpoint)==5 and all(k in endpoint for k in EXPECTED_EFFECTS),
         "MISSING_2015_EXPERIMENTAL_OUTCOMES")
    for label,p in EXPECTED_EFFECTS.items():
        need(math.isclose(endpoint[label]["p_exact"],p,
                          rel_tol=0,abs_tol=1e-12),
             "PUBLISHED_2015_NO_DETECTED_EFFECT_CELL_CHANGED")
    for label,chisq in (("final_seed_set_six_sites",16.710),
                       ("seed_predation_six_sites",92.808)):
        t=endpoint[label]
        need(math.isclose(t["chi2"],chisq,
                          rel_tol=0,abs_tol=1e-10) and
             t["p_exact"] is None and
             t["p_less_than"]==MAX_AUTHENTIC_REPORT_P,
             "PUBLISHED_2015_FINAL_OR_ENEMY_SUCCESS_CELL_CHANGED")
    sites=a["site_specific_seed_predation"]
    need(len(sites)==6 and {x["name"] for x in sites}==EXPECTED_SITES,
         "SOURCE_SIX_ORIGINAL_PREDATION_LOCATIONS_CHANGED")
    sorted_sites=sorted(sites,key=lambda x:x["name"])
    values=[x["p_upper"] for x in sorted_sites]
    expected_bh=conservative_bh(values)
    qmap={}
    for row,q in zip(sorted_sites,expected_bh):
        name=row["name"]
        need(math.isclose(row["bh_q_conservative_upper"],q,
                          rel_tol=0,abs_tol=1e-12),
             "ORIGINAL_SIX_SITE_BH_RECEIPT_DRIFT")
        need(row["passes_BH_conservative"]==(q<.05),
             "PUBLISHED_BH_SITE_CLASSIFICATION_CHANGED")
        need(row["p_upper"]==(row["p_exact"] if row["p_exact"] is not None
                              else row["p_lt"]),
             "SOURCE_P_LT_SYMBOL_WAS_FALSIFIED")
        qmap[name]={
            "p_exact":row["p_exact"],"p_reported_less_than":row["p_lt"],
            "BH_q_conservative_upper":q,
            "passes_conservative_BH_at_0p05":q<.05,
        }
    need(math.isclose(qmap["Deqin"]["BH_q_conservative_upper"],
                       .0552,rel_tol=0,abs_tol=1e-12),
         "DEQIN_P_0P046_IS_NOT_Q_0P05")
    need(qmap["Zhongdian"]["p_exact"]==.789,
         "NONRESPONDING_ZHONGDIAN_P_CHANGED")
    need(sum(x["passes_conservative_BH_at_0p05"]
             for x in qmap.values())==4,
         "PUBLISHED_2015_FOUR_CONSERVATIVE_SITE_DISCOVERIES_CHANGED")
    need(a["nominal_significant_sites"]==5 and
         a["conservative_BH_significant_sites"]==4,
         "AUTHOR_FIVE_NOMINAL_VS_FOUR_BH_SOURCE_MISMATCH")
    need(b["n_trait_populations"]==14 and
         b["n_pollen_plants"]==299 and
         b["n_pollen_flowers"]==598 and
         b["n_seed_populations"]==12 and
         b["n_linked_trait_predation_populations"]==7,
         "SOURCE_2016_DIFFERENT_POPULATION_DENOMINATORS_CHANGED")
    need(b["linked_populations"]==[1,3,5,8,9,10,11] and
         b["lost_morphology_labels_populations"]==[2,4,6,7,12],
         "FIVE_2016_SOURCE_MISSING_MORPHOLOGY_JOINS_NOT_PRESERVED")
    need("flower length -" in b["exsertion"] and
         "bract height" in b["exsertion"],
         "COMPOSITE_EXSERTION_NOT_A_SINGLE_GENETIC_MODULE")
    models={m["endpoint"]:m for m in b["source_table_1_models"]}
    need(set(models)=={"pollen_arrival","predispersal_seed_predation",
                       "initial_seed_set","final_viable_seed"},
         "2016_PUBLISHED_FOUR_STAGE_MODEL_LIST_CHANGED")
    need(models["pollen_arrival"]["exsertion_sign"]=="positive" and
         models["pollen_arrival"]["exsertion_p_less_than"]==.01 and
         models["predispersal_seed_predation"]["exsertion_sign"]=="positive" and
         models["predispersal_seed_predation"]["exsertion_p_less_than"]==.001,
         "2016_TRUE_CONFLICTING_EXSERTION_SELECTION_NOT_ANCHORED")
    need(models["final_viable_seed"]["best_model_has_population_interactions"] is True,
         "2016_FINAL_FITNESS_MOSAIC_INTERACTION_REMOVED")
    syn=evidence["synthesis"]
    need(syn["same_species_different_trials"] and
         not syn["not_same_individual_factorial_water_x_exsertion"] is False
         and syn["study_2016_original_raw_individual_reanalysis_not_performed"],
         "FALSE_SAME_HEAD_FACTORIAL_CAUSAL_LOSS")
    return {
        "version":"pedicularis_rex_q2_two_access_operators_audit_v1",
        "source_papers":{
            "wet_bract_intervention":a["doi"],
            "floral_exsertion_observational_selection":b["doi"],
        },
        "2015":{
            "nominal_seed_predation_significant_sites":5,
            "conservative_six_site_BH_significant_sites":4,
            "per_site_p_and_q_upper":qmap,
            "legitimate_pollinator_Shama_treatment_p":
                endpoint["legitimate_bumblebee_visitation_Shama_only"]["p_exact"],
            "nectar_robber_Shama_treatment_p":
                endpoint["nectar_robber_visitation_Shama_only"]["p_exact"],
            "initial_seed_set_treatment_p":
                endpoint["initial_seed_set_six_sites"]["p_exact"],
            "final_seed_set_treatment_p_upper":
                endpoint["final_seed_set_six_sites"]["p_less_than"],
            "seed_predation_treatment_p_upper":
                endpoint["seed_predation_six_sites"]["p_less_than"],
        },
        "2016":{
            "trait_pollen_populations":b["n_trait_populations"],
            "individual_pollen_plants":b["n_pollen_plants"],
            "individual_pollen_flowers":b["n_pollen_flowers"],
            "linked_predation_trait_populations":b["n_linked_trait_predation_populations"],
            "linked_population_codes":b["linked_populations"],
            "model_best_AICc":{m["endpoint"]:m["best_AICc"]
                               for m in b["source_table_1_models"]},
            "corolla_exsertion_positively_associated_with":{
                "stigmatic_pollen_receipt":True,
                "predispersal_seed_predation":True,
            },
        },
        "biology":{
            "operator_water_retention_plus_cut_changes_seed_enemies":
                True,
            "operator_exsertion_exposes_flower_to_mutualists_and_enemies":
                True,
            "not_jointly_manipulated_or_measured_in_same_organism":True,
            "separate_partner_access_axes_are_plausible_not_identified_genetic_modules":
                True,
            "new_Pedicularis_raw_plant_rows_analyzed":0,
            "new_Cirsium_field_events":0,
        },
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    data=json.loads(args.source.read_text(encoding="utf-8"))
    result=audit(data)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,
                                    allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
