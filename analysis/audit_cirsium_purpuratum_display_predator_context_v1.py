#!/usr/bin/env python3
"""Recompute SOURCE-PUBLISHED *Cirsium purpuratum* display/predation curves.

Ohashi & Yahara (2000) observed DIFFERENT populations in DIFFERENT years:
  Nikko 1993: 13 mothers, true damaged-seed fractions per mother.
  Kawamata 1994: 27 mothers, reconstructed predator losses from number of
    damaged heads using 451 florets/head and 52 seeds/larva from Nikko.
Their fitted P_damage(F)=1-exp(-a*F) coefficients are already published.

Ohashi & Yahara (1998) measured separate pollinator visit/display functions
with Latin squares. No original head-level seed/pollen physical interface
data or matched experimental fitness responses are in these two publications.

This is an honest numerical re-expression of AUTHOR-PUBLISHED curve
coefficients, NOT a refit to raw plants, confidence interval or species
effect of Cirsium true spine/orientation/guild identity.
"""
from __future__ import annotations
import argparse
import copy
import json
import math
from pathlib import Path

SOURCE=Path(__file__).resolve().parents[1]/(
    "data/evidence/cirsium_purpuratum_pollinator_predator_display_context_v1.json"
)
EXPECTED={
  "Nikko":{"year":1993,"n":13,"a":0.000063,"r_squared":0.44,
           "p":0.013,"larinus":0.44,"delia":0.56,
           "infested_pct":90.8,"mean_larvae":5.4,"max_heads":33,
           "assessment":"direct"},
  "Kawamata":{"year":1994,"n":27,"a":0.0000075,"r_squared":0.26,
           "p":0.007,"larinus":0.84,"delia":0.16,
           "infested_pct":31.3,"mean_larvae":0.26,"max_heads":65,
           "assessment":"inferred"},
}
DISPLAY_FLORETS=(451,1000,5000)
SOURCE_DOIS={"10.2307/2446309","10.1139/b99-182"}

def demand(condition,code):
    if not condition:raise ValueError(code)

def prediction(a,F):
    demand(isinstance(a,(int,float)) and math.isfinite(a) and a>0,
           "DISPLAY_PREDATOR_A_INVALID")
    demand(isinstance(F,(int,float)) and math.isfinite(F) and F>=0,
           "DISPLAY_FLORET_NUMBER_INVALID")
    return -math.expm1(-a*F)

def naive_optimum(a):
    # This is an author's MODEL endpoint under fixed S and no extra costs.
    demand(a>0 and math.isfinite(a),"INVALID_FITTED_DENSITY_PARAM")
    return 1/a

def check(source):
    demand(source["plant"]=="Cirsium purpuratum",
           "OTHER_PLANT_DATA_MASQUERADES_AS_PURPURATUM")
    poll=source["pollinator_behavior"]
    herb=source["herbivory"]
    demand(poll["doi"]=="10.2307/2446309"
           and herb["doi"]=="10.1139/b99-182",
           "POLLINATOR_OR_SEED_PREDATOR_PRIMARY_STUDY_REPLACED")
    demand("Latin square" in poll["design"],
           "POLLINATOR_SOURCE_DESIGN_CHANGED")
    demand(poll["pollination_outcomes_measured"].startswith(
        "visit rate; NOT independent"),
        "POLLINATOR_VISITS_FALSELY_DECLARED_EFFECTIVE_POLLEN")
    demand(herb["plant_scale_model"]==
           "P_damage(F)=1-exp(-a*F), F=total season florets per mother plant",
           "SEED_LOSS_MODEL_OR_UNIT_CHANGED")
    demand(set(herb["seed_feeder_species"])=={
        "Larinus meleagris","Delia platula"},
        "SEED_PREDATOR_GUILD_COMPOSITION_CHANGED")
    demand(len(herb["fit_sites"])==2,"UNSUPPORTED_THIRD_SAMPLE_GROUP")
    indexed={}
    for row in herb["fit_sites"]:
        site=row["name"]
        demand(site in EXPECTED and site not in indexed,
               "UNRECOGNIZED_OR_DUPLICATE_CIRSIUM_POPULATION")
        e=EXPECTED[site]
        for source_col,expected in (
            ("year",e["year"]),("n_plants_for_regression",e["n"]),
            ("a_display_damage",e["a"]),("r_squared",e["r_squared"]),
            ("reported_fit_p",e["p"]),
            ("seed_predator_Larinus_fraction",e["larinus"]),
            ("seed_predator_Delia_fraction",e["delia"]),
            ("infested_heads_pct_reported",e["infested_pct"]),
            ("mean_larvae_per_head",e["mean_larvae"]),
            ("observed_maximum_head_count_per_plant",e["max_heads"]),
        ):
            value=row[source_col]
            demand(isinstance(value,(int,float)) and
                   math.isclose(value,expected,abs_tol=1e-11,rel_tol=0),
                   "SOURCE_PUBLISHED_NUMBER_DRIFT_"+site+"_"+source_col)
        demand(abs(row["seed_predator_Larinus_fraction"]+
                   row["seed_predator_Delia_fraction"]-1)<1e-12,
               "SEED_FEEDER_SPECIES_COMPOSITION_SUM_INVALID")
        text=row["response"]
        if e["assessment"]=="direct":
            demand("directly assessed" in text and "INFERRED" not in text,
                   "NIKKO_TRUE_SEED_LOSS_FALSELY_RELABELED")
        else:
            demand("INFERRED" in text
                   and "451 florets/head" in text
                   and "52 seeds consumed/insect borrowed from Nikko" in text,
                   "KAWAMATA_RECONSTRUCTED_LOSS_FALSELY_DECLARED_DIRECT")
        indexed[site]=row
    demand(set(indexed)==set(EXPECTED),"SOURCE_SITE_MISSING")
    botanical=source["botanical_head_morphology_and_white_morph_context"]
    morph=botanical["white_flower_form"]
    demand(botanical["taxon"]=="Cirsium purpuratum"
           and botanical["japanese_name"]=="フジアザミ"
           and botanical["head_orientation"].startswith("normally downward")
           and morph["historical_name"]=="Cirsium purpuratum f. albiflorum (Kitam.) Kitam."
           and morph["status"]=="Kew synonym of Cirsium purpuratum, not independent evolutionary lineage",
           "PURPURATUM_BOTANICAL_IDENTITY_OR_FORMA_SYNOMYM_DRIFT")
    demand(not any(morph[k] for k in (
        "white_flower_occurence_voucher_collected_into_aza3",
        "within_same_population_white_purple_individuals_verified",
        "functional_anthocyanin_resource_retention_tested",
        "floral_RNA_and_pigment_pathway_measured",
        "white_flower_heritable_morph_origin_inferred",
        "reappearance_after_ancestral_loss_demonstrated"))
        and not botanical["original_Japan38_Comp1061_primary_taxon_sample_present"],
        "WHITE_FORMA_NAME_FALSELY_PROMOTED_TO_GENETIC_RETURN_OR_TREE_SAMPLE")
    demand(not botanical["not_a_single_species_proof_of_three_questions"] is False,
           "THREE_QUESTION_INFERENCE_SCOPE_EXPANDED_WITHOUT_SAMPLES")
    demand(not source["raw_individual_plant_records_loaded"] and
           not source["effect_of_changing_partner_identity_caused_by_randomization"]
           and not source["within_head_orientation_spine_trait_adaptive_fitness_effect_estimated"],
           "MODEL_WAS_FALSELY_PRESENTED_AS_PARTNER_CAUSAL_OR_SPINE_FITNESS")
    return indexed

def analyze(source):
    data=check(source)
    nic,kaw=data["Nikko"],data["Kawamata"]
    a_n=nic["a_display_damage"]
    a_k=kaw["a_display_damage"]
    sampled={}
    for site,row in sorted(data.items()):
        a=row["a_display_damage"]
        sampled[site]={
            "observed_population_year":row["year"],
            "n_real_mother_plants_used_by_original_fit":
                row["n_plants_for_regression"],
            "published_fit_coefficient_a_per_floret":a,
            "published_fit_r2":row["r_squared"],
            "published_fit_p":row["reported_fit_p"],
            "reported_Larinus_fraction":row["seed_predator_Larinus_fraction"],
            "reported_Delia_fraction":row["seed_predator_Delia_fraction"],
            "original_heads_with_detected_infestation_percent":
                row["infested_heads_pct_reported"],
            "original_mean_larvae_per_head":row["mean_larvae_per_head"],
            "original_max_heads_per_plant":row["observed_maximum_head_count_per_plant"],
            "damage_fraction_assessment":row["response"],
            "model_optimum_F_equal_1_over_a_rounded_coefficient":
                naive_optimum(a),
            "author_published_optimum_F_from_unrounded_fit":
                row["reported_f_star_florets_rounded"],
            "model_predictions_at_illustrative_F_NOT_OBSERVED": {
                str(f):prediction(a,f) for f in DISPLAY_FLORETS
            },
        }
    # Compare coefficients, NOT randomized contrasts or coefficient
    # confidence intervals, neither source level of replication permits.
    relation={
        "original_nikko_a_over_kawamata_a":a_n/a_k,
        "original_kawamata_optimum_F_over_nikko_optimum_F":
             (1/a_k)/(1/a_n),
        "published_site_infestation_percent_difference_nikko_minus_kawamata":
             nic["infested_heads_pct_reported"]-
             kaw["infested_heads_pct_reported"],
        "published_mean_larvae_per_head_ratio_nikko_over_kawamata":
             nic["mean_larvae_per_head"]/kaw["mean_larvae_per_head"],
        "seed_feeder_Larinus_fraction_kawamata_minus_nikko":
             kaw["seed_predator_Larinus_fraction"]-
             nic["seed_predator_Larinus_fraction"],
        "site_contrast_significance_test":"NOT_IDENTIFIABLE_FROM_TWO_SOURCE_FITS",
        "partner_causal_effect":"NOT_IDENTIFIABLE_YEAR_SITE_ASSESSMENT_CONFOUNDED",
    }
    demand(math.isclose(relation["original_nikko_a_over_kawamata_a"],8.4,
                        rel_tol=1e-12),
           "CONTRAST_ALGEBRA_NOT_REPRODUCIBLE")
    demand(math.isclose(relation["published_mean_larvae_per_head_ratio_nikko_over_kawamata"],
                        5.4/.26,rel_tol=1e-12),
           "ORIGINAL_LARVAE_INTENSITY_RATIO_MATH_ERROR")
    return {
      "status":"REAL_PUBLISHED_CIRSIUM_DISPLAY_PREDATION_CURVES_RECALCULATED_NOT_RAW_DATA",
      "original_study_2000_doi":source["herbivory"]["doi"],
      "pollinator_study_1998_doi":source["pollinator_behavior"]["doi"],
      "sources_share_species_but_not_study_individuals":True,
      "source_population_results":sampled,
      "numeric_context_differences_not_partner_randomized":relation,
      "max_claim":"SITE_AND_YEAR_CONTEXT_DIFFERENT_FITTED_SEED_PREDATION_DISPLAY_RESPONSE",
      "scientific_limitations":[
        "Nikko true seed-damage fraction; Kawamata seed-loss proportion reconstructed from insect counts under borrowed Nikko per-insect intake assumptions.",
        "Nikko 1993 and Kawamata 1994 are different sites/years/elevations, not within-study randomly exchanged guilds.",
        "Two different seed feeders occur at both sites in different percentages, but cause cannot be attributed to their identity from these data.",
        "1998 Latin-square bumblebee response and 2000 natural-predator study were not joint observations on the same individual heads.",
        "F=451,1000,5000 examples are arithmetic values of original published curves, not newly observed populations or data point estimates.",
        "F* derived from 1/a is an idealized per-season display optimum under original fixed S and no additional life-history tradeoffs; published rounded coefficients differ from author's unrounded fits.",
        "The floral display dimension is number of florets/heads, not measured physical head orientation or true botanical spine geometry.",
        "No confidence intervals for differences of population-specific a recovered; DO NOT infer treatment interaction p-value.",
      ],
      "causal_inference_denied":[
        "population_randomized_partner_composition",
        "same_head_pollinator_vs_seed_feeder_fitness_interaction",
        "genuine_spine_or_3D_angle_defense",
        "new_Cirsium_plant_individual_observations",
        "genetic_selection_or_EAzami_trait_memory",
      ],
    }

def negative_tests(source):
    source_report=analyze(source)
    demand(source_report["status"].startswith("REAL_PUBLISHED_CIRSIUM"),
           "REAL_PUBLISHED_SOURCE_MISSING")
    # Each fixture corrupts the source, never appended to real head ledgers.
    bad=[]
    x=copy.deepcopy(source)
    x["herbivory"]["fit_sites"][1]["response"]="directly assessed actual fraction"
    bad.append((x,"KAWAMATA_RECONSTRUCTED_LOSS_FALSELY_DECLARED_DIRECT"))
    x=copy.deepcopy(source)
    x["herbivory"]["fit_sites"][0]["a_display_damage"]=.000075
    bad.append((x,"SOURCE_PUBLISHED_NUMBER_DRIFT_Nikko_a_display_damage"))
    x=copy.deepcopy(source)
    x["raw_individual_plant_records_loaded"]=True
    bad.append((x,"MODEL_WAS_FALSELY_PRESENTED_AS_PARTNER_CAUSAL_OR_SPINE_FITNESS"))
    x=copy.deepcopy(source)
    x["effect_of_changing_partner_identity_caused_by_randomization"]=True
    bad.append((x,"MODEL_WAS_FALSELY_PRESENTED_AS_PARTNER_CAUSAL_OR_SPINE_FITNESS"))
    x=copy.deepcopy(source)
    x["herbivory"]["fit_sites"].append(copy.deepcopy(x["herbivory"]["fit_sites"][0]))
    bad.append((x,"UNSUPPORTED_THIRD_SAMPLE_GROUP"))
    x=copy.deepcopy(source)
    x["herbivory"]["fit_sites"][0]["seed_predator_Delia_fraction"]=.9
    bad.append((x,"SOURCE_PUBLISHED_NUMBER_DRIFT_Nikko_seed_predator_Delia_fraction"))
    x=copy.deepcopy(source)
    x["botanical_head_morphology_and_white_morph_context"]["white_flower_form"][
        "reappearance_after_ancestral_loss_demonstrated"]=True
    bad.append((x,"WHITE_FORMA_NAME_FALSELY_PROMOTED_TO_GENETIC_RETURN_OR_TREE_SAMPLE"))
    for wrong,error in bad:
        try:analyze(wrong)
        except ValueError as exc:
            demand(error in str(exc),"SCIENTIFIC_NEGATIVE_CONTROL_WRONG_REASON:"+str(exc))
        else:raise AssertionError("INVALID_SOURCE_FACT_PROMOTED_TO_NEW_CIRSIUM_EFFECT:"+error)
    demand(math.isclose(prediction(.000063,0),0,abs_tol=1e-15),
           "SATURATING_MODEL_INVALID_ORIGIN")
    demand(all(0<prediction(.000063,f)<1 for f in (100,1000,10000)),
           "SATURATING_CURVE_NOT_VALID_PROBABILITY")
    demand(math.isclose(naive_optimum(.000063),15873.015873015873,rel_tol=1e-12),
           "NAIVE_DISPLAY_FORMULA_WRONG")
    return {
      "status":"PASS_AUTHOR_PUBLISHED_CIRSIUM_CONTEXT_COEFFICIENTS_AND_CAUSAL_GATES",
      "source_corruption_cases_rejected":len(bad),
      "new_actual_Cirsium_head_field_rows":0,
      "true_author_published_population_fits_reexpressed":2,
      "raw_plants_reanalyzed":0,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,default=SOURCE)
    p.add_argument("--out",type=Path)
    x=p.parse_args()
    data=json.loads(x.manifest.read_text(encoding="utf-8"))
    receipt={
      "version":"cirsium_purpuratum_publish_model_separation_v1",
      "as_of":"2026-10-19",
      "original_published_parameter_audit":analyze(data),
      "synthetic_source_corruption_tests":negative_tests(data),
      "source_restrictions":"Published original population regression coefficients only; not the original per-plant data; different years; no partner factor experiment.",
    }
    if x.out:
        x.out.parent.mkdir(parents=True,exist_ok=True)
        x.out.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",
                         encoding="utf-8")
    a=receipt["original_published_parameter_audit"]["numeric_context_differences_not_partner_randomized"]
    print(json.dumps({
      "status":receipt["original_published_parameter_audit"]["status"],
      "coefficient_ratio_nikko_over_kawamata":a["original_nikko_a_over_kawamata_a"],
      "guild_Larinus_fraction_difference":a["seed_feeder_Larinus_fraction_kawamata_minus_nikko"],
      "raw_source_plants_reanalyzed":0,
      "guard":receipt["synthetic_source_corruption_tests"]["status"],
    },ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
