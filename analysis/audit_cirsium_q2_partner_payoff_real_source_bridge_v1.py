#!/usr/bin/env python3
"""Q2: Actual *Cirsium* partner-payoff evidence and missing same-head join.

All 4 sources are fetched at the EXACT original EAzami commit:
 c9870e91604ed1514344ddc8f765eae5237dbf27
- 9 seed-herbivory contrasts nested in 4 experiment articles (2 taxa)
- 6 pollination/reproductive assurance articles, ONE numeric bagged/open
- 15 earlier literature-map evidence entries + 3 appended later cases

What cannot be inferred:
* Herbivory suppression RR 2.67 across two taxa and C. erisithales
  open/bagged RR 3.92 are NOT matched same-species/payoff estimates,
  and cannot be compared as a selection tradeoff on one trait.
* 1983 C. discolor sticky effect positive, 1983 C. flodmani null,
  2003 C. discolor null: effect SIGN findings differ by site/time/species;
  these are distinct reports, NOT a factorial partner-switch experiment.
* 2007 C. discolor sticky guild contact routing is observed natural
  history, not measured fitness rescue by stickiness.
* Current original EAzami head/observation ledgers contain 0 field rows
  when pinned; PR #44 head-camera entry contracts are UNMERGED from PR50.

This script produces a source/data evidence map with strict original-grain
validation; it does not fake a combined meta-effect or new raw observations.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import re
from collections import Counter,defaultdict
from pathlib import Path

HERB_COLS={"study_cluster","effect_id","taxon","reduced_herbivory_mean",
    "reduced_herbivory_se","ambient_herbivory_mean","ambient_herbivory_se",
    "source_doi","source_locator"}
PA_COLS={"study_id","taxon_or_scope","source","design","pollinator_dependence",
    "pollen_limitation_open_conditions","exact_numeric_recovered",
    "quantitative_anchor"}
INT_COLS={"evidence_id","study_id","taxon","interaction_domain","partner_guild",
    "capitulum_module","design","response_chain","fitness_endpoint",
    "direction","doi"}
ORIGINAL_EAZAMI_REF="c9870e91604ed1514344ddc8f765eae5237dbf27"

def require(x,code):
    if not x: raise ValueError(code)

def read(path,cols):
    with open(path,encoding="utf-8-sig",newline="") as f:
        z=csv.DictReader(f)
        require(cols.issubset(set(z.fieldnames or [])),
                "EAzami_ORIGINAL_SOURCE_SCHEMA_DRIFT")
        return list(z)

def index_unique(rows,key,code):
    out={}
    for row in rows:
        val=row.get(key,"").strip()
        require(val and val not in out,code)
        out[val]=row
    return out

def parse_mean(anchor,term):
    match=re.search(r"\b"+re.escape(term)+r"\s+([0-9]+(?:\.[0-9]+)?)",anchor)
    require(match is not None,"SOURCE_NUMERIC_MEAN_NOT_RESOLVED_"+term)
    return float(match.group(1))

def analyze(herb,meta,pa,pa_meta,evidence,ext):
    require(len(herb)==9 and len(evidence)==15 and len(ext)==3
            and len(pa)==6,"Q2_SOURCE_EVIDENCE_ROW_COUNTS_CHANGED")
    hx=index_unique(herb,"effect_id","DUPLICATE_SOURCE_HERBIVORY_EFFECT")
    px=index_unique(pa,"study_id","DUPLICATE_POLLINATOR_ASSURANCE_STUDY")
    ex=index_unique(evidence+ext,"evidence_id",
                    "LITERATURE_EVIDENCE_ID_DOUBLE_COUNTED")
    require(meta["coverage"]["effect_rows"]==9
            and meta["coverage"]["independent_study_clusters"]==4
            and pa_meta["independent_study_count"]==6
            and pa_meta["exact_numeric_studies"]==1,
            "FROZEN_META_STUDY_GRAIN_MISMATCH")
    groups=defaultdict(list)
    for row in herb:
        groups[row["study_cluster"]].append(row)
        reduced=float(row["reduced_herbivory_mean"])
        ambient=float(row["ambient_herbivory_mean"])
        require(math.isfinite(reduced) and math.isfinite(ambient)
                and reduced>0 and ambient>0,
                "INVALID_SOURCE_TREATMENT_MEANS")
        effect=next((e for e in meta["effect_level"]
                     if e["effect_id"]==row["effect_id"]),None)
        require(effect is not None
                and effect["study"]==row["study_cluster"]
                and abs(effect["lnRR"]-math.log(reduced/ambient))<1e-10
                and abs(effect["response_ratio"]-reduced/ambient)<1e-10,
                "PUBLISHED_HERBIVORY_EFFECT_MATH_NOT_REPRODUCED")
    require(len(groups)==4
            and {x for r in herb for x in (r["taxon"],)}==
                {"Cirsium canescens","Cirsium occidentale var. occidentale"},
            "INDEPENDENT_4_STUDY_2_TAXON_SCOPE_DRIFT")
    require({x["study"] for x in meta["study_level"]}==set(groups),
            "SOURCE_META_DATA_GENERATION_STUDY_ORPHAN")
    reff=meta["random_effects"]
    require(reff["k"]==4
            and abs(math.exp(reff["pooled_lnRR"])-
                    reff["response_ratio"])<1e-10,
            "META_RESPONSE_RATIO_TRANSFORMATION_INVALID")
    require(reff["ci95_response_ratio"][0]<
            reff["response_ratio"]<reff["ci95_response_ratio"][1],
            "META_95CI_NOT_CONTAINING_ESTIMATE")
    # The EAzami source estimates already incorporate across-study weights.
    # Do NOT re-estimate variances or imply causal independent taxon replicates.

    exact_pollinator=[row for row in pa
                      if row["exact_numeric_recovered"]=="1"]
    require(len(exact_pollinator)==1
            and exact_pollinator[0]["study_id"]=="PA02"
            and exact_pollinator[0]["taxon_or_scope"]=="Cirsium erisithales",
            "EXACT_POLLINATION_ASSAY_SOURCE_CHANGED")
    pol=exact_pollinator[0]
    supplement=parse_mean(pol["quantitative_anchor"],"supplement")
    opened=parse_mean(pol["quantitative_anchor"],"control")
    bagged=parse_mean(pol["quantitative_anchor"],"bagged")
    require(abs(opened-47.88)<1e-9
            and abs(supplement-41.75)<1e-9
            and abs(bagged-12.2)<1e-9,
            "FROZEN_NUMERIC_POLLEN_TREATMENTS_CHANGED")
    require(pol["pollen_limitation_open_conditions"]=="absent"
            and pol["pollinator_dependence"]=="high",
            "POLLINATOR_DEPENDENCE_WRONGLY_CALLED_LIMITATION")

    by_species=defaultdict(list)
    for row in evidence+ext:by_species[row["taxon"]].append(row)
    # 1983 paper contains TWO species contrasts; do not count those as
    # two independent article-level studies, nor as one same-species contrast.
    target={name:ex[name] for name in
           ["INT006","INT016","INT017","INT018","INT001","INT002","INT015"]}
    require(target["INT006"]["taxon"]=="Cirsium discolor"
            and target["INT006"]["direction"]=="null"
            and target["INT006"]["capitulum_module"]=="stickiness"
            and target["INT016"]["taxon"]=="Cirsium discolor"
            and target["INT016"]["direction"]=="sticky_defence_benefit"
            and target["INT017"]["taxon"]=="Cirsium flodmani"
            and target["INT017"]["direction"]=="null"
            and target["INT018"]["taxon"]=="Cirsium discolor"
            and target["INT018"]["direction"]=="guild_selective",
            "TRUE_BRacteal_STICKINESS_CONTRADICTORY_CONTEXT_MISSING")
    require(all(target[k]["capitulum_module"]=="display_size"
                for k in ("INT001","INT002"))
            and target["INT015"]["capitulum_module"]=="flower_production"
            and all(target[k]["taxon"]=="Cirsium purpuratum"
                    for k in ("INT001","INT002","INT015")),
            "JAPANESE_CIRSIUM_DISPLAY_BOTH_GUILD_CONTEXT_NOT_PRESERVED")
    # Claim ranking not study heterogeneity tests:
    herb_taxa=set(row["taxon"] for row in herb)
    exact_poll_taxa={pol["taxon_or_scope"]}
    require(not herb_taxa & exact_poll_taxa,
            "ALREADY_MATCHED_HERBIVORE_AND_POLLINATOR_SAME_TAXON_DATA")
    # Critically, the six PA studies include several Cirsium congeners,
    # but only PA02 has exact compatible treatment mean values.
    # Do not infer that no pollination studies exist for herb_taxa.
    return {
      "version":"Cirsium_Q2_partner_payoff_real_source_coverage_v1",
      "provenance":{"repo":"zuizui0223/EAzami",
                    "commit":ORIGINAL_EAZAMI_REF,
                    "original_bounded_map_rows":len(evidence),
                    "independently_appended_source_rows":len(ext),
                    "unique_evidence_ids":len(ex),
                    "bounded_not_exhaustive_literature_review":True},
      "direct_herbivore_experimental_seed_cost":{
        "actual_published_effect_contrasts":len(hx),
        "independent_source_studies":len(groups),
        "taxa":sorted(herb_taxa),
        "experimental_seed_suppression_pooled_RR_existing_meta":
             reff["response_ratio"],
        "existing_meta_95CI":reff["ci95_response_ratio"],
        "existing_mixed_study_method_not_refitted_in_this_crosswalk":True,
        "not_a_true_spine_or_orientation_manipulation":True},
      "pollination_service":{
        "independent_studies_with_mixed_estimands":len(px),
        "recovered_exact_pollen_3arm_experiments":1,
        "numeric_source_taxon":pol["taxon_or_scope"],
        "mean_seed_output_supplement":supplement,
        "mean_seed_output_open":opened,
        "mean_seed_output_bagged":bagged,
        "open_vs_bagged_seed_ratio":opened/bagged,
        "supplement_vs_open_seed_ratio":supplement/opened,
        "same_study_open_pollen_limitation_detected":False,
        "bagging_changes_microclimate_and_visitation_not_a_trait_selection_experiment":True,
        "no_across_study_pollinator_RR_pooled":True},
      "key_actual_sticky_head_interventions":{
        "literature_rows":["INT016","INT017","INT006"],
        "distinct_experimental_articles":2,
        "independent_species_test_contexts":3,
        "species_positive_or_null_sign_evidence":{
          "Cirsium discolor":{
            "Willson_1983_Illinois_positive":"higher seed predation when sticky exudates occluded",
            "Thomas_2003_Illinois_null":"no increase seed predators or decrease seed output when traps neutralized",
            "natural_guild_observation_Thomas_2007":"pollinators do not touch traps; some seed predators bypass; ants and aphids deterred",
          },
          "Cirsium flodmani":{
            "Willson_1983_Minnesota_null":"occlusion without detectable head-seed-predation effect in 1980/1982"}
        },
        "species_environment_protocol_confounding_unresolved":True,
        "fitness_effect_size_incompatible_missing":True,
        "proven_sign_heterogeneity_does_not_prove_an_interacting_guild_caused_difference":True},
      "Japanese_same_species_separate_agent_sources":{
        "taxon":"Cirsium purpuratum",
        "pollinator_evidence_ids":["INT001","INT002"],
        "seed_predator_evidence_id":"INT015",
        "pollination_exposure_trait":"floral display size / flower density",
        "antagonist_trait":"seasonal flower production",
        "same_head_same_plant_both_fitness_outcomes_measured":False,
        "not_the_same_physical_intervention_axis":True},
      "joint_question_identifiability":{
        "explicit_taxa_with_exact_quantitative_antagonist_RR_and_exact_pollinator_3arm_seed_means":0,
        "independent_same_head_populations_with_both_guilds_and_true_spine_or_angle_manipulation":0,
        "trait_specific_partner_crossed_fitness_potential_outcomes_identified":False,
        "partner_by_trait_selection_gradient_identified":False,
        "population_or_lineage_fitness_gain_reversal_identified":False,
        "strongest_existing_real_biological_result":"head-level stickiness defensive benefit can be positive or undetectable in different source field contexts; pollinator/seed predator paths are not automatically the same",
        "still_unknown":"whether relative profitability switches with ACTUAL partner identity rather than site/species/time or trait state",
        "next_field_observation":"Independent same-stage insect taxon identification and genuine armature contact on prespecified Cirsium heads, followed by egg/pollen and viable-seed fate",
        "first_manipulation_if_permitted":"sham-controlled single capitulum axis, orthogonally varied interaction-agent exposure within one population; true head enrollment before outcomes",
        "full_Q2_causal_contrast":"W(z, pollinator, enemy) under crossed partner contexts; observe real pollen and seed endpoint; estimate interaction terms",
      },
      "CLAIM":"REAL_SOURCE_EVIDENCE_MAP_AND_CROSS_EXPERIMENTS_NOT_NEW_SAME_HEAD_PARTNER_PAYOFF_EFFECT",
    }

def negative_controls(herb,meta,pa,pa_meta,base,ext):
    import copy
    actual=analyze(herb,meta,pa,pa_meta,base,ext)
    assert actual["joint_question_identifiability"][
        "explicit_taxa_with_exact_quantitative_antagonist_RR_and_exact_pollinator_3arm_seed_means"]==0
    manipulations=[
      (lambda hh,mm,pp,pmm,bb,ee:hh.append(copy.deepcopy(hh[0])),
       "DUPLICATE_SOURCE_HERBIVORY_EFFECT"),
      (lambda hh,mm,pp,pmm,bb,ee:ee.append(copy.deepcopy(ee[0])),
       "LITERATURE_EVIDENCE_ID_DOUBLE_COUNTED"),
      (lambda hh,mm,pp,pmm,bb,ee:bb.append(copy.deepcopy(ee[0])),
       "LITERATURE_EVIDENCE_ID_DOUBLE_COUNTED"),
      (lambda hh,mm,pp,pmm,bb,ee:hh[0].update(
           ambient_herbivory_mean="-1"),
       "INVALID_SOURCE_TREATMENT_MEANS"),
      (lambda hh,mm,pp,pmm,bb,ee:pp[1].update(
           pollinator_dependence="low"),
       "POLLINATOR_DEPENDENCE_WRONGLY_CALLED_LIMITATION"),
      (lambda hh,mm,pp,pmm,bb,ee:next(
           r for r in ee if r["evidence_id"]=="INT016").update(direction="null"),
       "TRUE_BRacteal_STICKINESS_CONTRADICTORY_CONTEXT_MISSING"),
    ]
    for change,needle in manipulations:
        hh,mm,pp,pmm,bb,ee=copy.deepcopy((
            herb,meta,pa,pa_meta,base,ext))
        change(hh,mm,pp,pmm,bb,ee)
        try:analyze(hh,mm,pp,pmm,bb,ee)
        except ValueError as x:
            require(needle in str(x),"INCORRECT_Q2_NEGATIVE_TEST_REASON:"+str(x))
        else:raise AssertionError("INVALID_Q2_SOURCE_EVIDENCE_ACCEPTED:"+needle)
    return {
        "status":"PASS_SOURCE_BASED_Q2_PARTNER_EVIDENCE_AND_NO_FALSE_CROSSED_FITNESS",
        "corruption_tests_rejected":len(manipulations),
        "real_source_qualitative_sticky_trials":3,
        "matched_Cirsium_same_head_partner_payoff_records_created":0}

def main():
    p=argparse.ArgumentParser()
    for key in ("herbivory","herbivory-meta","pollination",
                "pollination-meta","interaction-base","interaction-extension"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    herb=read(a.herbivory,HERB_COLS)
    meta=json.loads(a.herbivory_meta.read_text(encoding="utf-8"))
    poll=read(a.pollination,PA_COLS)
    poll_meta=json.loads(a.pollination_meta.read_text(encoding="utf-8"))
    base=read(a.interaction_base,INT_COLS)
    ext=read(a.interaction_extension,INT_COLS)
    out={
        "as_of":"2026-10-19",
        "real_source_evidence":analyze(herb,meta,poll,poll_meta,base,ext),
        "synthetic_source_corruption_tests":negative_controls(
            herb,meta,poll,poll_meta,base,ext),
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",
                     encoding="utf-8")
    print(json.dumps({
        "status":out["real_source_evidence"]["CLAIM"],
        "herbivory_experiment_studies":out["real_source_evidence"][
            "direct_herbivore_experimental_seed_cost"]["independent_source_studies"],
        "pollination_numeric_experiments":out["real_source_evidence"][
            "pollination_service"]["recovered_exact_pollen_3arm_experiments"],
        "bracteal_sticky_source_studies":out["real_source_evidence"][
            "key_actual_sticky_head_interventions"]["distinct_experimental_articles"],
        "crossed_partner_fitness_data":out["real_source_evidence"][
            "joint_question_identifiability"]["trait_specific_partner_crossed_fitness_potential_outcomes_identified"],
        "source_integrity":out["synthetic_source_corruption_tests"]["status"],
    },indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
