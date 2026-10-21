#!/usr/bin/env python3
"""Source-verified published summaries: Cirsium pitcheri head hierarchy and 2 enemies.

Gijsman, Havens & Vitt 2020 Global Ecology and Conservation e00945
doi:10.1016/j.gecco.2020.e00945
The numbers are printed rounded descriptive proportions and original-author
model contrasts, not self-recovered plant, seed or insect individual records.

Head rank (terminal/secondary/tertiary) is not head *orientation*, and apparent
deer browsing of missing terminal capitula is the authors' HYPOTHESIS, not an
identified cause. Within-WFD weevil damages and between-year head absence must
not be pooled into a causal rank-by-enemy or heritable risk-spreading payoff.

The 2021 Dryad file is catalogued by an original exact SHA but download API
was 401 and protected HTML blocked binary access at 2026-10-21; it is NOT read.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path


def require(ok,code):
    if not ok:raise ValueError(code)


def numeric(x,name,lower=-math.inf,upper=math.inf):
    require(type(x) in (int,float) and math.isfinite(x) and lower<=x<=upper,
            "INVALID_NUMERIC_"+name)
    return float(x)


def analyze(s):
    require(s["source"]["doi"]=="10.1016/j.gecco.2020.e00945"
            and s["source"]["study_years"]==[2016,2017]
            and s["source"]["primary_published_article_transcription_not_raw_head_rows"] is True,
            "PUBLISHED_SOURCE_GRAIN_DRIFT")
    w=s["rank_weevil_2016_WFD"]
    f=w["infestation_evidence_fractions_approx"]
    require(list(f)==["terminal","secondary","tertiary"],
            "SOURCE_HEAD_RANKS_CHANGED")
    for key,expected in (("terminal",.27),("secondary",.50),("tertiary",.80)):
        require(abs(numeric(f[key],key,0,1)-expected)<1e-12,
                "ORIGINAL_WEEVIL_RANK_ESTIMATE_CHANGED")
    require(w["rank_causal_weevil_assignment"] is False,
            "WEEVIL_CONTRAST_NOT_RANDOMIZED")
    require(abs(w["mature_model_head_position_F"]-16.8)<1e-9
            and abs(w["mature_model_weevil_F"]-112.8)<1e-9,
            "ORIGINAL_AUTHOR_WEEVIL_GLM_DRIFT")
    require(abs(w["mature_seed_infested_to_uninfested_approx_ratio"]-.4)<1e-12
            and abs(w["unfilled_seed_infested_to_uninfested_approx_ratio"]-.6)<1e-12,
            "ORIGINAL_WEevil_SEED_COUNT_RATIOS_DRIFT")
    weight=s["rank_specific_seed_mass_2016_WFD"]
    terminal=weight["terminal"]
    secondary=weight["secondary"]
    require(terminal["infested_minus_uninfested_seed_mass_mg"]==-1.83
            and terminal["CI95_infested_minus_uninfested_mg"]==[-2.73,-.93],
            "ORIGINAL_TERMINAL_SEED_MASS_CONTRAST_CHANGED")
    require(secondary["infested_minus_uninfested_seed_mass_mg"] is None,
            "SECONDARY_SEED_WEIGHT_ESTIMATE_NOT_REPORTED")
    require(secondary["CI95_infested_minus_uninfested_mg"]==[-.76,.061]
            and abs(secondary["paper_p"]-.095)<1e-12
            and not weight["cross_rank_weevil_effect_difference_statistically_tested"],
            "SECONDARY_MASS_INTERVAL_OR_INTERACTION_DRIFT")
    yearly=s["field_seed_output_by_rank_year"]
    require(yearly["terminal_2017_mature_count_relative_to_2016_approx"]==.3
            and yearly["terminal_2017_reduction_observed_at_both_sites"] is True
            and yearly["tertiary_2017_mature_seeds_both_sites"]==0
            and abs(yearly["secondary_year_difference_p"]-.24)<1e-12,
            "ORIGINAL_2017_SEED_HEAD_OUTPUT_CONTEXT_DRIFT")
    br=s["terminal_head_absence_by_site_year"]
    require(br["specific_herbivore_causation_identified"] is False,
            "DEER_CAUSE_MUST_BE_UNCERTAIN")
    require(br["sampling_2017_parent_plants_per_site"]==15,
            "YEAR2017_PARENT_N_CHANGED")
    expected={"2016":{"SC":0,"WFD":.089},"2017":{"SC":.27,"WFD":.33}}
    for year in ("2016","2017"):
        for site in ("SC","WFD"):
            require(abs(numeric(br[year][site],"terminal_absence_"+year+"_"+site,0,1)
                        -expected[year][site])<1e-12,
                    "MISSING_TERMINAL_HEAD_ORIGINAL_SOURCE_DRIFT")
    site=s["site_level_mature_seed_attempted_2016"]
    require(abs(site["SC_no_weevils"]["mature_fraction_approx"]-.38)<1e-12
            and abs(site["WFD_weevils_present"]["mature_fraction_approx"]-.71)<1e-12
            and site["comparison_is_different_plant_level_estimand_than_infested_head_40percent_ratio"],
            "CROSS_SITE_SEED_DENOMINATOR_CHANGED")
    req=s["material_limit"]
    require(req["2020_and_2021_seeds_not_joined_by_original_head_id"],
            "DIFFERENT_STUDIES_HAVE_NO_SHARED_HEAD_IDS")
    require(req["terminal_missing_head_not_direct_deer_observation"]
            and req["no_valid_causal_rank_by_enemy_interaction"]
            and req["no_plant_level_total_net_fitness_estimated"]
            and req["rank_is_branch_position_not_head_inclination"]
            and req["no_measured_spines_phyllary_angle_or_mucilage"],
            "UNSUPPORTED_CAUSAL_HEAD_MODULE_CLAIM")
    require(req["new_field_observations_added"]==0
            and req["original_plant_and_head_rows_not_downloaded"],
            "FIELD_SAMPLE_COUNT_FABRICATED")
    dryad=s["dryad_blocked_separately"]
    require(dryad["api_returned_http"]==401
            and dryad["no_verified_raw_xlsx"] is True
            and dryad["individual_germination_file_id"]==483373,
            "DO_NOT_CLAIM_UNAVAILABLE_SEED_LEVEL_DATA")
    return {
      "version":"cirsium_pitcheri_2020_competing_risk_head_rank_published_reanalysis_v1",
      "source_doi":s["source"]["doi"],
      "source_is_published_rounded_data_not_raw_plant_heads":True,
      "within_WFD_2016_weevil_head_rank":{
        "fraction_of_heads_with_weevil_evidence":f,
        "secondary_vs_terminal_infestation_ratio":f["secondary"]/f["terminal"],
        "tertiary_vs_terminal_infestation_ratio":f["tertiary"]/f["terminal"],
        "rank_level_sample_sizes_and_confidence_intervals_unknown":True,
        "competing_effect_of_head_age_or_diameter_unresolved":True,
      },
      "WFD_2016_head_seed_quantity":{
        "infested_to_uninfested_mature_seed_count_ratio":.4,
        "infested_to_uninfested_unfilled_seed_count_ratio":.6,
        "author_mature_count_GLM_rank_F":16.8,
        "author_mature_count_GLM_weevil_F":112.8,
        "same_head_infestation_is_observed_not_randomized":True,
      },
      "WFD_2016_head_seed_weight":{
        "terminal_infested_minus_uninfested_mg":-1.83,
        "terminal_CI95_infested_minus_uninfested_mg":[-2.73,-.93],
        "secondary_CI95_infested_minus_uninfested_mg":[-.76,.061],
        "secondary_CI95_contains_zero":True,
        "secondary_exact_effect_point_estimate_missing":True,
        "difference_of_terminal_and_secondary_effects_significance_known":False,
      },
      "terminal_head_loss_cross_year_observations":{
        name:{
            "fraction_plants_missing_terminal_2016":br["2016"][name],
            "fraction_plants_missing_terminal_2017":br["2017"][name],
            "difference_in_missing_terminal_plant_fraction":
                br["2017"][name]-br["2016"][name],
            "original_cause_is_verified_deer_attack":False,
            "site_year_and_flowering_plant_cohorts_differ":True,
        } for name in ("SC","WFD")
      },
      "field_year_head_output":{
        "terminal_mature_seed_count_relative_to_2016_approx":.3,
        "secondary_year_comparison_original_p":.24,
        "secondary_year_effect_not_detected":True,
        "absence_of_detectable_secondary_year_difference_implies_causal_buffering":False,
        "tertiary_mature_seeds_reported_at_both_sites_in_2017":0,
      },
      "site_fraction_of_attempted_seeds_mature_2016":{
        "SC_weevil_free_site":.38,
        "WFD_weevil_present_site":.71,
        "WFD_minus_SC":.33,
        "WFD_site_difference_is_positive_effect_of_weevils":False,
        "do_not_confuse_whole_plant_fraction_with_head_level_count_ratio":True,
      },
      "inference_guards":{
        "cannot_compare_weevil_ranks_with_deer_caused_loss_as_one_same_cohort":True,
        "head_rank_is_not_phyllary_spine_or_gravity_angle":True,
        "risk_diversification_is_prospective_not_established_adaptation":True,
        "year_vs_enemy_observation_not_crossed_factorial":True,
        "no_fitted_same_head_seed_germination_with_weevil":True,
        "original_2021_Dryad_seed_file_access_still_blocked":True,
      },
      "prospective_test":"Pre-enrolled same head with nested parent plant and flowering time, rank, true phyllary/spine geometry; monitor independent large-herbivore browsing, adult weevil approaches/oviposition, actual mature seeds and source-linked germination before testing adaptive risk diversification.",
      "n_new_individual_head_field_observations":0,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    d=analyze(json.loads(a.source.read_text(encoding="utf-8")))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
       "status":"SOURCE_BOUNDED_ORIGINAL_2020_HEAD_RANK_AND_TWO_POTENTIAL_ENEMIES",
       "infestation_rank_relatives":d["within_WFD_2016_weevil_head_rank"],
       "missing_terminal_delta":d["terminal_head_loss_cross_year_observations"],
       "seed_mass_genuine_head_association":d["WFD_2016_head_seed_weight"],
       "new_field_samples":0,
    },indent=2))

if __name__=="__main__":main()
