#!/usr/bin/env python3
"""Cirsium pitcheri capitulum rank -> weevil evidence -> seed germination.

Independent, manually rechecked PUBLISHED NUMERIC SUMMARIES:
 Gijsman et al. (2020) 10.1016/j.gecco.2020.e00945, rank-specific
    approximate field infestation and observed mature seed association.
 Gijsman & Vitt (2021) 10.1002/ece3.7109, Table 1 (2016 seed origin),
    Table 2 (2017 seed origin), per site/rank/temp FITTED germination.

The two papers cannot be head-ID joined; germination was performed on
selected viable seeds. The author summary rows are NOT original individual
seeds or an intervention changing head rank or structural phyllaries.
No p-value or standard error for differences is invented from marginal SEs.
"""
from __future__ import annotations
import argparse
from collections import Counter
import copy
import hashlib
import json
import math
from pathlib import Path

SITES=("SC","WFD")
POSITIONS=("terminal","secondary")
DAY_TEMPS=(20,25,30)
YEARS=(2016,2017)
WE_PRESENT=("terminal","secondary","tertiary")
KEYS={"2016":"round3","2017":"total_fraction"}

def require(cond,msg):
    if not cond:raise ValueError(msg)

def number(x,name):
    require(isinstance(x,(float,int)) and not isinstance(x,bool) and math.isfinite(x), "INVALID_NUMERIC_"+name)
    return float(x)

def validate(records, strict=True):
    require(len(records)==24 if strict else bool(records),"EXPECTED_24_FITTED_SOURCE_ROWS")
    seen=set()
    for row in records:
        for col in ("collection_year","site","head_rank","day_temp_c","night_temp_c"):
            require(col in row,"SOURCE_CELL_COLUMN_MISSING_"+col)
        year=row["collection_year"]
        site=row["site"]
        pos=row["head_rank"]
        temp=row["day_temp_c"]
        require(year in YEARS and site in SITES and pos in POSITIONS and temp in DAY_TEMPS
                and row["night_temp_c"]==10,"NONORIGINAL_SITE_YEAR_RANK_TEMP")
        key=(year,site,pos,temp)
        require(key not in seen,"DUPLICATE_FITTED_SOURCE_CELL")
        seen.add(key)
        if year==2016:
            require(row["study_table"]=="Table 1" and
                    row["measurement"]=="fitted_cumulative_germination_three_cold_rounds",
                    "2016_SOURCE_ROUND3_DEFINITION_DRIFT")
            p=[]
            for k in ("round1","round2","round3"):
                require(isinstance(row.get(k),dict),"MISSING_COLD_ROUND")
                x=number(row[k].get("fraction"),k+"_p")
                se=number(row[k].get("se"),k+"_se")
                require(0<=x<=1 and 0<se<1,"INVALID_GERMINATION_FRACTION_OR_SE")
                p.append(x)
            require(p[0]<=p[1]+1e-10 and p[1]<=p[2]+1e-10,
                    "CUMULATIVE_GERMINATION_DECREASES_BETWEEN_ROUNDS")
            require("total_fraction" not in row,"2016_MIXED_2017_MEASUREMENT")
        else:
            require(row["study_table"]=="Table 2" and
                    row["measurement"]=="fitted_total_germination_after_three_rounds",
                    "2017_MODEL_ROW_DEFINITION_DRIFT")
            p=number(row.get("total_fraction"),"2017_p")
            se=number(row.get("total_se"),"2017_se")
            days=number(row.get("mean_days_to_germinate"),"mean_days")
            days_se=number(row.get("days_se"),"days_se")
            require(0<=p<=1 and 0<se<1 and days>0 and days_se>0,
                    "INVALID_2017_GERMINATION_MODEL_RESULT")
            require("round3" not in row,"2017_MIXED_2016_MEASUREMENT")
    if strict:
        expected={(y,s,p,t) for y in YEARS for s in SITES
                  for p in POSITIONS for t in DAY_TEMPS}
        require(seen==expected,"SOURCE_FACTORIAL_GRID_DRIFT")
    return {"n_source_model_cells":len(records),
            "n_per_collection_year":dict(Counter(z["collection_year"] for z in records)),
            "n_per_head_position":dict(Counter(z["head_rank"] for z in records)),
            "full_site_rank_temp_factorial":True if strict else None}


def germination_value(row):
    return (row["round3"]["fraction"] if row["collection_year"]==2016
            else row["total_fraction"])


def rank_contrasts(rows):
    by={(r["collection_year"],r["site"],r["head_rank"],r["day_temp_c"]):r
        for r in rows}
    output=[]
    for year in YEARS:
        for site in SITES:
            for day in DAY_TEMPS:
                a=by[(year,site,"terminal",day)]
                b=by[(year,site,"secondary",day)]
                ap=germination_value(a)
                bp=germination_value(b)
                require(bp>0,"MISSING_SECONDARY_DENOMINATOR")
                contrast={
                   "collection_year":year,"site":site,"day_temp_c":day,
                   "terminal_model_fitted_germination":ap,
                   "secondary_model_fitted_germination":bp,
                   "terminal_minus_secondary":ap-bp,
                   "terminal_divided_by_secondary":ap/bp,
                   "is_terminal_point_estimate_higher":ap>bp,
                   "is_secondary_point_estimate_higher":ap<bp,
                   "source_table":"Table 1 final cumulative round3" if year==2016
                                   else "Table 2 total model",
                   "independently_tested_rank_difference_from_original_se":False,
                }
                if year==2016:
                    contrast["by_cold_round_terminal_minus_secondary"]={
                        k:a[k]["fraction"]-b[k]["fraction"]
                        for k in ("round1","round2","round3")}
                else:
                    contrast["days_to_germination_terminal_minus_secondary"] = (
                        a["mean_days_to_germinate"]-b["mean_days_to_germinate"])
                output.append(contrast)
    return output


def describe_weather_rank(rows):
    output=rank_contrasts(rows)
    by={}
    for yr in YEARS:
        for site in SITES:
            g=[r for r in output if r["collection_year"]==yr and r["site"]==site]
            require(len(g)==3,"SOURCE_SITE_YEAR_TEMP_COVERAGE_DRIFT")
            by[f"{yr}_{site}"]={
                "n_temperature_conditions":len(g),
                "terminal_higher_point_estimates":sum(r["terminal_minus_secondary"]>0 for r in g),
                "secondary_higher_point_estimates":sum(r["terminal_minus_secondary"]<0 for r in g),
                "temperature_specific_rank_difference":{
                    str(r["day_temp_c"]):r["terminal_minus_secondary"] for r in g},
                "mean_difference_across_three_lab_temperatures_descriptive_only":
                    sum(r["terminal_minus_secondary"] for r in g)/3,
            }
    site_year_diff={}
    for site in SITES:
        x=by[f"2016_{site}"]["temperature_specific_rank_difference"]
        y=by[f"2017_{site}"]["temperature_specific_rank_difference"]
        site_year_diff[site]={
            temp:y[temp]-x[temp] for temp in ("20","25","30")
        }
    return {
        "by_site_collection_year":by,
        "between_study_year_method_point_contrast_changes_not_statistical_interactions":
            site_year_diff,
        "n_2016_terminal_better_out_of_six":
            sum(r["terminal_minus_secondary"]>0 for r in output if r["collection_year"]==2016),
        "n_2017_terminal_better_out_of_six":
            sum(r["terminal_minus_secondary"]>0 for r in output if r["collection_year"]==2017),
        "germination_sc20_2016_terminal_secondary_ratio":
            next(r["terminal_divided_by_secondary"] for r in output
                 if r["collection_year"]==2016 and r["site"]=="SC" and r["day_temp_c"]==20),
        "germination_wfd20_2017_terminal_secondary_ratio":
            next(r["terminal_divided_by_secondary"] for r in output
                 if r["collection_year"]==2017 and r["site"]=="WFD" and r["day_temp_c"]==20),
        "paired_model_rows":output,
        "not_an_inference_of_significant_between_year_site_by_rank_interaction":True,
    }


def weevil_summary(source):
    observed=source["rank_infestation_approximate_fraction"]
    require(set(observed)==set(WE_PRESENT),"WEEVIL_RANKS_NOT_TERM_SEC_TERT")
    for key in WE_PRESENT:
        p=number(observed[key],"weevil_"+key)
        require(0<p<1,"INVALID_APPROX_PUBLISHED_WEEVIL_PREVALENCE")
    if abs(observed["terminal"]-.27)>1e-12 or abs(observed["secondary"]-.50)>1e-12 or abs(observed["tertiary"]-.80)>1e-12:
        raise ValueError("APPROX_ORIGINAL_WEEVIL_RANK_VALUES_CHANGED")
    seed_ratio=source["mature_seed_count_infested_relative_to_uninfested"]
    require(abs(seed_ratio-.40)<1e-12,"PUBLISHED_MATURE_SEED_RATIO_DRIFT")
    return {
       "site":source["study_site_for_rank_infestation"],
       "study_context":"2016 WFD observed Larinus planus infestation evidence",
       "rank_infestation_approximate_fraction":observed,
       "secondary_to_terminal_approx_infested_fraction_ratio":observed["secondary"]/observed["terminal"],
       "tertiary_to_terminal_approx_infested_fraction_ratio":observed["tertiary"]/observed["terminal"],
       "mature_seeds_infested_vs_uninfested_observed_association_ratio":seed_ratio,
       "observational_head_rank_comparison_not_insect_randomization":True,
       "relative_risk_is_from_rounded_published_percentages_not_raw_binomial_estimator":True,
       "infestation_ratio_cannot_be_multipled_by_viable_seed_ratio_to_infer_per_rank_fitness":True,
       "rank_is_not_head_orientation_or_phyllary_armature":True,
    }


def assert_provenance(source):
    require(source["taxon"]=="Cirsium pitcheri","SOURCE_TAXON_ID_CHANGED")
    require(source["source2020"]["doi"]=="10.1016/j.gecco.2020.e00945" and
            source["source2021"]["doi"]=="10.1002/ece3.7109",
            "PUBLISHED_SOURCE_DOI_DRIFT")
    require(source["source2021"]["status"]==
        "MANUALLY_TRANSCRIBED_PUBLISHED_MODEL_FITTED_MEANS_AND_SE_NOT_ORIGINAL_INDIVIDUAL_ROWS",
        "DO_NOT_CLAIM_SEED_LEVEL_INPUT")
    require(source["source2020"]["status"]==
        "PUBLISHED_APPROXIMATE_POINT_SUMMARIES_NOT_ORIGINAL_HEAD_LEVEL_RAW_CSV",
        "DO_NOT_CLAIM_HEAD_LEVEL_DATA")
    c=source["study_population_sample_info"]
    require(c["source_collection_2016"]==
       {"SC_parent_plants":34,"WFD_parent_plants":43,"SC_germination_seeds":1217,"WFD_germination_seeds":1078},
       "STUDY_2016_ORIGINAL_COHORT_COUNT_CHANGED")
    require(c["source_collection_2017"]==
       {"SC_parent_plants":15,"WFD_parent_plants":15,"SC_germination_seeds":1073,"WFD_germination_seeds":1277},
       "STUDY_2017_ORIGINAL_COHORT_COUNT_CHANGED")
    require(all(source["ontology"][k] for k in
                ("head_rank_is_not_world_gravity_angle",
                 "head_rank_is_not_outer_phyllary_spread",
                 "germination_lab_temp_is_not_field_temperature")),
            "HEAD_STRUCTURE_AND_RANK_ONTOLOGY_CONFUSED")
    require(all(source["claim_ceiling"].values()),"INFERENCE_CEILING_CHANGED")


def run(source):
    assert_provenance(source)
    rows=source["published_2021_model_rows"]
    coverage=validate(rows,strict=True)
    contrasts=describe_weather_rank(rows)
    checks={
       (2016,"SC","terminal",20):.828,
       (2016,"SC","secondary",20):.197,
       (2016,"WFD","terminal",20):.750,
       (2016,"WFD","secondary",20):.436,
       (2017,"WFD","terminal",20):.360,
       (2017,"WFD","secondary",20):.393,
       (2017,"WFD","terminal",30):.071,
       (2017,"WFD","secondary",30):.081,
    }
    for r in rows:
        k=(r["collection_year"],r["site"],r["head_rank"],r["day_temp_c"])
        if k in checks:
            require(abs(germination_value(r)-checks[k])<1e-12,
                    "PUBLISHED_PROVENANCE_ANCHOR_CELL_CHANGED")
    return {
      "version":"cirsium_pitcheri_head_rank_weevil_germination_secondary_analysis_v1",
      "analysis_basis":source["version"],
      "original_authors_are_data_owners":True,
      "source_2020":source["source2020"],
      "source_2021":source["source2021"],
      "coverage":coverage,
      "reported_weevil_seed_stage":weevil_summary(source["source2020"]),
      "published_2016_and_2017_rank_germination":contrasts,
      "claim_ceiling":{
        "two_papers_do_not_yield_original_same_head_id_join":True,
        "rank_related_to_size_and_anthesis_order_not_independent_orientation_or_spine_effect":True,
        "weevil_infested_heads_are_not_randomized_to_weevils":True,
        "germination_fitted_means_not_raw_seed_data":True,
        "different_collection_years_different_germination_data_and_analyses":True,
        "sign_of_point_estimate_not_a_significant_rank_effect":True,
        "no_inferred_net_recruitment_reproductive_fitness_by_multiplication":True,
        "no_genetic_architecture_or_reexpression_measured":True,
      },
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--published",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    source=json.loads(a.published.read_text(encoding="utf-8"))
    result=run(source)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
       "published_table_rows":result["coverage"]["n_source_model_cells"],
       "2020_approx_infestation_ratios":{
         k:v for k,v in result["reported_weevil_seed_stage"].items()
         if "ratio" in k},
       "2016_terminal_higher_conditions":
         result["published_2016_and_2017_rank_germination"]["n_2016_terminal_better_out_of_six"],
       "2017_terminal_higher_conditions":
         result["published_2016_and_2017_rank_germination"]["n_2017_terminal_better_out_of_six"],
       "SC_2016_20C_ratio":
         result["published_2016_and_2017_rank_germination"]["germination_sc20_2016_terminal_secondary_ratio"],
       "WFD_2017_20C_ratio":
         result["published_2016_and_2017_rank_germination"]["germination_wfd20_2017_terminal_secondary_ratio"],
    },indent=2))


if __name__=="__main__":
    main()
