#!/usr/bin/env python3
"""Source-authentic ORIGINAL Cirsium vulgare life-stage herbivory and fate audit.

Klinkhamer, de Jong & van der Meijden (1988), Journal of Ecology 76:403-414,
DOI 10.2307/2260602. Author-uploaded original full-text text is available:
https://www.researchgate.net/publication/271759934_Production_Dispersal_and_Predation_of_Seeds_in_the_Biennial_Cirsium_Vulgare

Replays PUBLISHED Table 1, Table 2 and text values transcribed into the
source-registered JSON. The actual 117/67 original flowering plants and
268/97 collected heads are *different* observational grains. Five 1981/82
exclosure dates × scrub/herbaceous habitats and three mesh levels are ten
original aggregate table rows, not ten randomly assigned biological subjects.
No raw mother/head/cage individual records or author PDF binary verified.

1982 undamaged STEM comparison includes 50 ADDITIONAL NONRANDOM plants;
do not use a binomial randomized-treatment interpretation of 80% vs 15%.
1982/84 stem damage rate and predation within flowerheads are different
conditional seed-loss denominators: DO NOT sum, multiply or infer net
fitness with a causal path decomposition. Postdispersal bird-net 2 cm
permits rodents; 0.6 cm small-cage denies rodents and birds but admits
small insects. Author reports paired container comparisons p<0.05, not
a newly fitted exact p, and found bird-only net similar to open.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path

EXPECTED_ROWS=[
    ("1981-09-25","scrub",12,27,27),
    ("1981-10-22","scrub",9,22,19),
    ("1981-11-18","scrub",8,93,83),
    ("1982-01-20","scrub",7,67,71),
    ("1982-02-05","scrub",3,73,79),
    ("1981-09-25","herbaceous_vegetation",9,79,76),
    ("1981-10-22","herbaceous_vegetation",5,62,83),
    ("1981-11-18","herbaceous_vegetation",3,78,70),
    ("1982-01-20","herbaceous_vegetation",5,78,70),
    ("1982-02-05","herbaceous_vegetation",1,75,57),
]
EXPECTED_YEAR_COUNTS={
    "1982":(117,2,50,52,268,67,201),
    "1984":(67,51,0,51,97,49,48),
}

def check(cond,error):
    if not cond:
        raise ValueError(error)

def mean(values):
    check(len(values)>0,"EMPTY_SOURCE_AVERAGE")
    return sum(values)/len(values)

def audit_input(source):
    check(source["version"]=="cirsium_vulgare_original_1988_stage_heterogeneity_tables_v1",
          "ORIGINAL_PAPER_REGISTRY_VERSION_MISMATCH")
    meta=source["source"]
    check(meta["doi"]=="10.2307/2260602" and
          meta["paper_year"]==1988 and
          meta["study_locality"]=="Meijendel coastal sand-dunes, The Netherlands",
          "WRONG_CIRSIUM_VULGARE_PAPER_OR_STUDY_POPULATION")
    check(meta["reported_data_are_published_table_and_text_values_not_retrieved_individual_records"] is True
          and meta["original_publisher_pdf_bytes_verified"] is False,
          "FICTITIOUS_VERIFIED_INDIVIDUAL_SOURCE")
    check(set(source["study_years"])=={"1982","1984"},
          "FALSE_AUTHOR_YEAR_SERIES")
    for year,expected in EXPECTED_YEAR_COUNTS.items():
        cohort=source["study_years"][year]
        stem=cohort["stem_census"]
        head=cohort["before_dispersal_flowerhead_census"]
        observed=(
            stem["n_random_flowering_individuals"],
            stem["n_random_flowering_plants_undamaged_stems"],
            stem["n_additional_nonstratified_undamaged_stem_plants_nonnrandom"],
            stem["n_seed_sample_undamaged_stem_plants"],
            head["n_heads_with_visible_seed_predation"]+
                head["n_heads_without_visible_seed_predation"],
            head["n_heads_with_visible_seed_predation"],
            head["n_heads_without_visible_seed_predation"],
        )
        check(observed==expected,"SOURCE_ORIGINAL_SAMPLE_GRAIN_DRIFT_"+year)
        check(observed[1]+observed[2]==observed[3],
              "SOURCE_NONRANDOM_UNDAMAGED_STEM_SUPPLEMENT_LOST")
        check(stem["undamaged_plant_sample_not_random_population_sample"] is
              (year=="1982"),
              "MISSING_1982_NONRANDOM_ADDITIONAL_50_STEMS")
        expected_stem_attack=98 if year=="1982" else 24
        expected_head_attack=25 if year=="1982" else 51
        expected_stem_loss=80 if year=="1982" else 15
        expected_head_loss=3 if year=="1982" else 17
        check(stem["percent_stems_attacked_any_herbivore_reported"]==
              expected_stem_attack and
              head["pct_flowerheads_with_traces_of_seed_predation_reported"]==
              expected_head_attack and
              stem["author_stem_damage_seed_production_reduction_pct"]==
              expected_stem_loss and
              head["mean_seed_loss_to_predation_inside_flowerheads_pct_reported"]==
              expected_head_loss,
              "ORIGINAL_PAPER_TWO_STAGES_UNDER_TWO_YEARS_CHANGED")
        expected_seed={
            "1982":(50,246,37,43,258),
            "1984":(1470,1713,139,208,2080),
        }[year]
        check((
            stem["mean_good_seeds_per_random_flowering_plant"],
            stem["mean_good_seeds_per_undamaged_stem_plant"],
            head["mean_good_seeds_per_attacked_head"],
            head["mean_good_seeds_per_unattacked_head"],
            cohort["author_maximum_no_stem_or_head_predation_potential_seeds_per_plant"],
        )==expected_seed,"ORIGINAL_SEED_PRODUCTION_TABLE_CHANGED_"+year)
        check(head["head_records_linked_to_original_stem_census_plant_ids"] is False,
              "FALSE_ORIGINAL_STEM_AND_HEAD_PERSON_LINKS")
        check(not (year=="1982" and stem["n_random_flowering_plants_undamaged_stems"]>2),
              "SOURCE_CHERRY_PICKED_1982_UNDAMAGED_STEMS")

    exp=source["postdispersal_exclosure"]
    check(exp["table_number"]==2 and len(exp["rows"])==10 and
          exp["n_repeated_experiment_dates"]==5 and
          exp["n_containers_per_exposure_series_as_described_by_author"]==10
          and exp["n_seeds_per_container_as_described_by_author"]==10,
          "AUTHOR_RODENT_EXCLOSURE_TABLE_UNITS_CHANGED")
    observed=[
        (r["date"],r["habitat"],r["cage"],r["net"],r["open"])
        for r in exp["rows"]
    ]
    check(observed==EXPECTED_ROWS,
          "PUBLISHED_1981_1982_TABLE2_ORIGINAL_VALUES_DRIFT")
    check(exp["treatments"]["cage"]["mesh_cm"]==.6 and
          exp["treatments"]["net"]["mesh_cm"]==2 and
          exp["treatments"]["open"]["mesh_cm"] is None,
          "SOURCE_0P6CM_VS_2CM_EXCLOSURE_MESH_CONFUSED")
    check(exp["table2_values_are_group_level_percent_seed_loss_not_individual_container_rows"] is True
          and exp["separate_paired_container_level_raw_observations_retrieved"] is False,
          "FICTITIOUS_INDEPENDENT_SEED_CONTAINER_REPLICATES")
    check(exp["separate_seeds_sown_95_local_populations_1982"][
              "n_populations"]==95 and
          exp["separate_1981_3_week_sown_grid_experiment"][
              "n_grids"]==9,
          "POSTDISPERSAL_SEPARATE_SAMPLE_COHORTS_CONFLATED")
    limits=source["limits"]
    required=[
        "stem_plant_n_1982_undamaged_includes_50_nonnrandomly_added",
        "author_stem_reduction_percent_is_comparison_of_seed_production_not_a_randomized_stem_attack",
        "flowerhead_infestation_and_seeds_observational_not_randomized_infestation",
        "stem_and_head_denominator_nonidentical",
        "postdispersal_exclosure_data_from_1981_to_1982_not_same_year_as_1984_head_census",
        "postdispersal_habitat_date_cells_not_independent_organismal_replicates",
        "no_plant_head_oviposition_pollen_service_filled_fitness_link",
    ]
    check(all(limits[x] is True for x in required) and
          limits["no_new_own_field_samples"]==0,
          "ECOLOGICAL_DENOMINATOR_OR_FIELD_CLAIM_UNSAFE")

def analyze(original):
    audit_input(original)
    years={}
    for year in ("1982","1984"):
        y=original["study_years"][year]
        stem=y["stem_census"]
        head=y["before_dispersal_flowerhead_census"]
        n=head["n_heads_with_visible_seed_predation"]+head["n_heads_without_visible_seed_predation"]
        years[year]={
            "n_stem_random_plants":stem["n_random_flowering_individuals"],
            "n_head_samples":n,
            "head_damage_count":head["n_heads_with_visible_seed_predation"],
            "head_damage_fraction_from_counts":
                head["n_heads_with_visible_seed_predation"]/n,
            "head_good_seeds_attacked_minus_unattacked_mean":
                head["mean_good_seeds_per_attacked_head"]-
                head["mean_good_seeds_per_unattacked_head"],
            "head_good_seeds_attacked_over_unattacked_mean":
                head["mean_good_seeds_per_attacked_head"]/
                head["mean_good_seeds_per_unattacked_head"],
            "stem_damage_reported_pct":
                stem["percent_stems_attacked_any_herbivore_reported"],
            "head_damage_reported_pct":
                head["pct_flowerheads_with_traces_of_seed_predation_reported"],
            "stem_seed_reduction_pct_replayed":
                100*(1-stem["mean_good_seeds_per_random_flowering_plant"]/
                     stem["mean_good_seeds_per_undamaged_stem_plant"]),
            "within_head_predation_mean_seed_loss_pct_reported":
                head["mean_seed_loss_to_predation_inside_flowerheads_pct_reported"],
            "mean_seed_per_random_plant":
                stem["mean_good_seeds_per_random_flowering_plant"],
            "nondrandom_undamaged_added":
                stem["n_additional_nonstratified_undamaged_stem_plants_nonnrandom"],
        }
    change={
        "stem_damaged_prevalence_pp":
            years["1984"]["stem_damage_reported_pct"]-
            years["1982"]["stem_damage_reported_pct"],
        "head_visible_damaged_prevalence_pp":
            years["1984"]["head_damage_reported_pct"]-
            years["1982"]["head_damage_reported_pct"],
        "reported_stem_seed_reduction_pp":
            original["study_years"]["1984"]["stem_census"][
                "author_stem_damage_seed_production_reduction_pct"]-
            original["study_years"]["1982"]["stem_census"][
                "author_stem_damage_seed_production_reduction_pct"],
        "reported_predation_within_heads_seed_loss_pp":
            years["1984"]["within_head_predation_mean_seed_loss_pct_reported"]-
            years["1982"]["within_head_predation_mean_seed_loss_pct_reported"],
        "year_change_is_different_field_cohorts_not_random_experimental_assignment":
            True,
    }
    exp=original["postdispersal_exclosure"]
    rows=exp["rows"]
    treatments=("cage","net","open")
    avgs={k:mean([r[k] for r in rows]) for k in treatments}
    dated=[{
        "date":r["date"],"habitat":r["habitat"],
        "open_minus_cage":r["open"]-r["cage"],
        "net_minus_cage":r["net"]-r["cage"],
        "net_minus_open":r["net"]-r["open"],
    } for r in rows]
    by_habitat={}
    for h in exp["habitats"]:
        small=[r for r in rows if r["habitat"]==h]
        by_habitat[h]={
            "n_dates":len(small),
            "means":{k:mean([r[k] for r in small]) for k in treatments},
            "open_minus_cage":mean([r["open"]-r["cage"] for r in small]),
            "net_minus_open":mean([r["net"]-r["open"] for r in small]),
        }
    results={
        "schema":"aza3_original_cirsium_vulgare_year_stem_head_postdispersal_independent_JS_v1",
        "as_of":"2026-10-26",
        "source":{
            "doi":original["source"]["doi"],
            "original_author_year":1988,
            "source_fulltext_author_uploaded_url":
                original["source"]["public_author_uploaded_fulltext"],
            "source_input_type":
                "published tables/narrative, not individual raw field records",
            "source_table1_years":["1982","1984"],
            "source_table2_exclosure_dates":exp["study_dates"],
        },
        "year_biological_stage_contrasts":years,
        "year_change_1984_minus_1982":change,
        "postdispersal_exclosure":{
            "n_source_habitat_dates":len(rows),
            "n_distinct_dates":5,"n_habitats":2,
            "original_measured_table2_cells":len(rows)*3,
            "equal_habitat_date_mean_predation_pct":avgs,
            "source_original_date_habitat_rows":dated,
            "equal_habitat_date_mean_open_minus_cage_pp":
                mean([x["open_minus_cage"] for x in dated]),
            "equal_habitat_date_mean_net_minus_cage_pp":
                mean([x["net_minus_cage"] for x in dated]),
            "equal_habitat_date_mean_net_minus_open_pp":
                mean([x["net_minus_open"] for x in dated]),
            "n_habitat_dates_open_greater_cage":
                sum(x["open_minus_cage"]>0 for x in dated),
            "n_habitat_dates_net_greater_cage":
                sum(x["net_minus_cage"]>0 for x in dated),
            "per_habitat":by_habitat,
            "reported_author_bird_only_net_vs_open_no_detected_difference":
                True,
            "reported_author_rodent_and_bird_cage_below_both_at_source_p_lt":
                .05,
            "descriptive_10_cells_not_ten_independent_origin_or_plant_fitness_measures":
                True,
        },
        "warnings":{
            "not_legitimate_to_add_stem_seed_reduction_plus_flowerhead_predation_pct":
                True,
            "nonnrandom_1982_undamaged_stem_comparator":True,
            "different_study_years_and_postdispersal_year":True,
            "observed_head_groups_are_not_randomized_damage_treatment":True,
            "no_detailed_source_flowerhead_individual_ids":True,
            "mouse_exclosure_not_sticky_thistle_head_defense":True,
            "no_genetically_inherited_capitulum_LEGO_selection_or_trait_reexpression_inferred":
                True,
            "new_own_observations":0,
        },
    }
    check(years["1982"]["n_head_samples"]==268 and
          years["1984"]["n_head_samples"]==97 and
          change["stem_damaged_prevalence_pp"]==-74 and
          change["head_visible_damaged_prevalence_pp"]==26,
          "ORIGINAL_STAGE_CONTRAST_DRIFT")
    check(math.isclose(results["postdispersal_exclosure"][
          "equal_habitat_date_mean_open_minus_cage_pp"],57.3) and
          results["postdispersal_exclosure"][
          "n_habitat_dates_open_greater_cage"]==10,
          "ORIGINAL_SOURCE_EXCLOSURE_EXPERIMENT_DRIFT")
    return results

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    original=json.loads(args.source.read_text(encoding="utf-8"))
    result=analyze(original)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                                   allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "source":result["source"],
        "year":result["year_biological_stage_contrasts"],
        "change":result["year_change_1984_minus_1982"],
        "exclosure_means":result["postdispersal_exclosure"][
            "equal_habitat_date_mean_predation_pct"],
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
