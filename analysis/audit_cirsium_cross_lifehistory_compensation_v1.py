#!/usr/bin/env python3
"""Compare *published* monocarpic and *raw-plant* iterocarpic Cirsium 2x2 designs.

Source A: West & Louda, 2018, Oecologia DOI 10.1007/s00442-017-4027-9.
Original Table 1 means±SE, 134 C. canescens maternal plants across 2007/08.
Individual canescens plant rows have NOT been obtained.

Source B: West & Louda, 2021, AJB DOI 10.1002/ajb2.1744,
USDA DOI 10.15482/USDA.ADC/1522649 CC0, actual 106 C. undulatum mother
rows across 2007/08, source Git blob locked. Prior aza3 PR #73 audits
323 per-head rows and 2 conflicting maternal treatment identifiers.

The species-specific late-head seed-scoring methods and reduced-herbivory
strength have not been formally calibrated, and 1 monocarpic + 1 iterocarpic
species cannot identify *life-history type* as the causal moderator.
Observed means are descriptive; source authors originally reported the
biology. Do not assign canescens individual-level SEs or year-specific
means based on four published pooled group means.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math
from collections import Counter
from pathlib import Path

CAN_DOI="10.1007/s00442-017-4027-9"
UND_DOI="10.1002/ajb2.1744"
USDA_DOI="10.15482/USDA.ADC/1522649"
USDA_GIT_BLOB="458800272c3d0a03c7653702027a8f8578d9e88c"
C_GROUPS=("reduced_intact","reduced_damaged",
          "ambient_intact","ambient_damaged")
YEAR_ARMS={"2007":{"reduced_intact":10,"reduced_damaged":11,
                   "ambient_intact":9,"ambient_damaged":9},
           "2008":{"reduced_intact":21,"reduced_damaged":19,
                   "ambient_intact":15,"ambient_damaged":12}}
UND_ARM={"ndi":"reduced_intact","di":"reduced_damaged",
         "ndni":"ambient_intact","dni":"ambient_damaged"}


def require(ok,message):
    if not ok:raise ValueError(message)


def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def read_undulatum_source(path):
    raw=Path(path).read_bytes()
    require(git_blob(raw)==USDA_GIT_BLOB,"UNDULATUM_SOURCE_GIT_BLOB_DRIFT")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    require(len(rows)==106,"UNDULATUM_MOTHER_COUNT_DRIFT")
    expect=("plt","year","plt.gsd","dam","ins","trt","ap.hd","o.hds")
    require(all(set(expect).issubset(r) for r in rows),
            "UNDULATUM_SOURCE_SCHEMA_DRIFT")
    ids=[r["year"]+"|"+r["plt"] for r in rows]
    require(len(ids)==len(set(ids)),"DUPLICATE_UNDULATUM_MOTHER")
    for r in rows:
        require(r["year"] in YEAR_ARMS and
                r["trt"]==r["dam"]+r["ins"] and
                r["trt"] in UND_ARM,"INVALID_UNDULATUM_TREATMENT")
        for k in ("plt.gsd","ap.hd","o.hds"):
            try:
                x=float(r[k])
            except (ValueError,KeyError):
                raise ValueError("INVALID_UNDULATUM_GOOD_SEED_COUNT") from None
            require(math.isfinite(x) and x>=0,"INVALID_UNDULATUM_GOOD_SEED_COUNT")
    actual=Counter((r["year"],UND_ARM[r["trt"]]) for r in rows)
    require(actual==Counter({(y,a):n for y,d in YEAR_ARMS.items()
                            for a,n in d.items()}),
            "ORIGINAL_UNDULATUM_YEAR_TREATMENT_COUNTS_DRIFT")
    return rows


def validate_canescens_table(source):
    require(source["provenance"]["doi"]==CAN_DOI and
            source["provenance"]["original_raw_plant_rows_available_here"] is False
            and source["provenance"]["means_transcribed_from_published_article_not_independently_measured"] is True,
            "CANESCENS_PUBLISHED_NOT_RAW_SOURCE_CEILING_DRIFT")
    d=source["table1_plant_rows"]
    require(len(d)==4 and
            [r["treatment"] for r in d]==list(C_GROUPS),
            "CANESCENS_TABLE1_ARM_ORDER_DRIFT")
    require([r["n_plants"] for r in d]==[37,31,34,32] and
            sum(r["n_plants"] for r in d)==134,
            "CANESCENS_TABLE1_N_CHANGED")
    assert_cells=[
        ("undamaged_seeds_per_plant","mean",
            [355.9,477.8,186.3,112.0]),
        ("undamaged_seeds_per_plant","se",
            [49.0,81.8,29.6,27.9]),
        ("nonapical_heads_with_seed","mean",[1.9,3.5,1.5,1.3]),
        ("matured_seeds_per_nonapical_head","mean",[69.6,99.0,30.7,31.5]),
        ("nonapical_head_damage_percent","mean",[48.1,48.4,78.5,72.9])
    ]
    for field,part,expected in assert_cells:
        require(all(abs(r[field][part]-v)<1e-9
                    for r,v in zip(d,expected)),
                "CANESCENS_TABLE1_PUBLISHED_VALUE_CHANGED_"+field+"_"+part)
    require(source["provenance"]["author_reported_plant_level_apical_damage_by_herbivory_interaction_p"]==.04,
            "CANESCENS_AUTHORS_MODEL_P_CHANGED")
    return {r["treatment"]:r for r in d}


def means_from_undulatum(rows,exclude_id=()):
    excluded=set(exclude_id)
    records=[r for r in rows if
             r["year"]+"|"+r["plt"] not in excluded]
    groups={y:{a:[r for r in records if r["year"]==y and
                    UND_ARM[r["trt"]]==a] for a in C_GROUPS}
            for y in YEAR_ARMS}
    by_year={}
    pooled={}
    year_equal={}
    for yr,arm_groups in groups.items():
        by_year[yr]={}
        for arm,rs in arm_groups.items():
            require(rs,"EMPTY_UNDULATUM_YEAR_TREATMENT_CELL")
            seeds=[float(r["plt.gsd"]) for r in rs]
            sides=[float(r["o.hds"]) for r in rs]
            by_year[yr][arm]={
                "n_mothers":len(rs),
                "mean_good_seed_per_plant":sum(seeds)/len(rs),
                "mean_nonapical_good_seed_per_plant":sum(sides)/len(rs),
                "n_zero_good_seed_mothers":sum(s==0 for s in seeds),
            }
    for arm in C_GROUPS:
        rs=[r for y in groups for r in groups[y][arm]]
        pooled[arm]={
          "n_mothers":len(rs),
          "mean_good_seed_per_plant":
              sum(float(r["plt.gsd"]) for r in rs)/len(rs),
          "mean_nonapical_good_seed_per_plant":
              sum(float(r["o.hds"]) for r in rs)/len(rs),
        }
        year_equal[arm]={
          "n_year_cohorts":2,
          "mean_good_seed_per_plant":
              sum(by_year[y][arm]["mean_good_seed_per_plant"]
                  for y in YEAR_ARMS)/2,
        }
    return {"year_by_arm":by_year,
            "original_sample_weighted_pooled_by_arm":pooled,
            "equal_year_weighted_by_arm":year_equal,
            "excluded_mother_ids":sorted(excluded)}


def effect(g,field):
    means={k:v[field] for k,v in g.items()}
    require(all(means[k]>0 for k in ("reduced_intact","ambient_intact")),
            "ZERO_NO_DAMAGE_REFERENCE_MEAN")
    reduced_d=means["reduced_damaged"]-means["reduced_intact"]
    ambient_d=means["ambient_damaged"]-means["ambient_intact"]
    ratios={
       "reduced":means["reduced_damaged"]/means["reduced_intact"],
       "ambient":means["ambient_damaged"]/means["ambient_intact"],
    }
    return {
        "mean_by_arm":means,
        "damage_minus_intact_good_seeds_under_reduced_herbivory":reduced_d,
        "damage_minus_intact_good_seeds_under_ambient_herbivory":ambient_d,
        "damage_differences_difference_reduced_minus_ambient":reduced_d-ambient_d,
        "damage_to_intact_seed_mean_ratio_by_herbivory":ratios,
        "ratio_of_damage_ratios_reduced_to_ambient":
            ratios["reduced"]/ratios["ambient"],
        "no_source_specific_randomization_based_2x2_inference_from_means":True,
    }


def claim_limits():
    return {
       "canescens_134_original_mother_rows_recovered":False,
       "the_two_species_represent_two_independently_replicated_life_history_classes":False,
       "canescens_and_undulatum_seed_scoring_procedures_proven_identical":False,
       "published_canescens_2_year_means_are_year_specific":False,
       "between_species_life_history_effect_identified":False,
       "a_cross_species_significance_test_was_performed":False,
       "both_species_measured_true_phyllary_spine_or_head_inclination":False,
       "gene_level_reexpression_or_genetic_modules_measured":False,
       "new_canescens_field_samples_collected":False,
       "undulatum_323_head_rows_treated_as_independent_mothers":False,
       "undulatum_bad_original_treatment_labels_autocorrected":False,
    }


def analyze(can_source,und_rows):
    source=validate_canescens_table(can_source)
    original=means_from_undulatum(und_rows)
    excluded=means_from_undulatum(
        und_rows,exclude_id=("2007|323","2007|375","2007|327"))
    C={a:{
       "n_mothers":r["n_plants"],
       "mean_good_seed_per_plant":r["undamaged_seeds_per_plant"]["mean"],
       "author_published_mean_se":r["undamaged_seeds_per_plant"]["se"],
       "mean_nonapical_heads_with_seed":r["nonapical_heads_with_seed"]["mean"],
       "mean_mature_seed_per_nonapical_head":r["matured_seeds_per_nonapical_head"]["mean"],
       "mean_damage_pct_nonapical_head":r["nonapical_head_damage_percent"]["mean"],
    } for a,r in source.items()}
    c_eff=effect(C,"mean_good_seed_per_plant")
    u_eff=effect(original["original_sample_weighted_pooled_by_arm"],
                 "mean_good_seed_per_plant")
    equal_eff=effect(original["equal_year_weighted_by_arm"],
                     "mean_good_seed_per_plant")
    years={year:effect(group,"mean_good_seed_per_plant")
           for year,group in original["year_by_arm"].items()}
    excluded_eff={
        "by_year":{
             y:effect(groups,"mean_good_seed_per_plant")
             for y,groups in excluded["year_by_arm"].items()
        },
        "sample_weighted_pooled":effect(
             excluded["original_sample_weighted_pooled_by_arm"],
             "mean_good_seed_per_plant"),
    }
    return {
       "version":"aza3_20261021_comparable_cirsium_two_lifehistories_treatments_v1",
       "study_registry":{
          "canescens":{
             "life_history":"MONOCARPIC_ONE_REPRODUCTIVE_EPISODE",
             "source_doi":CAN_DOI,
             "author_original_plant_level_n":134,
             "data_grain":"PUBLISHED_TWO_YEAR_POOLED_GROUP_MEANS_AND_SE",
             "per_individual_plant_records_available":False,
             "original_author_apical_damage_by_herbivory_model_p":.04,
          },
          "undulatum":{
             "life_history":"ITEROCARPIC_REPEATED_REPRODUCTION_POSSIBLE",
             "source_doi":UND_DOI,
             "data_doi":USDA_DOI,
             "actual_plant_level_n":len(und_rows),
             "data_grain":"ORIGINAL_INDIVIDUAL_MOTHER_PLANT_CSV",
             "distinct_reproductive_year_cohorts":["2007","2008"],
             "can_infer_total_lifetime_fitness_from_2_seed_cohorts":False,
          }
       },
       "canescens_author_published_2_year_plant_treatment_stats":C,
       "undulatum_direct_original_mother_data":original,
       "compensation_contrasts":{
          "canescens_author_two_year_pooled":c_eff,
          "undulatum_original_mother_sample_weighted_two_year_pooled":u_eff,
          "undulatum_original_mother_equal_year_weighted":equal_eff,
          "undulatum_real_year_specific":years,
          "undulatum_exclude_3_conflicting_mothers":excluded_eff,
       },
       "cross_species_descriptive_ratio_pattern":{
          "ratio_of_damage_ratios_canescens_divided_by_undulatum_pooled":
             c_eff["ratio_of_damage_ratios_reduced_to_ambient"]/
             u_eff["ratio_of_damage_ratios_reduced_to_ambient"],
          "no_standardized_ecological_or_replicated_lifehistory_causal_effect":
             True,
          "within_source_raw_vs_published_and_pooling_protocols_differ":
             True,
          "conclusion":"C. canescens pooled published estimates are more herbivory-regime-sensitive in point means; C. undulatum individual mother outcomes vary strongly by year and source labels; no inference of a general semelparous/iteroparous divergence."
       },
       "no_canescens_raw_bootstrap_or_constructed_p_from_SEs":True,
       "all_scientific_inference_limits":claim_limits(),
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--canescens-published-json",type=Path,required=True)
    p.add_argument("--undulatum-mothers",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    z=analyze(
      json.loads(args.canescens_published_json.read_text(encoding="utf-8")),
      read_undulatum_source(args.undulatum_mothers))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(z,ensure_ascii=False,indent=2,
                                  allow_nan=False)+"\n",encoding="utf-8")
    for name,r in z["compensation_contrasts"].items():
        if name.startswith("undulatum_exclude"):continue
        if name=="undulatum_real_year_specific":
            for y,u in r.items():
                print("UND_YEAR",y,"damage ratios",
                      u["damage_to_intact_seed_mean_ratio_by_herbivory"])
        else:
            print(name, r["damage_to_intact_seed_mean_ratio_by_herbivory"],
                  r["ratio_of_damage_ratios_reduced_to_ambient"])
    print("NO_LIFE_HISTORY_CAUSAL_EFFECT_IDENTIFIED; 134 published summaries, 106 verified individual mothers.")

if __name__=="__main__":
    main()
