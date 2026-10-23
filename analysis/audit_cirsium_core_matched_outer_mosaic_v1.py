#!/usr/bin/env python3
"""Cirsium head-only paired phenotype mosaic: authentic measured plants.

This code asks whether plants with similar corolla + involucre core dimensions
can have different RELATIVE phyllary-tip envelope widths, conditionally within
source taxon and sex strata. It never declares mechanical spines, visitor
entry, heritable genetics, ancestry, trait regain or natural selection.

Source: Michalkova et al. (2023), Plant Systematics and Evolution,
DOI 10.1007/s00606-023-01854-2, original publisher Online Resource 6.
The published 129 individuals contain 114 source-vouchered nonhybrids; their
8 morphology cells each have been independently checked against original
publisher XLSX (114/114 rows, 912/912 cells) in the aza3 source audit.

Matched "pairs" are samples of one individual each and are greedy
descriptive sets, NOT randomized interventions, not plant population
prevalence estimates, and not independent Cirsium evolutionary origins.
No within-population assignment or habitat is verified from this derivative.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

SOURCE_DOI="10.1007/s00606-023-01854-2"
ORIGINAL_XLSX_SHA256="c06423cc6f9800ba4775a4428afa072ba1b3dc6fb48ccd05fdb7e42a2f03ee79"
SUBSET_GIT_BLOB_SHA="be95e6002ab45393786de05a7180eb77f38d941d"
MAIN_VALUES=("Cor_L","Inv_L","Inv_W","Ibr_W","Cap_W","Inv_Sp")
EXPECTED_STRATA={
    "Cirsium acaulon|f":10,
    "Cirsium acaulon|h":30,
    "Cirsium bertolonii - Monte Sagro|f":10,
    "Cirsium bertolonii - Monte Sagro|h":24,
    "Cirsium erisithales|f":10,
    "Cirsium erisithales|h":30,
}
CONFIGS={
    "primary":{"cor_tol":.10,"inv_tol":.10,"cap_tol":.15,"outer_min":.10},
    "strict_core":{"cor_tol":.05,"inv_tol":.05,"cap_tol":.10,"outer_min":.10},
    "relaxed_core":{"cor_tol":.15,"inv_tol":.15,"cap_tol":.15,"outer_min":.10},
    "modest_outer":{"cor_tol":.10,"inv_tol":.10,"cap_tol":.15,"outer_min":.05},
}


def demand(x,msg):
    if not x:
        raise ValueError(msg)


def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def read_original_subset(raw):
    demand(blob(raw)==SUBSET_GIT_BLOB_SHA,"PUBLISHER_VERIFIED_SUBSET_BLOB_CHANGED")
    rd=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    required={"excel_row","taxon","sex"}|set(MAIN_VALUES)
    demand(required.issubset(rd.fieldnames or []),"MORPHOLOGY_SOURCE_SCHEMA_CHANGED")
    rows=[]
    for r in rd:
        rec={
            "taxon":r["taxon"],"sex":r["sex"],
            "excel_row":int(r["excel_row"]),
        }
        demand(rec["sex"] in ("f","h") and rec["taxon"],
               "UNKNOWN_SOURCE_SEX_OR_TAXON")
        for k in MAIN_VALUES:
            v=float(r[k])
            demand(math.isfinite(v) and v>0,"MISSING_OR_NONPOSITIVE_"+k)
            rec[k]=v
        demand(math.isclose(rec["Inv_Sp"],rec["Ibr_W"]/rec["Inv_W"],
                            abs_tol=1e-10,rel_tol=0),
               "INV_SP_IS_NOT_ORIGINAL_OUTER_TIP_SPREAD_RATIO")
        rows.append(rec)
    demand(len(rows)==114 and len({r["excel_row"] for r in rows})==114,
           "ORIGINAL_114_SOURCE_PLANTS_DRIFT")
    demand(all(rows[i]["excel_row"]<rows[i+1]["excel_row"]
               for i in range(len(rows)-1)),
           "ORIGINAL_EXCEL_ROW_ORDER_CHANGED")
    strata=Counter(r["taxon"]+"|"+r["sex"] for r in rows)
    demand(dict(strata)==EXPECTED_STRATA,"ORIGINAL_SEX_TAXON_STRATA_DRIFT")
    return rows


def all_candidate_edges(rows,threshold):
    groups=defaultdict(list)
    for r in rows:
        groups[r["taxon"]+"|"+r["sex"]].append(r)
    edges=[]
    for group in sorted(groups):
        g=groups[group]
        for i in range(len(g)):
            for j in range(i+1,len(g)):
                a,b=g[i],g[j]
                co=abs(math.log(a["Cor_L"]/b["Cor_L"]))
                iv=abs(math.log(a["Inv_L"]/b["Inv_L"]))
                cp=abs(math.log(a["Cap_W"]/b["Cap_W"]))
                if not (co<=threshold["cor_tol"] and
                        iv<=threshold["inv_tol"] and
                        cp<=threshold["cap_tol"]):
                    continue
                ext=abs(math.log(a["Inv_Sp"]/b["Inv_Sp"]))
                score=math.sqrt((co/threshold["cor_tol"])**2+
                                (iv/threshold["inv_tol"])**2+
                                (cp/threshold["cap_tol"])**2)
                edges.append({
                    "taxon":a["taxon"],
                    "sex":a["sex"],
                    "excel_rows":[a["excel_row"],b["excel_row"]],
                    "core_distance_scaled":score,
                    "outer_log_contrast":ext,
                    "outer_spread_ratio_max_over_min":math.exp(ext),
                    "is_outer_divergent":ext>=threshold["outer_min"],
                    "core_absolute_log_differences":{
                        "Cor_L":co,"Inv_L":iv,"Cap_W":cp},
                })
    return edges


def greedy(edges):
    seen=set()
    selected=[]
    for edge in edges:
        a,b=edge["excel_rows"]
        if a in seen or b in seen:
            continue
        seen.update((a,b))
        selected.append(edge)
    demand(len({n for p in selected for n in p["excel_rows"]})==
           2*len(selected),"REUSED_MOTHER_PLANT_IN_INDEPENDENT_PAIR_SET")
    return selected


def matching(rows,config):
    edges=all_candidate_edges(rows,config)
    core=greedy(sorted(edges,key=lambda z:(
        z["core_distance_scaled"],z["excel_rows"][0],z["excel_rows"][1])))
    eligible=[z for z in edges if z["is_outer_divergent"]]
    outer=greedy(sorted(eligible,key=lambda z:(
        -z["outer_log_contrast"],z["core_distance_scaled"],
         z["excel_rows"][0],z["excel_rows"][1])))
    groups=sorted(EXPECTED_STRATA)
    def subset(arr,g):
        return [z for z in arr if z["taxon"]+"|"+z["sex"]==g]
    gc={}
    for g in groups:
        i=subset(edges,g)
        j=subset(eligible,g)
        h=subset(core,g)
        o=subset(outer,g)
        gc[g]={
            "n_individuals":EXPECTED_STRATA[g],
            "core_compatible_edges":len(i),
            "outer_divergent_edges":len(j),
            "core_first_disjoint_pairs":len(h),
            "core_first_pairs_outer_divergent":
                sum(z["is_outer_divergent"] for z in h),
            "outer_first_disjoint_pairs":len(o),
        }
    return {
        "thresholds":config,
        "n_core_compatible_edges":len(edges),
        "n_outer_divergent_edges":len(eligible),
        "n_disjoint_core_first_pairs":len(core),
        "n_disjoint_core_first_outer_divergent":
            sum(z["is_outer_divergent"] for z in core),
        "n_disjoint_outer_first_pairs":len(outer),
        "group_counts":gc,
        "core_first_pairs":core,
        "outer_first_pairs":outer,
    }


def analyze(raw):
    rows=read_original_subset(raw)
    analyses={name:matching(rows,threshold)
              for name,threshold in CONFIGS.items()}
    p=analyses["primary"]
    demand((p["n_core_compatible_edges"],p["n_outer_divergent_edges"],
            p["n_disjoint_core_first_pairs"],
            p["n_disjoint_core_first_outer_divergent"],
            p["n_disjoint_outer_first_pairs"])==(239,67,44,9,16),
           "ORIGINAL_PRIMARY_MATCHED_PHENOTYPE_MOSAIC_COUNTS_DRIFT")
    demand(p["group_counts"]["Cirsium acaulon|h"][
           "outer_divergent_edges"]==0 and
           p["group_counts"]["Cirsium acaulon|f"][
           "outer_divergent_edges"]==0,
           "SOURCE_ACAULON_TIP_RATIO_NEAR_FLOOR_CHANGED")
    demand(sum(v["outer_first_disjoint_pairs"]
               for v in p["group_counts"].values())==16,
           "SOURCE_DISJOINT_SET_STRATIFICATION_DRIFT")
    return {
        "version":"aza3_actual_cirsium_outer_head_mosaic_matched_source_JS_v1",
        "source":{
            "publication":"Michalkova et al 2023 Plant Systematics and Evolution",
            "doi":SOURCE_DOI,
            "source_original_xlsx_sha256":ORIGINAL_XLSX_SHA256,
            "aza3_publisher_cell_replay":
                "114/114 rows, 912/912 morphology cells verified",
            "subset_csv_git_blob_sha":blob(raw),
            "n_original_individuals":len(rows),
            "n_strata":len(EXPECTED_STRATA),
            "taxa_concept_count":3,
            "numerical_proxy":
                "Inv_Sp = Ibr_W / Inv_W; not physical spine length",
        },
        "estimand":(
            "Within identical taxon and reproductive-sex strata, are there "
            "same-individual reproductive core size-compatible pairs whose "
            "relative external phyllary-tip spread differs? Descriptive "
            "feasibility, not selection."),
        "matching":{
            "core_variables":["Cor_L","Inv_L","Cap_W"],
            "core_match_maximum_absolute_log_ratio":
                {"Cor_L":.1,"Inv_L":.1,"Cap_W":.15},
            "external_proxy":"log(Ibr_W/Inv_W) = log(Inv_Sp)",
            "external_contrast_absolute_log_ratio_at_least":.1,
            "core_first":(
                "greedy across eligible edges ranked by scaled normalized "
                "core distance ascending; row id tie break; does not select "
                "on outer contrast"),
            "outer_first":(
                "greedy across only outer-divergent eligible edges ranked by "
                "outer contrast descending; row id tie break"),
            "greedy_not_optimal_maximum_cardinality":True,
            "source_strata_individuals_not_assumed_independent_phylogenetic_lineages":
                True,
        },
        "analyses":analyses,
        "caveats":{
            "same_plant_trait_matrix_not_necessarily_same_exact_floret_timepoint":True,
            "external_ratio_wider_phyllary_tip_envelope_not_true_spine_geometry":True,
            "C_acaulon_spread_near_floor":True,
            "selection_of_pairs_not_population_prevalence_estimate":True,
            "within_taxon_sex_not_proven_same_field_population":True,
            "random_photo_measurement_error_not_reestimated":True,
            "head_access_events_or_pollen_or_seed_fitness_measured":0,
            "genetic_LEGO_ancestral_regain_or_selection_evidence":False,
            "new_own_biological_samples":0,
        },
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--subset",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    result=analyze(a.subset.read_bytes())
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                                allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "original_plant_count":result["source"]["n_original_individuals"],
        "primary_summary":{
            key:value for key,value in result["analyses"]["primary"].items()
            if key not in ("core_first_pairs","outer_first_pairs")
        },
        "sensitivity":{
            name:{
                "n_core_compatible_edges":r["n_core_compatible_edges"],
                "n_outer_divergent_edges":r["n_outer_divergent_edges"],
                "n_disjoint_core_first_pairs":r["n_disjoint_core_first_pairs"],
                "n_disjoint_core_first_outer_divergent":
                    r["n_disjoint_core_first_outer_divergent"],
                "n_disjoint_outer_first_pairs":r["n_disjoint_outer_first_pairs"],
            } for name,r in result["analyses"].items()},
    },ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
