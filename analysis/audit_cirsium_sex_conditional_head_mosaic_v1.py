#!/usr/bin/env python3
"""Original-plant Cirsium head phenotype contrasts by reproductive sex.

Michalkova et al (2023), DOI 10.1007/s00606-023-01854-2.
The subset is 114 originally recorded non-hybrid individuals (3 taxon labels,
female n=10 per taxon vs hermaphrodite n=24/30), 912 published numeric cells
previously source-verified. Source: aza3 data/evidence/...subset_v1.csv.

The measured endpoints here are relative source head architecture (corolla
length, outer phyllary-tip envelope Ibr_W/Inv_W); NOT calibrated spine geometry,
mechanical exclusion, a causal sex-change treatment, pollination, seed predation,
adaptive genetics or ancestral re-expression. Taxa (n=3) are not independent
evolutionary transitions or randomly sampled genus-wide replication.
Monte Carlo sex-label exchange is CONDITIONAL on the exact source taxon and
sample counts, NOT a random-assignment causal experiment.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

EXPECTED_BLOB = "be95e6002ab45393786de05a7180eb77f38d941d"
AUTHOR_ORIGINAL_XLSX_SHA256 = (
    "c06423cc6f9800ba4775a4428afa072ba1b3dc6fb48ccd05fdb7e42a2f03ee79"
)
SOURCE_DOI = "10.1007/s00606-023-01854-2"
SEED = 20261024
MASK = 0xffffffff
PERMUTATIONS = 20000
BOOTSTRAPS = 5000
PRIMARY = ["corolla_log", "outer_ratio_log"]
FEATURES = [
    "corolla_log", "outer_ratio_log", "relative_floret_log",
    "absolute_tip_width_log", "involucre_body_width_log",
    "involucre_length_log", "head_width_log"
]
TAXA = [
    "Cirsium acaulon", "Cirsium bertolonii - Monte Sagro",
    "Cirsium erisithales"
]
F_COUNTS = [10, 10, 10]
H_COUNTS = [30, 24, 30]


def check(expr, error):
    if not expr:
        raise ValueError(error)


def blob(raw):
    return hashlib.sha1(
        b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw
    ).hexdigest()


def mean(xs):
    check(len(xs) > 0, "EMPTY_OBSERVATION_GROUP")
    return sum(xs) / len(xs)


def quantile(xs, p):
    seq = sorted(xs)
    ix = p * (len(seq) - 1)
    lo = int(ix)
    return seq[lo] * (1 - (ix-lo)) + seq[min(lo+1, len(seq)-1)] * (ix-lo)


def xorshift32(seed=SEED):
    check(type(seed) is int and 0 < seed < 2**32,
          "INVALID_OR_ZERO_DETERMINISTIC_SEED")
    state = seed
    while True:
        state ^= (state << 13) & MASK
        state ^= (state >> 17)
        state ^= (state << 5) & MASK
        state &= MASK
        yield state


def source_features(r):
    pos = ["Cor_L", "Inv_L", "Inv_W", "Ibr_W", "Cap_W", "Inv_Sp"]
    values = {}
    for k in pos:
        try:
            values[k] = float(r[k])
        except (KeyError, ValueError, TypeError):
            raise ValueError("SOURCE_POSITIVE_HEAD_COLUMN_MISSING_"+k) from None
        check(math.isfinite(values[k]) and values[k]>0,
              "SOURCE_INVALID_HEAD_LENGTH_"+k)
    check(math.isclose(values["Inv_Sp"],
                       values["Ibr_W"] / values["Inv_W"],
                       abs_tol=1e-10, rel_tol=0),
          "ORIGINAL_OUTER_TIP_RATIO_INVALID")
    return {
        "corolla_log": math.log(values["Cor_L"]),
        "outer_ratio_log": math.log(values["Inv_Sp"]),
        "relative_floret_log": math.log(values["Cor_L"] / values["Inv_L"]),
        "absolute_tip_width_log": math.log(values["Ibr_W"]),
        "involucre_body_width_log": math.log(values["Inv_W"]),
        "involucre_length_log": math.log(values["Inv_L"]),
        "head_width_log": math.log(values["Cap_W"]),
    }


def read_source(raw, strict=True):
    sha=blob(raw)
    if strict:
        check(sha==EXPECTED_BLOB,"SOURCE_ORIGINAL_114_GIT_BLOB_DRIFT")
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    req = {"excel_row","taxon","sex","Cor_L","Inv_L","Inv_W",
           "Ibr_W","Cap_W","Inv_Sp"}
    check(req.issubset(reader.fieldnames or ()),
          "SOURCE_ORIGINAL_114_COLUMN_DRIFT")
    rows=list(reader)
    check(len(rows)==114,"SOURCE_NOT_114_ORIGINAL_PLANTS")
    ids=[]
    groups=[]
    for taxon,n_f,n_h in zip(TAXA,F_COUNTS,H_COUNTS):
        rr=[r for r in rows if r["taxon"]==taxon]
        check(len(rr)==n_f+n_h,
              "SOURCE_TAXON_PLANT_COUNTS_CHANGED")
        try:
            numeric_ids=[int(r["excel_row"]) for r in rr]
        except (ValueError,TypeError):
            raise ValueError("SOURCE_EXCEL_ROW_NOT_INTEGER") from None
        check(len(set(numeric_ids))==len(rr),
              "DUPLICATED_PLANT_ORIGINAL_EXCEL_ROW")
        ids+=numeric_ids
        female=[source_features(r) for r in rr if r["sex"]=="f"]
        herm=[source_features(r) for r in rr if r["sex"]=="h"]
        check(len(female)==n_f and len(herm)==n_h,
              "SOURCE_FEMALE_HERMAPHRODITE_COUNTS_CHANGED")
        groups.append({
            "taxon":taxon,"n_f":n_f,"n_h":n_h,
            "all":[{"sex":r["sex"],"values":source_features(r)} for r in rr],
            "female":female,"herm":herm,
        })
    check(len(set(ids))==len(rows) and
          set(r["taxon"] for r in rows)==set(TAXA) and
          all(r["sex"] in ("f","h") for r in rows),
          "SOURCE_PLANT_REPETITION_OR_UNSUPPORTED_SEX")
    return groups,sha


def contrast(female,hermaphrodite):
    return {
        k:mean([x[k] for x in female]) -
          mean([x[k] for x in hermaphrodite])
        for k in FEATURES
    }


def conditional_permutations(groups, observed, pooled, draws=PERMUTATIONS):
    check(draws==PERMUTATIONS,"SOURCE_MONTE_CARLO_REPLICATION_BUDGET_FIXED")
    counts={t+"|"+k:0 for t in TAXA for k in PRIMARY}
    counts.update({"pooled|"+k:0 for k in PRIMARY})
    rng=xorshift32(SEED)
    for _ in range(draws):
        diffs=[]
        for g in groups:
            perm=list(range(len(g["all"])))
            for j in range(len(perm)-1,0,-1):
                k=next(rng)%(j+1)
                perm[j],perm[k]=perm[k],perm[j]
            f=[g["all"][k]["values"] for k in perm[:g["n_f"]]]
            h=[g["all"][k]["values"] for k in perm[g["n_f"]:]]
            c={key:mean([v[key] for v in f]) -
                    mean([v[key] for v in h])
               for key in PRIMARY}
            diffs.append(c)
            for key in PRIMARY:
                if abs(c[key]) >= abs(observed[g["taxon"]][key])-1e-12:
                    counts[g["taxon"]+"|"+key]+=1
        for key in PRIMARY:
            if abs(mean([d[key] for d in diffs])) >= abs(pooled[key])-1e-12:
                counts["pooled|"+key]+=1
    ps={k:(v+1)/(draws+1) for k,v in counts.items()}
    keys=[t+"|"+k for t in TAXA for k in PRIMARY]
    keys.sort(key=lambda k:(ps[k],k))
    qvalues={}
    q=1.0
    for i in range(len(keys)-1,-1,-1):
        key=keys[i]
        q=min(q,ps[key]*len(keys)/(i+1))
        qvalues[key]=q
    return ps,qvalues


def within_stratum_bootstrap(groups,draws=BOOTSTRAPS):
    check(draws==BOOTSTRAPS,"SOURCE_BOOTSTRAP_REPLICATION_BUDGET_FIXED")
    distrib={t+"|"+k:[] for t in TAXA+["pooled"] for k in FEATURES}
    # Continuing RNG stream is intentionally NOT done: the original
    # independent JS starts a fresh PRNG before permutation then continues
    # the same stream through bootstrap. Reproduce that exact state by
    # consuming the prerequisite sex permutations in the outer caller.
    raise NotImplementedError("USE combined_replay; RNG must remain shared")


def combined_replay(groups):
    observed={g["taxon"]:contrast(g["female"],g["herm"]) for g in groups}
    pooled={k:mean([observed[t][k] for t in TAXA]) for k in FEATURES}
    rng=xorshift32(SEED)
    exceed={t+"|"+k:0 for t in TAXA for k in PRIMARY}
    exceed.update({"pooled|"+k:0 for k in PRIMARY})
    for _ in range(PERMUTATIONS):
        deltas=[]
        for g in groups:
            inds=list(range(len(g["all"])))
            for j in range(len(inds)-1,0,-1):
                k=next(rng)%(j+1)
                inds[j],inds[k]=inds[k],inds[j]
            ff=[g["all"][i]["values"] for i in inds[:g["n_f"]]]
            hh=[g["all"][i]["values"] for i in inds[g["n_f"]:]]
            d={k:mean([x[k] for x in ff])-mean([x[k] for x in hh])
               for k in PRIMARY}
            deltas.append(d)
            for k in PRIMARY:
                if abs(d[k])>=abs(observed[g["taxon"]][k])-1e-12:
                    exceed[g["taxon"]+"|"+k]+=1
        for k in PRIMARY:
            if abs(mean([x[k] for x in deltas]))>=abs(pooled[k])-1e-12:
                exceed["pooled|"+k]+=1
    ps={k:(v+1)/(PERMUTATIONS+1) for k,v in exceed.items()}
    keys=[t+"|"+k for t in TAXA for k in PRIMARY]
    keys.sort(key=lambda k:(ps[k],k))
    qvalues={}
    q=1.0
    for i in range(len(keys)-1,-1,-1):
        key=keys[i]
        q=min(q,ps[key]*len(keys)/(i+1))
        qvalues[key]=q
    distrib={t+"|"+k:[] for t in TAXA+["pooled"] for k in FEATURES}
    for _ in range(BOOTSTRAPS):
        deltas=[]
        for g in groups:
            ff=[g["female"][next(rng)%g["n_f"]] for j in range(g["n_f"])]
            hh=[g["herm"][next(rng)%g["n_h"]] for j in range(g["n_h"])]
            dd=contrast(ff,hh)
            deltas.append(dd)
            for k in FEATURES:
                distrib[g["taxon"]+"|"+k].append(dd[k])
        for k in FEATURES:
            distrib["pooled|"+k].append(mean([d[k] for d in deltas]))
    def metadata(label,effects):
        return {
            k:{
                "mean_log_female_minus_hermaphrodite":effects[k],
                "ratio_percent":100*math.expm1(effects[k]),
                "ci95":[quantile(distrib[label+"|"+k],.025),
                        quantile(distrib[label+"|"+k],.975)],
                "conditional_permutation_p_two_sided":(
                    ps[label+"|"+k] if k in PRIMARY else None),
                "within_taxa_six_test_BH_q":(
                    qvalues[label+"|"+k] if label!="pooled" and
                    k in PRIMARY else None),
            } for k in FEATURES
        }
    return observed,pooled,metadata,ps,qvalues


def analyze(raw):
    groups,sha=read_source(raw)
    observed,pooled,metadata,ps,qvals=combined_replay(groups)
    by_taxon={
        g["taxon"]:{
            "n_f":g["n_f"],"n_h":g["n_h"],
            "features":metadata(g["taxon"],observed[g["taxon"]])
        } for g in groups}
    left={t:{
        k:mean([observed[z][k] for z in TAXA if z!=t])
        for k in FEATURES
    } for t in TAXA}
    assert all(observed[t]["corolla_log"]<0 for t in TAXA)
    return {
        "schema":"aza3_original_Cirsium_head_sex_module_contrast_independent_JS_v1",
        "as_of":"2026-10-24",
        "source":{
            "paper_doi":SOURCE_DOI,
            "published_spreadsheet_sha256":AUTHOR_ORIGINAL_XLSX_SHA256,
            "source_derived_subset_git_blob":sha,
            "original_source_n_individuals":114,
            "verified_original_morphometry_cells":912,
            "taxa":list(TAXA),
            "source_strata":
                "sex f female; h hermaphroditic; one original plant per row",
        },
        "design":{
            "primary_outcomes":list(PRIMARY),
            "sex_contrast":
                "mean(log(feature)) among female - mean(log(feature)) among hermaphroditic in source taxon",
            "outer_spread_source_identity":
                "Inv_Sp = Ibr_W / Inv_W (not true spine)",
            "secondary_descriptive_features":
                [k for k in FEATURES if k not in PRIMARY],
            "pooled_effect":
                "unweighted mean of 3 taxon-specific original female-h log contrasts, not three independent phylogenetic origins",
            "randomization_reference":
                "Monte Carlo exchange original sex labels WITHIN each taxon conditional on sex counts, not sex-randomized treatment",
            "permutation_draws":PERMUTATIONS,
            "pseudocount_p":"+1 numerator / +1 denominator",
            "six_test_BH_family":
                "3 original taxon strata x 2 predesignated outcomes; exploratory discovered after original source inspection",
            "bootstrap_draws":BOOTSTRAPS,
            "bootstrap_unit":
                "original plant within taxon x sex, resample each stratum separately",
            "rng":"xorshift32",
            "rng_seed":SEED,
            "source_group_order":list(TAXA),
            "feature_order":list(FEATURES),
            "sort_rows_preserve_publisher_subset":True,
        },
        "by_taxon":by_taxon,
        "pooled_equal_taxon":metadata("pooled",pooled),
        "leave_one_original_taxon_out_equal_taxon":left,
        "limits":{
            "observational_sex_not_randomized":True,
            "original_taxa_three_not_independent_evolutionary_contrasts":True,
            "taxon_sex_sample_balance_not_population_prevalence":True,
            "corolla_shorter_is_not_direct_pollen_effect":True,
            "outer_ratio_mechanical_spines_not_measured":True,
            "outer_floor_effect_C_acaulon":True,
            "sampling_protocol_genetic_structure_not_adjusted":True,
            "within_sex_taxon_bootstrap_not_cross_genus_inference":True,
            "NO_same_plant_insect_access_or_viable_seed_fate":0,
            "new_own_field_samples":0,
            "adaptive_genetic_LEGO_or_ancestral_regain_proven":False,
        },
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--subset",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    result=analyze(args.subset.read_bytes())
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8")
    print(json.dumps({
        "plants":result["source"]["original_source_n_individuals"],
        "pooled":{k:result["pooled_equal_taxon"][k] for k in PRIMARY},
        "per_taxon":{
            t:{k:row["features"][k] for k in PRIMARY}
            for t,row in result["by_taxon"].items()
        },
    },ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
