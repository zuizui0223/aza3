#!/usr/bin/env python3
"""Original Ipomopsis sequential selection: source-conditional *higher-order* contrasts.

Uses the exact author 437-line public CSV from Campbell et al. (2022),
American Naturalist DOI 10.1086/716740; data DOI 10.7280/D1KM49,
Zenodo record 5083399, author-listed MD5:
12ed3dbe6d269c5f808f9c6b230d468c.

Source statuses:
- An ORIGINAL CSV record is NOT a known independent genotype (no plant ID).
- Type A/H/T = I. aggregata/hybrid/I. tenuituba, source Site agg/hyb/ten.
  The original author README states hyb+ten share the upper physical garden.
- Outcomes are seeds initiated and fraction of FRUITS escaping fly
  seed-predation, not measured pollen delivery and not a sequential per-seed
  survival product. Do not multiply initiation and escape to recreate seeds.
- width = corolla width; color = petal red/green reflectance; sepal = sepal
  width. They are NOT Cirsium phyllary spines, inclination or measured access.
- Contrast added quadratic width or two trait product term against original
  linear trait main effects within Site×Year×Type. The source labels are not
  randomized experimental trait assignments. Freedman–Lane null references
  and source-row bootstrap are EXPLORATORY, not causal/randomized inference.
- All nonlinear/cross-term analyses are selected after inspecting original
  data. BH is a limited multiplicity guard; source phylogeny and true
  independent plant identifiers are unavailable.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

AUTHOR_MD5="12ed3dbe6d269c5f808f9c6b230d468c"
HEADERS=(
    "Site","Year","Type","width","color","nectar","sepal",
    "totpinene","indoleN","seeds","seedsinit","flyescape"
)
METRICS=("seed_initiation_log1p","fruit_fly_escape_fraction")
MODELS=(
    {"id":"width_quadratic","fields":["width"],"extra":"square"},
    {"id":"width_by_color","fields":["width","color"],"extra":"interaction"},
    {"id":"width_by_sepal","fields":["width","sepal"],"extra":"interaction"},
)
MASK=0xffffffff
YEAR_RANGE=range(2010,2019)
PERMUTATIONS=5000
BOOTSTRAPS=1500

def require(ok,msg):
    if not ok:
        raise ValueError(msg)

def mean(vals):
    require(len(vals)>0,"INVALID_EMPTY_MEAN")
    return sum(vals)/len(vals)

def quantile(vals,q):
    require(vals and 0<=q<=1,"INVALID_SOURCE_QUANTILE")
    arr=sorted(vals)
    t=(len(arr)-1)*q
    i=math.floor(t)
    f=t-i
    return arr[i]*(1-f)+arr[min(i+1,len(arr)-1)]*f

def xorshift32(seed):
    require(isinstance(seed,int) and 0<seed<2**32,"INVALID_OR_ZERO_SEED")
    state=seed
    while True:
        state^=(state<<13)&MASK
        state^=state>>17
        state^=(state<<5)&MASK
        state&=MASK
        yield state

def read_author_csv(raw):
    require(hashlib.md5(raw).hexdigest()==AUTHOR_MD5,
            "ORIGINAL_AUTHOR_ZENODO_CSV_MD5_CHANGED")
    stream=io.StringIO(raw.decode("utf-8-sig"))
    r=csv.DictReader(stream)
    require(tuple(r.fieldnames or ())==HEADERS,"ORIGINAL_12_HEADER_FIELDS_CHANGED")
    source=[]
    original_codes=set()
    for lineno,row in enumerate(r,2):
        require(None not in row and all(v is not None for v in row.values()),
                "ORIGINAL_INCOMPLETE_SOURCE_RECORD")
        normalized={}
        for j,k in enumerate(HEADERS):
            value=row[k].strip()
            if value in (".",""):
                normalized[k]=None
            elif j in (0,2):
                normalized[k]=value
            else:
                v=float(value)
                require(math.isfinite(v),"SOURCE_NUMBER_NOT_FINITE")
                normalized[k]=v
        require(normalized["Site"] in ("agg","hyb","ten") and
                normalized["Type"] in ("A","H","T") and
                normalized["Year"] in YEAR_RANGE,
                "ORIGINAL_SOURCE_CONTEXT_DRIFT")
        require(normalized["flyescape"] is None or
                0<=normalized["flyescape"]<=1,
                "INVALID_FLY_ESCAPE_FRACTION")
        require(normalized["seedsinit"] is None or
                normalized["seedsinit"]>=0,"INVALID_SEEDS_INITIATED")
        normalized["source_original_csv_row_id"]=lineno
        serialized=tuple(row[k] for k in HEADERS)
        require(serialized not in original_codes,
                "DUPLICATED_WHOLE_SOURCE_ROW")
        original_codes.add(serialized)
        source.append(normalized)
    require(len(source)==437,"ORIGINAL_437_SOURCE_RECORDS_CHANGED")
    return source

def solve(A,b):
    """Ordinary full-rank Gaussian normal equations, fixed deterministic pivots."""
    n=len(b)
    matrix=[[*A[j],b[j]] for j in range(n)]
    for i in range(n):
        best=max(range(i,n),key=lambda j:abs(matrix[j][i]))
        if abs(matrix[best][i])<1e-10:
            return None
        matrix[i],matrix[best]=matrix[best],matrix[i]
        z=matrix[i][i]
        for k in range(i,n+1):
            matrix[i][k]/=z
        for j in range(n):
            if j==i:continue
            z=matrix[j][i]
            for k in range(i,n+1):
                matrix[j][k]-=z*matrix[i][k]
    return [matrix[i][n] for i in range(n)]

def fit(X,Y,k):
    require(len(X)==len(Y) and k>=1,"INVALID_MODEL_DIMENSIONS")
    A=[[0.0]*k for i in range(k)]
    b=[0.0]*k
    for x,y in zip(X,Y):
        for j in range(k):
            b[j]+=x[j]*y
            for m in range(k):
                A[j][m]+=x[j]*x[m]
    beta=solve(A,b)
    if beta is None:return None
    sse=0.0
    for x,y in zip(X,Y):
        residual=y-sum(beta[j]*x[j] for j in range(k))
        sse+=residual*residual
    return {"beta":beta,"sse":sse}

def assemble(source,model):
    fields=model["fields"]
    joint=[r for r in source if
           all(r[k] is not None for k in fields) and
           r["seedsinit"] is not None and r["flyescape"] is not None]
    groups={}
    for r in joint:
        key=r["Site"]+"|"+str(int(r["Year"]))+"|"+r["Type"]
        groups.setdefault(key,[]).append(r)
    kept=[(key,rs) for key,rs in groups.items() if
          len(rs)>=2 and
          all(len({r[f] for r in rs})>=2 for f in fields)]
    eligible=[r for key,rr in kept for r in rr]
    require(eligible,"NO_SOURCE_ELIGIBLE_TRAIT_GROUPS")
    averages=[mean([r[f] for r in eligible]) for f in fields]
    scales=[math.sqrt(mean([(r[f]-averages[j])**2 for r in eligible]))
            for j,f in enumerate(fields)]
    require(all(s>0 for s in scales),"NO_ORIGINAL_FLORAL_TRAIT_VARIATION")
    result=[]
    for key,rr in kept:
        vals=[]
        for r in rr:
            xf=[(r[f]-averages[j])/scales[j]
                for j,f in enumerate(fields)]
            extra=(xf[0]**2 if model["extra"]=="square"
                   else xf[0]*xf[1])
            vals.append({
                "id":r["source_original_csv_row_id"],
                "reg":[*xf,extra],
                "Y":[math.log1p(r["seedsinit"]),r["flyescape"]],
                "year":int(r["Year"]),
            })
        result.append({
            "key":key,"n":len(vals),"year":int(rr[0]["Year"]),
            "original_file_rows":[r["source_original_csv_row_id"] for r in rr],
            "v":vals,
        })
    return {
        "n_source_joint_records":len(joint),
        "n_eligible_records":len(eligible),
        "n_groups":len(result),
        "mean_center":averages,
        "pop_sd":scales,
        "groups":result,
    }

def center(groups):
    design=[]
    out=[[],[]]
    ranges=[]
    start=0
    for g in groups:
        rr=g["v"]
        n=len(rr)
        k=len(rr[0]["reg"])
        means=[mean([r["reg"][j] for r in rr]) for j in range(k)]
        outcomes=[mean([r["Y"][j] for r in rr]) for j in range(2)]
        for r in rr:
            design.append([r["reg"][j]-means[j] for j in range(k)])
            for j in range(2):
                out[j].append(r["Y"][j]-outcomes[j])
        ranges.append({"from":start,"to":start+n})
        start+=n
    return {"design":design,"outcome":out,"range":ranges}

def evaluate(source,model,index):
    frame=assemble(source,model)
    data=center(frame["groups"])
    X=data["design"]
    k=len(model["fields"])+1
    observed=[]
    for Y in data["outcome"]:
        base=fit(X,Y,k-1)
        full=fit(X,Y,k)
        require(base is not None and full is not None and base["sse"]>0,
                "SECOND_ORDER_MODEL_NOT_IDENTIFIABLE")
        partial=(base["sse"]-full["sse"])/base["sse"]
        observed.append({
            "base_sse":base["sse"],
            "full_sse":full["sse"],
            "partial_R2":partial,
            "extra_standardized_beta":full["beta"][-1],
            "linear_coefficients":base["beta"],
            "added_term_variation_explained_fraction":partial,
        })
    seed=20261025+index
    rng=xorshift32(seed)
    counts=[0,0]
    boot_beta=[[],[]]
    boot_partial=[[],[]]
    fitted=[]
    residual=[]
    for response,Y in enumerate(data["outcome"]):
        coeff=fit(X,Y,k-1)["beta"]
        yy=[sum(c*x[j] for j,c in enumerate(coeff)) for x in X]
        fitted.append(yy)
        residual.append([value-yy[i] for i,value in enumerate(Y)])
    for rep in range(PERMUTATIONS):
        ids=list(range(len(X)))
        for span in data["range"]:
            start,end=span["from"],span["to"]
            for j in range(end-1,start,-1):
                slot=start+(next(rng)%(j-start+1))
                ids[j],ids[slot]=ids[slot],ids[j]
        for response in range(2):
            Y=[fitted[response][i]+residual[response][ids[i]]
               for i in range(len(ids))]
            base=fit(X,Y,k-1)
            full=fit(X,Y,k)
            require(base is not None and full is not None,
                    "PERMUTATION_MODEL_SINGULAR")
            partial=(base["sse"]-full["sse"])/base["sse"]
            if partial>=observed[response]["partial_R2"]-1e-12:
                counts[response]+=1
    for rep in range(BOOTSTRAPS):
        sample=[]
        for group in frame["groups"]:
            n=group["n"]
            sample.append({
                "key":group["key"],"year":group["year"],
                "n":n,
                "v":[group["v"][next(rng)%n] for j in range(n)],
            })
        z=center(sample)
        for response in range(2):
            base=fit(z["design"],z["outcome"][response],k-1)
            full=fit(z["design"],z["outcome"][response],k)
            if base is None or full is None:
                boot_beta[response].append(None)
                boot_partial[response].append(None)
                continue
            boot_beta[response].append(full["beta"][-1])
            boot_partial[response].append(
                (base["sse"]-full["sse"])/base["sse"])
    endpoints={}
    for j,label in enumerate(METRICS):
        bb=[z for z in boot_beta[j] if z is not None]
        rr=[z for z in boot_partial[j] if z is not None]
        require(bb and rr,"ALL_BOOTSTRAP_REPLICATES_SINGULAR")
        endpoints[label]={
            **observed[j],
            "freedman_lane_within_source_group_p":
                (counts[j]+1)/(PERMUTATIONS+1),
            "bootstrap_beta_ci95":[quantile(bb,.025),quantile(bb,.975)],
            "bootstrap_partial_R2_ci95":[quantile(rr,.025),quantile(rr,.975)],
            "n_valid_bootstrap":len(bb),
            "note":(
                "The bootstrap interval for nonnegative partial R2 is not "
                "a test of model advantage"),
        }
    cross={}
    for year in YEAR_RANGE:
        training=[g for g in frame["groups"] if g["year"]!=year]
        testing=[g for g in frame["groups"] if g["year"]==year]
        if not testing:
            continue
        tr=center(training)
        te=center(testing)
        current={
            "n_test_rows":len(te["design"]),
            "n_test_groups":len(testing),
        }
        for j,label in enumerate(METRICS):
            b=fit(tr["design"],tr["outcome"][j],k-1)
            f=fit(tr["design"],tr["outcome"][j],k)
            if b is None or f is None:
                current[label]=None
                continue
            baseline=extended=0.0
            for x,y in zip(te["design"],te["outcome"][j]):
                err_b=y-sum(t*x[ii] for ii,t in enumerate(b["beta"]))
                err_f=y-sum(t*x[ii] for ii,t in enumerate(f["beta"]))
                baseline+=err_b*err_b
                extended+=err_f*err_f
            current[label]={
                "base_test_SSE":baseline,
                "extended_test_SSE":extended,
            }
        cross[str(year)]=current
    cv={}
    for label in METRICS:
        baseline=extended=0.0
        for item in cross.values():
            if item[label] is not None:
                baseline+=item[label]["base_test_SSE"]
                extended+=item[label]["extended_test_SSE"]
        cv[label]={
            "heldout_year_relative_SSE_gain":(baseline-extended)/baseline,
            "positive_means_improvement":True,
            "base_SSE":baseline,
            "extra_SSE":extended,
        }
    return {
        "model_id":model["id"],
        "shape":model["extra"],
        "original_fields":model["fields"],
        "n_original_joint_records":frame["n_source_joint_records"],
        "n_eligible_original_source_rows":frame["n_eligible_records"],
        "n_original_site_year_type_groups":frame["n_groups"],
        "model_source_group_rows":[{
            "key":g["key"],"n":g["n"],"source_file_rows":g["original_file_rows"],
        } for g in frame["groups"]],
        "scaling":{
            "source_trait_raw_means":frame["mean_center"],
            "source_trait_raw_sd":frame["pop_sd"],
            "global_standardization_then_within_original_stratum_centering":True,
        },
        "estimand":(
            "Incremental source within-garden/year/type second-order trait "
            "term conditional on linear traits; descriptive, no random "
            "floral trait allocation"),
        "inference":{
            "permutations":PERMUTATIONS,
            "bootstrap_draws":BOOTSTRAPS,
            "seed":seed,
            "algorithm":"xorshift32",
            "reference":(
                "Freedman-Lane permuted main-effects regression residuals "
                "within source Site|Year|Type strata, paired stage outcomes "
                "shuffled with shared assignment"),
            "bootstrap":(
                "Original source rows resampled within Site|Year|Type "
                "fixed strata; no verified plant identity"),
            "test_statistic":(
                "nonnegative partial R2, refit base and full"),
        },
        "endpoints":endpoints,
        "leave_one_source_year_out_within_test_stratum":cross,
        "leave_one_year_out_summary":cv,
    }

def analyze(raw):
    source=read_author_csv(raw)
    results=[evaluate(source,model,i)
             for i,model in enumerate(MODELS)]
    ps=[]
    for result in results:
        for stage,item in result["endpoints"].items():
            ps.append({
                "key":result["model_id"]+"|"+stage,
                "p":item["freedman_lane_within_source_group_p"],
            })
    ps.sort(key=lambda a:(a["p"],a["key"]))
    qs={}
    q=1.0
    for i in range(len(ps)-1,-1,-1):
        q=min(q,ps[i]["p"]*len(ps)/(i+1))
        qs[ps[i]["key"]]=q
    return {
        "schema":"aza3_original_ipomopsis_higher_order_selection_source_JS_v1",
        "as_of":"2026-10-24",
        "source":{
            "archive":"https://zenodo.org/records/5083399",
            "paper":"10.1086/716740",
            "file":"Selection_by_pollinators_and_predators.csv",
            "author_file_md5":AUTHOR_MD5,
            "source_original_records":len(source),
            "original_missing_code":".",
            "source_original_plant_ids_unavailable":True,
        },
        "design":{
            "models":[m["id"] for m in MODELS],
            "primary_stage_responses":list(METRICS),
            "original_site_labels_are_not_distinct_physical_sites":True,
            "year_levels":list(YEAR_RANGE),
            "conditional_exploratory_not_randomized":True,
            "interpret_permutation_with_source_exchangeability_caveat":True,
            "structural_bridge_to_Cirsium_only_not_original_Cirsium_phenotypic_fitness_data":True,
        },
        "models":results,
        "multiplicity":{
            "six_source_exploratory_tests_BH_q":qs,
            "n_tests":len(ps),
        },
        "limits":{
            "no_verified_unique_plant_ids":True,
            "different_organism_than_Cirsium":True,
            "original_authors_multivariate_selection_results_not_refuted":True,
            "whole_outcome_observational_and_exposure_not_randomized":True,
            "partial_R2_effect_sizes_subject_to_selection_bias":True,
            "bootstrap_partial_R2_interval_not_evidence_against_zero":True,
            "quadratic_extra_term_is_not_an_evolutionary_reversion":True,
            "interaction_geometry_not_proven_Pedicularis_or_Cirsium_mechanism":True,
            "no_new_field_observations":0,
        },
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--csv",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    summary=analyze(a.csv.read_bytes())
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(
        json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "source_original_records":summary["source"]["source_original_records"],
        "inference":{
            m["model_id"]:{
                "n_source":m["n_eligible_original_source_rows"],
                "stage_partial_R2":{
                    k:v["partial_R2"]
                    for k,v in m["endpoints"].items()},
                "cross_year_gains":{
                    k:v["heldout_year_relative_SSE_gain"]
                    for k,v in m["leave_one_year_out_summary"].items()},
            } for m in summary["models"]
        },
    },ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
