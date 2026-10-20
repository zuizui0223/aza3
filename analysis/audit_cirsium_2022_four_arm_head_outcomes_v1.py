#!/usr/bin/env python3
"""Independent source-authenticated FOUR-ARM Cirsium arvense 2022 field audit.

Unit: original recorded treatment plant proxy (site × block × arm), aggregated
across sampled heads, followed by within-site-block arm differences.
The publisher did not provide a verified unique plant code and described a
first/second/third/fourth plant arm ordering within blocks. Hence comparisons
are DESCRIPTIVE; never label as a randomized causal treatment effect.

Author data: Yogita Karpate et al., 2025 Research Square preprint
DOI 10.21203/rs.3.rs-6634172/v1, original 2022 field study.
Archive DOI 10.5281/zenodo.15364152; author GitHub commit pinned below.

Pollinator exclusion = organza mesh bag (compound physical intervention).
Predator exclusion = stem adhesive barriers + manual predator removal.
Predator model = green plastic spider + real predator exclusion.
NO treatment manipulates Cirsium true spine anatomy or head orientation.
"""
from __future__ import annotations
from collections import Counter, defaultdict
import argparse
import csv
import hashlib
import io
import json
import math
import random
from pathlib import Path

SOURCE_REPO="Yogita25/Multitrophic_Interactions"
SOURCE_COMMIT="8a151c62ed069f775a16b9bc5fbdce209cbb8691"
SOURCE_GIT_BLOB="ca99c762162d15f006f0c73b4603ca71fd320a29"
PAPER_DOI="10.21203/rs.3.rs-6634172/v1"
ARCHIVE_DOI="10.5281/zenodo.15364152"
SITES=("campus","construction","ditch","rapseed")
ARMS=("Control","Pollinator Exclusion","Predator Exclusion","Predator Model")
NUM=("Flowers","Viable_Seeds","NViable_Seeds","Damaged_Seeds","Missing_Seeds","Insects")
REQUIRED=set(NUM)|{"Site_Code","Block_Code","Treatment_Code","Vial_Code"}
PAIRS=(("Pollinator Exclusion","Control"),("Predator Exclusion","Control"),
       ("Predator Model","Control"),("Predator Model","Predator Exclusion"),
       ("Pollinator Exclusion","Predator Exclusion"))
OUTCOMES=("viable_fraction","nonviable_fraction","observed_damaged_fraction",
          "signed_missing_fraction","inferred_eaten_signed_fraction",
          "emerged_insects_per_head","viable_per_head")

def check(flag,message):
    if not flag:raise ValueError(message)

def mean(xs):
    check(len(xs)>0,"EMPTY_COMPARISON")
    return sum(xs)/len(xs)

def quantile(xs,p):
    a=sorted(xs);q=(len(a)-1)*p;i=int(q)
    return a[i]*(1-(q-i))+a[min(i+1,len(a)-1)]*(q-i)

def parse_rows(data:bytes, strict=True):
    rawhash=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if strict:check(rawhash==SOURCE_GIT_BLOB,"ORIGINAL_CSV_GIT_BLOB_CHANGED")
    lines=csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
    check(REQUIRED.issubset(set(lines.fieldnames or [])),"ORIGINAL_SCHEMA_MISSING")
    rows=list(lines)
    if strict:check(len(rows)==615,"SOURCE_615_HEAD_ROWS_CHANGED")
    for r in rows:
        for col in NUM:
            check(r[col] is not None and r[col].strip()!="","MISSING_NUMERIC_"+col)
            try:n=float(r[col])
            except ValueError:raise ValueError("NONNUMERIC_"+col) from None
            check(math.isfinite(n) and n.is_integer() and
                  (col=="Missing_Seeds" or n>=0),
                  "INVALID_HEAD_COUNT_"+col)
        f=float(r["Flowers"])
        check(f>0,"FLOWERS_DENOMINATOR_ZERO")
        check(abs(f-sum(float(r[n]) for n in
              ("Viable_Seeds","NViable_Seeds","Damaged_Seeds","Missing_Seeds")))<1e-9,
              "HEAD_SEED_ACCOUNTING_IDENTITY_BROKEN")
    return rows,rawhash

def group_estimate(rows):
    check(bool(rows),"EMPTY_PLANT_PROXY")
    s=lambda name:sum(int(float(r[name])) for r in rows)
    flowers=s("Flowers");head=len(rows)
    check(flowers>0,"NO_FLORETS")
    a={
      "n_sampled_heads":head,
      "viable_fraction":s("Viable_Seeds")/flowers,
      "nonviable_fraction":s("NViable_Seeds")/flowers,
      "observed_damaged_fraction":s("Damaged_Seeds")/flowers,
      "signed_missing_fraction":s("Missing_Seeds")/flowers,
      "inferred_eaten_signed_fraction":(s("Damaged_Seeds")+s("Missing_Seeds"))/flowers,
      "emerged_insects_per_head":s("Insects")/head,
      "viable_per_head":s("Viable_Seeds")/head,
    }
    check(math.isclose(sum(a[k] for k in ("viable_fraction","nonviable_fraction",
         "observed_damaged_fraction","signed_missing_fraction")),1,abs_tol=1e-10),
         "FRACTION_BOOKKEEPING_BROKEN")
    check(math.isclose(a["observed_damaged_fraction"]+a["signed_missing_fraction"],
        a["inferred_eaten_signed_fraction"],abs_tol=1e-10),
        "SIGNED_RESIDUAL_IDENTITY_BROKEN")
    return a

def sign_counts(vals):
    eps=1e-12
    return {"positive":sum(z>eps for z in vals),
            "negative":sum(z< -eps for z in vals),
            "zero":sum(abs(z)<=eps for z in vals)}

def prep(rows,strict):
    keep=[r for r in rows if r["Site_Code"] in SITES and r["Treatment_Code"] in ARMS]
    proxy=defaultdict(list)
    for r in keep:
        check(r["Block_Code"].strip()!="","MISSING_BLOCK")
        proxy[(r["Site_Code"],r["Block_Code"],r["Treatment_Code"])].append(r)
    blocks=defaultdict(dict)
    for (site,block,arm),rs in proxy.items():
        blocks[(site,block)][arm]=group_estimate(rs)
    complete=[{"site":s,"block":b,"arms":z} for (s,b),z in sorted(blocks.items())
              if all(a in z for a in ARMS)]
    if strict:
        check(len(keep)==612 and len(proxy)==164 and len(blocks)==45,
              "FIXED_SOURCE_UNIT_COUNT_DRIFT")
        check(len(complete)==36,"FIXED_COMPLETE_FOUR_ARM_BLOCKS_DRIFT")
        check(Counter(d["site"] for d in complete)=={
             "campus":12,"construction":9,"ditch":8,"rapseed":7},
             "FIXED_COMPLETE_BLOCK_SITE_DISTRIBUTION_DRIFT")
        check(Counter(r["Treatment_Code"] for r in keep)=={
             "Control":164,"Pollinator Exclusion":142,
             "Predator Exclusion":154,"Predator Model":152},
             "FIXED_TREATMENT_HEAD_COUNT_DRIFT")
    return keep,proxy,blocks,complete

def contrast(pair,arm_a,arm_b):
    x=pair["arms"]
    return {k:x[arm_a][k]-x[arm_b][k] for k in OUTCOMES}

def summarize_contrast(blocks, arm_a,arm_b):
    rows=[{"site":site,"block":block,"arms":units}
          for (site,block),units in sorted(blocks.items())
          if arm_a in units and arm_b in units]
    check(rows,"NO_MATCHED_BLOCKS")
    r={"arm_a":arm_a,"arm_b":arm_b,"n_matched_blocks":len(rows),
       "n_matched_blocks_per_site":{site:sum(v["site"]==site for v in rows)
                                    for site in SITES},
       "outcomes":{}}
    for metric in OUTCOMES:
        values=[contrast(z,arm_a,arm_b)[metric] for z in rows]
        rr={"mean_a_minus_b":mean(values),"sign":sign_counts(values)}
        rr["site_means"]={s:mean([z for z,v in zip(values,rows) if v["site"]==s])
                         for s in SITES if any(v["site"]==s for v in rows)}
        rr["leave_one_site_out"]={
            s:mean([z for z,v in zip(values,rows) if v["site"]!=s])
            for s in SITES if any(v["site"]!=s for v in rows)}
        # With one matched block, deletion leaves no comparable units. Keep
        # this sensitivity undefined; missing is not a biological zero.
        d=[mean(values[:i]+values[i+1:]) for i in range(len(values))] if len(values)>1 else []
        rr["leave_one_block_out_range"]=[min(d),max(d)] if d else None
        r["outcomes"][metric]=rr
    return r

def complete_block_summary(complete,bootstrap=5000,seed=20261020):
    check(bootstrap>=100,"BOOTSTRAP_TOO_FEW")
    siteblocks={s:[z for z in complete if z["site"]==s] for s in SITES}
    check(all(siteblocks.values()),"REQUIRE_ALL_FOUR_ORIGINAL_STUDY_SITES")
    d={
      "n_all_four_arms":len(complete),
      "site_n":{s:len(v) for s,v in siteblocks.items()},
      "contrasts":{},
      "bootstrap_replicates":bootstrap,
      "bootstrap_seed":seed,
      "sampling_scope":"Resampling site-matched blocks within four ORIGINAL fixed sites, descriptive only; not four independent studies or causal randomization."
    }
    randomgen=random.Random(seed)
    for arm in ARMS[1:]:
        diffs=[contrast(z,arm,"Control")["viable_fraction"] for z in complete]
        d["contrasts"][arm]={
          "mean_viable_fraction_minus_control":mean(diffs),
          "sign":sign_counts(diffs),
          "site_means":{s:mean([contrast(z,arm,"Control")["viable_fraction"]
                                for z in siteblocks[s]]) for s in SITES},
        }
    same=[contrast(z,"Pollinator Exclusion","Predator Exclusion")["viable_fraction"]
          for z in complete]
    d["contrasts"]["Pollinator Exclusion minus Predator Exclusion"]={
      "mean_viable_fraction_difference":mean(same),
      "sign":sign_counts(same),
      "site_means":{s:mean([contrast(z,"Pollinator Exclusion","Predator Exclusion")["viable_fraction"]
                          for z in siteblocks[s]]) for s in SITES},
    }
    distributions={k:[] for k in d["contrasts"]}
    for _ in range(bootstrap):
        draw=[xs[randomgen.randrange(len(xs))]
              for s in SITES for xs in [siteblocks[s]] for j in range(len(xs))]
        for arm in ARMS[1:]:
            distributions[arm].append(mean([contrast(z,arm,"Control")["viable_fraction"] for z in draw]))
        distributions["Pollinator Exclusion minus Predator Exclusion"].append(
          mean([contrast(z,"Pollinator Exclusion","Predator Exclusion")["viable_fraction"] for z in draw]))
    for k,ds in distributions.items():
        d["contrasts"][k]["conditional_site_stratified_bootstrap_interval_95"]=[
            quantile(ds,.025),quantile(ds,.975)]
    return d

def synthetic_controls():
    def row(site,block,arm,flower,viable,damaged=0,nonviable=0,insects=0):
        return {"Site_Code":site,"Block_Code":block,"Treatment_Code":arm,"Vial_Code":"1",
                "Flowers":str(flower),"Viable_Seeds":str(viable),
                "NViable_Seeds":str(nonviable),"Damaged_Seeds":str(damaged),
                "Missing_Seeds":str(flower-viable-damaged-nonviable),"Insects":str(insects)}
    a=[row("campus","A","Control",10,8,1),
       row("campus","A","Pollinator Exclusion",10,2,1),
       row("campus","A","Predator Exclusion",10,7,2),
       row("campus","A","Predator Model",10,9,0),
       row("campus","B","Control",10,8,0),
       row("campus","B","Pollinator Exclusion",10,4,0),
       row("campus","B","Predator Exclusion",10,7,0)]
    keep,proxy,blocks,complete=prep(a,False)
    assert len(complete)==1 and len(blocks)==2
    assert math.isclose(summarize_contrast(blocks,"Pollinator Exclusion","Control")[
      "outcomes"]["viable_fraction"]["mean_a_minus_b"],-.5)
    assert len(complete)==1 and complete[0]["block"]=="A"
    assert group_estimate([row("campus","C","Control",10,9,2)])[
      "signed_missing_fraction"]<0
    assert sign_counts([0,1,-1])=={"positive":1,"negative":1,"zero":1}
    return "PASS_SYNTHETIC_ACCOUNTING_AND_MISSING_BLOCK_CONTROLS"

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seeds",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=5000)
    args=p.parse_args()
    check(100<=args.bootstrap<=30000,"BOOTSTRAP_BUDGET_OUT_OF_RANGE")
    raw=args.seeds.read_bytes()
    rows,git_blob=parse_rows(raw,strict=True)
    keep,proxy,blocks,complete=prep(rows,True)
    r={
      "version":"cirsium_2022_four_arm_source_marginal_contrasts_v1",
      "date":"2026-10-20",
      "source":{"study_year":2022,"study_preprint_year":2025,"paper_doi":PAPER_DOI,
                "author_archive_doi":ARCHIVE_DOI,"author_repo":SOURCE_REPO,
                "author_commit":SOURCE_COMMIT,"verified_seeds_csv_git_blob":git_blob,
                "source_head_rows":len(rows),"eligible_head_rows":len(keep),
                "head_rows_excluded_due_to_missing_arm_or_site":len(rows)-len(keep),
                "inferred_proxy_plants":len(proxy),"original_distinct_site_blocks":len(blocks),
                "raw_negative_missing_rows":sum(int(float(z["Missing_Seeds"]))<0 for z in rows),
                "raw_negative_signed_eaten_rows":sum(int(z["Damaged_Seeds"])+int(z["Missing_Seeds"])<0 for z in rows)},
      "treated_head_rows":dict(Counter(r["Treatment_Code"] for r in keep)),
      "paired_block_contrasts":{
         arm_a+" MINUS "+arm_b:summarize_contrast(blocks,arm_a,arm_b)
         for arm_a,arm_b in PAIRS},
      "identical_complete_four_arm_blocks":complete_block_summary(complete,args.bootstrap),
      "synthetic_controls":synthetic_controls(),
      "claim_limit":{
         "source_preprint_not_peer_reviewed_verified":"Treat the versioned 2025 source as a preprint until journal status independently checked.",
         "proxy_plant_id":"Site x block x treatment, not a verified publication-level stable unique plant ID",
         "design_randomization_caution":"Source describes first control, second bagged, third predator-excluded, fourth predator-model plants; independently randomized assignment to plant identity is NOT verified.",
         "pollinator_exclusion":"Mesh bag also changes physical exposure, humidity and other insects; NOT isolated pollen effect.",
         "predator_exclusion":"Stem adhesive barriers and manual predator removal; NOT true botanical spine removal.",
         "predator_model":"Plastic spider plus predator exclusion, NOT predator addition or field trophic effect.",
         "head_units_not_randomized_replicates":True,
         "seed_count_measure":"Published Viable_Seeds author scoring, not a new germination test",
         "signed_missing":"Flowers - Viable - Nonviable - Damaged is accounting residual, may be negative; inferred Eaten is not independently observed consumption.",
         "not_measured":["head gravity orientation","true involucral spine length","pre-entry feeding insect route",
                         "pollen delivery assay","verified predator killing before seed damage",
                         "heritable evolutionary reexpression","same-head fully determined lifetime plant fitness"],
         "does_not_establish_zero_predator_effect":True,
         "new_Cirsium_individual_field_experiments_completed":0
       }}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":"ORIGINAL_2022_FOUR_ARM_REPLAY",
      "head_rows":len(rows),"matched_complete_blocks":len(complete),
      "observed":{k:v["mean_viable_fraction_minus_control"] for k,v
                  in r["identical_complete_four_arm_blocks"]["contrasts"].items()
                  if "mean_viable_fraction_minus_control" in v}},indent=2))
if __name__=="__main__":main()
