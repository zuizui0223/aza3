#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math, statistics
from pathlib import Path

def read_table(path: Path):
    out={}
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["locus"]]=r
    return out

def f(x):
    if x in (None,"","None","nan","NA"):
        return None
    try:
        return float(x)
    except Exception:
        return None

def median_abs_delta(rows, key):
    vals=[]
    for a,b in rows:
        x=f(a.get(key)); y=f(b.get(key))
        if x is not None and y is not None:
            vals.append(abs(x-y))
    return statistics.median(vals) if vals else None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--hap1", required=True, type=Path)
    ap.add_argument("--hap2", required=True, type=Path)
    ap.add_argument("--taxon", required=True)
    ap.add_argument("--hap1-id", required=True)
    ap.add_argument("--hap2-id", required=True)
    ap.add_argument("--output-json", required=True, type=Path)
    ap.add_argument("--output-tsv", required=True, type=Path)
    args=ap.parse_args()

    h1=read_table(args.hap1); h2=read_table(args.hap2)
    union=sorted(set(h1)|set(h2))
    common=sorted(set(h1)&set(h2))
    paired=[(h1[x],h2[x]) for x in common]

    rows=[]
    for locus in union:
        a=h1.get(locus); b=h2.get(locus)
        ca="NO_PAF_RECORD" if a is None else a["placement_class"]
        cb="NO_PAF_RECORD" if b is None else b["placement_class"]
        ua=ca=="UNIQUE_HIGH_CONF"; ub=cb=="UNIQUE_HIGH_CONF"
        rows.append({
            "locus":locus,
            "hap1_class":ca,
            "hap2_class":cb,
            "hap1_unique":ua,
            "hap2_unique":ub,
            "same_unique_status":ua==ub,
            "hap1_qcov":"" if a is None else a["query_coverage"],
            "hap2_qcov":"" if b is None else b["query_coverage"],
            "hap1_match_fraction":"" if a is None else a["match_fraction_of_query"],
            "hap2_match_fraction":"" if b is None else b["match_fraction_of_query"],
            "hap1_mapq":"" if a is None else a["mapq"],
            "hap2_mapq":"" if b is None else b["mapq"],
        })

    common_rows=[r for r in rows if r["locus"] in set(common)]
    n=len(common_rows)
    both_unique=sum(r["hap1_unique"] and r["hap2_unique"] for r in common_rows)
    h1_only=sum(r["hap1_unique"] and not r["hap2_unique"] for r in common_rows)
    h2_only=sum((not r["hap1_unique"]) and r["hap2_unique"] for r in common_rows)
    neither=sum((not r["hap1_unique"]) and (not r["hap2_unique"]) for r in common_rows)
    agreement=sum(r["same_unique_status"] for r in common_rows)/n if n else None

    # use original paired dictionaries for continuous deltas
    out={
        "result_version":"r1b0_same_individual_haplotype_control_v1",
        "taxon":args.taxon,
        "hap1_reference":args.hap1_id,
        "hap2_reference":args.hap2_id,
        "union_loci":len(union),
        "common_loci_with_paf_records_both":len(common),
        "unique_status_agreement_fraction":agreement,
        "both_unique_high_conf":both_unique,
        "hap1_unique_hap2_not":h1_only,
        "hap1_not_hap2_unique":h2_only,
        "neither_unique":neither,
        "both_unique_fraction_of_common":both_unique/n if n else None,
        "median_abs_query_coverage_delta":median_abs_delta(paired,"query_coverage"),
        "median_abs_match_fraction_delta":median_abs_delta(paired,"match_fraction_of_query"),
        "median_abs_mapq_delta":median_abs_delta(paired,"mapq"),
        "claim_boundary":(
            "Same-individual haplotype-reference placement-noise control only. "
            "Different scaffold/chromosome identifiers are not treated as coordinate discordance. "
            "No synteny, local ancestry, haplotype age, recombination, SV reuse, genomic modularity "
            "or genotype-phenotype inference is authorized."
        ),
    }

    args.output_tsv.parent.mkdir(parents=True,exist_ok=True)
    with args.output_tsv.open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["locus"],delimiter="\t")
        w.writeheader(); w.writerows(rows)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
