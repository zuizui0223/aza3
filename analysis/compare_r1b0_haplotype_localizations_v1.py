#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, statistics
from pathlib import Path

def read_table(path: Path):
    out={}
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["locus"]]=r
    return out

def fnum(x):
    try:
        return float(x)
    except Exception:
        return None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--hap1",required=True,type=Path)
    ap.add_argument("--hap2",required=True,type=Path)
    ap.add_argument("--taxon",required=True)
    ap.add_argument("--query-set-id",required=True)
    ap.add_argument("--output-json",required=True,type=Path)
    ap.add_argument("--output-tsv",required=True,type=Path)
    args=ap.parse_args()

    a=read_table(args.hap1); b=read_table(args.hap2)
    union=sorted(set(a)|set(b))
    common=sorted(set(a)&set(b))
    rows=[]
    for locus in union:
        ra=a.get(locus); rb=b.get(locus)
        ca="NO_PAF_RECORD" if ra is None else ra["placement_class"]
        cb="NO_PAF_RECORD" if rb is None else rb["placement_class"]
        ua=ca=="UNIQUE_HIGH_CONF"; ub=cb=="UNIQUE_HIGH_CONF"
        qa=None if ra is None else fnum(ra.get("query_coverage"))
        qb=None if rb is None else fnum(rb.get("query_coverage"))
        ma=None if ra is None else fnum(ra.get("match_fraction_of_query"))
        mb=None if rb is None else fnum(rb.get("match_fraction_of_query"))
        rows.append({
            "locus":locus,
            "hap1_class":ca,
            "hap2_class":cb,
            "hap1_unique":ua,
            "hap2_unique":ub,
            "unique_status_agree":ua==ub,
            "hap1_qcov":"" if qa is None else qa,
            "hap2_qcov":"" if qb is None else qb,
            "abs_delta_qcov":"" if qa is None or qb is None else abs(qa-qb),
            "hap1_match_fraction":"" if ma is None else ma,
            "hap2_match_fraction":"" if mb is None else mb,
            "abs_delta_match_fraction":"" if ma is None or mb is None else abs(ma-mb),
        })

    common_rows=[r for r in rows if r["locus"] in set(common)]
    if not common_rows:
        raise SystemExit("no common loci")

    ua=sum(r["hap1_unique"] for r in common_rows)
    ub=sum(r["hap2_unique"] for r in common_rows)
    agree=sum(r["unique_status_agree"] for r in common_rows)/len(common_rows)
    fa=ua/len(common_rows); fb=ub/len(common_rows)
    frac_diff=abs(fa-fb)
    qd=[float(r["abs_delta_qcov"]) for r in common_rows if r["abs_delta_qcov"]!=""]
    md=[float(r["abs_delta_match_fraction"]) for r in common_rows if r["abs_delta_match_fraction"]!=""]

    if agree>=0.95 and frac_diff<=0.03:
        cls="LOW_HAPLOTYPE_REFERENCE_NOISE"
    elif agree>=0.90 and frac_diff<=0.07:
        cls="MODERATE_HAPLOTYPE_REFERENCE_NOISE"
    else:
        cls="HIGH_HAPLOTYPE_REFERENCE_NOISE"

    args.output_tsv.parent.mkdir(parents=True,exist_ok=True)
    with args.output_tsv.open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader(); w.writerows(rows)

    out={
        "result_version":"r1b0_within_individual_haplotype_control_v1",
        "taxon":args.taxon,
        "query_set_id":args.query_set_id,
        "union_loci":len(union),
        "common_loci":len(common),
        "hap1_unique_fraction":fa,
        "hap2_unique_fraction":fb,
        "absolute_unique_fraction_difference":frac_diff,
        "unique_status_agreement_fraction":agree,
        "both_unique":sum(r["hap1_unique"] and r["hap2_unique"] for r in common_rows),
        "hap1_unique_hap2_not":sum(r["hap1_unique"] and not r["hap2_unique"] for r in common_rows),
        "hap1_not_hap2_unique":sum((not r["hap1_unique"]) and r["hap2_unique"] for r in common_rows),
        "median_absolute_query_coverage_difference":statistics.median(qd) if qd else None,
        "median_absolute_match_fraction_difference":statistics.median(md) if md else None,
        "classification":cls,
        "classification_rule":{
            "LOW_HAPLOTYPE_REFERENCE_NOISE":"unique-status agreement >=0.95 and absolute unique-fraction difference <=0.03",
            "MODERATE_HAPLOTYPE_REFERENCE_NOISE":"otherwise, agreement >=0.90 and absolute unique-fraction difference <=0.07",
            "HIGH_HAPLOTYPE_REFERENCE_NOISE":"otherwise"
        },
        "claim_boundary":"Same-individual reference-haplotype placement noise only; not synteny proof, local ancestry, recombination, haplotype age, SV reuse, genomic modularity or genotype-phenotype association."
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
