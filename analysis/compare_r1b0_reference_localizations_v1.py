#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
from itertools import combinations

REFS=("heterophyllum_hap1","dissectum_hap1","nipponicum")

def read_table(path: Path):
    out={}
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["locus"]]=r
    return out

def phi_binary(a,b):
    # simple agreement is more interpretable here because all loci are the same frozen set
    if not a:
        return None
    return sum(x==y for x,y in zip(a,b))/len(a)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--heterophyllum",required=True,type=Path)
    ap.add_argument("--dissectum",required=True,type=Path)
    ap.add_argument("--nipponicum",required=True,type=Path)
    ap.add_argument("--output-json",required=True,type=Path)
    ap.add_argument("--output-tsv",required=True,type=Path)
    ap.add_argument("--query-set-id",required=True)
    args=ap.parse_args()

    tables={
        "heterophyllum_hap1":read_table(args.heterophyllum),
        "dissectum_hap1":read_table(args.dissectum),
        "nipponicum":read_table(args.nipponicum),
    }
    union=sorted(set().union(*(set(t) for t in tables.values())))
    common=sorted(set.intersection(*(set(t) for t in tables.values()))) if tables else []

    rows=[]
    for locus in union:
        row={"locus":locus}
        n_unique=0
        n_useful=0
        for ref in REFS:
            r=tables[ref].get(locus)
            cls=None if r is None else r["placement_class"]
            row[f"{ref}_class"]=cls or "NO_PAF_RECORD"
            row[f"{ref}_target"]="" if r is None else r["top_target"]
            row[f"{ref}_qcov"]="" if r is None else r["query_coverage"]
            row[f"{ref}_mapq"]="" if r is None else r["mapq"]
            if cls=="UNIQUE_HIGH_CONF":
                n_unique += 1
            if cls in {"UNIQUE_HIGH_CONF","MULTIPLE_HIGH_CONF","LOW_CONF"}:
                n_useful += 1
        row["n_unique_high_conf"]=n_unique
        row["n_useful_refs"]=n_useful
        row["unique_all_three"]=n_unique==3
        row["unique_at_least_two"]=n_unique>=2
        rows.append(row)

    args.output_tsv.parent.mkdir(parents=True,exist_ok=True)
    fields=list(rows[0]) if rows else ["locus"]
    with args.output_tsv.open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t")
        w.writeheader(); w.writerows(rows)

    ref_summaries={}
    for ref,t in tables.items():
        classes=[r["placement_class"] for r in t.values()]
        ref_summaries[ref]={
            "loci_with_paf_record":len(t),
            "unique_high_conf":sum(c=="UNIQUE_HIGH_CONF" for c in classes),
            "multiple_high_conf":sum(c=="MULTIPLE_HIGH_CONF" for c in classes),
            "low_conf":sum(c=="LOW_CONF" for c in classes),
            "no_useful_hit":sum(c=="NO_USEFUL_HIT" for c in classes),
            "unique_high_conf_fraction":(sum(c=="UNIQUE_HIGH_CONF" for c in classes)/len(t) if t else None),
        }

    common_rows=[r for r in rows if r["locus"] in set(common)]
    pairwise={}
    for a,b in combinations(REFS,2):
        xa=[r[f"{a}_class"]=="UNIQUE_HIGH_CONF" for r in common_rows]
        xb=[r[f"{b}_class"]=="UNIQUE_HIGH_CONF" for r in common_rows]
        pairwise[f"{a}__{b}"]={
            "unique_status_agreement_fraction":phi_binary(xa,xb),
            "both_unique_high_conf":sum(x and y for x,y in zip(xa,xb)),
            "a_unique_b_not":sum(x and not y for x,y in zip(xa,xb)),
            "a_not_b_unique":sum((not x) and y for x,y in zip(xa,xb)),
        }

    out={
        "result_version":"r1b0_reference_localization_cross_reference_v1",
        "query_set_id":args.query_set_id,
        "reference_ids":list(REFS),
        "union_loci":len(union),
        "common_loci_with_paf_records_all_three":len(common),
        "common_unique_all_three":sum(r["unique_all_three"] for r in common_rows),
        "common_unique_all_three_fraction":(
            sum(r["unique_all_three"] for r in common_rows)/len(common_rows) if common_rows else None
        ),
        "common_unique_at_least_two":sum(r["unique_at_least_two"] for r in common_rows),
        "common_unique_at_least_two_fraction":(
            sum(r["unique_at_least_two"] for r in common_rows)/len(common_rows) if common_rows else None
        ),
        "reference_summaries":ref_summaries,
        "pairwise_unique_status":pairwise,
        "claim_boundary":(
            "Cross-reference spliced coding-locus placement viability only. "
            "This is not synteny proof, local ancestry, haplotype age, recombination, "
            "structural-variant reuse, genomic modularity or genotype-phenotype association."
        ),
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
