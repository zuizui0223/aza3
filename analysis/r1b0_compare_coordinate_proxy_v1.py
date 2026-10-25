#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
from collections import Counter,defaultdict

def read(path):
    return {r["locus"]:r for r in csv.DictReader(open(path))}

def truth(v):
    return str(v).lower() in {"true","1","yes"}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--heterophyllum",required=True)
    ap.add_argument("--dissectum",required=True)
    ap.add_argument("--nipponicum",required=True)
    ap.add_argument("--output",required=True)
    args=ap.parse_args()
    H,D,N=read(args.heterophyllum),read(args.dissectum),read(args.nipponicum)
    loci=sorted(set(H)&set(D)&set(N))
    # Learn coarse chromosome pairing H->D from loci with one consensus refseq in both.
    pairs=Counter()
    for g in loci:
        if H[g]["consensus_refseq"] and D[g]["consensus_refseq"] and truth(H[g]["all_high_conf_same_refseq"]) and truth(D[g]["all_high_conf_same_refseq"]):
            pairs[(H[g]["consensus_refseq"],D[g]["consensus_refseq"])]+=1
    htot=defaultdict(Counter); dtot=defaultdict(Counter)
    for (h,d),n in pairs.items():
        htot[h][d]+=n; dtot[d][h]+=n
    hmajor={h:c.most_common(1)[0] for h,c in htot.items()}
    dmajor={d:c.most_common(1)[0] for d,c in dtot.items()}
    reciprocal={}
    for h,(d,n) in hmajor.items():
        reciprocal[h]=d if d in dmajor and dmajor[d][0]==h else None

    rows=[]
    for g in loci:
        h=H[g];d=D[g];n=N[g]
        rec=bool(h["consensus_refseq"] and d["consensus_refseq"] and reciprocal.get(h["consensus_refseq"])==d["consensus_refseq"])
        rows.append({
            "locus":g,
            "heterophyllum_refseq":h["consensus_refseq"],
            "dissectum_refseq":d["consensus_refseq"],
            "nipponicum_refseq":n["consensus_refseq"],
            "heterophyllum_same_refseq":h["all_high_conf_same_refseq"],
            "dissectum_same_refseq":d["all_high_conf_same_refseq"],
            "nipponicum_same_refseq":n["all_high_conf_same_refseq"],
            "reciprocal_majority_H_D":rec
        })
    summary={
        "status":"ok",
        "loci":len(rows),
        "H_D_reciprocal_majority_loci":sum(r["reciprocal_majority_H_D"] for r in rows),
        "H_D_reciprocal_majority_fraction":sum(r["reciprocal_majority_H_D"] for r in rows)/len(rows) if rows else None,
        "H_distinct_refseqs_with_majority_pair":len(hmajor),
        "D_distinct_refseqs_with_majority_pair":len(dmajor),
        "reciprocal_chromosome_pairs":sum(v is not None for v in reciprocal.values()),
        "note":"This is a target-locus coordinate proxy only. It does not infer C. sieboldii local ancestry or genomic modularity."
    }
    Path(args.output).write_text(json.dumps({"summary":summary,"reciprocal_map":reciprocal,"rows":rows},indent=2)+"\n")
    print(json.dumps(summary,indent=2))
if __name__=="__main__":
    main()
