#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json
from collections import defaultdict
from pathlib import Path

def parse_paf(path: Path):
    byq=defaultdict(list)
    with path.open(encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            f=raw.rstrip("\n").split("\t")
            if len(f)<12:
                continue
            qname=f[0]; qlen=int(f[1]); qs=int(f[2]); qe=int(f[3])
            strand=f[4]; tname=f[5]; tlen=int(f[6]); ts=int(f[7]); te=int(f[8])
            nmatch=int(f[9]); alen=int(f[10]); mapq=int(f[11])
            tags={}
            for x in f[12:]:
                p=x.split(":",2)
                if len(p)==3:
                    tags[p[0]]=(p[1],p[2])
            byq[qname].append({
                "qname":qname,"qlen":qlen,"qstart":qs,"qend":qe,
                "strand":strand,"tname":tname,"tlen":tlen,"tstart":ts,"tend":te,
                "nmatch":nmatch,"alen":alen,"mapq":mapq,
                "query_coverage":(qe-qs)/qlen if qlen else 0.0,
                "match_fraction_of_query":nmatch/qlen if qlen else 0.0,
                "dv":float(tags["dv"][1]) if "dv" in tags and tags["dv"][0]=="f" else None,
                "de":float(tags["de"][1]) if "de" in tags and tags["de"][0]=="f" else None,
                "tp":tags.get("tp",(None,None))[1],
            })
    return byq

def locus_id(qname: str)->str:
    return qname.split("|",1)[0]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--paf",required=True,type=Path)
    ap.add_argument("--reference-id",required=True)
    ap.add_argument("--output-tsv",required=True,type=Path)
    ap.add_argument("--output-json",required=True,type=Path)
    ap.add_argument("--min-qcov",type=float,default=0.70)
    ap.add_argument("--min-mapq",type=int,default=20)
    args=ap.parse_args()

    byq=parse_paf(args.paf)
    rows=[]
    counts=defaultdict(int)
    for qname,hits in sorted(byq.items(),key=lambda kv:locus_id(kv[0])):
        hits=sorted(hits,key=lambda x:(x["mapq"],x["match_fraction_of_query"],x["nmatch"]),reverse=True)
        top=hits[0]
        high=[h for h in hits if h["query_coverage"]>=args.min_qcov and h["mapq"]>=args.min_mapq]
        if len(high)==1:
            cls="UNIQUE_HIGH_CONF"
        elif len(high)>1:
            cls="MULTIPLE_HIGH_CONF"
        elif top["query_coverage"]>=0.50 and top["mapq"]>=5:
            cls="LOW_CONF"
        else:
            cls="NO_USEFUL_HIT"
        counts[cls]+=1
        second=hits[1] if len(hits)>1 else None
        rows.append({
            "locus":locus_id(qname),
            "query_name":qname,
            "reference_id":args.reference_id,
            "placement_class":cls,
            "top_target":top["tname"],
            "top_start":top["tstart"],
            "top_end":top["tend"],
            "strand":top["strand"],
            "query_length":top["qlen"],
            "query_coverage":top["query_coverage"],
            "match_fraction_of_query":top["match_fraction_of_query"],
            "mapq":top["mapq"],
            "dv":top["dv"],
            "de":top["de"],
            "secondary_high_conf_count":max(0,len(high)-1),
            "second_target":None if second is None else second["tname"],
            "second_mapq":None if second is None else second["mapq"],
        })

    args.output_tsv.parent.mkdir(parents=True,exist_ok=True)
    fieldnames=list(rows[0]) if rows else [
        "locus","query_name","reference_id","placement_class","top_target","top_start","top_end",
        "strand","query_length","query_coverage","match_fraction_of_query","mapq","dv","de",
        "secondary_high_conf_count","second_target","second_mapq"
    ]
    with args.output_tsv.open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fieldnames,delimiter="\t")
        w.writeheader(); w.writerows(rows)

    total=len(rows)
    out={
        "status":"ok",
        "reference_id":args.reference_id,
        "paf":str(args.paf),
        "min_query_coverage":args.min_qcov,
        "min_mapq":args.min_mapq,
        "queries_with_any_paf_hit":total,
        "placement_counts":dict(counts),
        "unique_high_conf_fraction":counts["UNIQUE_HIGH_CONF"]/total if total else None,
        "multiple_high_conf_fraction":counts["MULTIPLE_HIGH_CONF"]/total if total else None,
        "claim_boundary":"Spliced coding-locus placement viability only. Not local ancestry, synteny proof, haplotype age, recombination, SV reuse, or genotype-phenotype association."
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
