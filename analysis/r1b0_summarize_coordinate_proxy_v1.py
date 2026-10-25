#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,statistics
from pathlib import Path
from collections import defaultdict,Counter

def read_fasta_headers(path):
    out=[]
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                q=line[1:].split()[0].strip()
                src,locus=q.split("__",1)
                out.append((q,src,locus))
    return out

def parse_paf(path):
    rows=[]
    with open(path) as f:
        for line in f:
            if not line.strip(): continue
            v=line.rstrip("\n").split("\t")
            if len(v)<12: continue
            q,qlen,qs,qe,strand,t,tlen,ts,te,nmatch,alen,mapq=v[:12]
            tags={}
            for z in v[12:]:
                p=z.split(":",2)
                if len(p)==3: tags[p[0]]=p[2]
            rows.append({
                "q":q,"qlen":int(qlen),"qs":int(qs),"qe":int(qe),"strand":strand,
                "t":t,"tlen":int(tlen),"ts":int(ts),"te":int(te),
                "nmatch":int(nmatch),"alen":int(alen),"mapq":int(mapq),
                "tp":tags.get("tp",""),"AS":int(tags["AS"]) if tags.get("AS","").lstrip("-").isdigit() else None
            })
    return rows

def med(x):
    return None if not x else statistics.median(x)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--queries",required=True)
    ap.add_argument("--paf",required=True)
    ap.add_argument("--reference-id",required=True)
    ap.add_argument("--out-csv",required=True)
    ap.add_argument("--out-json",required=True)
    args=ap.parse_args()

    qmeta={q:(src,locus) for q,src,locus in read_fasta_headers(args.queries)}
    alns=defaultdict(list)
    for r in parse_paf(args.paf):
        if r["q"] in qmeta:
            alns[r["q"]].append(r)

    best={}
    for q,rs in alns.items():
        # Prefer minimap2 primary; then MAPQ, alignment score/nmatch and span.
        rs2=sorted(rs,key=lambda r:(r["tp"]=="P",r["mapq"],r["AS"] if r["AS"] is not None else r["nmatch"],r["alen"]),reverse=True)
        best[q]=rs2[0]

    byloc=defaultdict(list)
    for q,(src,locus) in qmeta.items():
        r=best.get(q)
        if r is None:
            byloc[locus].append({"q":q,"src":src,"mapped":False})
        else:
            cov=(r["qe"]-r["qs"])/r["qlen"] if r["qlen"] else 0
            ident=r["nmatch"]/r["alen"] if r["alen"] else 0
            byloc[locus].append({
                "q":q,"src":src,"mapped":True,"t":r["t"],"ts":r["ts"],"te":r["te"],
                "mapq":r["mapq"],"cov":cov,"ident":ident,"strand":r["strand"]
            })

    out=[]
    for locus in sorted(byloc):
        items=byloc[locus]
        mapped=[x for x in items if x["mapped"]]
        high=[x for x in mapped if x["mapq"]>=20 and x["cov"]>=0.70]
        refs=Counter(x["t"] for x in high)
        consensus=refs.most_common(1)[0][0] if refs else ""
        same=bool(high) and len(refs)==1
        starts=[x["ts"] for x in high if x["t"]==consensus]
        ends=[x["te"] for x in high if x["t"]==consensus]
        out.append({
            "reference_id":args.reference_id,
            "locus":locus,
            "source_records":len(items),
            "mapped_sources":len(mapped),
            "high_conf_sources":len(high),
            "high_conf_fraction":len(high)/len(items) if items else 0,
            "distinct_high_conf_refseqs":len(refs),
            "consensus_refseq":consensus,
            "all_high_conf_same_refseq":same,
            "consensus_start_min":min(starts) if starts else "",
            "consensus_end_max":max(ends) if ends else "",
            "median_mapq":med([x["mapq"] for x in mapped]),
            "median_query_coverage":med([x["cov"] for x in mapped]),
            "median_identity":med([x["ident"] for x in mapped]),
        })

    Path(args.out_csv).parent.mkdir(parents=True,exist_ok=True)
    fields=list(out[0].keys())
    with open(args.out_csv,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
    summary={
        "status":"ok","reference_id":args.reference_id,"expected_loci":241,"loci_reported":len(out),
        "loci_any_mapped":sum(r["mapped_sources"]>0 for r in out),
        "loci_any_high_conf":sum(r["high_conf_sources"]>0 for r in out),
        "loci_all_source_records_high_conf":sum(r["high_conf_sources"]==r["source_records"] for r in out),
        "loci_all_high_conf_same_refseq":sum(bool(r["all_high_conf_same_refseq"]) for r in out),
        "median_locus_high_conf_fraction":med([r["high_conf_fraction"] for r in out]),
        "median_of_locus_median_identity":med([r["median_identity"] for r in out if r["median_identity"] is not None]),
        "median_of_locus_median_query_coverage":med([r["median_query_coverage"] for r in out if r["median_query_coverage"] is not None]),
    }
    Path(args.out_json).write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
if __name__=="__main__":
    main()
