#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path

def read_fasta(path: Path):
    name=None; seq=[]
    for raw in path.read_text(encoding="utf-8").splitlines():
        line=raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name=line[1:].split()[0]; seq=[]
        else:
            seq.append(line)
    if name is not None:
        yield name, "".join(seq)

def write_record(fh, name, seq, width=80):
    fh.write(f">{name}\n")
    for i in range(0,len(seq),width):
        fh.write(seq[i:i+width]+"\n")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--target",required=True,type=Path)
    ap.add_argument("--loci",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--receipt",required=True,type=Path)
    args=ap.parse_args()
    loci=[x.strip() for x in args.loci.read_text().splitlines() if x.strip()]
    want=set(loci)
    if len(want)!=241:
        raise SystemExit(f"expected 241 distinct loci, found {len(want)}")
    rows=[]
    seen=set()
    for h,s in read_fasta(args.target):
        if "-" not in h:
            continue
        src,locus=h.split("-",1)
        if locus in want:
            q=f"{src}__{locus}"
            rows.append((q,src,locus,s.upper()))
            seen.add(locus)
    missing=sorted(want-seen)
    if missing:
        raise SystemExit(f"missing target loci: {missing[:20]}")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("w") as fh:
        for q,src,locus,s in rows:
            write_record(fh,q,s)
    counts={}
    for _,src,locus,_ in rows:
        counts[locus]=counts.get(locus,0)+1
    receipt={
        "status":"ok",
        "target":str(args.target),
        "locus_manifest":str(args.loci),
        "distinct_loci":len(seen),
        "query_records":len(rows),
        "records_per_locus_min":min(counts.values()),
        "records_per_locus_max":max(counts.values()),
        "records_per_locus_hist":{str(k):list(counts.values()).count(k) for k in sorted(set(counts.values()))},
        "sources":sorted({src for _,src,_,_ in rows})
    }
    args.receipt.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))
if __name__=="__main__":
    main()
