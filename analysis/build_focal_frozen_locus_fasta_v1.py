#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path

def read_fasta(path: Path):
    name=None; seq=[]
    with path.open(encoding="utf-8", errors="replace") as fh:
        for raw in fh:
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
    ap.add_argument("--hybpiper-root",required=True,type=Path)
    ap.add_argument("--loci",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--report",required=True,type=Path)
    ap.add_argument("--sample-id",default="SRR30887308")
    args=ap.parse_args()

    loci=[x.strip() for x in args.loci.read_text(encoding="utf-8").splitlines() if x.strip()]
    if len(loci)!=len(set(loci)):
        raise SystemExit("duplicate locus IDs in frozen list")
    wanted=set(loci)

    found={}
    sources={}
    for p in args.hybpiper_root.rglob("*.FNA"):
        gene=p.stem
        if gene not in wanted:
            continue
        recs=[(h,s.upper()) for h,s in read_fasta(p) if s]
        if not recs:
            continue
        h,seq=max(recs,key=lambda x:len(x[1]))
        if gene not in found or len(seq)>len(found[gene]):
            found[gene]=seq
            sources[gene]=str(p)

    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("w",encoding="utf-8") as out:
        for gene in loci:
            if gene in found:
                write_record(out,f"{gene}|{args.sample_id}",found[gene])

    missing=[g for g in loci if g not in found]
    report={
        "status":"ok",
        "sample_id":args.sample_id,
        "requested_loci":len(loci),
        "recovered_loci":len(found),
        "recovery_fraction":len(found)/len(loci) if loci else None,
        "total_bp":sum(map(len,found.values())),
        "missing_loci":missing,
        "source_files":sources,
        "claim_boundary":"FASTA construction only; no reference-transferability conclusion."
    }
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k not in {"missing_loci","source_files"}},indent=2))

if __name__=="__main__":
    main()
