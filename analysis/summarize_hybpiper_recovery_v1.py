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
                name=line[1:].split()[0]
                seq=[]
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--expected-loci", type=int, default=1061)
    ap.add_argument("--output", required=True, type=Path)
    args=ap.parse_args()

    fna=list(args.root.rglob("*.FNA"))
    genes={}
    sources={}
    for p in fna:
        records=list(read_fasta(p))
        if not records:
            continue
        gene=p.stem
        # HybPiper normally yields one coding sequence per gene/sample.
        best=max((s for _,s in records), key=len)
        if gene not in genes or len(best)>len(genes[gene]):
            genes[gene]=best
            sources[gene]=str(p)

    paralog_files=[str(p) for p in args.root.rglob("*") if p.is_file() and "paralog" in p.name.casefold()]
    warning_lines=0
    warning_files=[]
    for s in paralog_files:
        p=Path(s)
        try:
            lines=[x.strip() for x in p.read_text(encoding="utf-8", errors="replace").splitlines() if x.strip()]
        except Exception:
            continue
        if lines:
            warning_files.append({"path":s,"nonempty_lines":len(lines),"preview":lines[:10]})
            warning_lines += len(lines)

    lengths=[len(s) for s in genes.values()]
    out={
        "status":"ok",
        "root":str(args.root),
        "expected_target_loci":args.expected_loci,
        "recovered_distinct_fna_loci":len(genes),
        "recovered_fraction":len(genes)/args.expected_loci if args.expected_loci else None,
        "total_recovered_bp":sum(lengths),
        "median_recovered_bp":None if not lengths else sorted(lengths)[len(lengths)//2],
        "min_recovered_bp":None if not lengths else min(lengths),
        "max_recovered_bp":None if not lengths else max(lengths),
        "paralog_named_files":len(paralog_files),
        "paralog_nonempty_line_total":warning_lines,
        "paralog_files_with_content":warning_files,
        "gene_source_files":sources,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in out.items() if k not in {"gene_source_files","paralog_files_with_content"}},indent=2))

if __name__=="__main__":
    main()
