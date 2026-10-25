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

def read_locus_list(path: Path) -> list[str]:
    vals=[x.strip() for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]
    if len(vals) != len(set(vals)):
        raise ValueError(f"duplicate loci in frozen list: {path}")
    return vals

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--expected-loci", type=int, default=1061,
                    help="Broad public named-locus universe; compatibility only, not historical author target count.")
    ap.add_argument("--frozen-loci", type=Path, default=None,
                    help="Optional frozen high-stringency locus list, e.g. Moreyra-compatible 241.")
    ap.add_argument("--frozen-loci-id", default="MOREYRA_COMPATIBILITY_241")
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
        "broad_locus_universe_id":"PUBLIC_1061_BROAD_RECOVERY",
        "expected_public_named_loci":args.expected_loci,
        "historical_author_target_count_reported":1064,
        "historical_target_identity_status":"EXACT_1064_TARGET_UNRECOVERED__THREE_LOCUS_DIFFERENCE_UNRESOLVED",
        "recovered_distinct_fna_loci":len(genes),
        "public_1061_recovery_fraction":len(genes)/args.expected_loci if args.expected_loci else None,
        "total_recovered_bp":sum(lengths),
        "median_recovered_bp":None if not lengths else sorted(lengths)[len(lengths)//2],
        "min_recovered_bp":None if not lengths else min(lengths),
        "max_recovered_bp":None if not lengths else max(lengths),
        "paralog_named_files":len(paralog_files),
        "paralog_nonempty_line_total":warning_lines,
        "paralog_files_with_content":warning_files,
        "gene_source_files":sources,
        "claim_boundary":"Compatibility recovery only. Public 1061 named loci are not the exact historical author 1064-target file or the published final 350-locus matrix."
    }

    if args.frozen_loci is not None:
        frozen=read_locus_list(args.frozen_loci)
        recovered=sorted(set(genes) & set(frozen))
        out["frozen_locus_layer_id"]=args.frozen_loci_id
        out["frozen_loci_path"]=str(args.frozen_loci)
        out["frozen_loci_count"]=len(frozen)
        out["frozen_loci_recovered"]=len(recovered)
        out["frozen_loci_recovery_fraction"]=len(recovered)/len(frozen) if frozen else None
        out["frozen_loci_missing"]=sorted(set(frozen)-set(genes))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in out.items()
                      if k not in {"gene_source_files","paralog_files_with_content","frozen_loci_missing"}},indent=2))

if __name__=="__main__":
    main()
