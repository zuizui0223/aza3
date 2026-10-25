#!/usr/bin/env python3
"""Summarize one HybPiper sample directory into a compact JSON receipt."""
from __future__ import annotations
import argparse, json, csv
from pathlib import Path


def count_nonempty_lines(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for x in path.read_text(encoding="utf-8", errors="ignore").splitlines() if x.strip())


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sample-dir", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--target-id", required=True)
    ap.add_argument("--output", required=True, type=Path)
    args=ap.parse_args()

    d=args.sample_dir
    genes=d/"genes_with_seqs.txt"
    recovered=[x.strip() for x in genes.read_text(encoding="utf-8").splitlines() if x.strip()] if genes.exists() else []

    fna=list(d.rglob("*.FNA"))
    paralog_candidates=[]
    for p in d.rglob("*paralog*"):
        if p.is_file():
            paralog_candidates.append(str(p.relative_to(d)))

    read_counts={}
    grc=d/"gene_read_counts.tsv"
    if grc.exists():
        with grc.open(encoding="utf-8", errors="ignore") as fh:
            reader=csv.reader(fh, delimiter="\t")
            rows=list(reader)
            if rows:
                header=rows[0]
                for row in rows[1:]:
                    if len(row)>=2:
                        try:
                            read_counts[row[0]]=int(float(row[1]))
                        except Exception:
                            pass

    receipt={
        "status":"ok" if recovered else "no_recovered_genes",
        "sample":args.sample,
        "target_id":args.target_id,
        "sample_dir":str(d),
        "genes_with_seqs_count":len(recovered),
        "genes_with_seqs":recovered,
        "fna_file_count":len(fna),
        "gene_read_counts_file_present":grc.exists(),
        "gene_read_counts_nonzero_genes":sum(v>0 for v in read_counts.values()),
        "gene_read_counts_total":sum(read_counts.values()) if read_counts else None,
        "paralog_related_files":paralog_candidates[:100],
        "paralog_related_file_count":len(paralog_candidates),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in receipt.items() if k!="genes_with_seqs"},indent=2))


if __name__=="__main__":
    main()
