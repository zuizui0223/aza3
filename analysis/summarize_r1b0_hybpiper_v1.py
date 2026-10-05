#!/usr/bin/env python3
"""Summarize a single HybPiper R1B-0 run without interpreting biology."""
from __future__ import annotations
import argparse, json
from pathlib import Path


def fasta_headers(path: Path):
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip().split()[0]


def target_locus(header: str) -> str:
    return header.split("-", 1)[1] if "-" in header else header


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sample-dir", type=Path, required=True)
    ap.add_argument("--target", type=Path, required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--taxon", required=True)
    ap.add_argument("--target-version", required=True)
    ap.add_argument("--output", type=Path, required=True)
    args=ap.parse_args()

    target_headers=list(fasta_headers(args.target))
    target_loci=sorted({target_locus(h) for h in target_headers})

    fna_files=sorted(args.sample_dir.rglob("*.FNA"))
    recovered={}
    for p in fna_files:
        records=list(fasta_headers(p))
        if not records:
            continue
        locus=p.stem
        recovered[locus]=str(p)

    # HybPiper output names should correspond to target loci. Report, do not silently coerce.
    matched=sorted(set(recovered) & set(target_loci))
    unexpected=sorted(set(recovered) - set(target_loci))

    # Search standard text outputs for paralog indicators without assuming one exact filename.
    paralog_files=sorted(
        p for p in args.sample_dir.rglob("*")
        if p.is_file() and "paralog" in p.name.casefold()
    )
    nonempty_paralog_files=[str(p) for p in paralog_files if p.stat().st_size > 0]

    receipt={
        "status":"ok",
        "run":args.run,
        "taxon":args.taxon,
        "target_version":args.target_version,
        "sample_dir":str(args.sample_dir),
        "target_sequence_records":len(target_headers),
        "target_distinct_loci":len(target_loci),
        "fna_files_found":len(fna_files),
        "recovered_target_loci":len(matched),
        "recovered_target_locus_fraction":(len(matched)/len(target_loci) if target_loci else None),
        "unexpected_fna_loci":unexpected,
        "nonempty_paralog_files":nonempty_paralog_files,
        "claim_boundary":"Captured-locus recovery only; no local ancestry, haplotype-age, recombination, structural-variant, genomic-modularity, or genotype-phenotype inference."
    }
    args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
