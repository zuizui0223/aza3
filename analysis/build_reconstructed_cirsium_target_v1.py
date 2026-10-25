#!/usr/bin/env python3
"""Reconstruct a Cirsium-augmented HybPiper nucleotide target file.

This does not claim byte identity with the intermediate target file used by
Moreyra et al. It reproduces the published logic: append homologous Cirsium
tioganum coding sequences, recovered with the original Compositae1061 target,
as additional source sequences for the same loci.

Input:
  --original-target  public comp1061_hybpiper_reference.fasta
  --tioganum-root    HybPiper output directory for SRR25265669
Output:
  --output           reconstructed target FASTA
  --report           JSON reconstruction receipt
"""
from __future__ import annotations
import argparse, json
from pathlib import Path


def read_fasta(path: Path):
    name = None
    seq = []
    with path.open(encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def write_record(fh, name: str, seq: str, width: int = 80):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i+width] + "\n")


def locus_from_target_header(header: str) -> str:
    if "-" not in header:
        raise ValueError(f"target header lacks source-locus separator: {header}")
    return header.split("-", 1)[1]


def find_tioganum_gene_sequences(root: Path):
    """Return gene->sequence from HybPiper per-gene FNA output.

    HybPiper stores extracted coding sequences under per-gene directories.
    We scan *.FNA recursively, use the filename stem as the gene identifier,
    and retain the longest non-empty sequence if duplicates are encountered.
    """
    out = {}
    source = {}
    for p in sorted(root.rglob("*.FNA")):
        records = [(h, s.upper()) for h, s in read_fasta(p) if s]
        if not records:
            continue
        _, seq = max(records, key=lambda x: len(x[1]))
        gene = p.stem
        if gene not in out or len(seq) > len(out[gene]):
            out[gene] = seq
            source[gene] = str(p)
    return out, source


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--original-target", required=True, type=Path)
    ap.add_argument("--tioganum-root", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--report", required=True, type=Path)
    args = ap.parse_args()

    original = list(read_fasta(args.original_target))
    if not original:
        raise SystemExit("original target FASTA is empty")

    original_loci = {}
    for header, seq in original:
        locus = locus_from_target_header(header)
        original_loci.setdefault(locus, 0)
        original_loci[locus] += 1

    tiog, tiog_sources = find_tioganum_gene_sequences(args.tioganum_root)
    appendable = {g: s for g, s in tiog.items() if g in original_loci and s}
    off_target = sorted(set(tiog) - set(original_loci))

    existing_headers = {h for h, _ in original}
    appended = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as out:
        for h, s in original:
            write_record(out, h, s.upper())
        for gene in sorted(appendable):
            header = f"tiog-{gene}"
            if header in existing_headers:
                continue
            write_record(out, header, appendable[gene])
            appended.append(gene)

    receipt = {
        "status": "ok",
        "method": "published-logic reconstruction; not byte-identical-author-file claim",
        "original_target": str(args.original_target),
        "tioganum_root": str(args.tioganum_root),
        "output": str(args.output),
        "original_sequence_records": len(original),
        "original_distinct_loci": len(original_loci),
        "tioganum_fna_genes_found": len(tiog),
        "tioganum_genes_matching_target_loci": len(appendable),
        "tioganum_genes_appended": len(appended),
        "target_loci_without_tioganum_sequence": len(set(original_loci) - set(appended)),
        "off_target_fna_gene_count": len(off_target),
        "off_target_fna_genes": off_target[:100],
        "source_files_for_appended_genes": {g: tiog_sources[g] for g in appended},
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in receipt.items() if k != "source_files_for_appended_genes"}, indent=2))


if __name__ == "__main__":
    main()
