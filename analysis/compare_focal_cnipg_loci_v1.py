#!/usr/bin/env python3
"""Compare focal C. sieboldii HybPiper loci with frozen C. nipponicum genome-derived loci.

This is a sequence-compatibility layer for aza3 R1B-0. It does not test chromosome
coordinates or local ancestry. The C. nipponicum pack is materialized from the frozen
EAzami durable payload and the focal sequences come from the independent aza3
HybPiper rerun.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import median

from Bio.Align import PairwiseAligner


def read_fasta(path: Path):
    name = None
    seq = []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq).upper()


def best_fna_by_locus(root: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for p in root.rglob("*.FNA"):
        records = [(h, s) for h, s in read_fasta(p) if s]
        if not records:
            continue
        seq = max((s for _, s in records), key=len)
        locus = p.stem
        if locus not in out or len(seq) > len(out[locus]):
            out[locus] = seq
    return out


def read_single_fasta(path: Path) -> str:
    recs = [(h, s) for h, s in read_fasta(path) if s]
    if len(recs) != 1:
        raise ValueError(f"{path}: expected exactly one non-empty FASTA record, found {len(recs)}")
    return recs[0][1]


def quantile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    x = sorted(values)
    pos = (len(x) - 1) * p
    lo = int(pos)
    hi = min(lo + 1, len(x) - 1)
    frac = pos - lo
    return x[lo] + (x[hi] - x[lo]) * frac


def local_pair_metrics(a: str, b: str, aligner: PairwiseAligner) -> dict[str, float | int]:
    aln = aligner.align(a, b)[0]
    # Alignment.aligned contains ungapped aligned blocks for target/query.
    ta, qa = aln.aligned
    aligned_pairs = 0
    matches = 0
    for (ts, te), (qs, qe) in zip(ta, qa):
        n = min(te - ts, qe - qs)
        if n <= 0:
            continue
        aligned_pairs += n
        matches += sum(x == y for x, y in zip(a[ts:ts+n], b[qs:qs+n]))
    mismatches = aligned_pairs - matches
    identity = matches / aligned_pairs if aligned_pairs else 0.0
    shorter = min(len(a), len(b))
    coverage_shorter = aligned_pairs / shorter if shorter else 0.0
    return {
        "focal_bp": len(a),
        "cnipg_bp": len(b),
        "aligned_ungapped_pairs": aligned_pairs,
        "matches": matches,
        "mismatches": mismatches,
        "ungapped_identity": identity,
        "coverage_of_shorter": coverage_shorter,
        "alignment_score": float(aln.score),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--focal-root", type=Path, required=True)
    ap.add_argument("--cnipg-root", type=Path, required=True)
    ap.add_argument("--frozen-loci", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    args = ap.parse_args()

    frozen = [x.strip() for x in args.frozen_loci.read_text(encoding="utf-8").splitlines() if x.strip()]
    if len(frozen) != 241 or len(set(frozen)) != 241:
        raise ValueError(f"expected 241 unique frozen loci, found {len(frozen)} / {len(set(frozen))}")

    focal = best_fna_by_locus(args.focal_root)
    cnipg_loci_dir = args.cnipg_root / "loci"
    cnipg = {}
    for locus in frozen:
        p = cnipg_loci_dir / f"{locus}.fasta"
        if p.exists():
            cnipg[locus] = read_single_fasta(p)

    common = sorted(set(frozen) & set(focal) & set(cnipg))

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5

    rows = []
    for locus in common:
        m = local_pair_metrics(focal[locus], cnipg[locus], aligner)
        rows.append({"locus": locus, **m})

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "locus", "focal_bp", "cnipg_bp", "aligned_ungapped_pairs", "matches",
        "mismatches", "ungapped_identity", "coverage_of_shorter", "alignment_score",
    ]
    with args.output_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    identities = [float(r["ungapped_identity"]) for r in rows]
    coverage = [float(r["coverage_of_shorter"]) for r in rows]
    focal_lengths = [int(r["focal_bp"]) for r in rows]
    cn_lengths = [int(r["cnipg_bp"]) for r in rows]

    summary = {
        "result_version": "r1b0_focal_cnipg_sequence_compatibility_v1",
        "status": "descriptive_sequence_compatibility_only",
        "frozen_locus_universe": 241,
        "focal_hybpiper_loci_found": len(set(focal) & set(frozen)),
        "cnipg_strict_loci_found": len(set(cnipg) & set(frozen)),
        "common_loci_compared": len(common),
        "common_fraction_of_cnipg_pack": len(common) / len(cnipg) if cnipg else None,
        "identity": {
            "q05": quantile(identities, 0.05),
            "q25": quantile(identities, 0.25),
            "median": quantile(identities, 0.50),
            "q75": quantile(identities, 0.75),
            "q95": quantile(identities, 0.95),
            "min": min(identities) if identities else None,
            "max": max(identities) if identities else None,
        },
        "coverage_of_shorter": {
            "q05": quantile(coverage, 0.05),
            "median": quantile(coverage, 0.50),
            "q95": quantile(coverage, 0.95),
        },
        "median_focal_bp": median(focal_lengths) if focal_lengths else None,
        "median_cnipg_bp": median(cn_lengths) if cn_lengths else None,
        "loci_identity_ge_0_90": sum(x >= 0.90 for x in identities),
        "loci_identity_ge_0_95": sum(x >= 0.95 for x in identities),
        "loci_identity_ge_0_98": sum(x >= 0.98 for x in identities),
        "claim_boundary": (
            "Sequence-level compatibility between recovered coding loci only. "
            "This does not establish chromosome-coordinate concordance, local ancestry, "
            "haplotype age, structural-variant reuse or genomic modularity."
        ),
    }
    args.output_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))

    if not rows:
        raise ValueError("no common focal/CNIPG loci were alignable")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
