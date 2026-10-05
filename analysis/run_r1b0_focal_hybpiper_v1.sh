#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-SRR30887308}"
OUT="${2:-r1b0_focal_hybpiper}"
TARGET="${3:-comp1061_hybpiper_reference.fasta}"
PREFIX="sieboldii_original_comp1061"

mkdir -p "$OUT"

echo "[1/6] fetch"
prefetch "$RUN" --max-size 20G
SRA_PATH="$(find "$PWD" \( -path "*/$RUN/$RUN.sra" -o -path "*/$RUN/$RUN" \) | head -n 1 || true)"
if [ -z "$SRA_PATH" ]; then SRA_PATH="$RUN"; fi
fasterq-dump --split-files --threads 4 --outdir "$OUT" "$SRA_PATH"
R1="$OUT/${RUN}_1.fastq"
R2="$OUT/${RUN}_2.fastq"
test -s "$R1"; test -s "$R2"

echo "[2/6] frozen trimming"
fastp \
  -i "$R1" -I "$R2" \
  -o "$OUT/${RUN}_R1.trim.fastq.gz" \
  -O "$OUT/${RUN}_R2.trim.fastq.gz" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 4 \
  --json "$OUT/fastp.json" \
  --html "$OUT/fastp.html"

echo "[3/6] HybPiper assemble"
cd "$OUT"
hybpiper assemble \
  -t_dna "../$TARGET" \
  -r "${RUN}_R1.trim.fastq.gz" "${RUN}_R2.trim.fastq.gz" \
  --bwa \
  --prefix "$PREFIX" \
  --cpu 4 \
  2>&1 | tee hybpiper_assemble.log

echo "[4/6] summarize recovered coding loci"
python - "$PREFIX" "../$TARGET" "summary.json" "summary.md" <<'PY'
import json, sys, statistics
from pathlib import Path

root=Path(sys.argv[1])
target=Path(sys.argv[2])
outj=Path(sys.argv[3]); outm=Path(sys.argv[4])

# target locus universe
headers=[]
for line in target.open():
    if line.startswith(">"):
        headers.append(line[1:].split()[0])
target_loci=sorted({h.split("-",1)[1] for h in headers if "-" in h})

# HybPiper FNA output; filename stem is gene/locus
fna=list(root.rglob("*.FNA"))
by={}
for p in fna:
    gene=p.stem
    seqs=[]
    name=None; seq=[]
    for raw in p.open(errors="ignore"):
        s=raw.strip()
        if not s: continue
        if s.startswith(">"):
            if name is not None: seqs.append("".join(seq))
            name=s[1:].split()[0]; seq=[]
        else: seq.append(s)
    if name is not None: seqs.append("".join(seq))
    if seqs:
        by[gene]=max(len(x) for x in seqs)

# Detect genes_with_seqs file if present, but don't depend on exact HybPiper path.
gws=[]
for p in root.rglob("*genes_with_seqs*.txt"):
    for x in p.read_text(errors="ignore").split():
        if x in target_loci:
            gws.append(x)
gws=sorted(set(gws))

# Paralog warning candidates.
paralog_lines=[]
paralog_files=[]
for p in root.rglob("*paralog*"):
    if p.is_file() and p.stat().st_size < 10_000_000:
        txt=p.read_text(errors="ignore")
        if txt.strip():
            paralog_files.append(str(p))
            paralog_lines += [x for x in txt.splitlines() if x.strip()]
warn_loci=sorted({x.split()[0] for x in paralog_lines if x.split() and x.split()[0] in target_loci})

recovered=sorted(set(by) & set(target_loci))
if gws:
    recovered=sorted(set(recovered) | set(gws))
lens=sorted(by[g] for g in by if g in target_loci)

def q(p):
    if not lens:return None
    pos=(len(lens)-1)*p
    lo=int(pos); hi=min(lo+1,len(lens)-1); f=pos-lo
    return lens[lo]*(1-f)+lens[hi]*f

d={
  "status":"complete",
  "scope":"focal C. sieboldii HybPiper recovery against public 1061-locus Compositae1061 target; no reference-transfer decision",
  "target_sequence_records":len(headers),
  "target_distinct_loci":len(target_loci),
  "fna_files_found":len(fna),
  "recovered_target_loci":len(recovered),
  "recovered_target_fraction":len(recovered)/len(target_loci),
  "zero_recovered_loci":len(set(target_loci)-set(recovered)),
  "genes_with_seqs_loci_detected":len(gws),
  "paralog_warning_loci_detected":len(warn_loci),
  "paralog_files_detected":paralog_files,
  "recovered_sequence_length_q25":q(.25),
  "recovered_sequence_length_median":q(.5),
  "recovered_sequence_length_q75":q(.75),
  "recovered_sequence_length_min":min(lens) if lens else None,
  "recovered_sequence_length_max":max(lens) if lens else None,
  "recovered_loci":recovered,
  "unrecovered_loci":sorted(set(target_loci)-set(recovered)),
  "paralog_warning_loci":warn_loci
}
outj.write_text(json.dumps(d,indent=2)+"\n")
lines=[
 "# R1B-0 focal HybPiper recovery",
 "",
 f"- **Target loci:** {d['target_distinct_loci']}",
 f"- **Recovered loci:** {d['recovered_target_loci']} ({100*d['recovered_target_fraction']:.2f}%)",
 f"- **Unrecovered loci:** {d['zero_recovered_loci']}",
 f"- **Detected paralog-warning loci:** {d['paralog_warning_loci_detected']}",
 f"- **Recovered sequence median length:** {d['recovered_sequence_length_median']}",
 "",
 "> Focal target-recovery result only. It does not authorize PUBLIC_TRANSFER_GREEN or local-genomic inference."
]
outm.write_text("\n".join(lines)+"\n")
print("\n".join(lines))
PY

echo "[5/6] inventory outputs"
find "$PREFIX" -maxdepth 4 -type f | sort > output_files.txt
hybpiper --version > hybpiper_version.txt 2>&1 || true
bwa 2>&1 | head -n 3 > bwa_version.txt || true
spades.py --version > spades_version.txt 2>&1 || true

echo "[6/6] cleanup large reads/intermediates"
rm -f "${RUN}_1.fastq" "${RUN}_2.fastq" "${RUN}_R1.trim.fastq.gz" "${RUN}_R2.trim.fastq.gz"
find "$PREFIX" -type f \( -name "*.bam" -o -name "*.sam" -o -name "*.fastq" -o -name "*.fastq.gz" \) -delete || true
echo done
