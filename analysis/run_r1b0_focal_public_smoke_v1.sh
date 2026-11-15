#!/usr/bin/env bash
set -euo pipefail

RUN="${1:-SRR30887308}"
OUT="${2:-r1b0_focal_smoke}"
TARGET="${3:-comp1061_hybpiper_reference.fasta}"

mkdir -p "$OUT"

echo "[1/7] versions"
{
  echo "run=$RUN"
  echo "date_utc=$(date -u +%FT%TZ)"
  command -v prefetch >/dev/null && prefetch --version | head -n 1 || true
  command -v fasterq-dump >/dev/null && fasterq-dump --version | head -n 1 || true
  bwa 2>&1 | head -n 3 || true
  samtools --version | head -n 1
} > "$OUT/tool_versions.txt"

echo "[2/7] target audit"
python - "$TARGET" "$OUT/target_audit.json" <<'PY'
import json, sys
from pathlib import Path
p=Path(sys.argv[1]); out=Path(sys.argv[2])
headers=[]
seq_bp=0
name=None
with p.open() as fh:
    for line in fh:
        line=line.strip()
        if not line: continue
        if line.startswith(">"):
            headers.append(line[1:].split()[0])
        else:
            seq_bp += len(line)
loci=sorted({h.split("-",1)[1] for h in headers if "-" in h})
prefixes=sorted({h.split("-",1)[0] for h in headers if "-" in h})
d={"target":str(p),"sequence_records":len(headers),"distinct_loci":len(loci),
   "reference_prefixes":prefixes,"total_target_bp":seq_bp}
out.write_text(json.dumps(d,indent=2)+"\n")
print(json.dumps(d,indent=2))
PY

echo "[3/7] fetch SRA"
prefetch "$RUN" --max-size 20G
SRA_PATH="$(find "$PWD" -path "*/$RUN/$RUN.sra" -o -path "*/$RUN/$RUN" | head -n 1 || true)"
if [ -z "$SRA_PATH" ]; then
  SRA_PATH="$RUN"
fi

echo "[4/7] FASTQ"
fasterq-dump --split-files --threads 4 --outdir "$OUT" "$SRA_PATH"
R1="$OUT/${RUN}_1.fastq"
R2="$OUT/${RUN}_2.fastq"
test -s "$R1"; test -s "$R2"

python - "$R1" "$R2" "$OUT/fastq_receipt.json" <<'PY'
import json, sys
from pathlib import Path
def stats(path):
    lines=0; bp=0; reads=0
    with open(path) as fh:
        while True:
            h=fh.readline()
            if not h: break
            seq=fh.readline().strip(); plus=fh.readline(); qual=fh.readline()
            reads += 1; bp += len(seq); lines += 4
    return {"path":path,"reads":reads,"bases":bp}
d={"R1":stats(sys.argv[1]),"R2":stats(sys.argv[2])}
Path(sys.argv[3]).write_text(json.dumps(d,indent=2)+"\n")
print(json.dumps(d,indent=2))
PY

echo "[5/7] direct BWA mapping to public Compositae1061 reference"
cp "$TARGET" "$OUT/target.fasta"
bwa index "$OUT/target.fasta"
bwa mem -t 4 "$OUT/target.fasta" "$R1" "$R2" \
  | samtools sort -@ 2 -o "$OUT/${RUN}.target.sorted.bam"
samtools index "$OUT/${RUN}.target.sorted.bam"
samtools flagstat "$OUT/${RUN}.target.sorted.bam" > "$OUT/flagstat.txt"
samtools idxstats "$OUT/${RUN}.target.sorted.bam" > "$OUT/idxstats.tsv"

echo "[6/7] locus-level smoke summary"
python - "$TARGET" "$OUT/idxstats.tsv" "$OUT/flagstat.txt" "$OUT/summary.json" <<'PY'
import json,re,sys
from pathlib import Path
target, idxp, flagp, outp = map(Path, sys.argv[1:])
record_to_locus={}
headers=[]
with target.open() as fh:
    for line in fh:
        if line.startswith(">"):
            h=line[1:].split()[0]
            headers.append(h)
            record_to_locus[h]=h.split("-",1)[1] if "-" in h else h
mapped_by_record={}
with idxp.open() as fh:
    for line in fh:
        ref,length,mapped,unmapped=line.rstrip("\n").split("\t")
        if ref=="*": continue
        mapped_by_record[ref]=int(mapped)
hit_records={r for r,n in mapped_by_record.items() if n>0}
hit_loci={record_to_locus[r] for r in hit_records if r in record_to_locus}
all_loci=set(record_to_locus.values())
flag=flagp.read_text()
m_total=re.search(r"^(\d+) \+ \d+ in total",flag,re.M)
m_mapped=re.search(r"^(\d+) \+ \d+ mapped \(([^%]+)%",flag,re.M)
m_proper=re.search(r"^(\d+) \+ \d+ properly paired \(([^%]+)%",flag,re.M)
d={
 "target_sequence_records":len(headers),
 "target_distinct_loci":len(all_loci),
 "target_records_with_any_mapped_read":len(hit_records),
 "target_loci_with_any_mapped_read":len(hit_loci),
 "target_locus_hit_fraction":len(hit_loci)/len(all_loci) if all_loci else None,
 "flagstat_total_reads":int(m_total.group(1)) if m_total else None,
 "flagstat_mapped_reads":int(m_mapped.group(1)) if m_mapped else None,
 "flagstat_mapped_percent":float(m_mapped.group(2)) if m_mapped else None,
 "flagstat_properly_paired_reads":int(m_proper.group(1)) if m_proper else None,
 "flagstat_properly_paired_percent":float(m_proper.group(2)) if m_proper else None,
 "scope":"technical smoke test only; direct mapping to original Compositae1061 target; not HybPiper recovery and not reference-transfer decision"
}
outp.write_text(json.dumps(d,indent=2)+"\n")
print(json.dumps(d,indent=2))
PY

echo "[7/7] compact markdown"
python - "$OUT/summary.json" "$OUT/summary.md" <<'PY'
import json,sys
from pathlib import Path
d=json.load(open(sys.argv[1]))
rows=[
("# R1B-0 focal technical smoke test",""),
("Target distinct loci",d["target_distinct_loci"]),
("Loci with >=1 mapped read",d["target_loci_with_any_mapped_read"]),
("Locus hit fraction",d["target_locus_hit_fraction"]),
("Mapped reads",d["flagstat_mapped_reads"]),
("Mapped percent",d["flagstat_mapped_percent"]),
("Properly paired percent",d["flagstat_properly_paired_percent"]),
]
lines=[]
for k,v in rows:
    if k.startswith("#"): lines.append(k)
    else: lines.append(f"- **{k}:** {v}")
lines += ["","> Technical smoke only. This does not replace HybPiper locus recovery or authorize any biological inference."]
Path(sys.argv[2]).write_text("\n".join(lines)+"\n")
PY

rm -f "$R1" "$R2"
echo "done"
