#!/usr/bin/env bash
set -euo pipefail

OUT="${1:-r1b0_panel_direct_smoke}"
TARGET="${2:-comp1061_hybpiper_reference.fasta}"
mkdir -p "$OUT/per_run" "$OUT/per_locus"

RUNS=(
"SRR30887308|Cirsium sieboldii"
"SRR30887271|Cirsium japonicum"
"SRR30887240|Cirsium lineare"
"SRR25265660|Cirsium nipponicum"
"SRR30887223|Cirsium nipponicum var. yoshinoi"
"SRR30887286|Cirsium nipponicum var. incomptum"
"SRR25265649|Cirsium pendulum"
"SRR30887259|Cirsium dipsacolepis"
"SRR30887291|Cirsium tanakae"
"SRR25265685|Cirsium japonicum var. maackii"
)

python - "$TARGET" "$OUT/target_audit.json" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]); out=Path(sys.argv[2])
headers=[]; bp=0
for line in p.open():
    s=line.strip()
    if not s: continue
    if s.startswith(">"): headers.append(s[1:].split()[0])
    else: bp+=len(s)
loci=sorted({h.split("-",1)[1] for h in headers if "-" in h})
prefixes=sorted({h.split("-",1)[0] for h in headers if "-" in h})
d={"sequence_records":len(headers),"distinct_loci":len(loci),"total_target_bp":bp,"reference_prefixes":prefixes}
assert d["sequence_records"]==2597, d
assert d["distinct_loci"]==1061, d
out.write_text(json.dumps(d,indent=2)+"\n")
PY

cp "$TARGET" "$OUT/target.fasta"
bwa index "$OUT/target.fasta"

echo "run,taxon,read_pairs,total_primary_reads,mapped_reads,mapped_percent,properly_paired_reads,properly_paired_percent,loci_ge_1,loci_ge_10,loci_ge_50,loci_ge_100,zero_hit_loci,median_mapped_reads_per_locus,q25,q75,max_mapped_reads_per_locus" > "$OUT/panel_summary.csv"
echo "run,taxon,locus,mapped_reads,target_records,hit_records" > "$OUT/per_locus_long.csv"

for item in "${RUNS[@]}"; do
  IFS='|' read -r RUN TAXON <<< "$item"
  echo "=== $RUN $TAXON ==="
  mkdir -p "$OUT/per_run/$RUN"

  prefetch "$RUN" --max-size 20G
  SRA_PATH="$(find "$PWD" \( -path "*/$RUN/$RUN.sra" -o -path "*/$RUN/$RUN" \) | head -n 1 || true)"
  if [ -z "$SRA_PATH" ]; then SRA_PATH="$RUN"; fi

  fasterq-dump --split-files --threads 4 --outdir "$OUT/per_run/$RUN" "$SRA_PATH"
  R1="$OUT/per_run/$RUN/${RUN}_1.fastq"
  R2="$OUT/per_run/$RUN/${RUN}_2.fastq"
  test -s "$R1"; test -s "$R2"

  bwa mem -t 4 "$OUT/target.fasta" "$R1" "$R2" \
    | samtools sort -@ 2 -o "$OUT/per_run/$RUN/${RUN}.bam"
  samtools flagstat "$OUT/per_run/$RUN/${RUN}.bam" > "$OUT/per_run/$RUN/flagstat.txt"
  samtools idxstats "$OUT/per_run/$RUN/${RUN}.bam" > "$OUT/per_run/$RUN/idxstats.tsv"

  python - "$RUN" "$TAXON" "$OUT/per_run/$RUN/idxstats.tsv" "$OUT/per_run/$RUN/flagstat.txt" "$OUT/panel_summary.csv" "$OUT/per_locus_long.csv" <<'PY'
import csv,re,sys,statistics
from collections import defaultdict
run,taxon,idxp,flagp,sumout,longout=sys.argv[1:]
records=defaultdict(list)
with open(idxp) as fh:
    for line in fh:
        ref,length,mapped,unmapped=line.rstrip("\n").split("\t")
        if ref=="*": continue
        locus=ref.split("-",1)[1] if "-" in ref else ref
        records[locus].append(int(mapped))
loci=sorted(records)
mapped_locus={l:sum(records[l]) for l in loci}
flag=open(flagp).read()
def grab(pattern, group=1, cast=int):
    m=re.search(pattern,flag,re.M)
    return cast(m.group(group).replace(",","")) if m else None
total_primary=grab(r"^(\d+) \+ \d+ primary$")
mapped=grab(r"^(\d+) \+ \d+ mapped \(([^%]+)%",1,int)
mapped_pct=grab(r"^(\d+) \+ \d+ mapped \(([^%]+)%",2,float)
proper=grab(r"^(\d+) \+ \d+ properly paired \(([^%]+)%",1,int)
proper_pct=grab(r"^(\d+) \+ \d+ properly paired \(([^%]+)%",2,float)
read_pairs=total_primary//2 if total_primary is not None else None
vals=sorted(mapped_locus.values())
def q(p):
    if not vals:return None
    pos=(len(vals)-1)*p
    lo=int(pos); hi=min(lo+1,len(vals)-1); frac=pos-lo
    return vals[lo]*(1-frac)+vals[hi]*frac
row=[
 run,taxon,read_pairs,total_primary,mapped,mapped_pct,proper,proper_pct,
 sum(v>=1 for v in vals),sum(v>=10 for v in vals),sum(v>=50 for v in vals),sum(v>=100 for v in vals),
 sum(v==0 for v in vals),q(.5),q(.25),q(.75),max(vals) if vals else None
]
with open(sumout,"a",newline="") as fh: csv.writer(fh).writerow(row)
with open(longout,"a",newline="") as fh:
    w=csv.writer(fh)
    for locus in loci:
        arr=records[locus]
        w.writerow([run,taxon,locus,sum(arr),len(arr),sum(x>0 for x in arr)])
PY

  rm -f "$R1" "$R2" "$OUT/per_run/$RUN/${RUN}.bam"
  rm -rf "$RUN"
done

python - "$OUT/panel_summary.csv" "$OUT/panel_comparison.json" "$OUT/panel_comparison.md" <<'PY'
import csv,json,statistics,sys
from pathlib import Path
rows=list(csv.DictReader(open(sys.argv[1])))
numcols=["mapped_percent","properly_paired_percent","loci_ge_1","loci_ge_10","loci_ge_50","loci_ge_100","zero_hit_loci","median_mapped_reads_per_locus"]
for r in rows:
    for c in numcols:r[c]=float(r[c])
focal=next(r for r in rows if r["run"]=="SRR30887308")
metrics={}
for c in ["mapped_percent","properly_paired_percent","loci_ge_1","loci_ge_10","loci_ge_50","loci_ge_100","median_mapped_reads_per_locus"]:
    vals=[r[c] for r in rows]
    asc=sorted(vals)
    metrics[c]={
      "focal":focal[c],
      "panel_min":min(vals),"panel_median":statistics.median(vals),"panel_max":max(vals),
      "focal_rank_ascending":1+sum(v<focal[c] for v in vals),
      "n":len(vals)
    }
d={
 "status":"complete",
 "scope":"direct mapping to original Compositae1061 target only; not HybPiper recovery and not R1B-0 final decision",
 "target_loci":1061,
 "focal_run":"SRR30887308",
 "metrics":metrics,
 "rows":rows
}
Path(sys.argv[2]).write_text(json.dumps(d,indent=2)+"\n")
lines=["# R1B-0 ten-run direct-mapping panel","",f"Focal: {focal['run']} {focal['taxon']}",""]
for c,m in metrics.items():
    lines.append(f"- **{c}:** focal {m['focal']}; panel median {m['panel_median']}; range {m['panel_min']}–{m['panel_max']}; ascending rank {m['focal_rank_ascending']}/{m['n']}")
lines += ["","> Direct-mapping technical comparison only. It does not replace HybPiper recovery, target-file sensitivity, or three-genome reference localization."]
Path(sys.argv[3]).write_text("\n".join(lines)+"\n")
PY

gzip -f "$OUT/per_locus_long.csv"
echo "done"
