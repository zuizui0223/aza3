#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path

POS={"TRANSCRIPT_POSITIVE","TRANSCRIPT_STRONG"}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    files=sorted(args.input_dir.rglob("*_result.json"))
    rows=[json.loads(p.read_text(encoding="utf-8")) for p in files]
    expected={"FC","WY","FB","TJ","NH","LT"}
    got={r["sample_code"] for r in rows}
    if got!=expected: raise SystemExit(f"expected six samples {sorted(expected)}, got {sorted(got)}")
    if {r["morph"] for r in rows}!={"W","BP"}: raise SystemExit("morph panel mismatch")

    counts={}
    for target in ("DFR","ANS"):
        counts[target]={}
        for morph in ("W","BP"):
            subset=[r for r in rows if r["morph"]==morph]
            pos=sum(r["targets"][target]["tier"] in POS for r in subset)
            strong=sum(r["targets"][target]["tier"]=="TRANSCRIPT_STRONG" for r in subset)
            counts[target][morph]={"positive_n":pos,"strong_n":strong,"total_n":len(subset)}

    both=all(counts[t][m]["positive_n"]>=2 for t in ("DFR","ANS") for m in ("W","BP"))
    one_across=any(all(counts[t][m]["positive_n"]>=2 for m in ("W","BP")) for t in ("DFR","ANS"))
    any_pos=sum(counts[t][m]["positive_n"] for t in ("DFR","ANS") for m in ("W","BP"))
    if both: decision="RETENTION_BOTH_TARGETS_ACROSS_BOTH_MORPHS"
    elif one_across: decision="RETENTION_PARTIAL"
    elif any_pos>0: decision="ANCESTRY_OR_TISSUE_LIMITED"
    else: decision="ANCESTRY_OR_TISSUE_LIMITED"

    out={
        "result_version":"takaoense_public_pathway_retention_aggregate_v1",
        "samples":rows,"counts":counts,"decision":decision,
        "reactivation_ladder_effect":{
            "machinery_persistence":"PUBLIC_TRANSCRIPT_SUPPORT" if decision in {"RETENTION_BOTH_TARGETS_ACROSS_BOTH_MORPHS","RETENTION_PARTIAL"} else "NOT_PASSED",
            "suppression":"NOT_TESTED",
            "floral_reactivation":"NOT_TESTED",
            "historical_regain":"SEPARATE_HISTORY_GATE"
        },
        "boundary":"This six-run young-leaf screen can support coding/transcript retention only. It cannot establish floral suppression, floral reactivation, colour causation, adaptation, or historical regain."
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__": main()
