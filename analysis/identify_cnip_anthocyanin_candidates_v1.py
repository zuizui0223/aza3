#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

QUERY_CLASS={
    "DFR_GERBERA_P51105":"DFR",
    "DFR_ARABIDOPSIS_P51102":"DFR",
    "ANS_PETUNIA_P51092":"ANS",
    "ANS_ARABIDOPSIS_Q96323":"ANS",
    "DECOY_ANR_ARABIDOPSIS_Q9SEV0":"ANR_DECOY",
    "DECOY_FLS1_ARABIDOPSIS_Q96330":"FLS_DECOY",
}

def read_fasta(path):
    name=None; seq=[]
    with path.open(encoding="utf-8",errors="ignore") as fh:
        for raw in fh:
            line=raw.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None: yield name,"".join(seq)
                name=line[1:].split()[0]; seq=[]
            else: seq.append(line)
    if name is not None: yield name,"".join(seq)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--hits",required=True,type=Path)
    ap.add_argument("--proteins",required=True,type=Path)
    ap.add_argument("--output-json",required=True,type=Path)
    ap.add_argument("--output-fasta",required=True,type=Path)
    args=ap.parse_args()

    rows=[]
    with args.hits.open(encoding="utf-8") as fh:
        for r in csv.reader(fh,delimiter="\t"):
            if not r: continue
            rows.append({
                "query":r[0],"subject":r[1],"pident":float(r[2]),"aln_len":int(r[3]),
                "qlen":int(r[4]),"slen":int(r[5]),"evalue":float(r[6]),
                "bitscore":float(r[7]),"qcov":float(r[8]),"scov":float(r[9])
            })

    byq=defaultdict(list)
    for r in rows: byq[r["query"]].append(r)
    for q in byq: byq[q].sort(key=lambda x:(-x["bitscore"],x["evalue"]))

    missing=[q for q in QUERY_CLASS if not byq.get(q)]
    if missing: raise SystemExit(f"no DIAMOND hit for required queries: {missing}")

    top={q:byq[q][0] for q in QUERY_CLASS}
    dfr_subjects=[top["DFR_GERBERA_P51105"]["subject"],top["DFR_ARABIDOPSIS_P51102"]["subject"]]
    ans_subjects=[top["ANS_PETUNIA_P51092"]["subject"],top["ANS_ARABIDOPSIS_Q96323"]["subject"]]
    dfr_consensus=(dfr_subjects[0] if len(set(dfr_subjects))==1 else None)
    ans_consensus=(ans_subjects[0] if len(set(ans_subjects))==1 else None)
    anr_decoy=top["DECOY_ANR_ARABIDOPSIS_Q9SEV0"]["subject"]
    fls_decoy=top["DECOY_FLS1_ARABIDOPSIS_Q96330"]["subject"]

    def subject_hit(query,subject):
        for h in byq[query]:
            if h["subject"]==subject:
                return h
        return None

    def top_second_ratio(query):
        if len(byq[query])<2 or byq[query][1]["bitscore"]<=0:
            return None
        return byq[query][0]["bitscore"]/byq[query][1]["bitscore"]

    def target_decoy_ratio(target_queries,decoy_query,subject):
        th=[subject_hit(q,subject) for q in target_queries]
        dh=subject_hit(decoy_query,subject)
        if any(x is None for x in th) or dh is None or dh["bitscore"]<=0:
            return None
        return min(x["bitscore"] for x in th)/dh["bitscore"]

    dfr_queries=["DFR_GERBERA_P51105","DFR_ARABIDOPSIS_P51102"]
    ans_queries=["ANS_PETUNIA_P51092","ANS_ARABIDOPSIS_Q96323"]
    dfr_top2={q:top_second_ratio(q) for q in dfr_queries}
    ans_top2={q:top_second_ratio(q) for q in ans_queries}
    dfr_decoy_ratio=target_decoy_ratio(dfr_queries,"DECOY_ANR_ARABIDOPSIS_Q9SEV0",dfr_consensus) if dfr_consensus else None
    ans_decoy_ratio=target_decoy_ratio(ans_queries,"DECOY_FLS1_ARABIDOPSIS_Q96330",ans_consensus) if ans_consensus else None

    dfr_quality=bool(dfr_consensus and all(
        top[q]["pident"]>=60 and top[q]["qcov"]>=75 and (dfr_top2[q] or 0)>=1.5
        for q in dfr_queries
    ))
    ans_quality=bool(ans_consensus and all(
        top[q]["pident"]>=60 and top[q]["qcov"]>=75 and (ans_top2[q] or 0)>=1.5
        for q in ans_queries
    ))
    dfr_strong=bool(dfr_quality and dfr_decoy_ratio is not None and dfr_decoy_ratio>=1.5)
    ans_strong=bool(ans_quality and ans_decoy_ratio is not None and ans_decoy_ratio>=1.5)

    seqs=dict(read_fasta(args.proteins))
    chosen={}
    if dfr_consensus and dfr_consensus in seqs: chosen["DFR"]=dfr_consensus
    if ans_consensus and ans_consensus in seqs: chosen["ANS"]=ans_consensus

    args.output_fasta.parent.mkdir(parents=True,exist_ok=True)
    with args.output_fasta.open("w",encoding="utf-8") as out:
        for cls,sid in chosen.items():
            out.write(f">{cls}|{sid}\n")
            s=seqs[sid]
            for i in range(0,len(s),80): out.write(s[i:i+80]+"\n")

    def compact(x):
        return {k:x[k] for k in ("subject","pident","aln_len","qlen","slen","evalue","bitscore","qcov","scov")}

    receipt={
        "result_version":"cnip_anthocyanin_candidate_identification_v1",
        "query_top_hits":{q:compact(top[q]) for q in QUERY_CLASS},
        "DFR":{
            "reference_top_subjects":dfr_subjects,
            "consensus_subject":dfr_consensus,
            "ANR_decoy_top_subject":anr_decoy,
            "target_top_to_second_bitscore_ratio":dfr_top2,
            "minimum_target_to_same_subject_ANR_bitscore_ratio":dfr_decoy_ratio,
            "sequence_thresholds":{"pident_min":60,"qcov_min":75,"top_to_second_bitscore_ratio_min":1.5,"target_to_decoy_bitscore_ratio_min":1.5},
            "classification":"SEQUENCE_ORTHOLOGY_CANDIDATE_STRONG" if dfr_strong else "AMBIGUOUS"
        },
        "ANS":{
            "reference_top_subjects":ans_subjects,
            "consensus_subject":ans_consensus,
            "FLS_decoy_top_subject":fls_decoy,
            "target_top_to_second_bitscore_ratio":ans_top2,
            "minimum_target_to_same_subject_FLS_bitscore_ratio":ans_decoy_ratio,
            "sequence_thresholds":{"pident_min":60,"qcov_min":75,"top_to_second_bitscore_ratio_min":1.5,"target_to_decoy_bitscore_ratio_min":1.5},
            "classification":"SEQUENCE_ORTHOLOGY_CANDIDATE_STRONG" if ans_strong else "AMBIGUOUS"
        },
        "candidate_fasta_records":chosen,
        "boundary":"This screen nominates C. nipponicum protein candidates. Reference convergence, top-hit separation and target-vs-decoy score margins reduce family-level misassignment but do not by themselves prove biochemical function in Cirsium."
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__": main()
