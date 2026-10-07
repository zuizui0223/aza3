#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,statistics
from collections import defaultdict
from pathlib import Path

def read_fasta(path):
    name=None; seq=[]
    with path.open(encoding="utf-8") as fh:
        for raw in fh:
            line=raw.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None: yield name,"".join(seq).replace("*","")
                name=line[1:].split()[0]; seq=[]
            else: seq.append(line)
    if name is not None: yield name,"".join(seq).replace("*","")

def union_coverage(intervals,length):
    if not intervals or length<=0: return 0.0
    xs=sorted((min(a,b),max(a,b)) for a,b in intervals)
    merged=[]
    for a,b in xs:
        if not merged or a>merged[-1][1]+1: merged.append([a,b])
        else: merged[-1][1]=max(merged[-1][1],b)
    covered=sum(b-a+1 for a,b in merged)
    return covered/length

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--hits",required=True,type=Path)
    ap.add_argument("--candidates",required=True,type=Path)
    ap.add_argument("--ena-metadata",required=True,type=Path)
    ap.add_argument("--sample-code",required=True)
    ap.add_argument("--run",required=True)
    ap.add_argument("--morph",required=True)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()

    seqs=dict(read_fasta(args.candidates))
    target_map={}
    for sid,seq in seqs.items():
        target=sid.split("|",1)[0]
        target_map[sid]=(target,len(seq))

    accepted=defaultdict(list)
    if args.hits.exists():
        with args.hits.open(encoding="utf-8") as fh:
            for r in csv.reader(fh,delimiter="\t"):
                if not r: continue
                q,sid=r[0],r[1]
                if sid not in target_map: continue
                pident=float(r[2]); aln=int(r[3]); evalue=float(r[4]); bits=float(r[5])
                sstart=int(r[6]); send=int(r[7]); slen=int(r[8])
                if pident<75 or aln<30 or evalue>1e-10: continue
                target,_=target_map[sid]
                accepted[target].append({
                    "query":q,"pident":pident,"aln_aa":aln,"evalue":evalue,
                    "bitscore":bits,"sstart":sstart,"send":send,"slen":slen
                })

    result={}
    for target in ("DFR","ANS"):
        hs=accepted.get(target,[])
        reads=sorted({h["query"] for h in hs})
        if hs:
            slen=max(h["slen"] for h in hs)
            cov=union_coverage([(h["sstart"],h["send"]) for h in hs],slen)
            med_id=statistics.median(h["pident"] for h in hs)
            max_bits=max(h["bitscore"] for h in hs)
        else:
            slen=len(next((seq for sid,seq in seqs.items() if sid.startswith(target+"|")),""))
            cov=0.0; med_id=None; max_bits=None
        n=len(reads)
        if n>=20 and cov>=0.50: tier="TRANSCRIPT_STRONG"
        elif n>=5 and cov>=0.20: tier="TRANSCRIPT_POSITIVE"
        else: tier="NO_POSITIVE_EVIDENCE"
        result[target]={
            "tier":tier,"unique_read_hits":n,"subject_union_coverage":cov,
            "median_identity":med_id,"max_bitscore":max_bits,"candidate_length_aa":slen
        }

    meta=json.loads(args.ena_metadata.read_text(encoding="utf-8"))
    out={
        "result_version":"takaoense_public_pathway_retention_sample_v1",
        "sample_code":args.sample_code,"run":args.run,"morph":args.morph,
        "tissue":"young_leaf","targets":result,"ena_fastq_metadata":meta,
        "positive_filter":{"pident_min":75,"aligned_aa_min":30,"evalue_max":1e-10},
        "boundary":"Positive transcript evidence supports presence of a homologous expressed coding sequence in this young-leaf library. NO_POSITIVE_EVIDENCE is not gene deletion, pathway loss, or floral non-expression."
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__": main()
