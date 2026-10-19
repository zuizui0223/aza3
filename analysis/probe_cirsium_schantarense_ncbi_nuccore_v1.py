#!/usr/bin/env python3
"""Attempt a public NCBI Nucleotide accession discovery without guessing IDs.

The scientific target is 2015 Bae Cirsium shantarense (schantarense)
NUCLEAR rDNA ITS. Cirsium SHANSIENSE plastome MN871982 must not be
treated as the same accession/taxon. Discovery via NCBI E-utilities
cannot independently validate Bae's specimen unless source crosswalk
and voucher/publication metadata agree.

This script explicitly distinguishes:
  NOT_QUERIED / NCBI_UNAVAILABLE / ZERO_INDEXED_HITS_FOR_QUERY /
  PUBLIC_CANDIDATE_NOT_LINKED_TO_BAE_2015.
No output ever claims a verified 2015 accession without a verified
publication/biosample/voucher join (not yet implemented).
"""
from __future__ import annotations
import argparse
import json
import os
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

NCBI="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
QUERIES={
    "schantarense_exact_text":"\"Cirsium schantarense\"[All Fields]",
    "shantarense_paper_spelling":"\"Cirsium shantarense\"[All Fields]",
    "schantarense_and_ITS":"\"Cirsium schantarense\"[All Fields] AND (\"ITS\"[All Fields] OR \"internal transcribed spacer\"[All Fields])",
}
MAX_PER_QUERY=30
TIMEOUT=18

def fetch_json(path,params):
    url=NCBI+path+"?"+urllib.parse.urlencode(params)
    req=urllib.request.Request(url,headers={
        "User-Agent":"CirsiumSourceProvenanceResearch/1.0 (read-only)",
        "Accept":"application/json",
    })
    with urllib.request.urlopen(req,timeout=TIMEOUT) as response:
        raw=response.read(1000000)
    return json.loads(raw.decode("utf-8"))

def live_query():
    result={
      "status":"NCBI_PUBLIC_EUTILS_READ_ONLY_SOURCE_DISCOVERY",
      "query_scope":"Public NCBI Nucleotide E-utilities (NOT raw Bae lab files)",
      "query_results":{},
      "discovery_not_a_verified_Bae_2015_accession":True,
      "known_wrong_candidate_excluded":"MN871982 is original Cirsium shansiense chloroplast genome",
      "NCBI_record_Bae_voucher_verified":False,
      "comparably_placed_nuclear_phylogeny_obtained":False,
    }
    ids=set()
    for name,q in QUERIES.items():
        try:
            item=fetch_json("esearch.fcgi",{
               "db":"nuccore","term":q,"retmode":"json",
               "retmax":MAX_PER_QUERY,
            })["esearchresult"]
            these=[str(x) for x in item.get("idlist",[])]
            ids.update(these)
            result["query_results"][name]={
                "status":("PUBLIC_INDEX_CANDIDATES_ONLY"
                          if item.get("count","0")!="0"
                          else "ZERO_INDEXED_HITS_FOR_QUERY_NOT_NO_DNA"),
                "source_total_count":int(item.get("count","0")),
                "returned_uids":these,
                "all_available_results_collected":(
                    int(item.get("count","0"))<=len(these)),
            }
        except (KeyError,ValueError,TimeoutError,urllib.error.URLError,
                OSError,json.JSONDecodeError) as err:
            result["query_results"][name]={
              "status":"NCBI_UNAVAILABLE_NOT_ZERO",
              "error_type":type(err).__name__,
              "error_detail":str(err)[:300],
              "source_total_count":None,
              "returned_uids":[],
            }
    result["n_unique_candidate_nuccore_uids"]=len(ids)
    if ids:
        try:
            info=fetch_json("esummary.fcgi",{
               "db":"nuccore","id":",".join(sorted(ids)),
               "retmode":"json",
            })["result"]
            candidates=[]
            for uid in info.get("uids",[]):
                item=info.get(uid,{})
                candidates.append({
                    "uid":uid,
                    "accession_version":item.get("accessionversion") or
                       item.get("caption"),
                    "record_title":item.get("title"),
                    "taxid":item.get("taxid"),
                    "length":item.get("slen"),
                    "is_proven_Bae_2015_material":False,
                })
            result["indexed_NCBI_candidates"]=candidates
            result["candidate_resolution_status"]="PUBLIC_RECORDS_UNVERIFIED_BAE_2015_SAMPLE_IDENTITY"
        except (ValueError,KeyError,urllib.error.URLError,
                TimeoutError,OSError,json.JSONDecodeError) as err:
            result["indexed_NCBI_candidates"]=[]
            result["candidate_resolution_status"]="NCBI_ESUMMARY_UNAVAILABLE"
            result["candidate_detail_error"]=type(err).__name__
    else:
        result["indexed_NCBI_candidates"]=[]
        result["candidate_resolution_status"]="NO_VISIBLE_CANDIDATES_OR_API_ACCESS"
    result["direct_verified_Bae_2015_schantarense_ITS_accessions"]=[]
    result["phylogenetic_memory_estimate"]=None
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--live",action="store_true")
    args=ap.parse_args()
    receipt=live_query() if args.live else {
        "status":"NOT_QUERIED",
        "query_results":{},
        "indexed_NCBI_candidates":[],
        "direct_verified_Bae_2015_schantarense_ITS_accessions":[],
        "phylogenetic_memory_estimate":None,
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",
                        encoding="utf-8")
    print(json.dumps({
       "status":receipt["status"],
       "searches":{k:{"status":v["status"],
                       "count":v["source_total_count"]}
                   for k,v in receipt["query_results"].items()},
       "candidate_titles":receipt["indexed_NCBI_candidates"],
       "Bae_2015_accession_directly_verified":0,
       "memory_estimate":None,
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
