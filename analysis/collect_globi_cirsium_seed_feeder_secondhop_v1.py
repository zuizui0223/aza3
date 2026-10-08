#!/usr/bin/env python3
"""Second-hop parasitoid discovery for head-feeding Cirsium insects.

Species-insect edges and insect-parasitoid edges must not be treated as if
the parasitoid attacks that insect in the same Cirsium head without a primary
study confirming location and host association. Not a population network.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import io
import json
import sys
import time
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urlencode

URL="https://api.globalbioticinteractions.org/interaction.csv"
TYPES={"parasitoidof", "hasparasitoid"}
COLS=["head_seed_feeder","parasitoid","source_relation","source_taxon",
      "target_taxon","dataset_id","source_reference","study_title","study_doi",
      "study_external_id","event_date","source_page_digest","provenance_key",
      "same_cirsium_head_confirmed","ecological_conclusion"]

def txt(row,*names):
    for name in names:
        z=str(row.get(name,"") or "").strip()
        if z:
            return " ".join(z.split())
    return ""

def normalize(t):
    return "".join(ch for ch in str(t).lower() if ch.isalnum())

def exact_species_prefix(name,taxon):
    """Only same binomial, permitting authorship suffix, not a congeneric host."""
    parts=name.strip().casefold().split()
    host=taxon.strip().casefold().split()
    return len(parts)>=2 and len(host)==2 and parts[:2]==host

def parse_record(row,host,page_sha):
    src=txt(row,"source_taxon_name","sourceTaxonName")
    tgt=txt(row,"target_taxon_name","targetTaxonName")
    relation=txt(row,"interaction_type","interactionTypeName")
    rel=normalize(relation)
    if rel not in TYPES:
        return None
    if rel=="parasitoidof":
        parasitoid, insect=src,tgt
    else:
        parasitoid, insect=tgt,src
    if not exact_species_prefix(insect,host):
        return None
    provenance=txt(row,"study_source_id","studySourceId","study_external_id",
                   "studyExternalId","source_namespace","sourceNamespace")
    study_title=txt(row,"study_title","studyTitle")
    reference=txt(row,"study_source_citation","studySourceCitation",
                  "study_citation","studyCitation","study_url","studyUrl") or study_title
    row_arg=txt(row,"argument_type","argumentTypeName",
                "argument_type_name","argumentType","argument_type_id","argumentTypeId")
    if "refut" in " ".join((provenance,reference,row_arg)).lower():
        return None
    doi=txt(row,"study_doi","studyDoi","study_source_doi","studySourceDoi")
    event=txt(row,"event_date","eventDate")
    if not parasitoid:
        return None
    # A host-parasitoid edge alone does not establish occurrence in a given plant.
    key=hashlib.sha256("|".join((host,parasitoid,rel,provenance,reference,
                                doi,event)).encode()).hexdigest()[:24]
    return {
      "head_seed_feeder":host,"parasitoid":parasitoid,
      "source_relation":relation, "source_taxon":src,"target_taxon":tgt,
      "dataset_id":provenance,"source_reference":reference,
      "study_title":study_title, "study_doi":doi,
      "study_external_id":txt(row,"study_external_id","studyExternalId"),
      "event_date":event,"source_page_digest":page_sha,
      "provenance_key":key,
      "same_cirsium_head_confirmed":"no",
      "ecological_conclusion":"HOST_PARASITOID_RECORD_ONLY__HEAD_CONTEXT_REQUIRES_PRIMARY_SOURCE"
    }

def fetch(host,side,limit,offset,timeout):
    key="sourceTaxon" if side=="source" else "targetTaxon"
    url=URL+"?"+urlencode({key:host,"includeObservations":"true",
                           "limit":limit,"offset":offset})
    last=None
    for i in range(3):
        try:
            request=Request(url,headers={"User-Agent":"aza3-Cirsium-secondhop/1.0",
                                         "Accept":"text/csv"})
            with urlopen(request,timeout=timeout) as response:
                raw=response.read()
            reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")))
            fields=set(reader.fieldnames or [])
            if not (({"sourceTaxonName","targetTaxonName","interactionTypeName"} <= fields) or
                    ({"source_taxon_name","target_taxon_name","interaction_type"} <= fields)):
                raise ValueError("No valid GloBI CSV source/target/type header")
            return list(reader),url,hashlib.sha256(raw).hexdigest()
        except Exception as e:
            last=str(e)
            if i<2:time.sleep(1.5*(i+1))
    raise RuntimeError(f"{host} {side} offset={offset}: {last}")

def run(a):
    anchor=Path(a.anchors)
    original=list(csv.DictReader(anchor.open(encoding="utf-8-sig")))
    hosts=sorted({r["partner"] for r in original
                 if r.get("role")=="head_seed_feeder"
                 and r.get("organ")=="capitulum"
                 and len(r["partner"].split())==2})
    if not hosts:
        raise SystemExit("No two-token head seed-feeders in frozen source anchors")
    out=Path(a.output_dir)
    out.mkdir(parents=True,exist_ok=True)
    records={}
    errors=[]
    ledger=[]
    ends={}
    for host in hosts:
        for side in ("source","target"):
            exhausted=False
            for page in range(a.max_pages):
                offset=page*a.page_size
                try:
                    rows,url,digest=fetch(host,side,a.page_size,offset,a.timeout)
                except RuntimeError as exc:
                    errors.append(str(exc));break
                valid=0
                for raw in rows:
                    candidate=parse_record(raw,host,digest)
                    if candidate is None:continue
                    records[candidate["provenance_key"]]=candidate
                    valid+=1
                ledger.append({"host":host,"side":side,"offset":offset,"rows":len(rows),
                               "candidate_rows":valid,"page_sha256":digest,"url":url})
                print(f"host={host} {side} offset={offset} rows={len(rows)} hit={valid}",flush=True)
                if len(rows)<a.page_size:
                    exhausted=True
                    break
                time.sleep(a.pause)
            ends[f"{host}|{side}"]="EXHAUSTED" if exhausted else "NOT_EXHAUSTED"
    output=out/"seed_feeder_parasitoid_candidate_edges.csv"
    with output.open("w",encoding="utf-8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=COLS)
        wr.writeheader()
        wr.writerows(sorted(records.values(),key=lambda r:(r["head_seed_feeder"],
                                                  r["parasitoid"],r["provenance_key"])))
    (out/"secondhop_query_ledger.json").write_text(json.dumps(ledger,indent=2)+"\n")
    status="QUERY_EXHAUSTED" if all(x=="EXHAUSTED" for x in ends.values()) and not errors else "PARTIAL_OR_FAILED"
    summary={"status":status,"hosts":hosts,"n_hosts":len(hosts),
       "n_host_parasitoid_candidate_rows":len(records),
       "distinct_host_parasitoid_name_pairs":len({(r["head_seed_feeder"],r["parasitoid"])
                                                  for r in records.values()}),
       "errors":errors,"terminal":ends,
       "file_sha256":hashlib.sha256(output.read_bytes()).hexdigest(),
       "decision_ceiling":"Parasitoid–seed-feeder association, no Cirsium head-localized triad without common primary-source evidence; no ecological effect, defence or adaptation inferred",
       "negative_inference":"No API record is not proof of absent parasitism",
       "head_context_verified_rows":0}
    (out/"secondhop_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)
    if errors or not records:raise SystemExit(2)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--anchors",default="data/evidence/cirsium_capitulum_interaction_source_anchors_v1.csv")
    p.add_argument("--output-dir",default="outputs/globi_cirsium_secondhop_v1")
    p.add_argument("--page-size",type=int,default=200)
    p.add_argument("--max-pages",type=int,default=3)
    p.add_argument("--pause",type=float,default=0.35)
    p.add_argument("--timeout",type=int,default=40)
    a=p.parse_args()
    if not (1<=a.page_size<=1024 and 1<=a.max_pages<=100 and a.pause>=0):
        p.error("invalid request size")
    run(a)

if __name__=="__main__":
    main()
