#!/usr/bin/env python3
"""Read-only public metadata probe for the Karpate et al. Cirsium dataset.

No secrets or access bypass; two fixed Zenodo public URLs, bounded retries.
Not an analysis of plant-level data; no zero values inferred from errors.

Preprint DOI: 10.21203/rs.3.rs-6634172/v1 (posted June 2025)
Archived manuscript data/code DOI claim: 10.5281/zenodo.15364152
"""
from __future__ import annotations
import argparse
import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

RECORD_ID=15364152
URLS=(
  f"https://zenodo.org/api/records/{RECORD_ID}",
  f"https://zenodo.org/records/{RECORD_ID}/export/json",
)

def json_probe(url):
    req=urllib.request.Request(url,headers={
      "Accept":"application/json",
      "User-Agent":"CirsiumPublicEvidenceAudit/1.0 (public metadata; read only)",
    })
    for retry in range(2):
        try:
            with urllib.request.urlopen(req,timeout=17) as response:
                raw=response.read(1000000)
                code=response.status
            data=json.loads(raw.decode("utf-8"))
            return {"status":"PUBLIC_METADATA_RETRIEVED","http_status":code,
                    "payload":data}
        except urllib.error.HTTPError as ex:
            if ex.code in (429,500,502,503,504) and retry==0:
                time.sleep(1.6)
                continue
            return {"status":"METADATA_UNAVAILABLE_NOT_ZERO_EXPERIMENT",
                    "http_status":ex.code,"failure_type":"HTTPError"}
        except (OSError,TimeoutError,UnicodeError,ValueError,
                json.JSONDecodeError) as ex:
            return {"status":"METADATA_UNAVAILABLE_NOT_ZERO_EXPERIMENT",
                    "http_status":None,"failure_type":type(ex).__name__}
    return {"status":"METADATA_UNAVAILABLE_NOT_ZERO_EXPERIMENT",
            "http_status":None,"failure_type":"RETRIES_EXHAUSTED"}

def normalized(data):
    if not isinstance(data,dict):
        raise ValueError("ZENODO_PUBLIC_RECORD_NOT_OBJECT")
    root=data.get("data",data)
    if not isinstance(root,dict) or root.get("id")!=RECORD_ID:
        raise ValueError("ZENODO_RECORD_ID_DOES_NOT_MATCH_DECLARED_SOURCE")
    meta=root.get("metadata",{})
    if not isinstance(meta,dict):
        raise ValueError("ZENODO_METADATA_IS_NOT_OBJECT")
    dois={str(root.get("doi") or "").lower(),str(meta.get("doi") or "").lower()}
    if "10.5281/zenodo.15364152" not in dois:
        raise ValueError("ZENODO_DOI_NOT_MATCHING_MANUSCRIPT_REFERENCE")
    files=root.get("files",[])
    if not isinstance(files,list):
        raise ValueError("ZENODO_FILES_PAYLOAD_INVALID")
    result=[]
    for f in files:
        if not isinstance(f,dict):
            raise ValueError("ZENODO_FILE_METADATA_INVALID")
        key=f.get("key") or f.get("filename") or ""
        if not isinstance(key,str):
            raise ValueError("ZENODO_FILENAME_INVALID")
        links=f.get("links") or {}
        result.append({
          "name":key,
          "bytes":f.get("size") or f.get("filesize"),
          "checksum":f.get("checksum"),
          "metadata_is_public":True,
          "direct_file_download_tested":False,
          "download_url_status":"NOT_TESTED",
          "public_file_url":links.get("self") if isinstance(links,dict) else None,
        })
    return {
        "record_id":root["id"],
        "doi":"10.5281/zenodo.15364152",
        "title":meta.get("title"),
        "license":meta.get("license"),
        "access_right":meta.get("access_right") or
            ((root.get("access") or {}).get("record") or
             (root.get("access") or {}).get("files")),
        "publication_date":meta.get("publication_date"),
        "n_files_in_metadata":len(result),
        "file_manifest":result,
        "real_analyzed_Cirsium_head_rows":0,
        "original_recomputed_experimental_contrasts":0,
        "head_spine_orientation_effect":"NOT_MEASURED_BY_PREPRINT_DESIGN",
    }

def test_synthetic():
    example={"id":RECORD_ID,"doi":"10.5281/zenodo.15364152",
             "metadata":{"title":"synthetic test only"},
             "files":[{"key":"synthetic.csv","size":12,
                       "checksum":"md5:synthetic"}]}
    o=normalized(example)
    assert o["record_id"]==RECORD_ID
    assert o["file_manifest"][0]["name"]=="synthetic.csv"
    assert o["real_analyzed_Cirsium_head_rows"]==0
    for poison in [
        {**example,"id":17},
        {**example,"doi":"10.5281/zenodo.OTHER"},
        {**example,"files":"not-an-array"},
    ]:
        try:normalized(poison)
        except ValueError: pass
        else:raise AssertionError("UNVERIFIED_OR_INVALID_ZENODO_SOURCE_ACCEPTED")
    return {"status":"PASS_SYNTHETIC_ZENODO_METADATA_SOURCE_CHECKS",
            "invalid_source_records_rejected":3,"actual_data_rows_inserted":0}

def run(live):
    out={
       "version":"cirsium_arvense_2025_zenodo_public_access_audit_v1",
       "as_of_UTC":datetime.now(timezone.utc).isoformat(),
       "source_preprint_doi":"10.21203/rs.3.rs-6634172/v1",
       "source_claimed_zenodo_doi":"10.5281/zenodo.15364152",
       "source_type":"2025 RESEARCH SQUARE PREPRINT; PEER REVIEW NOT VERIFIED",
       "synthetic_checks":test_synthetic(),
       "source_access_status":"NOT_QUERIED" if not live else "METADATA_UNAVAILABLE_NOT_ZERO_EXPERIMENT",
       "attempts":[],
       "verified_source_metadata":None,
       "no_data_downloaded_or_plant_reanalysis_claimed":True,
       "biological_preprint_reports":{
           "experiment_sites":4,
           "treatment_names":["control","predator_exclusion",
                              "pollinator_exclusion","predator_model"],
           "qualitative_predator_treatment_fitness_effect":
                "NO_DETECTABLE_DIFFERENCE_REPORTED",
           "pollinator_exclusion_effect_on_fertilization":
                "FERTILIZED_PROPORTION_REDUCED",
           "pollinator_exclusion_effect_on_emerged_florivores":
                "FEWER_EMERGED_INSECTS_REPORTED",
           "pollinator_exclusion_effect_on_fraction_predated":
                "HIGHER_FRACTION_REPORTED",
           "no_direct_head_armature_orientation_treatment":True,
       }
    }
    if not live:return out
    for u in URLS:
        ans=json_probe(u)
        audit={k:v for k,v in ans.items() if k!="payload"}
        audit["url"]=u
        out["attempts"].append(audit)
        if ans["status"]=="PUBLIC_METADATA_RETRIEVED":
            try:meta=normalized(ans["payload"])
            except ValueError as err:
                out["source_access_status"]="PUBLIC_METADATA_SOURCE_MISMATCH"
                out["validation_error"]=str(err)
                break
            out["source_access_status"]="ORIGINAL_ARCHIVE_METADATA_VERIFIED_NO_RAW_HEAD_REANALYSIS"
            out["verified_source_metadata"]=meta
            break
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--live",action="store_true")
    a=ap.parse_args()
    x=run(a.live)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,indent=2,ensure_ascii=False)+"\n",
                      encoding="utf-8")
    print(json.dumps({
        "status":x["source_access_status"],
        "queries":x["attempts"],
        "metadata":x["verified_source_metadata"],
        "real_field_rows_recomputed":0,
        "not_a_physical_spine_treatment":True
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
