#!/usr/bin/env python3
"""Download immutable publisher-authorized supplementary leaf/capitulum morphometrics.

This is an independent source-acquisition step. Data are retained only in
a GitHub Actions artifact, not committed, and are not yet analyzed.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen
from zipfile import ZipFile
from io import BytesIO
from datetime import datetime, timezone

URL = "https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs00606-023-01854-2/MediaObjects/606_2023_1854_MOESM6_ESM.xlsx"
DOI = "10.1007/s00606-023-01854-2"
DATA_FILENAME = "michalkova_2023_resource6_leaf_capitulum.xlsx"

def main():
    dest=Path("work/leaf_capitulum_published_source")
    dest.mkdir(parents=True,exist_ok=True)
    req=Request(URL,headers={"User-Agent":"Mozilla/5.0 (research data replication)","Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"})
    with urlopen(req,timeout=65) as rsp:
        content=rsp.read()
        final_url=rsp.geturl()
        status=getattr(rsp,"status",200)
    if status!=200 or not content.startswith(b"PK\x03\x04"):
        raise ValueError(f"Supplement is not a valid XLSX transport, status={status}; n={len(content)}")
    with ZipFile(BytesIO(content)) as z:
        paths=z.namelist()
        if "[Content_Types].xml" not in paths or not any(x.startswith("xl/worksheets/sheet") for x in paths):
            raise ValueError("Received ZIP is not an XLSX workbook")
        sheet_paths=[x for x in paths if x.startswith("xl/worksheets/sheet")]
    file=dest/DATA_FILENAME
    file.write_bytes(content)
    sha=hashlib.sha256(content).hexdigest()
    summary={"download_status":"PASS_VERIFIED_XLSX","source_title":"Hybridization may endanger the rare North Apennine endemic Cirsium bertolonii","source_doi":DOI,"supplement":"Online Resource 6: primary morphometric character measurements except achenes","source_url":URL,"resolved_url":final_url,"retrieved_at_utc":datetime.now(timezone.utc).isoformat(),"filename":DATA_FILENAME,"size_bytes":len(content),"sha256":sha,"excel_sheets":len(sheet_paths),"input_provenance_only":True,"claim_ceiling":"Not yet a statistically validated leaf vs capitulum variability test"}
    (dest/"source_manifest.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)
if __name__=="__main__":
    main()
