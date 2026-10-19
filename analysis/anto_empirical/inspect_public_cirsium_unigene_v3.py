# Empirical source recovery adopted from zuizui0223/anto/scripts/inspect_cirsium_japonicum_unigene_annotation_v6.py; canonical ownership aza3.
"""Extract verifiable *Cirsium japonicum* unigene ANNOTATION candidates, NOT ortholog calls.

Park et al. 2020 Supporting Information includes:
  si_002.txt: 104,890 unigene FASTA headers preceded by author cover text
  si_003.xlsx: 104,890 unigene information / BLAST annotations
The public Europe PMC supplemental archive can be downloaded as a bounded ZIP.

Data logic: do not infer any Ryukyu brevicaule/irumtiense sequence, expression,
flower colour, DFR activity or gene age. This script only links 2020 source
unigene FASTA IDs to publisher-annotated candidate rows when exact IDs match.
All hits are ORTHOLOGY_UNVERIFIED. No executable workbook formulas/macros run.
"""
from collections import Counter
from io import BytesIO
import argparse
import hashlib
import json
import re
import sys
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7675961/supplementaryFiles"
SI_TEXT = "ao0c04001_si_002.txt"
SI_XLSX = "ao0c04001_si_003.xlsx"
NS = {"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
EXACT = re.compile(r"\b(?:dihydroflavonol[-\s]?4[-\s]?reductase|anthocyanidin synthase|leucoanthocyanidin dioxygenase|(?<!\w)DFR(?!\w)|(?<!\w)ANS(?!\w))\b", re.I)
MAX_BYTES = 50_000_000
MAX_EXPANDED = 180_000_000


def _find_text(element,tag):
    x=element.find(tag,NS)
    return x.text if x is not None else None


def _string_table(z):
    """Read bounded shared strings; preserve original cell content."""
    with z.open("xl/sharedStrings.xml") as stream:
        strings=[]
        for _,el in ET.iterparse(stream,events=("end",)):
            if el.tag.endswith("}si"):
                text="".join(t.text or "" for t in el.iter() if t.tag.endswith("}t"))
                strings.append(text)
                el.clear()
                if len(strings)>600000:
                    raise ValueError("Too many distinct string-table entries")
    return strings


def _decode_cell(c,shared):
    t=c.attrib.get("t")
    v=_find_text(c,"m:v")
    if t=="s" and v is not None:
        n=int(v)
        if n<0 or n>=len(shared):
            raise ValueError("Shared string index is not valid")
        return shared[n]
    if t=="inlineStr":
        return "".join(x.text or "" for x in c.iter() if x.tag.endswith("}t"))
    return v or ""


def _rows(xlsx):
    """Yield row number and {cell reference: text}, one row at a time."""
    with zipfile.ZipFile(BytesIO(xlsx),"r") as z:
        if "xl/sharedStrings.xml" not in z.namelist() or "xl/worksheets/sheet1.xml" not in z.namelist():
            raise ValueError("Published supplemental workbook format changed")
        total=sum(x.file_size for x in z.infolist())
        if total>MAX_EXPANDED:
            raise ValueError("Expanded workbook size beyond bound")
        shared=_string_table(z)
        with z.open("xl/worksheets/sheet1.xml") as stream:
            for event,el in ET.iterparse(stream,events=("end",)):
                if el.tag.endswith("}row"):
                    yield int(el.attrib.get("r","0")), {
                        c.attrib.get("r",""): _decode_cell(c,shared)
                        for c in el.findall("m:c",NS)
                    }
                    el.clear()


def worksheet_report(xlsx):
    """Inspect labels first; annotation hit≠gene function or causal orthology."""
    first_rows=[]
    hits=[]
    row_count=0
    for row_number,cells in _rows(xlsx):
        row_count+=1
        if len(first_rows)<7:
            first_rows.append({
                "row":row_number,
                "cells":{k:v[:145] for k,v in list(cells.items())[:15]}
            })
        positives=[(k,v) for k,v in cells.items() if EXACT.search(v)]
        if positives and len(hits)<150:
            left=sorted(cells.items(),key=lambda item: (
                len(re.match(r"[A-Z]+",item[0]).group(0)) if re.match(r"[A-Z]+",item[0]) else 100,
                item[0]))[:8]
            hits.append({
                "row":row_number,
                "annotation_cells":{k:v[:260] for k,v in positives[:5]},
                "adjacent_cells":{k:v[:180] for k,v in left},
                "no_validated_orthology":True
            })
    return {"workbook_row_count":row_count,
            "first_rows":first_rows,
            "candidate_annotation_hit_rows_capped":hits,
            "candidate_annotation_row_count_at_least":len(hits),
            "hit_list_truncated_at_150":len(hits)==150}


def fasta_report(text_bytes):
    """Count records and sample header IDs without materializing 104k sequences."""
    records=0
    length_total=0
    headers=[]
    active=False
    for line in BytesIO(text_bytes):
        item=line.strip()
        if item.startswith(b">"):
            active=True
            records+=1
            if len(headers)<7:
                headers.append(item.decode("utf-8",errors="replace")[:200])
            continue
        if active and item and not item.startswith(b"#"):
            # Source FASTA was wrapped; count literal nucleotide content only.
            s=item.decode("ascii",errors="ignore").upper()
            if all(ch in "ACGTUNRYKMSWBDHV-" for ch in s):
                length_total+=len(s)
    return {"n_fasta_headers":records,"first_fasta_headers":headers,
            "total_length_nucleotide_compatible_lines":length_total,
            "note":"Sequence count and format only: source FASTA gene IDs still require exact workbook linkage"}


def inspect_public_archive(data):
    if len(data)>MAX_BYTES or not zipfile.is_zipfile(BytesIO(data)):
        raise ValueError("Source is not a bounded ZIP file")
    with zipfile.ZipFile(BytesIO(data),"r") as z:
        names={p.filename:p for p in z.infolist()}
        if sum(p.file_size for p in z.infolist())>MAX_EXPANDED:
            raise ValueError("Public supplementary archive outside decompression budget")
        text=next((p for p in z.infolist() if p.filename.endswith(SI_TEXT)),None)
        xlsx=next((p for p in z.infolist() if p.filename.endswith(SI_XLSX)),None)
        if text is None or xlsx is None:
            raise ValueError("Expected 2020 Cirsium unigene supplementary sources missing")
        fa=fasta_report(z.read(text))
        ann=worksheet_report(z.read(xlsx))
        return {"publisher_archived_study":"10.1021/acsomega.0c04001",
                "focal_ryukyu_data":False,
                "source_zip_sha256":hashlib.sha256(data).hexdigest(),
                "source_text_uncompressed_bytes":text.file_size,
                "source_xlsx_uncompressed_bytes":xlsx.file_size,
                "fasta_evidence":fa,
                "workbook_annotation_evidence":ann,
                "genuine_orthologs_identified":False,
                "reexpression_confirmed":False,
                "claim_ceiling":"Author source unigene/BLAST candidate IDs only, if exact corroboration succeeds."}


def download(max_bytes=MAX_BYTES,timeout=40):
    req=urllib.request.Request(URL,headers={
        "User-Agent":"anto-literature-audit/1.0",
        "Accept":"application/zip,application/octet-stream"})
    with urllib.request.urlopen(req,timeout=timeout) as f:
        expected=f.headers.get("Content-Length")
        if expected and int(expected)>max_bytes:
            raise ValueError("Source download exceeds predeclared budget")
        content=f.read(max_bytes+1)
        if len(content)>max_bytes:
            raise ValueError("Source exceeds predeclared download budget")
        return content


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fetch-public",action="store_true")
    ap.add_argument("--archive",help="local copy of legally available Europe PMC supplementary ZIP")
    ap.add_argument("--output")
    args=ap.parse_args()
    if bool(args.archive)==bool(args.fetch_public):
        raise ValueError("Choose exactly one of --fetch-public or --archive")
    try:
        b=open(args.archive,"rb").read() if args.archive else download()
        result={"status":"INSPECTED", **inspect_public_archive(b)}
    except (urllib.error.URLError,urllib.error.HTTPError,OSError,ValueError,
            zipfile.BadZipFile,ET.ParseError) as exc:
        result={"status":"NOT_INSPECTED", "error_type":type(exc).__name__,
                "detail":str(exc)[:300],
                "not_evidence_of_gene_absence":True,
                "no_Ryukyu_gene_claims":True}
    result["analysis_version"]="anto_v6_source_unigene_introspection"
    content=json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n"
    if args.output:
        from pathlib import Path
        path=Path(args.output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(content,encoding="utf-8")
    print(content)


if __name__=="__main__":
    main()
