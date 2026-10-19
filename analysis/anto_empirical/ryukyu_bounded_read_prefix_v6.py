"""Strictly capped *public* Ryukyu Cirsium 2026 young-leaf RNA pilot.

Reads ONLY the first N compressed bytes of each ENA R1 FASTQ.gz, not whole
libraries, not R2, and not an unbiased subsample of biological transcripts.

Aligns the fully decoded FASTQ-prefix records to 26 *different species* 2020
Cirsium japonicum BLAST-labelled pigment-family transcripts, and to 26
sequence-shuffled negatives, with minimap2 short-read preset.

A mapped read is only a putative similarity to one of 26 *source
transcripts*, NOT a gene-orthology call, intact toolkit, measured abundance,
floral anthocyanin, trait history or Ryukyu white-to-coloured re-expression.

No read data are committed to GitHub; only summary, download provenance and
false-positive controls. The optional live run must be run with --fetch-public.
Standard-library code apart from the external minimap2 executable.
"""
from collections import Counter,defaultdict
from hashlib import sha256
from pathlib import Path
from random import Random
import argparse
import gzip
import json
import re
import shutil
import subprocess
import tempfile
import urllib.error
import urllib.request
import zlib

ROOT=Path(__file__).resolve().parents[2]
MANIFEST=ROOT/"data/anto_source_handoff/public_sra/ryukyu_six_R1_prefix_source_v6.json"
REF=ROOT/"data/anto_source_handoff/molecular/cirsium_japonicum_2020_26_source_sequences_v6.fasta"
REF_META=ROOT/"data/anto_source_handoff/molecular/cirsium_japonicum_2020_26_sequence_provenance_v6.json"
DEFAULT_CAP=8*1024*1024
MAX_CAP=32*1024*1024
UNCOMPRESSED_CAP=160*1024*1024
STRICT_MIN_MAPQ=20
STRICT_MIN_ALIGNED=85
STRICT_MIN_IDENTITY=.85
WEAK_MIN_ALIGNED=65
WEAK_MIN_IDENTITY=.75
CIGAR=re.compile(r"(\d+)([MIDNSHP=X])")


def source_manifest(manifest=MANIFEST):
    data=json.loads(Path(manifest).read_text(encoding="utf-8"))
    if data.get("version")!="aza3_cirsium_Ryukyu_public_R1_prefix_manifest_v6":
        raise ValueError("Wrong public source manifest")
    runs=data.get("runs",[])
    if len(runs)!=6 or len({r["run_accession"] for r in runs})!=6:
        raise ValueError("Six distinct Ryukyu run accessions are required")
    if sum(x["taxon_as_experiment_titled"]=="Cirsium brevicaule" for x in runs)!=3:
        raise ValueError("The Cirsium brevicaule source count is not three")
    if sum(x["taxon_as_experiment_titled"]=="Cirsium irumtiense" for x in runs)!=3:
        raise ValueError("The Cirsium irumtiense source count is not three")
    for r in runs:
        run=r["run_accession"]
        if not re.fullmatch(r"SRR\d+",run):
            raise ValueError("Unexpected SRA Run ID")
        url=r["R1_https_url"]
        if (not url.startswith("https://ftp.sra.ebi.ac.uk/vol1/fastq/")
                or not url.endswith(f"/{run}_1.fastq.gz")):
            raise ValueError("URL does not match official ENA first mate")
        if (not re.fullmatch("[a-fA-F0-9]{32}",r["R1_provider_md5"]) or
                r["R1_full_file_bytes"]<1000000 or
                r["whole_run_complete"] or not r["mate2_not_retrieved"] or
                r["per_individual_corolla_pigmentation_quantified"] or
                r["tissue"]!="young leaf"):
            raise ValueError("Read provenance and experimental tissue inconsistent")
    return data


def read_fasta(path):
    rows={}
    name=None
    seq=[]
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                if name in rows:raise ValueError("Duplicate transcript FASTA ID")
                rows[name]="".join(seq).upper()
            name=line[1:].split()[0]
            seq=[]
        elif line.strip() and name is not None:
            seq.append(line.strip())
    if name is not None:
        if name in rows:raise ValueError("Duplicate transcript FASTA ID")
        rows[name]="".join(seq).upper()
    return rows


def prepare_reference(reference=REF, provenance=REF_META,seed=20261019):
    meta=json.loads(Path(provenance).read_text(encoding="utf-8"))
    if meta.get("contract")!="anto_park2020_26_source_transcripts_v6":
        raise ValueError("Wrong 2020 comparison transcript source")
    ref=read_fasta(reference)
    records=meta.get("records",[])
    if len(ref)!=26 or len(records)!=26 or set(ref)!={r["id"] for r in records}:
        raise ValueError("Mismatch of 26 independently published source transcript IDs")
    index={}
    lines=[]
    rng=Random(seed)
    for i,record in enumerate(sorted(records,key=lambda x:x["id"]),1):
        orig=ref[record["id"]].upper().replace("U","T").replace("-","N")
        if len(orig)!=record["source_fasta_length"]:
            raise ValueError("Published transcript nucleotide lengths no longer match")
        orig_unmodified=ref[record["id"]]
        if sha256(orig_unmodified.encode("ascii")).hexdigest()!=record["sha256_nucleotide"]:
            raise ValueError("Published source sequence checksum changed")
        if any(x not in "ACGTNRYKMSWBDHV" for x in orig):
            raise ValueError("Unsupported source reference nucleotide")
        rid=f"SRC{i:03d}"
        index[rid]={"kind":"source_transcript",
                    "source_transcript_id":record["id"],
                    "source_family_annotation":record["family_by_publisher_annotation_only"],
                    "source_sequence_nt":len(orig)}
        lines.append(f">{rid}\n{orig}\n")
        # The decoy has the same exact mono-nucleotide composition and length;
        # it is a sequence-mapping null, not a calibrated biological FDR.
        scrambled=list(orig)
        rng.shuffle(scrambled)
        did=f"DEC{i:03d}"
        index[did]={"kind":"composition_shuffled_decoy",
                    "source_transcript_id":record["id"],
                    "source_family_annotation":"NULL_NOT_A_GENE",
                    "source_sequence_nt":len(scrambled)}
        lines.append(f">{did}\n{''.join(scrambled)}\n")
    return "".join(lines),index


def fetch_first_compressed_bytes(url,limit,timeout=90):
    """Strict transfer read limit regardless of HTTP Range support.

    Only allow source URLs already validated by the frozen manifest;
    redirects after source lookup must still be HTTPS.
    """
    if not 1024<=limit<=MAX_CAP:
        raise ValueError("Only bounded public compressed-prefix retrieval is allowed")
    request=urllib.request.Request(
        url,headers={"Range":f"bytes=0-{limit-1}",
                     "User-Agent":"aza3-bounded-public-rna-pilot/1.0"})
    with urllib.request.urlopen(request,timeout=timeout) as response:
        if response.geturl().split(":",1)[0]!="https":
            raise ValueError("Download redirected to unencrypted resource")
        status=response.status
        if status not in (200,206):
            raise ValueError(f"Unexpected HTTP status {status}")
        content_range=response.headers.get("Content-Range","")
        pieces=[]
        n=0
        while n<limit:
            chunk=response.read(min(64*1024,limit-n))
            if not chunk:break
            pieces.append(chunk)
            n+=len(chunk)
        payload=b"".join(pieces)
        if len(payload)<4096 or payload[:2]!=b"\x1f\x8b":
            raise ValueError(f"Public ENA payload is not a valid FASTQ.gz prefix ({len(payload)} bytes)")
        if status==206 and not content_range.startswith("bytes 0-"):
            raise ValueError("Range server returned an unexpected nonzero-offset prefix")
        return payload,{
            "http_status":status,
            "range_content_range":content_range,
            "range_honored":status==206,
            "compressed_prefix_downloaded_bytes":len(payload),
            "compressed_prefix_sha256":sha256(payload).hexdigest(),
            "full_original_gzip_integrity_checked":False
        }


def decode_complete_fastq_prefix(gzip_prefix,decompressed_cap=UNCOMPRESSED_CAP):
    if gzip_prefix[:2]!=b"\x1f\x8b":
        raise ValueError("Source bytes do not begin with gzip magic")
    stream=zlib.decompressobj(zlib.MAX_WBITS|16)
    data=bytearray()
    for i in range(0,len(gzip_prefix),65536):
        b=gzip_prefix[i:i+65536]
        try:
            out=stream.decompress(b)
        except zlib.error as exc:
            raise ValueError("Public FASTQ gzip prefix is corrupt") from exc
        if len(data)+len(out)>decompressed_cap:
            raise ValueError("Inflated FASTQ exceeds strict memory budget")
        data.extend(out)
    try:
        text=data.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ValueError("Decoded FASTQ is not ASCII") from exc
    records=text.split("\n")
    nread=0
    ndiscard=0
    valid=[]
    read_len=Counter()
    qual_total=qual_sites=ambig=0
    # Keep only fully terminated 4-line records. Final truncated record
    # is expected with a compressed-byte-prefix, never interpreted as error.
    complete_line_count=len(records)-1
    ncomplete_lines=(complete_line_count//4)*4
    for i in range(0,ncomplete_lines,4):
        name,seq,plus,qual=[records[i+j].rstrip("\r") for j in range(4)]
        if not name.startswith("@") or not plus.startswith("+"):
            raise ValueError(f"FASTQ structure invalid near 4-line record {i//4}")
        seq=seq.upper()
        if not seq or len(seq)!=len(qual) or any(x not in "ACGTN" for x in seq):
            raise ValueError(f"FASTQ bases/quality mismatch at record {i//4}")
        if any(ord(x)<33 or ord(x)>126 for x in qual):
            raise ValueError("Invalid FASTQ quality encoding")
        read_len[len(seq)]+=1
        ambig+=seq.count("N")
        qual_total+=sum(ord(ch)-33 for ch in qual)
        qual_sites+=len(qual)
        valid.append(f"{name}\n{seq}\n+\n{qual}\n")
        nread+=1
    ndiscard=int(ncomplete_lines!=complete_line_count or len(records[-1])>0)
    if nread<1000:
        raise ValueError("Too few complete FASTQ reads under bounded prefix; source pilot not interpretable")
    return "".join(valid),{
        "read1_prefix_fully_completed_records":nread,
        "discarded_truncated_terminal_record_possible":bool(ndiscard),
        "decoded_uncompressed_bytes":len(data),
        "gzip_member_reached_true_eof":bool(stream.eof),
        "read_length_histogram":dict(sorted(read_len.items())),
        "mean_quality_phred33_descriptive":round(qual_total/qual_sites,3),
        "ambiguous_base_fraction_descriptive":round(ambig/qual_sites,6),
        "independent_biological_replicates_added":0,
    }


def _cigar_alignment_span(cigar):
    parts=CIGAR.findall(cigar)
    if "".join(a+b for a,b in parts)!=cigar:
        raise ValueError("SAM CIGAR includes invalid operators")
    aligned=sum(int(a) for a,b in parts if b in ("M","I","D","=","X"))
    aligned_read=sum(int(a) for a,b in parts if b in ("M","I","=","X"))
    return aligned,aligned_read


def parse_minimap2_sam(samfile,refids):
    """Triage read-to-published-transcript mappings with negative controls.

    One primary alignment per short read; neither MAPQ nor NM establishes
    molecular gene orthology. Decoy false-hit summaries are retained.
    """
    totals=Counter()
    raw_source=Counter()
    strict_source=Counter()
    weak_source=Counter()
    strict_decoy=Counter()
    mapq_hist=Counter()
    for line in Path(samfile).read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("@"):
            continue
        fields=line.split("\t")
        if len(fields)<11:raise ValueError("Invalid SAM row")
        flags=int(fields[1])
        if flags&0x100 or flags&0x800:
            totals["secondary_or_supplementary_alignments"]+=1
            continue
        totals["primary_alignment_rows"]+=1
        if flags&0x4 or fields[2]=="*":
            totals["primary_unmapped"]+=1
            continue
        ref=fields[2]
        if ref not in refids:raise ValueError("SAM target reference was not generated by this experiment")
        nm_field=next((s for s in fields[11:] if s.startswith("NM:i:")),None)
        if nm_field is None:raise ValueError("A mapped read lacks minimap2 NM edit count")
        nm=int(nm_field[5:])
        aligned,aligned_read=_cigar_alignment_span(fields[5])
        if aligned<=0:raise ValueError("Zero-length alignment")
        identity=max(0.,1.-nm/aligned)
        mapq=int(fields[4])
        mapq_hist[mapq]+=1
        name=refids[ref]["source_transcript_id"]
        category=refids[ref]["kind"]
        totals[f"primary_mapped_{category}"]+=1
        if category=="source_transcript":raw_source[name]+=1
        weak=aligned_read>=WEAK_MIN_ALIGNED and identity>=WEAK_MIN_IDENTITY
        if weak and category=="source_transcript":weak_source[name]+=1
        strict=(
            mapq>=STRICT_MIN_MAPQ and aligned_read>=STRICT_MIN_ALIGNED
            and identity>=STRICT_MIN_IDENTITY
        )
        if strict:
            totals[f"strict_mapped_{category}"]+=1
            if category=="source_transcript":strict_source[name]+=1
            else:strict_decoy[name]+=1
    return {
        "n_primary_records":totals["primary_alignment_rows"],
        "n_primary_unmapped":totals["primary_unmapped"],
        "n_primary_mapped_source":totals["primary_mapped_source_transcript"],
        "n_primary_mapped_shuffled_decoy":totals["primary_mapped_composition_shuffled_decoy"],
        "n_secondary_or_supplementary_alignments":totals["secondary_or_supplementary_alignments"],
        "n_strict_source_read_alignments":totals["strict_mapped_source_transcript"],
        "n_strict_shuffled_decoy_alignments":totals["strict_mapped_composition_shuffled_decoy"],
        "unique_source_transcripts_with_strict_match":len(strict_source),
        "per_source_primary_counts":dict(sorted(raw_source.items())),
        "per_source_weak_coverage_identity_counts":dict(sorted(weak_source.items())),
        "per_source_strict_counts":dict(sorted(strict_source.items())),
        "per_decoy_strict_counts":dict(sorted(strict_decoy.items())),
        "primary_mapq_histogram":dict(sorted(mapq_hist.items())),
        "sequence_similarity_not_orthology":True,
        "absence_cannot_establish_gene_loss":True,
    }


def align_short_reads(source_fasta,short_fastq,samfile):
    exe=shutil.which("minimap2")
    if not exe:raise RuntimeError("minimap2 executable required for short-read alignment")
    with Path(samfile).open("w",encoding="utf-8") as handle:
        process=subprocess.run(
            [exe,"-ax","sr","-k","15","-w","5","-t","2",str(source_fasta),str(short_fastq)],
            stdout=handle,stderr=subprocess.PIPE,text=True,timeout=150)
    if process.returncode:
        raise RuntimeError(f"minimap2 mapping failed with return code {process.returncode}: {process.stderr[-350:]}")


def analyze_public_prefixes(limit=DEFAULT_CAP,manifest=MANIFEST,reference=REF,
                            reference_meta=REF_META,only_run=None,strict_transfer=True):
    src=source_manifest(manifest)
    if limit<1024 or limit>MAX_CAP:raise ValueError("Prefix cannot exceed 32 MiB per run")
    reference_text,refs=prepare_reference(reference,reference_meta)
    jobs=src["runs"]
    if only_run is not None:
        jobs=[x for x in jobs if x["run_accession"]==only_run]
        if len(jobs)!=1:raise ValueError("Requested SRR is not among approved six accessions")
    out=[]
    total_bytes=0
    with tempfile.TemporaryDirectory() as work:
        folder=Path(work)
        ref=folder/"26_Cirsium_2020_sources_and_26_shuffles.fasta"
        ref.write_text(reference_text,encoding="ascii")
        for run in jobs:
            rec={"run_accession":run["run_accession"],
                 "source_specimen_code":run["source_specimen_code"],
                 "taxon_as_experiment_titled":run["taxon_as_experiment_titled"],
                 "biosample_accession":run["biosample_accession"],
                 "tissue":"young leaf",
                 "full_original_mate1_provider_md5_not_checked":run["R1_provider_md5"],
                 "R2_mate_not_retrieved":True,
                 "full_run_not_analysed":True}
            try:
                gz,meta=fetch_first_compressed_bytes(run["R1_https_url"],limit)
                total_bytes+=len(gz)
                fastq,qc=decode_complete_fastq_prefix(gz)
                fq=folder/(run["run_accession"]+"_read1_prefix.fastq")
                sam=folder/(run["run_accession"]+".sam")
                fq.write_text(fastq,encoding="ascii")
                align_short_reads(ref,fq,sam)
                map_report=parse_minimap2_sam(sam,refs)
                if map_report["n_primary_records"]!=qc["read1_prefix_fully_completed_records"]:
                    raise ValueError("FASTQ read count and SAM primary-read count differ")
                rec.update({"download_status":"PREFIX_ANALYZED",
                            "retrieval":meta,"read_quality":qc,
                            "mapping":map_report,
                            "interpretation":"SOURCE_SEQUENCE_SIMILARITY_IN_LEAF_PREFIX_ONLY" if map_report["n_strict_source_read_alignments"] else "NO_STRICT_MATCH_IN_UNREPRESENTATIVE_PREFIX"})
            except (ValueError,RuntimeError,OSError,urllib.error.URLError,
                    urllib.error.HTTPError,TimeoutError) as exc:
                rec.update({"download_status":"INCOMPLETE_OR_FAILED",
                            "error_type":type(exc).__name__,
                            "error_message":str(exc)[:300],
                            "no_biological_absence_inference":True})
                if strict_transfer:
                    out.append(rec)
                    raise
            out.append(rec)
    return {
        "version":"aza3_cirsium_Ryukyu_R1_bounded_pilot_v6",
        "data_type":"actual_public_R1_compressed_prefix" if any(x["download_status"]=="PREFIX_ANALYZED" for x in out)
                    else "NO_PUBLIC_READ_PREFIX_SUCCESSFULLY_ANALYZED",
        "source_manifest_url":src["source_url"],
        "source_project":"PRJNA1311153",
        "scope":"six source-identified 2026 Ryukyu Cirsium young-leaf samples, NOT flowers",
        "per_run_compressed_cap_bytes":limit,
        "total_compressed_bytes_transferred":total_bytes,
        "max_allowed_total_bytes":len(jobs)*limit,
        "n_runs_attempted":len(out),
        "n_successfully_analyzed":sum(x["download_status"]=="PREFIX_ANALYZED" for x in out),
        "n_downloads_failed":sum(x["download_status"]!="PREFIX_ANALYZED" for x in out),
        "reference_dataset":"26 source-verified 2020 Cirsium japonicum candidate transcripts plus 26 composition-shuffled decoys",
        "thresholds":{"strict_primary_min_mapq":STRICT_MIN_MAPQ,
                      "strict_min_aligned_read_bases":STRICT_MIN_ALIGNED,
                      "strict_min_identity":STRICT_MIN_IDENTITY,
                      "weak_min_aligned_read_bases":WEAK_MIN_ALIGNED,
                      "weak_min_identity":WEAK_MIN_IDENTITY},
        "original_2026_fastq_whole_run_MD5_verified":False,
        "phenotype_or_orthology_estimated":False,
        "can_compare_between_species_expression":False,
        "floral_reexpression_demonstrated":False,
        "early_gzip_prefix_read_order_unrepresentative":True,
        "samples":out,
        "claim_ceiling":"Putative 2020 Cirsium source transcript nucleotide similarity in bounded single-mate young-leaf FASTQ prefixes, not a verified functional gene, expression abundance, genomic allele, gene loss or floral trait re-evolution."
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--fetch-public",action="store_true",
                   help="Explicit bounded public ENA read retrieval; omit for offline validation only")
    p.add_argument("--manifest",type=Path,default=MANIFEST)
    p.add_argument("--source-fasta",type=Path,default=REF)
    p.add_argument("--source-provenance",type=Path,default=REF_META)
    p.add_argument("--cap-mib",type=int,default=8)
    p.add_argument("--only-run")
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    if args.cap_mib not in (1,8,32):
        raise ValueError("Allowed exact caps: 1, 8 or 32 MiB")
    if args.fetch_public:
        report=analyze_public_prefixes(args.cap_mib*1024*1024,
            args.manifest,args.source_fasta,args.source_provenance,
            args.only_run,strict_transfer=False)
    else:
        data=source_manifest(args.manifest)
        ref,ind=prepare_reference(args.source_fasta,args.source_provenance)
        report={
            "version":"aza3_cirsium_Ryukyu_R1_bounded_pilot_v6",
            "retrieval_status":"NOT_REQUESTED",
            "n_focal_source_runs":len(data["runs"]),
            "n_2020_reference_transcripts":26,
            "n_sequence_shuffled_decoys":26,
            "planned_cap_mib":args.cap_mib,
            "no_read_based_evidence_claimed":True,
            "source_reference_fasta_sha256":sha256(ref.encode("ascii")).hexdigest()
        }
    data=json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(data,encoding="utf-8")
    else:print(data)
    print("RYUKYU_PREFIX_PILOT_SUMMARY_JSON:"+json.dumps({
        k:report.get(k) for k in (
            "version","data_type","n_runs_attempted","n_successfully_analyzed",
            "total_compressed_bytes_transferred","floral_reexpression_demonstrated")
    },separators=(",",":")))
    for rec in report.get("samples",[]):
        mapping=rec.get("mapping",{})
        print("RYUKYU_PREFIX_SAMPLE_JSON:"+json.dumps({
            "run_accession":rec["run_accession"],
            "specimen":rec["source_specimen_code"],
            "taxon":rec["taxon_as_experiment_titled"],
            "status":rec["download_status"],
            "error":rec.get("error_message"),
            "n_reads":rec.get("read_quality",{}).get("read1_prefix_fully_completed_records"),
            "mapped_to_source":mapping.get("n_primary_mapped_source"),
            "strict_source":mapping.get("n_strict_source_read_alignments"),
            "strict_decoy":mapping.get("n_strict_shuffled_decoy_alignments"),
            "per_source_strict":mapping.get("per_source_strict_counts"),
            "retrieved_bytes":rec.get("retrieval",{}).get("compressed_prefix_downloaded_bytes")
        },separators=(",",":")))
    return report


if __name__=="__main__":
    main()
