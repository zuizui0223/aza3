"""Resource-capped, fail-closed *real* Ryukyu Cirsium leaf-RNA prefix screen.

Data: official ENA PRJNA1311153 per-run R1 gzip URLs, 3 C. brevicaule
and 3 C. irumtiense young-leaf samples. We fetch **only the first
specified byte range of R1**, never R2, full RNA-seq, or a randomized
representative subset. Prefix selection may be lane/library-order biased.

Candidates: source-verified Park 2020 *C. japonicum* assembled
DFR/ANS-like transcript sequences, 26 in 23 assembly groups.
They are **not** confirmed orthologs or anthocyanin functional genes.

Mapping: minimap2 short-read mapping to all 26 original source transcripts
PLUS same-base-composition shuffled decoys; samtools flags, MAPQ,
aligned length, edit rate. A high-confidence local read match means
merely SOURCE_CIRSIUM_LIKE_TRANSCRIPT_FRAGMENT_IN_LEAF_R1_PREFIX.

No hit => NOT_OBSERVED_IN_CAPPED_PREFIX, NEVER gene absence / pathway
loss. A hit => candidate sequence homolog, NEVER floral regulatory
competence, gene function, adaptation or evolutionary re-expression.
No genotype variants, gene copy numbers or differential expression are
inferred. Paired read evidence is deliberately unavailable.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import urlparse
import argparse
import csv
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
import zlib

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"data/anto_source_handoff"
SOURCE_MANIFEST=DATA/"public_sra/ryukyu_six_R1_prefix_source_v6.json"
RUNS=DATA/"public_sra/prjna1311153_run_manifest_20261017.csv"
SOURCE_FASTA=DATA/"molecular/cirsium_japonicum_2020_26_source_sequences_v6.fasta"
SOURCE_PROVENANCE=DATA/"molecular/cirsium_japonicum_2020_26_sequence_provenance_v6.json"
DEFAULT_BYTES=8*1024*1024
MAX_BYTES=32*1024*1024
MIN_READ_LEN=75
MIN_ALIGNED_NT=90
MIN_MAPQ=20
MAX_EDIT_FRAC=.2
MAX_READS_PER_RUN=450000
MAX_UNCOMPRESSED_PER_RUN=250*1024*1024
ENA_ALLOWED_HOST="ftp.sra.ebi.ac.uk"


def read_fasta(path):
    data={}
    header=None
    chunks=[]
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line=line.strip()
        if line.startswith(">"):
            if header is not None:
                if header in data:raise ValueError("Duplicate FASTA entry")
                data[header]="".join(chunks).upper()
            header=line[1:].split()[0]
            chunks=[]
        elif line and header is not None:
            chunks.append(line)
    if header is not None:
        if header in data:raise ValueError("Duplicate FASTA entry")
        data[header]="".join(chunks).upper()
    return data


def validate_manifest(path=SOURCE_MANIFEST, csv_path=RUNS,
                      fasta_path=SOURCE_FASTA, provenance_path=SOURCE_PROVENANCE):
    registry=json.loads(Path(path).read_text(encoding="utf-8"))
    if registry.get("version")!="aza3_cirsium_Ryukyu_public_R1_prefix_manifest_v6":
        raise ValueError("Wrong official ENA manifest revision")
    runs=registry.get("runs",[])
    if len(runs)!=6 or len({r["run_accession"] for r in runs})!=6:
        raise ValueError("Require six independently identified Ryukyu public R1 records")
    with Path(csv_path).open(newline="",encoding="utf-8") as f:
        source={r["run_accession"]:r for r in csv.DictReader(f)}
    counts=Counter()
    for r in runs:
        run=r["run_accession"]
        record=source.get(run)
        if (record is None or r["biosample_accession"]!=record["biosample_accession"]
            or r["source_specimen_code"]!=record["specimen_code_from_experiment_title"]
            or r["taxon_as_experiment_titled"]!=record["taxon_as_experiment_titled"]
            or r["tissue"]!="young leaf"
            or r["per_individual_corolla_pigmentation_quantified"] is not False):
            raise ValueError(f"ENA source-to-Ryukyu specimen mismatch: {run}")
        url=urlparse(r["R1_https_url"])
        if (url.scheme!="https" or url.hostname!=ENA_ALLOWED_HOST
                or not re.fullmatch(
                    rf"/vol1/fastq/SRR[0-9]+/[0-9]{{3}}/{run}/{run}_1\.fastq\.gz",
                    url.path)):
            raise ValueError("Unverified or nonofficial R1 source URL")
        if (r["R1_full_file_bytes"]<10**6
            or not re.fullmatch(r"[0-9a-f]{32}",r["R1_provider_md5"])
            or r["whole_run_complete"] is not False or
            r["mate2_not_retrieved"] is not True):
            raise ValueError("Public raw source size/provenance not verified")
        counts[r["taxon_as_experiment_titled"]]+=1
    if counts!={"Cirsium brevicaule":3,"Cirsium irumtiense":3}:
        raise ValueError("Unbalanced/changed Ryukyu biological source panel")

    source_info=json.loads(Path(provenance_path).read_text(encoding="utf-8"))
    dna=read_fasta(fasta_path)
    original=source_info.get("records",[])
    if (source_info.get("contract")!="anto_park2020_26_source_transcripts_v6"
        or len(dna)!=26 or len(original)!=26
        or {r["id"] for r in original}!=set(dna)):
        raise ValueError("2020 source-derived read targets no longer reproducible")
    for r in original:
        sequence=dna[r["id"]]
        if (len(sequence)!=r["source_fasta_length"]
            or sha256(sequence.encode("ascii")).hexdigest()!=r["sha256_nucleotide"]):
            raise ValueError("2020 source transcript checksum mismatch")
    return runs,dna,original


def download_official_r1_prefix(url,declared_full_file_bytes,requested_bytes,
                                output_path,timeout=50):
    """Download exactly Range bytes; never silently accept full 1-2GB file.

    The server must respond 206 Partial Content with accurate offset and
    declared original total size. The partial gzip cannot be checked with
    the full-run provider MD5.
    """
    if (type(requested_bytes) is not int or requested_bytes<1024
        or requested_bytes>MAX_BYTES or requested_bytes>=declared_full_file_bytes):
        raise ValueError("R1 requested byte budget invalid")
    parsed=urlparse(url)
    if parsed.scheme!="https" or parsed.hostname!=ENA_ALLOWED_HOST:
        raise ValueError("Only official HTTPS ENA data host permitted")
    request=urllib.request.Request(
        url,headers={"User-Agent":"aza3-capped-public-research-pilot-v6/1.0",
                     "Range":f"bytes=0-{requested_bytes-1}",
                     "Accept-Encoding":"identity"})
    with urllib.request.urlopen(request,timeout=timeout) as response:
        status=response.status
        cr=response.headers.get("Content-Range","")
        expected=f"bytes 0-{requested_bytes-1}/{declared_full_file_bytes}"
        if status!=206 or cr!=expected:
            raise ValueError(
                f"ENA did not honor a verified bounded byte range: HTTP {status}, "
                f"Content-Range {cr!r}; no full file will be downloaded"
            )
        content=response.read(requested_bytes+1)
    if len(content)!=requested_bytes or content[:2]!=b"\x1f\x8b":
        raise ValueError("Truncated/corrupt R1 prefix gzip or server exceeded byte budget")
    Path(output_path).write_bytes(content)
    return {"HTTP_status":status,"Content_Range":cr,
            "compressed_bytes_retrieved":len(content),
            "prefix_sha256":sha256(content).hexdigest(),
            "complete_full_run_file_md5_verified":False}


def partial_gzip_to_fastq(gzip_path,fastq_path,
                          max_uncompressed=MAX_UNCOMPRESSED_PER_RUN,
                          max_reads=MAX_READS_PER_RUN):
    """Extract only FULL four-line FASTQ reads from incomplete gzip prefix.

    The gzip EOF may legitimately be missing. Truncated trailing read is
    dropped and recorded, not silently treated as a valid 150nt read.
    """
    if max_uncompressed<1000 or max_reads<1:
        raise ValueError("Positive bounded read/decompression limits required")
    inflated=0
    pending=b""
    queue=[]
    valid=0
    invalid=0
    bases=0
    members=0
    compressor=zlib.decompressobj(31)
    finished_prefix=False
    with Path(gzip_path).open("rb") as source,Path(fastq_path).open("wb") as dest:
        def digest_bytes(data):
            nonlocal pending,valid,invalid,bases
            # Avoid repeatedly copying the remaining megabyte buffer for
            # each FASTQ line. Split once per bounded compressed chunk.
            lines=(pending+data).split(b"\n")
            pending=lines.pop()  # last line is always incomplete until LF
            for line in lines:
                queue.append(line.rstrip(b"\r"))
                if len(queue)<4:
                    continue
                name,seq,plus,quality=queue
                queue.clear()
                if (not name.startswith(b"@") or not plus.startswith(b"+")
                    or len(seq)!=len(quality) or len(seq)<MIN_READ_LEN
                    or re.fullmatch(rb"[ACGTNacgtn]+",seq) is None):
                    invalid+=1
                    if invalid>4:
                        raise ValueError("Too many malformed complete reads in ENA prefix")
                    continue
                if valid>=max_reads:
                    return True
                dest.write(b"\n".join((name,seq,plus,quality))+b"\n")
                valid+=1
                bases+=len(seq)
            return False
        for chunk in iter(lambda:source.read(256*1024),b""):
            remaining=chunk
            while remaining:
                decompressed=compressor.decompress(remaining)
                inflated+=len(decompressed)
                if inflated>max_uncompressed:
                    raise ValueError("R1 uncompressed prefix exceeds processing budget")
                if digest_bytes(decompressed):
                    finished_prefix=True
                    break
                if compressor.eof:
                    members+=1
                    remaining=compressor.unused_data
                    compressor=zlib.decompressobj(31)
                else:
                    remaining=b""
            if finished_prefix:break
    if valid<10:
        raise ValueError("Fewer than 10 complete reads in source prefix")
    return {
        "complete_valid_reads":valid,
        "total_nucleotide_bases":bases,
        "invalid_complete_record_count":invalid,
        "trailing_incomplete_lines_dropped":len(queue),
        "trailing_incomplete_line_bytes_dropped":len(pending),
        "full_fastq_file_validated":False,
        "gzip_full_stream_eof_verified":False,
        "last_completed_gzip_members":members,
        "read_count_capped":finished_prefix,
        "uncompressed_bytes_seen":inflated,
        "source_read_pairs_exhaustively_sampled":False,
        "representative_random_sample":False,
    }


def prepare_target_panel(dna,source_rows,out_fasta):
    """Source transcript reference panel plus mono-nucleotide shuffled decoys."""
    rng=random.Random(20261019)
    entries=[]
    byid={r["id"]:r for r in source_rows}
    for name in sorted(dna):
        seq=dna[name]
        if len(seq)<MIN_READ_LEN or not set(seq).issubset(set("ACGTUNRYKMSWBDHV-")):
            raise ValueError("Unusable source published nucleotide record")
        seq=seq.replace("U","T").replace("-","N")
        entries.append((name,seq))
        chars=list(seq)
        rng.shuffle(chars)
        entries.append(("DECOY_"+name,"".join(chars)))
    with Path(out_fasta).open("w",encoding="utf-8") as f:
        for name,seq in entries:
            f.write(">"+name+"\n")
            for i in range(0,len(seq),75):f.write(seq[i:i+75]+"\n")
    return {
        "n_real_publisher_transcripts":len(dna),
        "n_sequence_shuffle_decoys":len(dna),
        "reference_map":{"transcript_count":len(byid),"source_study":"Park 2020 Cirsium japonicum"},
        "decoy_method":"fixed-seed independent shuffled nucleotide composition, not a phylogenetic negative control"
    }


def cigar_aligned_span(cigar):
    """Read-consuming aligned span M/= /X only; excludes soft clips."""
    if not cigar or cigar=="*":
        return 0
    parts=re.findall(r"(\d+)([MIDNSHP=X])",cigar)
    if "".join(n+c for n,c in parts)!=cigar:raise ValueError("Malformed SAM CIGAR")
    return sum(int(n) for n,c in parts if c in "M=X")


def parse_mapped_sam(lines,run_id):
    """Read only SAM mapped primary hits. Inference at transcript level ONLY."""
    per_transcript=Counter()
    raw_primary=Counter()
    low_mapq=Counter()
    decoys=Counter()
    observed=set()
    high_conf_keys=set()
    aligned_lengths=[]
    identities=[]
    for line in lines:
        if not line.strip() or line.startswith("@"):continue
        fields=line.rstrip("\n").split("\t")
        if len(fields)<11:raise ValueError("Malformed samtools SAM output")
        flag=int(fields[1])
        if (flag&4) or (flag&256) or (flag&2048):
            continue
        name,ref,mpq,cigar=fields[0],fields[2],int(fields[4]),fields[5]
        raw_primary[ref]+=1
        span=cigar_aligned_span(cigar)
        nm=next((int(f[5:]) for f in fields[11:] if f.startswith("NM:i:")),None)
        if nm is None or nm<0:raise ValueError("Mapped read requires NM edit distance")
        edit_fraction=nm/max(1,span)
        if mpq<MIN_MAPQ or span<MIN_ALIGNED_NT or edit_fraction>MAX_EDIT_FRAC:
            low_mapq[ref]+=1
            continue
        if (name,ref) in high_conf_keys:continue
        high_conf_keys.add((name,ref))
        if ref.startswith("DECOY_"):
            decoys[ref]+=1
        else:
            per_transcript[ref]+=1
            aligned_lengths.append(span)
            identities.append(1-edit_fraction)
        observed.add(name)
    return {
        "study_run":run_id,
        "primary_mapped_alignments":sum(raw_primary.values()),
        "passing_high_confidence_unique_read_names":len(observed),
        "passing_real_transcript_reads":sum(per_transcript.values()),
        "passing_decoy_reads":sum(decoys.values()),
        "source_transcript_counts":dict(sorted(per_transcript.items())),
        "decoy_counts":dict(sorted(decoys.items())),
        "low_or_ambiguous_alignment_counts":dict(sorted(low_mapq.items())),
        "raw_primary_mapped_counts":dict(sorted(raw_primary.items())),
        "alignments_below_mapq_or_length_or_identity_are_not_evidence_of_gene_absence":True,
        "high_confidence_span_threshold":MIN_ALIGNED_NT,
        "high_confidence_mapq_threshold":MIN_MAPQ,
        "high_confidence_max_NM_fraction":MAX_EDIT_FRAC,
        "alignment_span_range":(
            [min(aligned_lengths),max(aligned_lengths)] if aligned_lengths else None),
        "mean_approx_alignment_identity":(
            sum(identities)/len(identities) if identities else None)
    }


def map_prefix_fastq(reference_fasta,reads_fastq,run_id,threads=2,timeout=180):
    """Capped external alignment: minimap2 -> samtools view mapped alignments."""
    if shutil.which("minimap2") is None or shutil.which("samtools") is None:
        raise RuntimeError("Require locally installed minimap2 and samtools")
    cmd=["minimap2","-t",str(threads),"-ax","sr","--secondary=no",
         str(reference_fasta),str(reads_fastq)]
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    b=subprocess.Popen(["samtools","view","-F","4","-"],
                       stdin=a.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                       text=True)
    a.stdout.close()
    try:
        out,err=b.communicate(timeout=timeout)
        aerr=a.stderr.read().decode("utf-8",errors="replace")
        a.returncode=a.wait(timeout=10)
    except (subprocess.TimeoutExpired,KeyboardInterrupt):
        a.kill()
        b.kill()
        raise
    if a.returncode!=0 or b.returncode!=0:
        raise RuntimeError(f"Alignment failed: minimap={aerr[-500:]}, samtools={err[-500:]}")
    return parse_mapped_sam(out.splitlines(),run_id)


def run_pilot(manifest=SOURCE_MANIFEST,run_manifest=RUNS,
              source_fasta=SOURCE_FASTA,source_provenance=SOURCE_PROVENANCE,
              requested_bytes=DEFAULT_BYTES,
              selected_runs=None):
    runs,dna,records=validate_manifest(manifest,run_manifest,source_fasta,source_provenance)
    if type(requested_bytes) is not int or not 1024<=requested_bytes<=MAX_BYTES:
        raise ValueError("Only capped 1KiB..32MiB per-run R1 prefixes permitted")
    if selected_runs is not None:
        selection=set(selected_runs)
        if not selection or not selection<={r["run_accession"] for r in runs}:
            raise ValueError("Selected source run absent or empty")
        runs=[r for r in runs if r["run_accession"] in selection]
    result={
        "version":"aza3_ryukyu_R1_prefix_gene_family_screen_v6",
        "source_project":"PRJNA1311153",
        "biological_data_source":"Actual public 2026 Ryukyu young-leaf RNA, restricted to gzip beginning of mate R1",
        "per_run_requested_compressed_prefix_bytes":requested_bytes,
        "project_declared_full_R1_MD5s_are_not_validated_by_prefix":True,
        "genotype_or_allele_state_inferred":False,
        "gene_copy_numbers_inferred":False,
        "RNA_expression_quantified_or_compared_between_species":False,
        "floral_gene_expression_evaluated":False,
        "molecular_function_or_orthology_verified":False,
        "floral_reexpression_demonstrated":False,
        "sampling_bias":"FASTQ file order may correlate with lane/technical effects; only first R1 gzip bytes were read",
        "reference_bias":"Park 2020 C. japonicum assembled BLAST DFR/ANS-like candidates, many not known functional genes",
        "full_run_gene_absence_evidence":False,
        "runs":[]
    }
    with TemporaryDirectory(prefix="aza3_ryukyu_source_") as folder:
        root=Path(folder)
        panel=prepare_target_panel(dna,records,root/"all_source_targets_plus_decoy.fasta")
        result["reference_panel"]=panel
        for sample in runs:
            run=sample["run_accession"]
            prefix=root/(run+".gz")
            fastq=root/(run+".fq")
            download=download_official_r1_prefix(
                sample["R1_https_url"],sample["R1_full_file_bytes"],
                requested_bytes,prefix)
            parsed=partial_gzip_to_fastq(prefix,fastq)
            mapped=map_prefix_fastq(root/"all_source_targets_plus_decoy.fasta",fastq,run)
            result["runs"].append({
                "run":run,
                "taxon":sample["taxon_as_experiment_titled"],
                "specimen_code":sample["source_specimen_code"],
                "biosample":sample["biosample_accession"],
                "original_declared_R1_total_bytes":sample["R1_full_file_bytes"],
                "original_declared_read_pairs":sample["original_read_pairs"],
                "original_provider_R1_md5":sample["R1_provider_md5"],
                "download":download,
                "parsed_R1_reads":parsed,
                "mapping":mapped,
                "interpretation":"SOURCE_TRANSCRIPT_FRAGMENT_MAPPING_ONLY" if mapped["passing_real_transcript_reads"] else "NO_HIGH_CONFIDENCE_SOURCE_MATCH_IN_BIASED_PREFIX_NOT_GENE_ABSENCE"
            })
            prefix.unlink(missing_ok=True)
            fastq.unlink(missing_ok=True)
    result["n_Ryukyu_source_specimens_queried"]=len(result["runs"])
    result["n_taxon_level_R1_source_positive_prefixes"]=sum(x["mapping"]["passing_real_transcript_reads"]>0 for x in result["runs"])
    result["total_mapped_nondecoy_read_fragments"]=sum(x["mapping"]["passing_real_transcript_reads"] for x in result["runs"])
    result["total_mapped_decoy_read_fragments"]=sum(x["mapping"]["passing_decoy_reads"] for x in result["runs"])
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--bytes-per-run",type=int,default=DEFAULT_BYTES)
    ap.add_argument("--runs",nargs="*",default=None)
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    res=run_pilot(requested_bytes=args.bytes_per_run,selected_runs=args.runs)
    text=json.dumps(res,ensure_ascii=False,sort_keys=True,indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding="utf-8")
    else:print(text)
    print("RYUKYU_PREFIX_SUMMARY_JSON:"+json.dumps({
        "n":res["n_Ryukyu_source_specimens_queried"],
        "prefix_bytes":res["per_run_requested_compressed_prefix_bytes"],
        "source_mapping_positive_R1_prefixes":res["n_taxon_level_R1_source_positive_prefixes"],
        "matching_source_read_fragments":res["total_mapped_nondecoy_read_fragments"],
        "decoy_read_fragments":res["total_mapped_decoy_read_fragments"],
        "floral_reexpression":res["floral_reexpression_demonstrated"]
    },separators=(",",":")))
    for r in res["runs"]:
        print("RYUKYU_PREFIX_RUN_JSON:"+json.dumps({
            "run":r["run"],"taxon":r["taxon"],"specimen":r["specimen_code"],
            "prefix_bytes":r["download"]["compressed_bytes_retrieved"],
            "reads":r["parsed_R1_reads"]["complete_valid_reads"],
            "mapping":r["mapping"]["source_transcript_counts"],
            "decoy":r["mapping"]["passing_decoy_reads"],
            "raw_primary":r["mapping"]["primary_mapped_alignments"],
            "interpretation":r["interpretation"],
        },separators=(",",":")))
    return res


if __name__=="__main__":
    main()
