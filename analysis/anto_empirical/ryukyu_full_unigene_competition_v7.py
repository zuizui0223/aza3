"""Genome-free specificity check using ALL 104,890 Park 2020 Cirsium transcripts.

Live public source:
- Park (2020) annotated unigene FASTA from Europe PMC / PMC7675961
- 2026 PRJNA1311153 SIX Ryukyu Cirsium young-leaf R1 FASTQ prefixes.
- Previously chosen 26 Park transcript BLAST DFR/ANS-LIKE candidates.

This tests whether earlier read-to-26-candidate hits remain UNIQUE after
competing against the author's ENTIRE assembled 104,890-transcript reference.
It is NOT a complete Cirsium genome, an enzyme activity validation, a
genotype call, a full-RNA expression comparison, or evolutionary re-expression.

Both public files are downloaded only within strict byte caps:
- 2020 source archive compressed <= 50 MB, audited uncompressed member SHA.
- 2026 R1 <= 8 MiB per run, must return HTTP 206 before reading bytes.
No raw 2026 reads or 2020 104k-source sequences enter Git history.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from zipfile import ZipFile
from io import BytesIO

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"data/anto_source_handoff"
sys.path.insert(0,str(ROOT/"analysis/anto_empirical"))

from ryukyu_ena_prefix_homology_v6 import (
    validate_manifest,download_official_r1_prefix,partial_gzip_to_fastq,
    cigar_aligned_span,
)
from cirsium_source_homology_v3 import read_fasta

SOURCE_URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7675961/supplementaryFiles"
SOURCE_UNIGENE="ao0c04001_si_002.txt"
EXPECTED_SOURCE_SHA="0c7b7d61a8e2ad5e89b9199347a50545be707329715cdbe0a8a351b0eef2815f"
EXPECTED_SOURCE_BYTES=90263821
EXPECTED_SOURCE_COUNT=104890
ARCHIVE_MAX=50*1000*1000
UNCOMPRESSED_ARCHIVE_MAX=180*1000*1000
MAX_PER_RUN=8*1024*1024
MAPQ=20
MIN_ALIGNED=90
MAX_EDIT_FRACTION=.20
EXPECTED_RUNS={"SRR35152724","SRR35152725","SRR35152729",
               "SRR35152730","SRR35152731","SRR35152732"}
TRANSCRIPT_HEADER=re.compile(rb"^>(TR\d+\|c\d+_g\d+_i\d+)\s+len=(\d+)(?:\s|$)")


def _bounded_public_archive(maximum=ARCHIVE_MAX,timeout=50):
    req=urllib.request.Request(SOURCE_URL,headers={
        "User-Agent":"aza3-cirsium-2020-reference-background-check/1.0",
        "Accept":"application/zip,application/octet-stream"})
    with urllib.request.urlopen(req,timeout=timeout) as r:
        declared=r.headers.get("Content-Length")
        if declared and int(declared)>maximum:
            raise ValueError("2020 public source archive exceeds upload budget")
        content=r.read(maximum+1)
    if len(content)>maximum:
        raise ValueError("2020 public archive exceeded compressed source cap")
    if not content.startswith(b"PK"):
        raise ValueError("2020 source archive is not a ZIP")
    return content


def _nucleotides(seq):
    s=seq.upper().replace("U","T")
    if any(c not in "ACGTNRYKMSWBDHV-" for c in s):
        raise ValueError("Unexpected published 2020 unigene base alphabet")
    return s.replace("-","N")


def prepare_competing_reference(archive, source_subset_fasta, output_fasta, expected_count=EXPECTED_SOURCE_COUNT, expected_sha=EXPECTED_SOURCE_SHA, expected_bytes=EXPECTED_SOURCE_BYTES, subset_n=26):
    """Stream-author-verified 104,890 unigenes to temporary reference.

    The 26 candidates must exactly match the author's full-source fasta.
    Preserve source unigene IDs, and add 26 mononucleotide-shuffled decoys.
    No fake copies or approximated homologs are admitted.
    """
    source_reference=read_fasta(source_subset_fasta)
    if len(source_reference)!=subset_n:
        raise ValueError("Unexpected number of audited source candidates")
    with ZipFile(BytesIO(archive),"r") as z:
        if sum(m.file_size for m in z.infolist())>UNCOMPRESSED_ARCHIVE_MAX:
            raise ValueError("Source archive decompression expansion too large")
        item=next((m for m in z.infolist()
                   if Path(m.filename).name==SOURCE_UNIGENE),None)
        if item is None or item.file_size!=expected_bytes:
            raise ValueError("Publisher supplementary 104k FASTA member missing or changed")
        seen=set()
        found={}
        count=0
        digest=sha256()
        current=None
        n=0
        chunks=[]
        def finalize(writer):
            nonlocal current,n,chunks
            if current is None:return
            seq="".join(chunks)
            if len(seq)!=n:raise ValueError(f"Original 2020 author source FASTA len differs for {current}")
            if current in source_reference:
                expected=_nucleotides(source_reference[current])
                if _nucleotides(seq)!=expected:
                    raise ValueError("Source candidate sequence mismatch against published 2020 subset: "+current)
                found[current]=sha256(seq.encode("ascii")).hexdigest()
            writer.write(">"+current+"\n")
            for offset in range(0,len(seq),75):
                writer.write(seq[offset:offset+75]+"\n")
            current=None;n=0;chunks=[]
        with z.open(item) as stream,Path(output_fasta).open("w",encoding="utf-8") as writer:
            for bline in stream:
                digest.update(bline)
                line=bline.strip()
                if line.startswith(b">"):
                    finalize(writer)
                    m=TRANSCRIPT_HEADER.match(line)
                    if m is None:
                        raise ValueError("Unexpected author unigene FASTA header: "+line[:130].decode("ascii",errors="replace"))
                    current=m.group(1).decode("ascii")
                    n=int(m.group(2))
                    if current in seen:raise ValueError("Repeated publisher unigene ID")
                    seen.add(current)
                    count+=1
                elif current is not None and line:
                    chunks.append(_nucleotides(line.decode("ascii")))
            finalize(writer)
            if digest.hexdigest()!=expected_sha or count!=expected_count:
                raise ValueError(f"Source FASTA hash/count changed: {digest.hexdigest()} / {count}")
            if set(found)!=set(source_reference):
                raise ValueError("Not all frozen source candidates were independently found in the original author FASTA")
            # Negative controls are 26 independent sequence-shuffled
            # counterparts; small source composition-matched comparison only.
            randomizer=random.Random(20261019)
            for gene,original in sorted(source_reference.items()):
                bases=list(_nucleotides(original))
                randomizer.shuffle(bases)
                seq="".join(bases)
                label="DECOY_"+gene
                writer.write(">"+label+"\n")
                for offset in range(0,len(seq),75):
                    writer.write(seq[offset:offset+75]+"\n")
    return {
        "full_publisher_unigene_count":count,
        "candidate_sequence_ids_verified_in_full_reference":len(found),
        "full_2020_author_fasta_uncompressed_sha256":digest.hexdigest(),
        "added_composition_shuffled_candidate_decoys":26,
        "background_reference_is_complete_genome":False
    }


def build_index(reference_fasta,output_index,timeout=360):
    if shutil.which("minimap2") is None:
        raise RuntimeError("minimap2 binary required")
    cmd=["minimap2","-x","sr","-k","15","-w","5","-I","1G",
         "-d",str(output_index),str(reference_fasta)]
    run=subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,
                       timeout=timeout,text=True)
    if run.returncode or not Path(output_index).is_file():
        raise RuntimeError("Whole-author-transcriptome mapping index construction failed: "+run.stderr[-250:])
    return {"index_created":True,"index_disk_bytes":Path(output_index).stat().st_size}


def parse_competitor_sam(lines,selected_ids):
    """Count primary source matches under full 104k reference competition."""
    counters=Counter()
    target=Counter()
    background=Counter()
    decoy=Counter()
    for line in lines:
        if not line or line.startswith("@"):continue
        fields=line.split("\t")
        if len(fields)<11:raise ValueError("Incomplete minimap2 alignment SAM record")
        flags=int(fields[1])
        if flags&4 or flags&256 or flags&2048:continue
        name,ref=fields[0],fields[2]
        counters["mapped_primary_to_any_104k_or_decoy"]+=1
        if ref.startswith("DECOY_"):
            cls="composition_shuffled_decoy"
        elif ref in selected_ids:
            cls="original_26_candidate"
        else:
            cls="other_104864_public_transcript"
        counters["raw_"+cls]+=1
        mapq=int(fields[4])
        aligned=cigar_aligned_span(fields[5])
        nm=next((int(s[5:]) for s in fields[11:] if s.startswith("NM:i:")),None)
        if nm is None:
            raise ValueError("Mapped SAM row has no edit-distance tag")
        ok=(mapq>=MAPQ and aligned>=MIN_ALIGNED
            and nm/max(aligned,1)<=MAX_EDIT_FRACTION)
        if not ok:
            counters["low_mapq_or_coverage_or_identity_"+cls]+=1
            continue
        counters["strict_"+cls]+=1
        if cls=="original_26_candidate":target[ref]+=1
        elif cls=="composition_shuffled_decoy":decoy[ref]+=1
        else:background[ref]+=1
    return {
        "primary_counts":dict(counters),
        "source_candidate_strict_26_counts":dict(sorted(target.items())),
        "background_104864_strict_target_groups":len(background),
        "composition_shuffled_decoy_strict_counts":dict(sorted(decoy.items())),
        "scoring_rules":{
            "min_mapq":MAPQ,
            "min_aligned_nt":MIN_ALIGNED,
            "max_edit_fraction":MAX_EDIT_FRACTION
        },
        "observed_not_gene_function":True
    }


def align_to_full_reference(index,fastq,timeout=350):
    if not shutil.which("samtools"):
        raise RuntimeError("samtools is required for output streaming")
    a=subprocess.Popen(["minimap2","-ax","sr","--secondary=no",
                        "-t","2",str(index),str(fastq)],
                        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    b=subprocess.Popen(["samtools","view","-F","4","-"],
                        stdin=a.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                        text=True)
    a.stdout.close()
    selected=read_fasta(DATA/"molecular/cirsium_japonicum_2020_26_source_sequences_v6.fasta")
    result=parse_competitor_sam(b.stdout,selected)
    b.stdout.close()
    be=b.stderr.read().decode("utf-8",errors="replace")[-500:]
    ae=a.stderr.read().decode("utf-8",errors="replace")[-500:]
    b_code=b.wait(timeout=timeout)
    a_code=a.wait(timeout=timeout)
    if a_code or b_code:
        raise RuntimeError(f"Reference competition mapping failed minimap2={a_code}, samtools={b_code}, {ae}, {be}")
    return result


def run(runs=None,cap_bytes=MAX_PER_RUN):
    if type(cap_bytes) is not int or cap_bytes<256*1024 or cap_bytes>MAX_PER_RUN:
        raise ValueError("Total R1 source budget limited to <=8 MiB per run")
    selected=EXPECTED_RUNS if runs is None else set(runs)
    if not selected or not selected<=EXPECTED_RUNS:
        raise ValueError("Select only source-audited Ryukyu run identifiers")
    registry,_,_=validate_manifest()
    registry={x["run_accession"]:x for x in registry}
    source=DATA/"molecular/cirsium_japonicum_2020_26_source_sequences_v6.fasta"
    result={
        "analysis_version":"aza3_Ryukyu_whole_2020_unigene_competition_v7",
        "biological_source":"2026 PRJNA1311153 Ryukyu young-leaf R1 prefix; Park 2020 complete author unigene reference",
        "preselected_2020_reference_subset_size":26,
        "source_2020_individual_coloured_flower_phenotypes_assayed":False,
        "Ryukyu_2026_individual_flower_RNA_assayed":False,
        "complete_run_or_R2_downloaded":False,
        "R1_prefix_per_sample_requested_bytes":cap_bytes,
        "specimens":[],
        "any_complete_DFR_ANS_orthology_established":False,
        "any_ancestral_white_to_colour_reexpression_established":False,
        "differential_expression_estimate_permitted":False,
    }
    with tempfile.TemporaryDirectory(prefix="aza3-public-competition-") as work:
        folder=Path(work)
        archive=_bounded_public_archive()
        source_fasta=folder/"2020_all_104890_unigenes_plus_decoys.fa"
        result["source_validation"]=prepare_competing_reference(
            archive,source,source_fasta)
        index=folder/"2020_all_104890_unigenes.mmi"
        result["index_build"]=build_index(source_fasta,index)
        for run_id in sorted(selected):
            r=registry[run_id]
            gz=folder/(run_id+".R1.prefix.fastq.gz")
            fq=folder/(run_id+".R1.complete.fastq")
            download=download_official_r1_prefix(
                r["R1_https_url"],r["R1_full_file_bytes"],
                cap_bytes,gz)
            qc=partial_gzip_to_fastq(gz,fq)
            mapping=align_to_full_reference(index,fq)
            row={
                "run_accession":run_id,
                "biosample_accession":r["biosample_accession"],
                "source_specimen_code":r["source_specimen_code"],
                "taxon_from_experiment_title":r["taxon_as_experiment_titled"],
                "tissue":"young leaf",
                "compressed_prefix_provenance":download,
                "complete_valid_R1_reads_from_prefix":qc["complete_valid_reads"],
                "source_competition_mappings":mapping,
                "no_full_file_MD5_verified":True,
                "plant_function_or_flower_regain_assayed":False,
            }
            result["specimens"].append(row)
            gz.unlink(missing_ok=True)
            fq.unlink(missing_ok=True)
    result["n_successful_bounded_specimens"]=len(result["specimens"])
    result["all_target_26_strict_reads_under_full_background"]=sum(
        sum(x["source_competition_mappings"]["source_candidate_strict_26_counts"].values())
        for x in result["specimens"])
    result["all_candidate_TR39766_R1_reads_under_full_background"]=sum(
        x["source_competition_mappings"]["source_candidate_strict_26_counts"].get("TR39766|c0_g1_i2",0)
        for x in result["specimens"])
    result["all_candidate_TR41106_R1_reads_under_full_background"]=sum(
        x["source_competition_mappings"]["source_candidate_strict_26_counts"].get("TR41106|c0_g1_i5",0)
        for x in result["specimens"])
    result["interpretation"]="Testing genomic family specificity under competing 104890 author transcript models, not functional gene presence or pigment competence. No hit is NOT_DETECTED_IN_PREFIX, never absent locus."
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--runs",nargs="*",default=None)
    p.add_argument("--limit-bytes",type=int,default=MAX_PER_RUN)
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    result=run(args.runs,args.limit_bytes)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    else:
        print(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True))
    print("WHOLE_SOURCE_COMPETITION_SUMMARY_JSON:"+json.dumps({
        "n_source_samples":result["n_successful_bounded_specimens"],
        "real_2020_source_unigenes":result["source_validation"]["full_publisher_unigene_count"],
        "strict_preselected_26_fragments":result["all_target_26_strict_reads_under_full_background"],
        "TR39766":result["all_candidate_TR39766_R1_reads_under_full_background"],
        "TR41106":result["all_candidate_TR41106_R1_reads_under_full_background"],
        "flower_regain":result["any_ancestral_white_to_colour_reexpression_established"]
    },separators=(",",":")))
    for item in result["specimens"]:
        a=item["source_competition_mappings"]
        print("WHOLE_SOURCE_COMPETITION_RUN_JSON:"+json.dumps({
            "run":item["run_accession"],
            "taxon":item["taxon_from_experiment_title"],
            "reads_scanned":item["complete_valid_R1_reads_from_prefix"],
            "strict26":sum(a["source_candidate_strict_26_counts"].values()),
            "TR39766":a["source_candidate_strict_26_counts"].get("TR39766|c0_g1_i2",0),
            "TR41106":a["source_candidate_strict_26_counts"].get("TR41106|c0_g1_i5",0),
            "strict_decoy":sum(a["composition_shuffled_decoy_strict_counts"].values())
        },separators=(",",":")))
    return result

if __name__=="__main__":
    main()
