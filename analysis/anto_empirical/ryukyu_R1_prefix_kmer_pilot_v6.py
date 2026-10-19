"""Bounded real-data pilot for six 2026 Ryukyu *Cirsium* young-leaf R1 runs.

SOURCE: ENA PRJNA1311153 first compressed R1 bytes; 2020 *C. japonicum*
publisher-authenticated, 26 DFR/ANS-LIKE (BLAST-only) transcript sequences.

A positive is a **cross-species exact 31-mer homology candidate** only.
A negative means NOT_DETECTED_IN_THIS_ORDERED_R1_PREFIX. NO inference
about functional orthology, anthocyanin loss, flower-stage expression,
colour recurrence, expression levels or evolutionary history.

All FASTQ bytes remain memory-only and are not committed. Network fetch
*requires* HTTP 206 with a zero-based Content-Range. 200/full-file
responses fail closed before downloading multi-GB FASTQs.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import argparse
import gzip
import json
import math
import random
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
import zlib

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"data/anto_source_handoff"
MANIFEST=DATA/"public_sra/ryukyu_six_R1_prefix_source_v6.json"
TARGETS=DATA/"molecular/cirsium_japonicum_2020_26_source_sequences_v6.fasta"
PROVENANCE=DATA/"molecular/cirsium_japonicum_2020_26_sequence_provenance_v6.json"

CAP_MAX=8*1024*1024
EXPANDED_MAX=96*1024*1024
K=31
Q_MIN=20
LOW_COMPLEXITY_ENTROPY=1.5
EXPECTED={"SRR35152724","SRR35152725","SRR35152729",
          "SRR35152730","SRR35152731","SRR35152732"}
A2B=str.maketrans("ACGT","TGCA")
DNA=set("ACGT")


def reverse_complement(seq):
    return seq.translate(A2B)[::-1]


def canonical(seed):
    other=reverse_complement(seed)
    return seed if seed<=other else other


def sequence_entropy(seed):
    counts=Counter(seed)
    n=len(seed)
    return -sum((v/n)*math.log2(v/n) for v in counts.values())


def allowed_seed(seed):
    if len(seed)!=K or any(c not in DNA for c in seed):
        return False
    if len(set(seed))<3 or sequence_entropy(seed)<LOW_COMPLEXITY_ENTROPY:
        return False
    if any(c*8 in seed for c in "ACGT"):
        return False
    return True


def dinucleotide_preserving_shuffle(seq, rng):
    """Eulerian shuffle, conserving 16 adjacent dinucleotide multiplicities.

    Not a phylogenetic null. Used only as a base-composition and local
    dinucleotide control against nonrandom low-complexity exact matches.
    """
    seq=seq.upper().replace("U","T")
    if not seq or not set(seq)<=DNA:
        raise ValueError("Expected unambiguous ACGT sequence for this shuffle")
    outgoing={c:[] for c in "ACGT"}
    for a,b in zip(seq,seq[1:]):
        outgoing[a].append(b)
    for edges in outgoing.values():
        rng.shuffle(edges)
    path,stack=[],[seq[0]]
    while stack:
        c=stack[-1]
        if outgoing[c]:
            stack.append(outgoing[c].pop())
        else:
            path.append(stack.pop())
    result="".join(reversed(path))
    if len(result)!=len(seq):
        raise ArithmeticError("Dinucleotide shuffled transcript length drifted")
    return result


def parse_source_fasta(path=TARGETS):
    records={}
    current=None
    parts=[]
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line=line.strip()
        if line.startswith(">"):
            if current:
                records[current]="".join(parts).upper()
            current=line[1:].split()[0]
            if not current or current in records:
                raise ValueError("Source transcript header missing or duplicate")
            parts=[]
        elif line:
            if current is None:
                raise ValueError("FASTA sequence before identifier")
            parts.append(line)
    if current:
        records[current]="".join(parts).upper()
    return records


def checked_targets(fasta=TARGETS, provenance=PROVENANCE):
    registry=json.loads(Path(provenance).read_text(encoding="utf-8"))
    if registry.get("contract")!="anto_park2020_26_source_transcripts_v6":
        raise ValueError("Unexpected publisher source provenance")
    records=parse_source_fasta(fasta)
    documented=registry.get("records",[])
    if len(records)!=26 or len(documented)!=26 or set(records)!={x["id"] for x in documented}:
        raise ValueError("Missing one or more of the 26 source transcript records")
    for r in documented:
        seq=records[r["id"]]
        if (len(seq)!=r["source_fasta_length"] or
                sha256(seq.encode("ascii")).hexdigest()!=r["sha256_nucleotide"]):
            raise ValueError("Source sequence does not match publisher-derived checksum")
    return records


def seed_index(transcripts, shuffled=False):
    """Canonical 31-mer lookup from 23 source *assembly groups* not 26 genes."""
    result=defaultdict(set)
    rng=random.Random(20261019)
    eligible=0
    for ident,raw in sorted(transcripts.items()):
        seq=raw.upper().replace("U","T")
        grp=ident.rsplit("_i",1)[0]
        # Keep all source IUPAC-containing transcripts, but do not add an
        # ambiguous 31-mer to the exact-match index or its decoy.
        if shuffled:
            seq="".join(
                dinucleotide_preserving_shuffle(part,rng)
                if part and set(part)<=DNA else part
                for part in re.split(r"([^ACGT]+)",seq) if part
            )
        for i in range(len(seq)-K+1):
            seed=seq[i:i+K]
            if allowed_seed(seed):
                result[canonical(seed)].add(grp)
                eligible+=1
    return {s:tuple(sorted(groups)) for s,groups in result.items()},eligible


def indexes(fasta=TARGETS,provenance=PROVENANCE):
    transcripts=checked_targets(fasta,provenance)
    real,n_real=seed_index(transcripts)
    control,n_control=seed_index(transcripts,shuffled=True)
    # Shared background motifs invalidate the 'negative' control; remove
    # identical reference k-mers rather than calling them independent.
    overlap=set(real)&set(control)
    for kmer in overlap:
        control.pop(kmer)
    return {
        "real":real,"decoy":control,
        "eligible_original_windows":n_real,
        "eligible_dinuc_shuffled_windows":n_control,
        "decoy_windows_removed_for_source_overlap":len(overlap),
        "source_transcripts":len(transcripts),
        "source_assembly_groups":len({x.rsplit("_i",1)[0] for x in transcripts})
    }


def validate_manifest(path=MANIFEST, max_per_run=CAP_MAX):
    v=json.loads(Path(path).read_text(encoding="utf-8"))
    if v.get("version")!="aza3_cirsium_Ryukyu_public_R1_prefix_manifest_v6":
        raise ValueError("Wrong frozen public-source manifest")
    runs=v.get("runs",[])
    if len(runs)!=6 or {r.get("run_accession") for r in runs}!=EXPECTED:
        raise ValueError("All six Ryukyu public runs must be tracked")
    if not isinstance(max_per_run,int) or not 0<max_per_run<=CAP_MAX:
        raise ValueError("Maximum per-run compressed read budget is 8 MiB")
    for x in runs:
        url=x.get("R1_https_url","")
        parsed=urllib.parse.urlparse(url)
        accession=x["run_accession"]
        if (parsed.scheme!="https" or parsed.hostname!="ftp.sra.ebi.ac.uk" or
                not parsed.path.endswith(f"/{accession}_1.fastq.gz") or parsed.query):
            raise ValueError("FASTQ URL not a valid allowlisted ENA R1 file")
        if (x["tissue"]!="young leaf" or x["mate2_not_retrieved"] is not True
                or x["whole_run_complete"] is not False
                or x["per_individual_corolla_pigmentation_quantified"] is not False):
            raise ValueError("Public source must not be recoded as a flower dataset")
        if not isinstance(x["R1_full_file_bytes"],int) or x["R1_full_file_bytes"]<=max_per_run:
            raise ValueError("Unexpected unbounded or fully contained source file")
        if not re.fullmatch("[0-9a-fA-F]{32}",x["R1_provider_md5"]):
            raise ValueError("No declared source full-file MD5")
    return sorted(runs,key=lambda r:r["run_accession"])


def fetch_range(r,cap=CAP_MAX,timeout=35):
    """Never read a byte from an HTTP 200 multi-GB full-file response."""
    req=urllib.request.Request(r["R1_https_url"],
        headers={
            "Range":f"bytes=0-{cap-1}",
            "Accept-Encoding":"identity",
            "User-Agent":"aza3-cirsium-source-prefix-validation/1.0"
        })
    with urllib.request.urlopen(req,timeout=timeout) as res:
        status=getattr(res,"status",res.getcode())
        content_range=res.headers.get("Content-Range","")
        if status!=206:
            raise ValueError(f"REMOTE_RANGE_UNAVAILABLE_HTTP_{status}")
        match=re.fullmatch(r"bytes\s+0-(\d+)/(\d+)",content_range)
        if not match or int(match.group(1))!=cap-1 or int(match.group(2))!=r["R1_full_file_bytes"]:
            raise ValueError("REMOTE_RANGE_HEADER_INVALID:"+content_range[:110])
        data=bytearray()
        while len(data)<cap:
            part=res.read(min(256*1024,cap-len(data)))
            if not part:break
            data.extend(part)
        if len(data)!=cap:
            raise ValueError("REMOTE_RANGE_TRUNCATED")
    if not data.startswith(b"\x1f\x8b"):
        raise ValueError("REMOTE_RANGE_NOT_GZIP")
    return bytes(data)


def partial_gzip(compressed):
    """Return decompressed prefix, no gzip EOF claim or full-run MD5 check."""
    if not compressed.startswith(b"\x1f\x8b"):
        raise ValueError("Invalid gzip beginning")
    d=zlib.decompressobj(16+zlib.MAX_WBITS)
    blocks=[]
    expanded=0
    for i in range(0,len(compressed),128*1024):
        part=d.decompress(compressed[i:i+128*1024],EXPANDED_MAX-expanded+1)
        expanded+=len(part)
        if expanded>EXPANDED_MAX or d.unconsumed_tail:
            raise ValueError("DECOMPRESSED_PREFIX_EXCEEDS_CAP")
        blocks.append(part)
    return b"".join(blocks),d.eof


def fastq_records(prefix,max_records=None):
    """Yield COMPLETE 4-line R1 FASTQ records; final cut line is discarded."""
    complete_lines=prefix.split(b"\n")[:-1]
    i=0
    n=0
    while i+3<len(complete_lines):
        h,s,plus,q=[line.rstrip(b"\r") for line in complete_lines[i:i+4]]
        i+=4
        if (not h.startswith(b"@") or not plus.startswith(b"+")
                or not s or len(s)!=len(q)
                or any(b not in b"ACGTNacgtn" for b in s)
                or any(c<33 or c>126 for c in q)):
            raise ValueError(f"Invalid complete FASTQ quartet at read index {n}")
        yield s.decode("ascii").upper(),q
        n+=1
        if max_records is not None and n>=max_records:
            break


def supporting_groups(sequence,quality,index):
    """Require two non-overlapping reference-unique 31-mers on the SAME read.

    The 2020 reference is from a DIFFERENT species, so this is a
    homology-compatible sequence signal only. No mapper/MAPQ, depth
    normalization, replicate or allele call is inferred.
    """
    if len(sequence)!=len(quality):
        raise ValueError("Read quality and bases lengths differ")
    low=[0]
    for c in quality:
        low.append(low[-1]+(c-33<Q_MIN))
    hits=defaultdict(list)
    multi_reference_seed_hits=0
    for pos in range(len(sequence)-K+1):
        if low[pos+K]!=low[pos]:
            continue
        seed=sequence[pos:pos+K]
        groups=index.get(canonical(seed))
        if groups is None:continue
        if len(groups)>1:
            multi_reference_seed_hits+=1
            continue
        hits[groups[0]].append(pos)
    positives=[]
    weak=[]
    for grp,positions in hits.items():
        if positions:
            weak.append(grp)
        if len(positions)>=2 and positions[-1]-positions[0]>=K:
            positives.append(grp)
    return {
        "unique_confirmed_group":positives[0] if len(positives)==1 else None,
        "ambiguous_multiple_confirmed_groups":len(positives)>1,
        "single_seed_or_overlapping_groups":len(weak)>0 and not positives,
        "multi_reference_seed_hits":multi_reference_seed_hits,
    }


def scan_prefix(compressed, real_index, decoy_index, max_records=None):
    expanded,is_complete_gzip=partial_gzip(compressed)
    values=Counter()
    groups=Counter()
    decoy_groups=Counter()
    n=0
    for seq,qual in fastq_records(expanded,max_records):
        n+=1
        value=supporting_groups(seq,qual,real_index)
        mock=supporting_groups(seq,qual,decoy_index)
        if value["unique_confirmed_group"]:
            groups[value["unique_confirmed_group"]]+=1
            values["reads_with_single_group_2nonoverlap_31mers"]+=1
        if mock["unique_confirmed_group"]:
            decoy_groups[mock["unique_confirmed_group"]]+=1
            values["decoy_reads_with_single_group_2nonoverlap_31mers"]+=1
        if value["ambiguous_multiple_confirmed_groups"]:
            values["reads_with_conflicting_multigroup_support"]+=1
        if value["single_seed_or_overlapping_groups"]:
            values["reads_with_unconfirmed_seed_support"]+=1
        values["real_multigroup_reference_seed_positions"]+=value["multi_reference_seed_hits"]
    if n<=0:
        raise ValueError("No complete FASTQ records inside bounded range")
    return {
        "reads_scanned":n,
        "decompressed_prefix_bytes":len(expanded),
        "compressed_prefix_bytes":len(compressed),
        "compressed_prefix_sha256":sha256(compressed).hexdigest(),
        "gzip_end_of_stream_observed":is_complete_gzip,
        "full_file_md5_verified":False,
        "counts":dict(sorted(values.items())),
        "source_assembly_groups_supported":dict(sorted(groups.items())),
        "dinucleotide_preserving_decoy_groups_supported":dict(sorted(decoy_groups.items())),
        "biological_interpretation":"SAME-READ 2 nonoverlapping exact 31-mer high-quality matches to 2020 Cirsium source assembly groups; not orthology, genotype, enzyme function or absence when no hits."
    }


def run_six(limit_bytes=CAP_MAX,manifest=MANIFEST,targets=TARGETS,
            provenance=PROVENANCE,max_records=None):
    runs=validate_manifest(manifest,limit_bytes)
    prepared=indexes(targets,provenance)
    report={
        "analysis_version":"aza3_Ryukyu_R1_prefix_probe_v6",
        "source_project":"PRJNA1311153",
        "genome_tissue":"young leaf",
        "source_reference":"Park 2020 Cirsium japonicum, 26 BLAST-labelled transcripts; NOT focal pigment alleles",
        "downloaded_six_entire_fastqs":False,
        "downloaded_second_mate":False,
        "download_budget_compressed_bytes_per_R1":limit_bytes,
        "planned_compressed_download_max_six":limit_bytes*len(runs),
        "reference_source_n_transcripts":prepared["source_transcripts"],
        "reference_source_assembly_groups":prepared["source_assembly_groups"],
        "real_seed_index_unique_kmers":len(prepared["real"]),
        "dinuc_shuffled_decoy_index_unique_kmers":len(prepared["decoy"]),
        "source_decoy_kmer_overlap_excluded":prepared["decoy_windows_removed_for_source_overlap"],
        "reference_marker_kmer_length":K,
        "min_phred_per_seed":Q_MIN,
        "minimum_per_read_support":"two nonoverlapping quality-filtered, group-unique 31-mers",
        "runs":[],
        "number_of_successfully_read_R1_prefixes":0,
        "valid_Ryukyu_orthology_calls":0,
        "corolla_expression_observations":0,
        "ancestral_pigment_reexpression_confirmed":False,
        "source_2020_transcripts_are_not_focal_2026_alleles":True,
    }
    for source in runs:
        row={
            "run_accession":source["run_accession"],
            "biosample_accession":source["biosample_accession"],
            "specimen_code":source["source_specimen_code"],
            "taxon_as_submitted":source["taxon_as_experiment_titled"],
            "tissue":"young leaf",
            "R1_original_full_size_bytes":source["R1_full_file_bytes"],
            "original_paired_read_count_reported":source["original_read_pairs"],
            "provider_full_R1_md5_not_reevaluated":source["R1_provider_md5"],
        }
        try:
            compressed=fetch_range(source,limit_bytes)
            score=scan_prefix(compressed,prepared["real"],prepared["decoy"],max_records)
            row["status"]="REAL_R1_PREFIX_SCANNED"
            row.update(score)
            row["candidate_signal_status"]=(
                "CROSS_SPECIES_HOMOLOGY_COMPATIBLE_SIGNAL_IN_PREFIX"
                if score["counts"].get("reads_with_single_group_2nonoverlap_31mers",0)
                else "NOT_DETECTED_IN_ORDERED_R1_PREFIX"
            )
            report["number_of_successfully_read_R1_prefixes"]+=1
        except (urllib.error.URLError,TimeoutError,OSError,ValueError,
                EOFError,zlib.error) as exc:
            row["status"]="SOURCE_PREFIX_UNAVAILABLE_OR_INVALID"
            row["retrieval_error_type"]=type(exc).__name__
            row["retrieval_error_message"]=str(exc)[:220]
            row["candidate_signal_status"]="NOT_ASSESSED"
            row["not_a_gene_absence_result"]=True
        report["runs"].append(row)
    report["total_reads_scanned"]=sum(r.get("reads_scanned",0) for r in report["runs"])
    report["total_strict_candidate_reads"]=sum(
        r.get("counts",{}).get("reads_with_single_group_2nonoverlap_31mers",0)
        for r in report["runs"])
    report["total_decoy_support_reads"]=sum(
        r.get("counts",{}).get("decoy_reads_with_single_group_2nonoverlap_31mers",0)
        for r in report["runs"])
    report["claim_ceiling"]=(
        "Partial ordered single-mate young-leaf RNA read matches to another species' unverified assembled transcript, NEVER DFR/ANS functional orthology, loss, regulation in a floret, historical reversal or expression comparison among species."
    )
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--limit-bytes",type=int,default=CAP_MAX)
    p.add_argument("--max-records",type=int,default=None)
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    report=run_six(args.limit_bytes,max_records=args.max_records)
    output=json.dumps(report,sort_keys=True,ensure_ascii=False,indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(output,encoding="utf-8")
    else:
        print(output)
    print("RYUKYU_PREFIX_SUMMARY_JSON:"+json.dumps({
        k:report[k] for k in ("number_of_successfully_read_R1_prefixes",
           "total_reads_scanned","total_strict_candidate_reads",
           "total_decoy_support_reads","reference_source_n_transcripts",
           "source_decoy_kmer_overlap_excluded",
           "ancestral_pigment_reexpression_confirmed")},separators=(",",":")))
    for row in report["runs"]:
        print("RYUKYU_PREFIX_RUN_JSON:"+json.dumps(row,separators=(",",":")))
    return report


if __name__=="__main__":
    main()
