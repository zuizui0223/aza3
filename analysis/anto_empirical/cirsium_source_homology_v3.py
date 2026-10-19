# Source migrated from zuizui0223/anto/src/anto_2020_homology_triage_v6.py (anto v4-v6); empirical ownership now aza3.
"""Protein-level homology TRIAGE on 26 real, author-published Cirsium transcripts.

Two Asteraceae reference proteins are derived from verified annotated CDS:
*Saussurea medusa* DFR EF600682.1 and ANS AY547342.1.

For each *Cirsium japonicum* 2020 transcript, enumerate translated segments
from six reading frames and compare local BLOSUM62 alignment to both anchors.
This tests **reference similarity**, NOT single-copy orthology, enzyme
function, floral expression, chromosomal ancestry, evolutionary reversal,
or any allele of the six 2026 Ryukyu *Cirsium* individuals.

Biopython alignment scores and thresholding are exploratory triage and
MUST NOT be interpreted as biochemical significance or gene count.
"""
import hashlib
import json
import re
from pathlib import Path
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices

REFS={"DFR":"EF600682.1","ANS":"AY547342.1"}
MIN_ORF_AA=35
REFERENCE_COVERAGE_SCREEN=.50
PERCENT_IDENT_SCREEN=.30
ALIGNED_AA_SCREEN=100


def read_fasta(path):
    rows={}
    current=None
    parts=[]
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line=line.strip()
        if line.startswith(">"):
            if current is not None:
                if current in rows: raise ValueError("Duplicate transcript ID")
                rows[current]="".join(parts).upper()
            current=line[1:].split()[0]
            if not current: raise ValueError("Missing sequence identifier")
            parts=[]
        elif line:
            if current is None: raise ValueError("Sequence without FASTA identifier")
            parts.append(line)
    if current is not None:
        if current in rows:raise ValueError("Duplicate FASTA identifier")
        rows[current]="".join(parts).upper()
    return rows


def ref_proteins(fasta_path):
    nt=read_fasta(fasta_path)
    out={}
    for gene, accession in REFS.items():
        hits=[s for k,s in nt.items() if k.startswith(accession+"|"+gene+"|")]
        if len(hits)!=1 or len(hits[0])%3:
            raise ValueError("Reference protein CDS absent or invalid")
        p=str(Seq(hits[0]).translate())
        if not p.endswith("*") or "*" in p[:-1]:
            raise ValueError("Invalid reference protein translation")
        out[gene]=p[:-1]
    return out


def six_frame_segments(nucleotides, min_aa=MIN_ORF_AA):
    """Yield longest stop-delimited translated segments in both orientations.

    RNA/U and IUPAC ambiguity are translated with X; literal alignment gaps
    are replaced with N. Partial segments are allowed, not full CDS claims.
    """
    if not nucleotides:
        return
    seq=nucleotides.upper().replace("U","T").replace("-","N")
    if not set(seq).issubset(set("ACGTNRYKMSWBDHV")):
        raise ValueError("Transcript contains unsupported base")
    for strand, oriented in (("+",seq),("-",str(Seq(seq).reverse_complement()))):
        for frame in range(3):
            n=((len(oriented)-frame)//3)*3
            if n < min_aa*3:
                continue
            aa=str(Seq(oriented[frame:frame+n]).translate())
            for match in re.finditer(r"[^*]+",aa):
                segment=match.group(0)
                if len(segment)>=min_aa:
                    yield {"strand":strand,"reading_frame":frame,
                           "aa_segment":segment,
                           "aa_offset":match.start(),
                           "segment_length_aa":len(segment)}


def _aligner():
    a=PairwiseAligner()
    a.mode="local"
    a.substitution_matrix=substitution_matrices.load("BLOSUM62")
    a.open_gap_score=-10
    a.extend_gap_score=-.5
    return a


def alignment_summary(reference,candidate,aligner=None):
    """Local protein alignment and precise identity / reference coverage."""
    a=_aligner() if aligner is None else aligner
    alignment=a.align(reference,candidate)[0]
    blocks_ref,blocks_query=alignment.aligned
    matches=aligned=0
    for (rs,re),(qs,qe) in zip(blocks_ref,blocks_query):
        if int(re-rs)!=int(qe-qs):
            raise ArithmeticError("Paired aligned blocks must have same length")
        aa1=reference[int(rs):int(re)]
        aa2=candidate[int(qs):int(qe)]
        aligned+=len(aa1)
        matches+=sum(x==y for x,y in zip(aa1,aa2))
    return {"protein_alignment_score":round(float(alignment.score),4),
            "aligned_aa":int(aligned),
            "exact_matches":int(matches),
            "aligned_identity":round(matches/aligned,6) if aligned else 0.0,
            "reference_coverage":round(aligned/len(reference),6),
            "normalized_alignment_score":round(float(alignment.score)/len(reference),6)}


def triage_sequence(seq, proteins):
    segments=list(six_frame_segments(seq))
    a=_aligner()
    by_gene={}
    for gene,reference in proteins.items():
        best=None
        for segment in segments:
            summ=alignment_summary(reference,segment["aa_segment"],a)
            row={**summ,
                 "orientation":segment["strand"],
                 "reading_frame":segment["reading_frame"],
                 "translated_segment_length":segment["segment_length_aa"],
                 "aa_offset":segment["aa_offset"]}
            if best is None or row["protein_alignment_score"]>best["protein_alignment_score"]:
                best=row
        by_gene[gene]=best
    if all(v is None for v in by_gene.values()):
        return {"family_screen":"NO_ORF_LONG_ENOUGH",
                "by_reference":by_gene,"orthology_confirmed":False}
    winner=max((g for g in by_gene if by_gene[g] is not None),
               key=lambda g:by_gene[g]["normalized_alignment_score"])
    w=by_gene[winner]
    passes=(w["reference_coverage"]>=REFERENCE_COVERAGE_SCREEN
            and w["aligned_identity"]>=PERCENT_IDENT_SCREEN
            and w["aligned_aa"]>=ALIGNED_AA_SCREEN)
    return {"family_screen":"STRONG_REFERENCE_SIMILARITY_CANDIDATE" if passes else "WEAK_OR_PARTIAL_REFERENCE_SIMILARITY",
            "best_reference_family":winner,
            "highest_score_features":w,
            "by_reference":by_gene,
            "thresholds":{"min_reference_coverage":REFERENCE_COVERAGE_SCREEN,
                          "min_aligned_identity":PERCENT_IDENT_SCREEN,
                          "min_aligned_aa":ALIGNED_AA_SCREEN},
            "orthology_confirmed":False,
            "enzyme_function_confirmed":False}


def run_study(source_fasta,source_metadata,reference_fasta):
    data=json.loads(Path(source_metadata).read_text(encoding="utf-8"))
    if data.get("contract")!="anto_park2020_26_source_transcripts_v6":
        raise ValueError("Required verified Cirsium 2020 transcript provenance")
    records=data.get("records")
    fasta=read_fasta(source_fasta)
    if len(records)!=26 or len(fasta)!=26:
        raise ValueError("Exactly 26 vetted published transcripts required")
    if set(fasta)!={r["id"] for r in records}:
        raise ValueError("FASTA identifiers differ from published-source metadata")
    proteins=ref_proteins(reference_fasta)
    result=[]
    for record in records:
        key=record["id"]
        seq=fasta[key]
        if len(seq)!=record["source_fasta_length"] or (
            hashlib.sha256(seq.encode("ascii")).hexdigest()!=record["sha256_nucleotide"]
        ):
            raise ValueError(f"Source transcript missing, corrupted or swapped: {key}")
        res=triage_sequence(seq,proteins)
        labelled=record["family_by_publisher_annotation_only"].startswith("DFR")
        predicted=res.get("best_reference_family")
        concordance=(
            predicted==("DFR" if labelled else "ANS")
            if predicted else None
        )
        result.append({
            "source_transcript_id":key,
            "published_best_blast_hit":record["published_blast_annotation"],
            "publisher_annotated_family":record["family_by_publisher_annotation_only"],
            "sequence_length_nt":len(seq),
            "recorded_source_gene_group":key.rsplit("_i",1)[0],
            "blast_annotation_and_reference_similarity_concordant":concordance,
            "potential_nonplant_best_hit_warning":"Leptospira" in record["published_blast_annotation"],
            **res
        })
    return {
        "analysis_version":"anto_2020_cirsium_homology_triage_v6",
        "source":"Park 2020, DOI 10.1021/acsomega.0c04001; real author transcript sequences",
        "reference_species":"Saussurea medusa (DFR, ANS)",
        "n_transcripts":len(result),
        "n_distinct_assembly_groups":len({r["recorded_source_gene_group"] for r in result}),
        "n_reference_screen_positive":sum(r["family_screen"]=="STRONG_REFERENCE_SIMILARITY_CANDIDATE" for r in result),
        "n_publisher_and_reference_family_concordant":sum(r["blast_annotation_and_reference_similarity_concordant"] is True for r in result),
        "n_publisher_and_reference_family_discordant":sum(r["blast_annotation_and_reference_similarity_concordant"] is False for r in result),
        "n_nonplant_best_hit_warnings":sum(r["potential_nonplant_best_hit_warning"] for r in result),
        "records":result,
        "evidence_scope":"protein sequence similarity among published Cirsium japonicum transcript candidates",
        "orthology_resolved":False,
        "Ryukyu_2026_sequences_included":False,
        "floral_reexpression_proven":False
    }
