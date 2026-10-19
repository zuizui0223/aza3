"""Cirsium 2020 source transcript family specificity against seven real proteins.

Biological question: do 2020 publisher BLAST-labelled DFR/ANS-like
transcripts align to DFR/ANS more strongly than to RELATED proteins:
ANR/BANYULS (DFR-like reductase) and FLS1 (ANS-like oxygenase)?

This is a source-backed *paralog ambiguity screen*, NOT an inference of
enzyme substrate specificity, orthology, copy number, floral gene expression
or historical anthocyanin reacquisition in the Ryukyu focal pair.

Reference protein sequences are accession-specific NCBI records. Use the
published 26 Cirsium japonicum transcript sequences and their frozen source
metadata within aza3. All threshold outcomes remain provisional.
"""
import hashlib
import json
from pathlib import Path

from cirsium_source_homology_v3 import (
    read_fasta, six_frame_segments, alignment_summary, _aligner
)

GROUPS=("DFR_reference","ANS_reference",
        "ANR_negative_control","FLS_negative_control")
POSITIVE_GROUPS=("DFR_reference","ANS_reference")
CAUTION_GROUP={
    "DFR_reference":"ANR_negative_control",
    "ANS_reference":"FLS_negative_control",
}
MIN_ALIGNED_AA=100
MIN_REFERENCE_COVERAGE=.50
MIN_IDENTITY=.30


def load_panel(reference_fasta, reference_provenance):
    raw=json.loads(Path(reference_provenance).read_text(encoding="utf-8"))
    if raw.get("version")!="aza3_cirsium_DFR_ANS_related_family_references_v4":
        raise ValueError("Unexpected accession-verified reference panel version")
    entries=raw.get("references")
    if not isinstance(entries,list) or len(entries)!=7:
        raise ValueError("Expected exactly seven positive and negative reference proteins")
    fasta=read_fasta(reference_fasta)
    if len(fasta)!=7:
        raise ValueError("Reference FASTA must have 7 protein records")
    ids=set()
    proteins=[]
    for r in entries:
        if r["group"] not in GROUPS or r["accession"] in ids:
            raise ValueError("Unsupported/duplicate reference class or accession")
        ids.add(r["accession"])
        header=f'{r["accession"]}|{r["group"]}|{r["taxon"].replace(" ","_")}'
        if header not in fasta:
            raise ValueError(f"Reference accession or taxon header mismatch: {header}")
        prot=fasta[header]
        if not prot or not set(prot)<=set("ACDEFGHIKLMNPQRSTVWXY") or len(prot)!=r["length_aa"]:
            raise ValueError("Reference protein sequence invalid")
        proteins.append({
            "accession":r["accession"],
            "group":r["group"],
            "gene":r["gene"],
            "taxon":r["taxon"],
            "protein":prot,
        })
    count={g:sum(r["group"]==g for r in proteins) for g in GROUPS}
    if count!={"DFR_reference":3,"ANS_reference":2,
               "ANR_negative_control":1,"FLS_negative_control":1}:
        raise ValueError("Panel lacks necessary DFR/ANS or related-family negatives")
    return proteins


def match_protein_segments(segments, panel, aligner=None):
    """Local BLOSUM62 score after searching both strands and three frames.

    Normalization by reference AA length makes raw scores more comparable,
    but it is not a likelihood ratio or an orthology posterior. Scores may
    differ by protein completeness and domain architecture.
    """
    a=_aligner() if aligner is None else aligner
    groups={}
    for r in panel:
        best=None
        for seg in segments:
            z=alignment_summary(r["protein"],seg["aa_segment"],a)
            z.update({
                "subject_accession":r["accession"],
                "subject_taxon":r["taxon"],
                "subject_gene":r["gene"],
                "group":r["group"],
                "orientation":seg["strand"],
                "frame":seg["reading_frame"],
                "segment_length_aa":seg["segment_length_aa"],
            })
            if best is None or (
                z["normalized_alignment_score"],z["aligned_aa"]
            ) > (best["normalized_alignment_score"],best["aligned_aa"]):
                best=z
        if best is None:
            continue
        old=groups.get(r["group"])
        if old is None or (
            best["normalized_alignment_score"],best["aligned_aa"]
        )>(old["normalized_alignment_score"],old["aligned_aa"]):
            groups[r["group"]]=best
    return groups


def analyze_transcript(source_sequence,panel,aligner=None):
    segments=list(six_frame_segments(source_sequence))
    scores=match_protein_segments(segments,panel,aligner)
    if not scores:
        return {
            "screen":"NO_TRANSLATED_SEGMENT_LONG_ENOUGH",
            "strong_reference_similarity":False,
            "best_family":"UNRESOLVED",
            "all_group_best":{},
            "orthology_verified":False,
            "enzyme_function_verified":False,
        }
    top_group=max(scores,key=lambda g:scores[g]["normalized_alignment_score"])
    top=scores[top_group]
    above=(top["aligned_aa"]>=MIN_ALIGNED_AA
           and top["reference_coverage"]>=MIN_REFERENCE_COVERAGE
           and top["aligned_identity"]>=MIN_IDENTITY)
    top_positive=max(POSITIVE_GROUPS,key=lambda g:
                     scores.get(g,{"normalized_alignment_score":-1})["normalized_alignment_score"])
    pos=scores.get(top_positive)
    neg=scores.get(CAUTION_GROUP[top_positive])
    gap=None if pos is None or neg is None else round(
        pos["normalized_alignment_score"]-neg["normalized_alignment_score"],6)
    # All classifications are *reference evidence*, never enzyme certainty.
    if above and top_group in POSITIVE_GROUPS:
        state="POSITIVE_FAMILY_REFERENCE_SUPPORTED_PARALOG_STILL_POSSIBLE"
    elif above and top_group not in POSITIVE_GROUPS:
        state="RELATED_NEGATIVE_FAMILY_REFERENCE_PREFERRED"
    else:
        state="WEAK_PARTIAL_OR_AMBIGUOUS_FAMILY"
    return {
        "screen":state,
        "best_family":top_group,
        "best_hit":top,
        "strong_reference_similarity":above,
        "best_positive_family":top_positive,
        "best_positive_vs_related_negative_normalized_gap":gap,
        "all_group_best":scores,
        "orthology_verified":False,
        "enzyme_function_verified":False,
        "floral_pathway_function_verified":False,
        "diagnostic_limit":"Local protein identity/coverage, not enzyme specificity, gene tree, dosage, floral expression or ancestry",
    }


def run_panel(source_fasta,source_metadata,reference_fasta,reference_metadata):
    pdata=json.loads(Path(source_metadata).read_text(encoding="utf-8"))
    if pdata.get("contract")!="anto_park2020_26_source_transcripts_v6":
        raise ValueError("Expected frozen 2020 Cirsium published transcripts")
    source=read_fasta(source_fasta)
    records=pdata.get("records",[])
    if len(records)!=26 or len(source)!=26 or set(source)!={r["id"] for r in records}:
        raise ValueError("Wrong source candidate record count or gene IDs")
    panel=load_panel(reference_fasta,reference_metadata)
    a=_aligner()
    results=[]
    for record in records:
        nt=source[record["id"]]
        if len(nt)!=record["source_fasta_length"] or (
            hashlib.sha256(nt.encode("ascii")).hexdigest()!=record["sha256_nucleotide"]
        ):
            raise ValueError(f"Corrupted source transcript: {record['id']}")
        score=analyze_transcript(nt,panel,a)
        p_family="DFR_reference" if record["family_by_publisher_annotation_only"].startswith("DFR") else "ANS_reference"
        results.append({
            "source_transcript_id":record["id"],
            "publisher_family_annotation":record["family_by_publisher_annotation_only"],
            "publisher_best_hit":record["published_blast_annotation"],
            "publisher_gene_group":record["id"].rsplit("_i",1)[0],
            "publisher_group_prefers_same_family_under_extended_panel":
                score.get("best_family")==p_family,
            "source_sequence_nt":len(nt),
            **score,
        })
    counts={g:sum(x["best_family"]==g for x in results) for g in GROUPS}
    strong=[x for x in results if x["strong_reference_similarity"]
            and x["best_family"] in POSITIVE_GROUPS]
    negative=[x for x in results if x["best_family"] not in POSITIVE_GROUPS
              and x["strong_reference_similarity"]]
    return {
        "analysis_version":"aza3_2020_cirsium_paralog_control_v4",
        "study":"Park et al. 2020 Cirsium japonicum published transcript set",
        "n_publisher_cirsium_transcripts":len(results),
        "n_source_assembly_groups":len({x["publisher_gene_group"] for x in results}),
        "n_reference_proteins":len(panel),
        "protein_reference_counts":{
            g:sum(x["group"]==g for x in panel) for g in GROUPS},
        "family_top_hit_counts":counts,
        "positive_family_supported_by_extended_panel":len(strong),
        "negative_related_family_strong_matches":len(negative),
        "strong_positive_candidate_ids":[x["source_transcript_id"] for x in strong],
        "negative_related_family_candidate_ids":[x["source_transcript_id"] for x in negative],
        "records":results,
        "functional_orthology_demonstrated":False,
        "2026_Ryukyu_focal_genotypes_included":False,
        "flower_colour_regain_demonstrated":False,
        "uncertainty":"Protein family reference similarity may depend on segment, isoform, partial coding region and taxonomic sampling; gene tree + biochemical assays required"
    }
