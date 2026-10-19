"""Fail-closed validation of first *real*, six-source Ryukyu leaf R1 homology data.

This rechecks immutable run/BioSample/source transcript joins and the source
snapshot against the 26 original Park 2020 BLAST-annotated transcript IDs.
The empirical conclusion is restricted to *capped R1 prefix fragment
similarity* in leaf RNA, not gene function/orthology/abundance, head colour,
introgression, adaptive fitness or ancestral trait re-expression.

The 16 MiB R1 prefix is nonrandom; R2 and >99% of full run FASTQs were not
used. An unobserved candidate can be expressed elsewhere and is NOT absent.
"""
from collections import Counter
from pathlib import Path
import csv
import json

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/anto_source_handoff"
RAW=DATA/"public_sra"
GENES=DATA/"molecular"
EVIDENCE=RAW/"ryukyu_six_leaf_R1_fragment_hits_v6.json"
TABLE=RAW/"ryukyu_six_R1_candidate_fragment_matrix_v6.csv"
SOURCE=RAW/"ryukyu_six_R1_prefix_source_v6.json"
CANDIDATES=GENES/"cirsium_japonicum_2020_26_sequence_provenance_v6.json"


def _json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def validate(evidence=EVIDENCE,table=TABLE,source=SOURCE,candidates=CANDIDATES):
    d=_json(evidence)
    registry=_json(source)
    original=_json(candidates)
    if d.get("version")!="aza3_Ryukyu_six_public_leaf_R1_candidate_fragment_results_v6":
        raise ValueError("Unexpected six-sample empirical result version")
    if registry.get("version")!="aza3_cirsium_Ryukyu_public_R1_prefix_manifest_v6":
        raise ValueError("Unexpected independent source manifest")
    if original.get("contract")!="anto_park2020_26_source_transcripts_v6":
        raise ValueError("Reference source provenance missing")
    source_by_run={r["run_accession"]:r for r in registry["runs"]}
    if len(source_by_run)!=6:
        raise ValueError("Must have exactly six source R1 accession registrations")
    candidate_ids={r["id"] for r in original["records"]}
    if len(candidate_ids)!=26:
        raise ValueError("Original source candidate family reference IDs changed")
    rows=d.get("rows",[])
    if len(rows)!=6 or len({r["run_accession"] for r in rows})!=6:
        raise ValueError("Six disjoint biological accession rows required")
    if d["summary"]["C_brevicaule_n"]!=3 or d["summary"]["C_irumtiense_n"]!=3:
        raise ValueError("One-to-one Ryukyu source/taxon join not checked")
    if d["read_processing"]["compressed_prefix_bytes_each"]!=16777216:
        raise ValueError("Frozen per-run budget differs from audited source execution")
    if d["read_processing"]["total_compressed_bytes_retrieved"]!=6*16777216:
        raise ValueError("Frozen compressed bytes exceed or differ from advertised cap")
    if d["read_processing"]["read_order_random_sample"] or not d["read_processing"]["R1_only"]:
        raise ValueError("Misrepresented R1 prefix as random complete RNA-seq")
    for fact in ("gene_copy_numbers_estimated",
                 "differential_expression_conclusion_allowed",
                 "confirmed_orthology",
                 "pigment_structural_gene_function_verified",
                 "individual_corolla_pigment_and_expression_measured",
                 "historical_ancestral_reexpression_event_inferred"):
        if d["claim_ceiling"][fact] is not False:
            raise ValueError(f"Unsupported scientific claim: {fact}")
    if d["claim_ceiling"]["congeneric_source_candidate_fragment_signal_in_both_taxa"] is not True:
        raise ValueError("Cross-taxon observation scope changed without new evidence")

    counts=Counter()
    total_reads=total_matched=total_decoys=total_39766=0
    observed_ids=set()
    for r in rows:
        orig=source_by_run.get(r["run_accession"])
        if (orig is None or r["sample_accession"]!=orig["biosample_accession"]
            or r["specimen_code"]!=orig["source_specimen_code"]
            or r["taxon_as_study_title"]!=orig["taxon_as_experiment_titled"]
            or r["tissue"]!="young leaf"):
            raise ValueError("Run, source taxon, BioSample or organ mismatch")
        if r["R1_compressed_prefix_bytes"]!=16777216:
            raise ValueError("Sampled FASTQ prefix length not equal across six runs")
        source_hits=r["matches_by_2020_Cirsium_transcript_ID"]
        if not set(source_hits)<=candidate_ids:
            raise ValueError("Source homology result contains unregistered target")
        if any(type(v)is not int or v<=0 for v in source_hits.values()):
            raise ValueError("Individual matching source transcripts must carry positive integer counts")
        if (sum(source_hits.values()) != r["high_confidence_source_candidate_mapped_fragments"]
            or len(source_hits)!=r["n_matching_selected_source_transcripts"]):
            raise ValueError("Per-sample transcript count accounting inconsistent")
        if r["R1_fragments_to_2020_TR39766"]!=source_hits.get("TR39766|c0_g1_i2",0):
            raise ValueError("TR39766 candidate evidence modified")
        if not (r["nuc_reference_panel_not_complete_transcriptome"] and
                not r["transcript_copy_number_known"] and
                not r["full_gene_sequence_known"] and
                not r["flower_specific_expression_known"] and
                not r["ancestral_reexpression_demonstrated"]):
            raise ValueError("Scope claims were silently promoted beyond source reads")
        counts[r["taxon_as_study_title"]]+=1
        total_reads+=r["complete_R1_reads_passing_QC"]
        total_matched+=r["high_confidence_source_candidate_mapped_fragments"]
        total_decoys+=r["shuffled_reference_decoy_hits"]
        total_39766+=r["R1_fragments_to_2020_TR39766"]
        observed_ids.update(source_hits)
    if counts!={"Cirsium brevicaule":3,"Cirsium irumtiense":3}:
        raise ValueError("Group sample counts changed")
    if (total_reads != d["summary"]["accepted_complete_R1_reads_total"]
        or total_matched!=d["summary"]["high_confidence_mapped_source_candidate_fragments"]
        or total_decoys!=d["summary"]["shuffled_sequence_decoy_fragment_matches"]
        or total_39766!=d["summary"]["candidate_TR39766_fragment_matches_in_all_six"]
        or len(observed_ids)!=d["summary"]["tested_2020_source_candidate_ids_with_any_hit"]):
        raise ValueError("Aggregate source evidence does not equal per-run observations")
    # Version-frozen result from actual first 6-run source execution; a
    # different prefix budget or updated source needs a NEW version.
    if (total_reads,total_matched,total_decoys,total_39766,len(observed_ids))!=(
        1499192,636,0,161,18):
        raise ValueError("First validated 6-source R1 result fingerprint differs")
    if any(x["R1_fragments_to_2020_TR39766"]<1 for x in rows):
        raise ValueError("Original source TR39766 fragment observations changed")

    with Path(table).open(newline="",encoding="utf-8") as f:
        tab=list(csv.DictReader(f))
    if len(tab)!=6 or {r["run_accession"] for r in tab}!={r["run_accession"] for r in rows}:
        raise ValueError("Companion source fragment matrix missing/swapped")
    by_run={r["run_accession"]:r for r in rows}
    for row in tab:
        expected=by_run[row["run_accession"]]
        for k in ("sample_accession","specimen_code","taxon_as_study_title"):
            if row[k]!=expected[k]:
                raise ValueError("Source biological ID drift in report CSV")
        for k in ("complete_R1_reads_passing_QC",
                  "high_confidence_source_candidate_mapped_fragments",
                  "n_matching_selected_source_transcripts",
                  "R1_fragments_to_2020_TR39766",
                  "shuffled_reference_decoy_hits"):
            if int(row[k])!=expected[k]:
                raise ValueError("Frozen six-source report CSV count mismatch")
        for candidate in candidate_ids:
            if candidate in row and int(row[candidate])!=expected[
                "matches_by_2020_Cirsium_transcript_ID"].get(candidate,0):
                raise ValueError("Transcript-by-sample table conflicts with source evidence")
    return {
        "version":d["version"],
        "cross_source_6_of_6_verified":True,
        "independent_source_accessions":6,
        "real_leaf_R1_reads_processed":total_reads,
        "candidate_source_transcript_R1_fragment_alignments":total_matched,
        "original_2020_candidate_ids_detected_in_selected_reference_panel":len(observed_ids),
        "TR39766_mapping_reads_across_all_six":total_39766,
        "shuffled_decoy_matches":total_decoys,
        "not_an_entire_RNA_library":True,
        "not_a_random_expression_sample":True,
        "no_gene_function_or_flower_reexpression_inference":True,
        "complete_2020_transcriptome_reference_competition_done":False,
    }


if __name__=="__main__":
    print(json.dumps(validate(),sort_keys=True,indent=2,ensure_ascii=False))
