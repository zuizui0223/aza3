"""Fail-closed audit of genuine six-Ryukyu leaf-RNA full 2020-source competition.

Key positive finding: Short fragments from the same six source-confirmed leaf
libraries can align to 2020 Cirsium DFR/ANS-like sequences when ALL 104890
published Park source transcripts compete, not only the selected 26.

Key negative/uncertainty: short read homology in a *young leaf* is NOT
functional pigmentation enzyme evidence, whole genomic locus retention,
floret competence, colour regain, independent biological time-course,
or differential expression across species.

This validator runs OFFLINE; the source archive and raw RNA files are NOT
redownloaded on every GitHub commit.
"""
from pathlib import Path
from collections import Counter
import json

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"data/anto_source_handoff"
DEFAULT=DIR/"public_sra/ryukyu_six_full_104890_unigene_competitive_mapping_v7.json"
RUN_SOURCE=DIR/"public_sra/ryukyu_six_R1_prefix_source_v6.json"
CANDIDATES=DIR/"molecular/cirsium_japonicum_2020_26_sequence_provenance_v6.json"


def audit(report=DEFAULT,manifest=RUN_SOURCE,refs=CANDIDATES):
    full=json.loads(Path(report).read_text(encoding="utf-8"))
    src=json.loads(Path(manifest).read_text(encoding="utf-8"))
    ref=json.loads(Path(refs).read_text(encoding="utf-8"))
    if full.get("version")!="aza3_Ryukyu_six_public_leaf_R1_whole_2020_source_competition_v7":
        raise ValueError("Full-reference result version is not verified")
    if src.get("version")!="aza3_cirsium_Ryukyu_public_R1_prefix_manifest_v6" or len(src["runs"])!=6:
        raise ValueError("Wrong six-source Ryukyu run manifest")
    if ref.get("contract")!="anto_park2020_26_source_transcripts_v6" or len(ref["records"])!=26:
        raise ValueError("Wrong 2020 full-reference candidate provenance")
    reference_ids={r["id"] for r in ref["records"]}
    base=full["reference_source_validation"]
    if (base["n_public_2020_unigene_entries"]!=104890
            or base["n_original_pigment_family_candidates_in_full_reference"]!=26
            or base["source_uncompressed_fasta_sha256"]!=
              "0c7b7d61a8e2ad5e89b9199347a50545be707329715cdbe0a8a351b0eef2815f"
            or base["not_complete_Ryukyu_genome"] is not True):
        raise ValueError("Source reference count/hash did not match published FASTA")
    process=full["read_processing"]
    if (process["compressed_prefix_bytes_per_run"]!=8388608
            or process["compressed_prefix_total_bytes"]!=6*8388608
            or process["original_R2_full_library_downloaded"] is not False
            or process["full_original_R1_file_MD5_verified"] is not False
            or process["is_random_representative_read_sample"] is not False):
        raise ValueError("Claim or budget exceeds bounded, biased source RNA data")
    if len(full.get("source_specimens",[]))!=6:
        raise ValueError("Expected six independent public source libraries")
    byid={r["run_accession"]:r for r in src["runs"]}
    observed=Counter()
    total_reads=total_strict=total_dfr=total_ans=total_null=0
    taxa=Counter()
    for r in full["source_specimens"]:
        run=r["run_accession"]
        if run in observed:
            raise ValueError("Repeated biological source RNA library")
        observed[run]+=1
        original=byid.get(run)
        if (original is None or
                r["biosample_accession"]!=original["biosample_accession"] or
                r["specimen_code"]!=original["source_specimen_code"] or
                r["taxon_as_source_title"]!=original["taxon_as_experiment_titled"] or
                r["analysis_tissue"]!="young leaf"):
            raise ValueError("Source run/BioSample/taxon identity changed")
        if (r["compressed_R1_prefix_bytes_requested"]!=8388608 or
                r["strict_shuffled_decoy_read_fragments"]!=0 or
                r["strict_source_candidates_from_original_26"]<
                   r["strict_read_fragments_assigned_TR39766"]+
                   r["strict_read_fragments_assigned_TR41106"]):
            raise ValueError("Per-source strict mapping totals or decoy data changed")
        if r["strict_read_fragments_assigned_TR39766"]<=0 or r["strict_read_fragments_assigned_TR41106"]<=0:
            raise ValueError("Either of two published candidate IDs no longer supported in all six source libraries")
        if r["is_2020_transcript_assay"] or r["is_individual_flower_colour_measurement"]:
            raise ValueError("Young leaf R1 cannot become a flower/2020 individual assay")
        taxa[r["taxon_as_source_title"]]+=1
        total_reads+=r["complete_R1_reads_analysed"]
        total_strict+=r["strict_source_candidates_from_original_26"]
        total_dfr+=r["strict_read_fragments_assigned_TR39766"]
        total_ans+=r["strict_read_fragments_assigned_TR41106"]
        total_null+=r["strict_shuffled_decoy_read_fragments"]
    if (set(observed)!=set(byid) or
            taxa!={"Cirsium brevicaule":3,"Cirsium irumtiense":3}):
        raise ValueError("Original six source samples or taxon balance drifted")
    stats=full["observed_results"]
    if (total_reads,total_strict,total_dfr,total_ans,total_null)!=(747234,215,66,106,0):
        raise ValueError("Results no longer match original source-affiliated workflow")
    if (total_reads!=process["n_complete_R1_prefix_reads"] or
            total_strict!=stats["all_six_candidate_source_transcript_strict_read_fragments"] or
            total_dfr!=stats["TR39766_DFR_related_candidate_fragments"] or
            total_ans!=stats["TR41106_ANS_Blast_annotated_candidate_fragments"] or
            total_null!=stats["strict_shuffled_decoy_fragments"] or
            stats["n_other_selected_candidate_fragments"]!=43 or
            stats["n_taxon_libraries_TR39766_detected"]!=6 or
            stats["n_taxon_libraries_TR41106_detected"]!=6):
        raise ValueError("Frozen source totals and per-individual read evidence disagree")
    if {"TR39766|c0_g1_i2","TR41106|c0_g1_i5"}-reference_ids:
        raise ValueError("Positive source transcript not in the 26 original publisher records")
    for label,status in full["claim_ceiling"].items():
        if status is not False:
            raise ValueError("Unsupported phenotype/genotype/orthology elevation: "+label)
    if not full["importantly_not_independent_biological_validation_of_26_transcript_run"]:
        raise ValueError("The source read prefix is not independent biological replication")
    return {
        "evidence_version":full["version"],
        "source_Cirsium_2020_full_transcripts_validated":104890,
        "focal_Ryukyu_2026_R1_biological_libraries":6,
        "young_leaf_R1_prefix_reads":total_reads,
        "full_104890_reference_competition_selected_26_candidate_matches":total_strict,
        "full_104890_reference_competition_TR39766_fragments":total_dfr,
        "full_104890_reference_competition_TR41106_fragments":total_ans,
        "full_104890_reference_competition_shuffled_null_fragments":total_null,
        "matching_TR39766_and_TR41106_each_in_all_six_R1":True,
        "focal_gene_orthology_verified":False,
        "functional_flower_anthocyanin_proven":False,
        "reexpression_re_evolution_supported":False,
        "source_run_url":full["source_run"]
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,ensure_ascii=False,sort_keys=True))
