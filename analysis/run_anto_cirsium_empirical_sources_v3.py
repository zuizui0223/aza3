"""aza3 is the canonical executable source-audit home for all Cirsium cases.

Public 2026 Ryukyu young-leaf data: *only accession provenance*, not
flower-expression assays or evidence of pigment reactivation.

Public 2020 Cirsium japonicum transcripts: 26 real transcript sequences
have published best-hit DFR/ANS-like names. Six-frame protein local
alignment against Saussurea reference proteins is a TRIAGE, never
phylogenetic orthology, flower developmental capability, Cirsium
brevicaule/irumtiense genotype evidence, or a genetic LEGO finding.

General history-dependent accessibility theory remains in anto, but all
empirical Cirsium source reads, alignment checks, samples and experiments
are now managed here in aza3.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
MODULES=ROOT/"analysis/anto_empirical"
sys.path.insert(0,str(MODULES))

from cirsium_run_intake_v3 import audit as sra_audit
from cirsium_molecular_reference_v3 import screen as reference_screen
from cirsium_organ_evidence_v3 import report as organ_report
from cirsium_source_homology_v3 import run_study as homology_study

SOURCES=ROOT/"data/anto_source_handoff"


def main(output=None):
    public=sra_audit(
        SOURCES/"public_sra/prjna1311153_run_manifest_20261017.csv",
        SOURCES/"public_sra/takaoense_aza3_morph_crosswalk_20261017.csv",
        SOURCES/"public_sra/prjna1311153_ena_platform_fingerprints_20261017.json"
    )
    preflight=reference_screen(
        SOURCES/"molecular/asteraceae_reference_cds_v5.fasta",
        SOURCES/"molecular/asteraceae_reference_provenance_v5.json",
        SOURCES/"molecular/ryukyu_gene_screen_register_v5.csv",
        SOURCES/"public_sra/prjna1311153_run_manifest_20261017.csv"
    )
    organs=organ_report(
        SOURCES/"molecular/cirsium_japonicum_organ_evidence_v6.json"
    )
    proteins=homology_study(
        SOURCES/"molecular/cirsium_japonicum_2020_26_source_sequences_v6.fasta",
        SOURCES/"molecular/cirsium_japonicum_2020_26_sequence_provenance_v6.json",
        SOURCES/"molecular/asteraceae_reference_cds_v5.fasta"
    )
    if proteins["n_transcripts"]!=26 or public["manifest_n_runs"]!=25:
        raise ValueError("Unexpected Cirsium biological source count")
    if proteins["Ryukyu_2026_sequences_included"] or proteins["floral_reexpression_proven"]:
        raise ValueError("Cirsium evidence claimed beyond admitted genotype/tissue scopes")
    if preflight["n_candidate_sequences_supplied"]!=0 or preflight["focal_gene_state"]!="UNMEASURED":
        raise ValueError("Cannot assume Ryukyu focal DFR/ANS alleles from 2020 source")
    if organs["ryukyu_competence_assessed"]:
        raise ValueError("Organ evidence unrelated to focal Ryukyu genomic sample")
    # Frozen source-derived shortlist is a reproduction target, not a
    # declaration of 5 enzyme loci or biological flower-colour function.
    short_path=SOURCES/"molecular/cirsium_japonicum_2020_source_similarity_shortlist_v3.json"
    short=json.loads(short_path.read_text(encoding="utf-8"))
    if short.get("version")!="aza3_cirsium_japonicum_source_similarity_shortlist_v3":
        raise ValueError("Unexpected source similarity evidence contract")
    current={r["source_transcript_id"]:r for r in proteins["records"]}
    frozen={r["id"]:r for r in short["records"]}
    passing={r["source_transcript_id"] for r in proteins["records"]
             if r["family_screen"]=="STRONG_REFERENCE_SIMILARITY_CANDIDATE"}
    if set(frozen)!=passing or len(frozen)!=5:
        raise ValueError("New source run no longer matches the five accepted reference-similarity candidates")
    for name,row in frozen.items():
        current_row=current.get(name)
        if (current_row is None or
            current_row["best_reference_family"]!=row["best_family"] or
            current_row["publisher_annotated_family"]!=row["annotated_family"] or
            current_row["highest_score_features"]["aligned_aa"]!=row["aa_aligned"] or
            abs(current_row["highest_score_features"]["aligned_identity"]-row["identity"])>1e-6 or
            abs(current_row["highest_score_features"]["reference_coverage"]-row["reference_coverage"])>1e-6):
            raise ValueError(f"Frozen source-protein homology metric changed: {name}")
    original_discordant={r["id"] for r in short["weak_or_partial_annotation_disagreements"]}
    current_discordant={r["source_transcript_id"] for r in proteins["records"]
                        if r["blast_annotation_and_reference_similarity_concordant"] is False}
    if current_discordant!=original_discordant or len(original_discordant)!=5:
        raise ValueError("Source BLAST family disagreements changed from frozen audit")

    result={
        "analysis_version":"aza3_cirsium_empirical_sources_v3",
        "source_owner":"zuizui0223/aza3",
        "upstream_theory_owner":"zuizui0223/anto",
        "empirical_source_archive_frozen":True,
        "cirsium_public_2026_sra":public,
        "cirsium_public_2020_organ_contrasts":organs,
        "public_asteraceae_reference_preflight":preflight,
        "cirsium_public_2020_protein_homology_triage":proteins,
        "frozen_2020_strong_similarity_shortlist":{
            "n_transcripts":len(frozen),
            "distinct_assembly_groups":short["distinct_assembly_groups_in_shortlist"],
            "family_counts":short["family_counts"],
            "validation":"SOURCE_PROTEIN_LOCAL_ALIGNMENT_REPRODUCED_FUNCTION_UNVERIFIED",
            "file":"data/anto_source_handoff/molecular/cirsium_japonicum_2020_source_similarity_shortlist_v3.json"
        },
        "interpretation_gates":{
            "focal_2026_Ryukyu_sequences_extracted":False,
            "2026_flower_stage_expression_assayed":False,
            "2020_trinucleotide_or_gene_identity_assigned_by_gene_tree":False,
            "documented_ancestor_white_to_coloured_transition":False,
            "retained_pathway_functionally_demonstrated_in_white_Ryukyu_florets":False,
            "genomic_independence_of_capitulum_LEGO_modules":False,
            "fitness_adaptation_demonstrated":False
        },
        "scientific_result":"Published Cirsium 2020 full transcript sequences can be screened for Saussurea DFR/ANS-like protein similarity. Gene-family annotation ≠ independently functional inherited pigment toolkit or allele history.",
    }
    text=json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n"
    if output is not None:
        Path(output).parent.mkdir(parents=True,exist_ok=True)
        Path(output).write_text(text,encoding="utf-8")
    else:
        print(text)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path)
    args=p.parse_args()
    main(args.out)
