"""Fail-closed verification of actual public Ryukyu six-leaf R1 prefix findings.

The two stored result files come from GitHub Actions run 39228212947
(16 MiB, established source script) and run 39228426774 (8 MiB,
alternative read-mapping filter). The 8 MiB reads are nested within
the same six real libraries and are NOT independent biological samples.

Allowed biological claim: source-read fragments resembling selected 2020
Cirsium japonicum flavonoid-related transcript sequences occur in *young
leaves* of both Ryukyu taxa. No full gene, anthocyanin enzyme activity,
floral expression, trait polarity, resumption of expression or allele origin
has been established.
"""
import csv
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/anto_source_handoff"
MANIFEST=DATA/"public_sra/ryukyu_six_R1_prefix_source_v6.json"
PRIMARY=DATA/"public_sra/ryukyu_six_source_leaf_R1_mapping_16MiB_v6.json"
SECONDARY=DATA/"public_sra/ryukyu_six_source_leaf_R1_mapping_8MiB_alternate_v6.json"
REF_PROVENANCE=DATA/"molecular/cirsium_japonicum_2020_26_sequence_provenance_v6.json"

REQUIRED_IDS={"TR39766|c0_g1_i2","TR41106|c0_g1_i5"}


def audit(primary=PRIMARY,secondary=SECONDARY,manifest=MANIFEST,
          reference_provenance=REF_PROVENANCE):
    a=json.loads(Path(primary).read_text(encoding="utf-8"))
    b=json.loads(Path(secondary).read_text(encoding="utf-8"))
    source=json.loads(Path(manifest).read_text(encoding="utf-8"))
    ref=json.loads(Path(reference_provenance).read_text(encoding="utf-8"))
    if (a.get("version")!="aza3_Ryukyu_six_leaf_R1_fragment_screen_v6" or
        b.get("version")!="aza3_Ryukyu_six_leaf_R1_fragment_8MiB_sensitivity_v6" or
        source.get("version")!="aza3_cirsium_Ryukyu_public_R1_prefix_manifest_v6"):
        raise ValueError("Wrong source dataset/analysis contract")
    if ref.get("contract")!="anto_park2020_26_source_transcripts_v6":
        raise ValueError("Wrong frozen 2020 reference transcript origin")
    accepted_ref={r["id"] for r in ref.get("records",[])}
    if len(accepted_ref)!=26:
        raise ValueError("Expected 26 source-verified 2020 transcripts, not a genome-wide reference")
    registry={r["run_accession"]:r for r in source["runs"]}
    if len(registry)!=6 or len(a.get("per_run",[]))!=6 or len(b.get("per_run",[]))!=6:
        raise ValueError("Must cover six original biological young-leaf libraries twice")
    nread=nstrict=nprimary=0
    source_by_gene=Counter()
    occurrence=Counter()
    for r in a["per_run"]:
        id=r["run_accession"]
        if id not in registry or r["biosample_accession"]!=registry[id]["biosample_accession"]:
            raise ValueError("BioSample/source accession mismatch")
        if (r["taxon_title"]!=registry[id]["taxon_as_experiment_titled"] or
            r["specimen_code"]!=registry[id]["source_specimen_code"] or
            r["tissue"]!="young leaf" or
            r["strict_shuffled_decoy_fragment_count"]!=0):
            raise ValueError("Sample identity/tissue/decoy result changed")
        if r["compressed_R1_prefix_bytes"]!=16777216:
            raise ValueError("16MiB planned primary prefix changed")
        reads=r["read1_records"]
        mapped=r["strict_source_fragment_count"]
        if not isinstance(reads,int) or not isinstance(mapped,int) or reads<1000 or mapped<0:
            raise ValueError("Impossible primary read counts")
        hits=r["per_published_2020_Cirsium_source_transcript_strict_counts"]
        if not set(hits)<=accepted_ref or any(not isinstance(n,int) or n<=0 for n in hits.values()):
            raise ValueError("A candidate mapping is not an accepted public source transcript")
        if sum(hits.values())!=mapped:
            raise ValueError("Per-sample source counts not conserved")
        nread+=reads
        nstrict+=mapped
        nprimary+=r["primary_raw_alignment_rows"]
        for k,v in hits.items():
            source_by_gene[k]+=v
            occurrence[k]+=1
    if len({x["run_accession"] for x in a["per_run"]})!=6:
        raise ValueError("Repeated same biological source accession")
    cum=a["n_observed_R1_records"],a["n_strict_source_like_fragments"],a["n_raw_primary_reference_hits"]
    if (nread,nstrict,nprimary)!=cum or cum!=(1499192,636,1144):
        raise ValueError("Observed real-data totals changed without review")
    if (a["n_strict_composition_shuffled_decoy_fragments"]!=0 or
        a["total_compressed_download_bytes"]!=6*16777216):
        raise ValueError("Original source budget/decoy result changed")
    if source_by_gene!=Counter(a["candidate_aggregate_strict_fragments"]):
        raise ValueError("Aggregated candidate totals inconsistent with original source run observations")
    if occurrence!=Counter(a["candidate_occurrence_n_of_six"]):
        raise ValueError("Candidate recurrence inconsistent with real samples")
    if any(occurrence[k]!=6 for k in REQUIRED_IDS):
        raise ValueError("The previously verified all-six source candidates changed")
    if source_by_gene["TR39766|c0_g1_i2"]!=161 or source_by_gene["TR41106|c0_g1_i5"]!=270:
        raise ValueError("Source-backed candidate mapped-fragment totals changed")

    q8=Counter()
    rd8=nd8=0
    for r in b["per_run"]:
        if r["run_accession"] not in registry or r["retrieval_status"]!="PREFIX_ANALYZED":
            raise ValueError("Secondary source run missing or not really retrieved")
        if r["specimen"]!=registry[r["run_accession"]]["source_specimen_code"]:
            raise ValueError("Secondary sample accession relabeled")
        if (r["strict_decoy_fragments"]!=0 or
            not isinstance(r["strict_source_fragments"],int) or r["strict_source_fragments"]<0):
            raise ValueError("Secondary mapping/negative control changed")
        if sum(r["per_source"].values())!=r["strict_source_fragments"]:
            raise ValueError("Secondary source target count not conserved")
        if not set(r["per_source"])<=accepted_ref:
            raise ValueError("Unexpected target in secondary method")
        rd8+=r["complete_R1_prefix_reads"]
        nd8+=r["strict_source_fragments"]
        for key,n in r["per_source"].items():
            if n>0:q8[key]+=1
    if {r["run_accession"] for r in b["per_run"]}!=set(registry):
        raise ValueError("Alternative 8MiB references different biological source material")
    if (rd8,nd8)!=(747995,378) or b["total_prefix_gzip_bytes"]!=6*8388608:
        raise ValueError("Alternative method compressed budget or read count changed")
    if any(q8[k]!=6 for k in REQUIRED_IDS):
        raise ValueError("All-six target recurrence not stable under technical prefix/mapping sensitivity")
    for field in ("functional_orthology_assayed","floral_stage_phenotype_or_RNA_assayed",
                  "flower_colour_reexpression_demonstrated",
                  "any_species_differential_expression_test_valid",
                  "genotype_phenotype_ancestral_polarity_test_completed",
                  "genome_wide_reference_competition_completed"):
        if a["biological_claim_ceiling"][field] is not False:
            raise ValueError("Biological inference exceeds original read-level evidence: "+field)
    return {
        "analysis_version":"aza3_ryukyu_real_R1_prefix_evidence_audit_v6",
        "n_independent_source_libraries":6,
        "primary_16MiB_total_R1_reads":nread,
        "primary_16MiB_stringent_2020_source_like_fragments":nstrict,
        "primary_16MiB_shuffled_decoy_fragments":0,
        "primary_16MiB_distinct_source_transcript_IDs_with_matches":len(source_by_gene),
        "supporting_8MiB_total_R1_reads":rd8,
        "supporting_8MiB_stringent_source_fragments":nd8,
        "supporting_8MiB_shuffled_decoy_fragments":0,
        "same_biological_libraries_not_independent_experiments":True,
        "TR39766_DFR_related_source_in_six_leaf_libraries":True,
        "TR41106_ANS_Blast_annotated_source_in_six_leaf_libraries":True,
        "TR39766_16MiB_fragments":source_by_gene["TR39766|c0_g1_i2"],
        "TR41106_16MiB_fragments":source_by_gene["TR41106|c0_g1_i5"],
        "functional_orthology_assayed":False,
        "floral_expression_assayed":False,
        "genetic_evolutionary_reacquisition_demonstrated":False,
        "between_species_differential_expression_permitted":False,
        "source_reference_panel_not_genome_wide":True,
        "interpretation":"Qualitative source RNA candidate-fragment detection in both taxa; no inference of floral enzyme functionality or historical regained flower phenotype."
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True,ensure_ascii=False))
