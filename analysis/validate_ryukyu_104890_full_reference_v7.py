"""Fail-closed intake of REAL six-Ryukyu 8MiB R1 whole-source competition.

Factual result: 104,890 verified 2020 Cirsium japonicum transcripts
competed with 26 selected source homology targets for each of six public
2026 Ryukyu Cirsium YOUNG-LEAF R1 sequence prefixes. This supports
cross-congeneric fragment homology only, not coding gene orthology,
protein catalysis, floret pigmentation, allele age or evolutionary regain.
"""
from collections import Counter
from pathlib import Path
import argparse
import json

ROOT=Path(__file__).resolve().parents[1]
SOURCES=ROOT/"data/anto_source_handoff"
PRIMARY=SOURCES/"public_sra/ryukyu_full_104890_reference_competition_8MiB_v7.json"
REGISTER=SOURCES/"public_sra/ryukyu_six_R1_prefix_source_v6.json"
OLD16=SOURCES/"public_sra/ryukyu_six_source_leaf_R1_mapping_16MiB_v6.json"
OLD8=SOURCES/"public_sra/ryukyu_six_source_leaf_R1_mapping_8MiB_alternate_v6.json"
MOLEC=SOURCES/"molecular/cirsium_japonicum_2020_26_sequence_provenance_v6.json"

EXPECTED_RUNS={
    "SRR35152724","SRR35152725","SRR35152729",
    "SRR35152730","SRR35152731","SRR35152732",
}
KEYS={
    "TR39766|c0_g1_i2":"TR39766_cirsium_japonicum_like_strict_fragments",
    "TR41106|c0_g1_i5":"TR41106_cirsium_japonicum_like_strict_fragments",
}


def _json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def validate(primary=PRIMARY,source=REGISTER,old16=OLD16,old8=OLD8,molecular=MOLEC):
    d=_json(primary)
    src=_json(source)
    original=_json(molecular)
    prior16=_json(old16)
    prior8=_json(old8)

    if (d.get("version")!="aza3_Ryukyu_whole_Cirsium_japonicum_reference_competition_v7"
        or d.get("source_project")!="PRJNA1311153"
        or d.get("reference_transcripts")!=104890):
        raise ValueError("Whole-source empirical result version/provenance changed")
    if (original["decompressed_fasta_sha256"]!=d["background_reference_member_sha256"]
        or original["source_fasta_entry_count"]!=104890
        or original["contract"]!="anto_park2020_26_source_transcripts_v6"
        or not d["selected_26_prior_blast_candidates_appear_in_full_original_source"]):
        raise ValueError("Source-author 104890 unigene reference checksum mismatch")
    if (src["version"]!="aza3_cirsium_Ryukyu_public_R1_prefix_manifest_v6"
        or len(src["runs"])!=6):
        raise ValueError("Original public Ryukyu source manifest incomplete")
    provider={r["run_accession"]:r for r in src["runs"]}
    if set(provider)!=EXPECTED_RUNS:
        raise ValueError("Unexpected run ID or duplicate BioSample registration")
    if d["read_prefix_compressed_bytes_per_run"]!=8388608 or not d["full_compressed_fastqs_not_downloaded"]:
        raise ValueError("Capped source bytes or full library status changed")
    if (not d["original_library_R1_full_MD5_not_verified"]
        or d["original_source_read_pair_mates_downloaded"]!="R1 prefixes only"
        or not d["assay"]["additional_104864_original_author_transcripts_compete_for_each_read"]):
        raise ValueError("Unverified sequencing integrity or incomplete background")
    if (d["assay"]["min_MAPQ"],d["assay"]["min_aligned_nt"],
        d["assay"]["max_NM_over_aligned_nt"])!=(20,90,.2):
        raise ValueError("Mapping QC rules changed")
    if d["assay"]["source_2020_candidates_functionally_assayed_for_DFR_ANS"]:
        raise ValueError("Molecular function was never experimentally measured")
    if d["assay"]["differential_expression_analysis_permitted"]:
        raise ValueError("An ordered R1 prefix cannot support expression comparisons")
    for statement in ("Not whole-locus functional DFR/ANS orthology or intact coding gene",
                      "Not a flower anthocyanin measurement for any plant"):
        if statement not in d["interpretation"]["unsupported"]:
            raise ValueError("Scientific claim ceiling was silently relaxed")

    rows=d.get("samples",[])
    if len(rows)!=6 or {r["run_accession"] for r in rows}!=EXPECTED_RUNS:
        raise ValueError("Expected exactly six public biological library records")
    if {r["biosample_accession"] for r in rows}!={x["biosample_accession"] for x in src["runs"]}:
        raise ValueError("BioSample cardinality error")
    observed=Counter()
    total_reads=0
    total_strict=0
    n_positive=Counter()
    for r in rows:
        run=r["run_accession"]
        parent=provider[run]
        if (r["biosample_accession"]!=parent["biosample_accession"]
            or r["specimen_code"]!=parent["source_specimen_code"]
            or r["taxon_as_source_record"]!=parent["taxon_as_experiment_titled"]
            or r["source_tissue"]!="young leaf"):
            raise ValueError("Run accession, taxon, original specimen or tissue wrongly joined")
        if (r["compressed_first_R1_bytes_observed"]!=8388608
            or r["original_provider_R1_full_MD5_not_reevaluated"] is not True
            or r["strict_composition_shuffled_decoy_fragment_matches"]!=0):
            raise ValueError("Source byte and computational-negative control mismatch")
        n=r["complete_valid_R1_records_analyzed"]
        nc=r["source_2020_26_candidate_strict_unique_fragments"]
        if not(type(n) is int and n>10000 and type(nc) is int and nc>=0):
            raise ValueError("Malformed measured read or candidate record count")
        if nc>n:raise ValueError("Mapped reads exceed complete FASTQ records")
        total_reads+=n
        total_strict+=nc
        observed[r["taxon_as_source_record"]]+=1
        for gene,key in KEYS.items():
            v=r[key]
            if not isinstance(v,int) or v<1 or v>nc:
                raise ValueError(f"Source transcript {gene} read count invalid")
            n_positive[gene]+=v
    if observed!={"Cirsium brevicaule":3,"Cirsium irumtiense":3}:
        raise ValueError("Missing or renamed original submitted taxon source")
    if (total_reads!=d["reads_observed_sum"] or total_reads!=747234 or
        total_strict!=d["total_source_26_strict_fragment_hits"] or total_strict!=215):
        raise ValueError("Public full-source result totals do not reconcile")
    if d["strict_26_candidate_composition_shuffled_decoy_matches"]!=0:
        raise ValueError("Frozen decoy background changed")
    for id,count in (("TR39766|c0_g1_i2",66),("TR41106|c0_g1_i5",106)):
        if n_positive[id]!=count or d["strict_source_fragment_by_pub2020_candidate"][id]!=count:
            raise ValueError("Source candidate whole-reference evidence changed")
    if (d["TR39766_source_like_fragments_detected_in_n_of_6"]!=6
        or d["TR41106_source_like_fragments_detected_in_n_of_6"]!=6):
        raise ValueError("All-six leaf-source candidate occurrence invalid")

    # Earlier 8 and 16 MiB source studies used DIFFERENT reference/sample filters,
    # but both contain the two same candidate IDs in all six source runs.
    if (prior16["n_distinct_source_transcripts_strict_mapped"]!=18
        or prior16["n_distinct_biological_library_sources"]!=6
        or prior8["n_source_runs"]!=6):
        raise ValueError("Historic restricted-reference audit missing")
    occurrence16=prior16["candidate_occurrence_n_of_six"]
    occurrence8=prior8["source_gene_occupancy"]
    for candidate in KEYS:
        if (occurrence16.get(candidate)!=6
            or occurrence8.get(candidate,{}).get("samples")!=6):
            raise ValueError("Original limited-reference source candidate occurrence no longer reproducible")

    return {
        "analysis_version":"aza3_Ryukyu_104890_reference_competition_evidence_gate_v7",
        "full_source_2020_unigene_background_verified":104890,
        "six_public_young_leaf_R1_specimens":6,
        "white_taxon_libraries":3,
        "coloured_taxon_libraries":3,
        "read_prefix_only_R1_bytes_per_library":8388608,
        "full_source_competition_validated_complete_reads":total_reads,
        "strict_26_source_candidate_reads_under_full_background":total_strict,
        "fragment_occurrences_in_all_six_sources":{
            "TR39766|c0_g1_i2":6,
            "TR41106|c0_g1_i5":6,
        },
        "candidate_counts_under_full_reference":{
            key:n_positive[key] for key in KEYS
        },
        "source_8MiB_16MiB_limited_reference_presence_agree":True,
        "source_prefixes_are_nested_not_independent_replication":True,
        "filter_thresholds_differ_so_counts_not_expression_fold_changes":True,
        "confirmed_functional_DFR_ANS_orthology":False,
        "confirmed_floret_anthocyanin_expression":False,
        "confirmed_ancestral_phenotype_reexpression":False,
        "claim_ceiling":"Public source-derived 2026 leaf RNA fragments homologous to preselected 2020 thistle sequences; none is a verified catalytic or developmental floral state.",
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-result",type=Path,default=PRIMARY)
    args=p.parse_args()
    print(json.dumps(validate(primary=args.source_result),ensure_ascii=False,indent=2,sort_keys=True))
