"""Validate source-derived NJ vs FastTree result without promoting topology to enzyme function.

The frozen comparison is from real Park (2020) *Cirsium japonicum*
transcripts. It does not include the 2026 Ryukyu *Cirsium* samples.
A "BETWEEN_FAMILY_ANCHORS" outcome is UNRESOLVED; method concordance
does NOT make the placement a confident protein-function discovery.

Pass optional --analysis to validate fresh FastTree results against
the frozen output with exactly the same independent source records.
"""
import argparse
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DEFAULT=ROOT/"data/anto_source_handoff/molecular/cirsium_2020_nj_fasttree_method_comparison_v5.json"
EXPECTED={
    ("DFR_ANR_primary","TR35526|c0_g1_i1"),
    ("DFR_ANR_primary","TR38875|c0_g4_i5"),
    ("DFR_ANR_primary","TR39766|c0_g1_i2"),
    ("DFR_ANR_alt_isoform","TR35526|c0_g1_i1"),
    ("DFR_ANR_alt_isoform","TR38875|c0_g4_i7"),
    ("DFR_ANR_alt_isoform","TR39766|c0_g1_i2"),
    ("ANS_FLS","TR58302|c0_g1_i1"),
}
UNRESOLVED="BETWEEN_FAMILY_ANCHORS"


def audit_frozen(path=DEFAULT):
    data=json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("version")!="aza3_cirsium_NJ_FastTree_family_concordance_v5":
        raise ValueError("Unrecognized gene family source audit")
    rows=data.get("records",[])
    if len(rows)!=7:
        raise ValueError("Expected seven, not seven independent biological genes")
    ids={(r["case"],r["source_transcript"]) for r in rows}
    if ids!=EXPECTED:
        raise ValueError("Frozen transcript-case identities changed, duplicate or missing")
    if data["n_tree_configurations"]!=3 or data["n_query_configurations"]!=7:
        raise ValueError("Frozen analysis count changed")
    if data["distinct_source_transcripts"]!=5 or data["n_same_family_placement_label"]!=7:
        raise ValueError("Biological source or method result count incorrect")
    if data["n_candidates_confidently_identified_as_functional_dfr_or_ans"]!=0:
        raise ValueError("Unsupported functional enzyme claim")
    if not (data["no_orthology_certified"] and data["no_gene_family_function_certified"] and
            data["no_ryukyu_regain_evidence"] and data["no_source_gene_evolution_date_inferred"]):
        raise ValueError("Claim ceiling silently removed")
    for r in rows:
        if (r["NJ_unrooted_reference_family_placement"] != UNRESOLVED or
                r["FastTree_WAG_unrooted_reference_family_placement"] != UNRESOLVED):
            raise ValueError("Family placement label changed; inspect biological source")
        if (r["classification_stable_across_inference_method"] is not True or
                r["functional_orthology_verified"] is not False):
            raise ValueError("Inconsistent inference or invalid orthology call")
    return data


def compare_analysis(frozen,analysis):
    if analysis.get("analysis_version")!="aza3_cirsium_source_family_method_comparison_v5":
        raise ValueError("Unexpected FastTree output version")
    if (analysis.get("n_inference_cases")!=3 or
            analysis.get("n_query_case_comparisons")!=7 or
            analysis.get("n_identical_NJ_vs_approxML_placements")!=7):
        raise ValueError("New phylogeny disagrees with frozen source report")
    if analysis["allele_orthology_certified"] or analysis["enzyme_function_verified"] or analysis["Ryukyu_white_to_pigmented_regain_demonstrated"]:
        raise ValueError("Phylogenetic support was promoted beyond available data")
    previous={(r["case"],r["source_transcript"]):r for r in frozen["records"]}
    actual={}
    for case,data in analysis["cases"].items():
        for transcript,value in data["queries"].items():
            actual[(case,transcript)]=value
    if set(previous)!=set(actual):
        raise ValueError("Incomplete or swapped source case and transcript ID")
    for key,prior in previous.items():
        current=actual[key]
        for value in ("NJ_unrooted_reference_family_placement",
                      "FastTree_WAG_unrooted_reference_family_placement",
                      "classification_stable_across_inference_method",
                      "functional_orthology_verified"):
            if current[value]!=prior[value]:
                raise ValueError(f"Family inference no longer matches frozen source result: {key} {value}")
    return True


def main(frozen_path=DEFAULT,analysis_path=None):
    frozen=audit_frozen(frozen_path)
    replayed=False
    if analysis_path:
        actual=json.loads(Path(analysis_path).read_text(encoding="utf-8"))
        compare_analysis(frozen,actual)
        replayed=True
    return {
        "version":frozen["version"],
        "n_source_transcripts":frozen["distinct_source_transcripts"],
        "query_evaluations":len(frozen["records"]),
        "NJ_and_FastTree_method_concordance_labels":frozen["n_same_family_placement_label"],
        "confident_DFR_ANS_family_labels":0,
        "read_from_actual_published_sequence":True,
        "biochemical_activity_assayed":False,
        "Ryukyu_coloured_flower_regain_tested":False,
        "rerun_against_frozen_report":replayed,
        "interpretation":"Method robustness of unresolved family classification, not orthology or enzyme confirmation."
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--report",type=Path,default=DEFAULT)
    parser.add_argument("--analysis",type=Path,default=None)
    args=parser.parse_args()
    print(json.dumps(main(args.report,args.analysis),ensure_ascii=False,indent=2,sort_keys=True))
