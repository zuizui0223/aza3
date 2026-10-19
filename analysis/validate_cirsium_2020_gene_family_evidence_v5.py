"""Fail-closed audit of genuine 2020 Cirsium gene-family phylogeny results.

Source manuscript: Park et al. 2020, not the 2026 Ryukyu pair.
Validates numbers from independent GitHub Actions and, if --analysis provided,
compares those against new MAFFT/NJ/bootstrap source calculations.

An unresolved family split does not mean a unique intermediate gene family;
a candidate residue substitution does NOT prove loss of enzyme activity.
"""
import argparse
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
LEDGER=ROOT/"data/anto_source_handoff/molecular/cirsium_japonicum_2020_gene_family_topology_evidence_v5.json"

FROZEN_IDS={
    ("DFR_ANR_primary","TR35526|c0_g1_i1"),
    ("DFR_ANR_primary","TR38875|c0_g4_i5"),
    ("DFR_ANR_primary","TR39766|c0_g1_i2"),
    ("DFR_ANR_alt_isoform","TR35526|c0_g1_i1"),
    ("DFR_ANR_alt_isoform","TR38875|c0_g4_i7"),
    ("DFR_ANR_alt_isoform","TR39766|c0_g1_i2"),
    ("ANS_FLS","TR58302|c0_g1_i1"),
}


def check_source_record(path=LEDGER):
    data=json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("version")!="aza3_cirsium_2020_gene_family_topology_report_v5":
        raise ValueError("Unexpected source version")
    records=data.get("reported_results",[])
    if len(records)!=7 or {(r["case"],r["source_transcript"]) for r in records}!=FROZEN_IDS:
        raise ValueError("Incomplete or duplicated seven observed gene-family evaluations")
    if data["number_of_distinct_query_transcripts"]!=5 or data["number_of_tree_configurations"]!=3:
        raise ValueError("Incorrect biological/analysis replication grain")
    if data["source_claim_ceiling"]["confirmed_ancestral_white_to_coloured_reexpression"]:
        raise ValueError("Published gene tree is not a floral colour reversal experiment")
    for r in records:
        t=r["bootstrap_tree_counts"]
        if r["valid_site_bootstraps"]!=200 or r["invalid_site_bootstraps"]!=0 or sum(t.values())!=200:
            raise ValueError("Bootstrap counts or invalid site resamples not conserved")
        if t!={"BETWEEN_FAMILY_ANCHORS":200} or r["observed_unrooted_NJ_reference_split"]!="BETWEEN_FAMILY_ANCHORS":
            raise ValueError("Source topology changed: re-review original family assignment")
        lo,hi=r["bootstrap_95pct_interval"]
        if not lo<=0<=hi or not r["interval_contains_zero"]:
            raise ValueError("Reference identity CI no longer contains zero; do not silently change inference")
        if r["site_replicates"]!=1000:
            raise ValueError("Identity uncertainty sample size mismatch")
        m=r["biochemical_motif_evidence"]
        if m["covered"]!=m["matching"]+m["substituted"]:
            raise ValueError("Motif substitution/match accounting inconsistent")
        if any(s["evidence_state"]=="NO_ALIGNMENT_COVERAGE" and s["source_aligned_residue"] is not None for s in m["sites"]):
            raise ValueError("Missing aligned site improperly called a residue")
        if r["enzyme_function_certified"]:
            raise ValueError("Cannot upgrade author assembly to enzymatic function by tree topology")
    return data


def compare_tree_analysis(ledger,analysis):
    """New tree run must match original NJ source result, independent of FastTree."""
    if analysis.get("analysis_version")!="aza3_cirsium_2020_protein_family_phylogeny_v5":
        raise ValueError("Wrong tree run version")
    rows={}
    for name,group in analysis.get("trees",{}).items():
        for taxon,entry in group["source_translation"].items():
            source=entry["source_id"]
            row={
                "topology":group["queries"][taxon],
                "identity":group["relative_identity_site_bootstrap"][taxon],
                "motif":group["reference_site_alignment_audit"][taxon]
            }
            rows[(name,source)]=row
    if set(rows)!=FROZEN_IDS:
        raise ValueError("Replayed family trees lost or duplicated expected source candidates")
    for r in ledger["reported_results"]:
        x=rows[(r["case"],r["source_transcript"])]
        tree=x["topology"]
        if (tree["observed_unrooted_NJ_family_split"]!=r["observed_unrooted_NJ_reference_split"]
            or tree["bootstrap_counts"]!=r["bootstrap_tree_counts"]
            or tree["bootstrap_valid_replicates"]!=r["valid_site_bootstraps"]):
            raise ValueError("Frozen NJ gene-family placement is not reproduced")
        identity=x["identity"]
        if abs(identity["observed_positive_minus_related_identity_margin"]
                -r["observed_positive_minus_related_identity_margin"])>1e-9:
            raise ValueError("Frozen relative amino acid identity margin differs")
        if any(abs(a-b)>1e-9 for a,b in zip(
                identity["site_bootstrap_identity_margin_95pct_interval"],
                r["bootstrap_95pct_interval"])):
            raise ValueError("Source identity bootstrap interval differs")
        if identity["reference_identity_diff_interval_includes_zero"]!=r["interval_contains_zero"]:
            raise ValueError("Source identity interval classification changed")
        motif=x["motif"]
        if (motif["n_sites_covered"]!=r["biochemical_motif_evidence"]["covered"]
            or motif["n_matching_residues"]!=r["biochemical_motif_evidence"]["matching"]
            or motif["n_aligned_substitutions"]!=r["biochemical_motif_evidence"]["substituted"]):
            raise ValueError("Source-motif reference mapping changed")
    return True


def main(report=LEDGER,analysis_path=None):
    data=check_source_record(report)
    if analysis_path is not None:
        actual=json.loads(Path(analysis_path).read_text(encoding="utf-8"))
        compare_tree_analysis(data,actual)
    q=data["reported_results"]
    return {
        "version":data["version"],
        "source_evaluations":len(q),
        "distinct_source_transcripts":len({r["source_transcript"] for r in q}),
        "complete_NJ_family_classifications_supported":0,
        "noisy_relative_aa_identity_intervals_contain_zero":sum(r["interval_contains_zero"] for r in q),
        "parent_motif_TR39766_DFR_covered":next(r for r in q if r["case"]=="DFR_ANR_primary" and
                                               r["source_transcript"]=="TR39766|c0_g1_i2")["biochemical_motif_evidence"],
        "source_replay_tested":analysis_path is not None,
        "gene_function_verified":False,
        "floral_phenotype_regain_demonstrated":False,
        "taxon":"Cirsium japonicum Park 2020, not Ryukyu Cirsium 2026",
    }


if __name__=="__main__":
    args=argparse.ArgumentParser()
    args.add_argument("--report",type=Path,default=LEDGER)
    args.add_argument("--analysis",type=Path)
    opts=args.parse_args()
    print(json.dumps(main(opts.report,opts.analysis),ensure_ascii=False,indent=2))
