"""Positive-reference leave-taxon-out sensitivity for Cirsium 2020 DFR/ANS-like transcripts.

Investigates whether the apparent *Cirsium japonicum* family similarity
disappears when the closest related positive reference (Saussurea/Cardueae)
is omitted. This is still homology sensitivity, NOT phylogenetic orthology
certification or a functional assay in Cirsium brevicaule/irumtiense.

All empirical code and sequence data are in aza3, while anto retains theory.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis/anto_empirical"))
from cirsium_source_homology_v3 import read_fasta
from cirsium_paralog_specificity_v4 import load_panel,analyze_transcript

DATA=ROOT/"data/anto_source_handoff/molecular"
SCENARIOS={
 "all_seven":[],
 "leave_out_Saussurea_medusa":["Saussurea medusa"],
 "leave_out_Gerbera_hybrida":["Gerbera hybrida"],
 "leave_out_Arabidopsis_positive_only":["BAD95233.1","NP_194019.1"],
}


def main(output=None):
    raw=json.loads((DATA/"cirsium_japonicum_2020_source_similarity_shortlist_v3.json").read_text(encoding="utf-8"))
    selected=raw.get("records",[])
    if len(selected)!=5 or {r["id"].rsplit("_i",1)[0] for r in selected}.__len__()!=4:
        raise ValueError("Expected frozen 5 transcripts in four source assembly groups")
    panel=load_panel(
        DATA/"extended_DFR_ANS_ANR_FLS_reference_proteins_v4.fasta",
        DATA/"extended_reference_protein_provenance_v4.json"
    )
    sources=read_fasta(DATA/"cirsium_japonicum_2020_26_source_sequences_v6.fasta")
    if len(sources)!=26:
        raise ValueError("Only frozen original published 26 transcripts admissible")

    evidence=[]
    for r in selected:
        key=r["id"]
        if key not in sources:raise ValueError("Missing publisher source transcript: "+key)
        baseline_family=r["best_family"]+"_reference"
        comparisons={}
        for name, omitted in SCENARIOS.items():
            subset=[entry for entry in panel
                    if entry["taxon"] not in omitted and entry["accession"] not in omitted]
            if not any(z["group"]==baseline_family for z in subset):
                raise ValueError("Sensitivity removes all references from focal family")
            p=analyze_transcript(sources[key],subset)
            comparisons[name]={
                "reference_count":len(subset),
                "best_group":p.get("best_family"),
                "strong_similarity":p.get("strong_reference_similarity",False),
                "match_to_frozen_family":p.get("best_family")==baseline_family,
                "strong_and_same_family":p.get("strong_reference_similarity",False)
                                      and p.get("best_family")==baseline_family,
                "aligned_aa":p.get("best_hit",{}).get("aligned_aa"),
                "reference_coverage":p.get("best_hit",{}).get("reference_coverage"),
                "aligned_identity":p.get("best_hit",{}).get("aligned_identity"),
                "best_reference_id":p.get("best_hit",{}).get("subject_accession"),
                "normalized_positive_minus_related_paralog_score":
                    p.get("best_positive_vs_related_negative_normalized_gap"),
            }
        evidence.append({
            "transcript_id":key,
            "source_assembly_group":key.rsplit("_i",1)[0],
            "initial_2020_source_family":r["best_family"],
            "comparisons":comparisons,
            "stable_across_all_reference_omissions":
                all(x["strong_and_same_family"] for x in comparisons.values()),
            "gene_function_verified":False,
        })
    result={
        "version":"aza3_cirsium_2020_reference_taxon_jackknife_v4",
        "source_study":"Park et al. 2020 DOI 10.1021/acsomega.0c04001",
        "number_of_screened_public_transcripts":len(evidence),
        "reference_scenarios":list(SCENARIOS),
        "stable_family_and_strong_in_all_scenarios":
            sum(r["stable_across_all_reference_omissions"] for r in evidence),
        "strong_when_Saussurea_removed":sum(
            r["comparisons"]["leave_out_Saussurea_medusa"]["strong_and_same_family"]
            for r in evidence),
        "records":evidence,
        "claim_limits":{
            "these_are_2026_Ryukyu_species_sequences":False,
            "orthology_certified":False,
            "enzyme_function_demonstrated":False,
            "floral_pigment_reexpression_demonstrated":False,
            "no_plant_histories_dated":True
        },
    }
    frozen=json.loads((DATA/"cirsium_2020_reference_omission_sensitivity_v4.json").read_text(encoding="utf-8"))
    if frozen.get("version")!="aza3_cirsium_reference_taxon_sensitivity_v4":
        raise ValueError("Frozen external-reference robustness report has unexpected version")
    if (result["stable_family_and_strong_in_all_scenarios"]!=frozen["stable_family_and_strong_in_all_scenarios"] or
        result["strong_when_Saussurea_removed"]!=frozen["strong_when_Saussurea_removed"]):
        raise ValueError("Strong-family reference omission sensitivity changed")
    frozen_by_id={row["transcript_id"]:row for row in frozen["scenarios_by_transcript"]}
    if set(frozen_by_id)!={row["transcript_id"] for row in evidence}:
        raise ValueError("Frozen original five transcribed IDs changed")
    for row in evidence:
        ref=frozen_by_id[row["transcript_id"]]
        if ref["stable_across_all_reference_omissions"]!=row["stable_across_all_reference_omissions"]:
            raise ValueError("Reference stability classification changed for "+row["transcript_id"])
        for case,meas in row["comparisons"].items():
            prior=ref["comparisons"][case]
            for field in ("best_group","strong_similarity","strong_and_same_family",
                          "best_reference_id","aligned_aa"):
                if meas[field]!=prior[field]:
                    raise ValueError(f"Frozen alignment reference sensitivity changed: {row['transcript_id']} {case} {field}")
            for field in ("reference_coverage","aligned_identity"):
                if abs(meas[field]-prior[field])>1e-6:
                    raise ValueError(f"Frozen protein alignment evidence changed: {row['transcript_id']} {case} {field}")
    result["frozen_published_source_reproducibility_verified"]=True
    result["frozen_record"]="data/anto_source_handoff/molecular/cirsium_2020_reference_omission_sensitivity_v4.json"
    text=json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n"
    if output is None:print(text)
    else:
        output=Path(output)
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(text,encoding="utf-8")
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path)
    args=p.parse_args()
    main(args.out)
