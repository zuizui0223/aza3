"""Source-locked evidence audit for prospective lost-flower-trait competence work.

Only already-published repository-backed Cirsium records are admitted.
No published leaf RNA run, taxon label, morphological code or gene-family
placement is promoted to flower-specific functional competence, genotype
restoration, ancestor polarity, viable seed fitness or a new field sample.

Stdlib-only; fails closed on a drifted canonical source Git blob.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = {
    "manifest": (
        "data/anto_source_handoff/public_sra/prjna1311153_run_manifest_20261017.csv",
        "3e044b4e3d4d1c3a278bd79872f84c1e9e0722bd",
    ),
    "qc": (
        "data/anto_source_handoff/public_sra/prjna1311153_qc_20261017.json",
        "137db65e483ee252d3ba58da4e7561db77aaaa10",
    ),
    "takaoense": (
        "data/anto_source_handoff/public_sra/takaoense_aza3_morph_crosswalk_20261017.csv",
        "1f92769dfa4876f5d4b080ce2faccc7a0081eed4",
    ),
    "family": (
        "data/anto_source_handoff/molecular/cirsium_2020_nj_fasttree_method_comparison_v5.json",
        "427abcbb94dd66f8248441c42c75ceef8a072744",
    ),
    "panel": (
        "data/planning/aza3_nature_genomic_sampling_panel_v1.csv",
        "e3fca68010629661e9e17e1915f09e9a28ab17",
    ),
}

UNMEASURED = "NOT_MEASURED_IN_REGISTERED_SOURCE"
UNKNOWN = "NOT_IDENTIFIED"
PLANNED = "PROSPECTIVE_ONLY"

def git_blob_sha(b: bytes) -> str:
    """Git's blob identity, not a SHA256 of just the payload."""
    return hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()

def load_source(name, root=ROOT):
    path, expected_sha = SOURCE[name]
    data = (Path(root)/path).read_bytes()
    observed_sha = git_blob_sha(data)
    if observed_sha != expected_sha:
        raise ValueError(f"{name}: source blob identity changed: {observed_sha}")
    return data

def csv_rows(name, root=ROOT):
    return list(csv.DictReader(io.StringIO(load_source(name,root).decode("utf-8-sig"))))

def json_data(name, root=ROOT):
    return json.loads(load_source(name,root))

def group_readiness(manifest, qc, morph, family, panel):
    if len(manifest)!=25 or qc.get("n_runs")!=25 or qc.get("n_biosamples")!=25:
        raise ValueError("The original source 25-run denominator changed")
    required=("run_accession","biosample_accession","specimen_code_from_experiment_title",
              "taxon_as_experiment_titled","tissue_from_project_description",
              "individual_flower_colour","library_strategy")
    for key in required:
        if any(not row.get(key) for row in manifest):
            raise ValueError(f"Missing required source field: {key}")
    ids={field:[r[field] for r in manifest] for field in ("run_accession","biosample_accession")}
    if any(len(vals)!=len(set(vals)) for vals in ids.values()):
        raise ValueError("Repeated source accession: do not pseudo-replicate")
    if any(r["tissue_from_project_description"].strip().lower()!="young leaf" for r in manifest):
        raise ValueError("Tissue not uniformly young leaf: requalify biological assays")
    if any(r["individual_flower_colour"] != "NOT_MEASURED_IN_THIS_DATASET" for r in manifest):
        raise ValueError("Individual-level floral phenotype status changed")
    if any(r["library_strategy"]!="RNA-Seq" for r in manifest):
        raise ValueError("Unexpected strategy in source manifest")
    counts=Counter(r["taxon_as_experiment_titled"] for r in manifest)
    if dict(counts)!=qc.get("counts_by_taxon") or sum(counts.values())!=25:
        raise ValueError("Taxon counts disagree with separately preserved QC")
    if counts.get("Cirsium brevicaule")!=3 or counts.get("Cirsium irumtiense")!=3:
        raise ValueError("Ryukyu six-plant source cohort drift")
    run_map={r["run_accession"]:r for r in manifest}
    if len(morph)!=6:
        raise ValueError("Six published takaoense library provenance entries expected")
    color_counts=Counter()
    for r in morph:
        if r["run_accession"] not in run_map:
            raise ValueError("Morph entry does not link to source library")
        source=run_map[r["run_accession"]]
        if (r["biosample_accession"]!=source["biosample_accession"]
            or r["specimen_code_from_experiment_title"]!=source["specimen_code_from_experiment_title"]
            or source["taxon_as_experiment_titled"]!="Cirsium japonicum var. takaoense"):
            raise ValueError("Crosswalk would join wrong biological specimen")
        if (r["laboratory_tissue"]!="young leaf" or
            r["floral_rna_measured"].lower()!="false" or
            r["floral_pigment_chemistry_measured"].lower()!="false" or
            r["independent_population_replication"].lower()!="false"):
            raise ValueError("Takaoense evidence status changed without requalification")
        color_counts[r["morph_code_from_aza3"]]+=1
    if dict(color_counts)!={"W":3,"BP":3}:
        raise ValueError("Morph count or code drift")

    records=family.get("records",[])
    if len(records)!=7 or len({r["source_transcript"] for r in records})!=5:
        raise ValueError("Family screen unit count drift")
    if (family.get("n_candidates_confidently_identified_as_functional_dfr_or_ans")!=0 or
        family.get("no_orthology_certified") is not True):
        raise ValueError("Molecular evidence promoted, freeze must be re-audited")
    if any(r.get("NJ_unrooted_reference_family_placement")!="BETWEEN_FAMILY_ANCHORS"
           or r.get("FastTree_WAG_unrooted_reference_family_placement")!="BETWEEN_FAMILY_ANCHORS"
           or r.get("functional_orthology_verified") is not False for r in records):
        raise ValueError("Different family placement / no longer unidentifiable")

    if any("recommended_n" not in r or "phase" not in r for r in panel):
        raise ValueError("Broken prospective plan")
    plan={r["taxon"]:r for r in panel if r["phase"]=="A"}
    if not all(r["recommended_n"].isdigit() for r in plan.values()):
        raise ValueError("Phase-A denominator is not numerical")
    if sum(int(r["recommended_n"]) for r in plan.values())!=298 or len(plan)!=6:
        raise ValueError("Must not silently change registered Phase-A design")
    # Surface observed source taxa and planned but unobserved taxa separately.
    all_taxa=sorted(set(counts)|set(plan))
    table=[]
    for taxon in all_taxa:
        n=counts.get(taxon,0)
        planned=int(plan[taxon]["recommended_n"]) if taxon in plan else None
        table.append({
            "taxon_as_published_or_planned":taxon,
            "public_young_leaf_RNA_seq_runs":n,
            "phase_A_target_n":planned,
            "actual_phase_A_individuals_collected_by_this_work":0,
            "public_flower_stage_RNA_seq_runs":0,
            "public_flower_pigment_metabolite_assays":0,
            "same_plant_verified_functional_toolkit":UNMEASURED,
            "induced_floral_reexpression":UNMEASURED,
            "heritable_repeated_reexpression":UNMEASURED,
            "reexpression_origin_history":UNKNOWN,
            "viable_achene_fitness_contrast":UNMEASURED,
            "readiness_label":"LEAF_RNA_METADATA_ONLY" if n>0 else PLANNED,
        })
    return {
        "schema":"aza3_p04_afterlife_public_source_readiness_v1",
        "scientific_status":"SOURCE_ONLY_PROVENANCE_AUDIT_NOT_FIELD_OR_FUNCTIONAL_RESULT",
        "source_git_blob_pins":{name:{"path":p,"sha":s} for name,(p,s) in SOURCE.items()},
        "source_summary":{
            "public_run_count":len(manifest),
            "public_source_taxon_names":len(counts),
            "young_leaf_runs":sum(counts.values()),
            "floral_RNA_runs_in_this_source":0,
            "Ryukyu_brevicaule_leaf_runs":counts["Cirsium brevicaule"],
            "Ryukyu_irumtiense_leaf_runs":counts["Cirsium irumtiense"],
            "takaoense_W_leaf_libraries":color_counts["W"],
            "takaoense_BP_leaf_libraries":color_counts["BP"],
            "takaoense_replication_confounded_by_geography":True,
            "published_2020_C_japonicum_family_evaluations":len(records),
            "distinct_2020_candidate_transcripts":5,
            "unresolved_2020_NJ_and_FastTree_evaluations":len(records),
            "verified_functional_pigmentation_enzymes":0,
            "distinct_target_phaseA_taxa":len(plan),
            "phaseA_target_individuals":298,
            "phaseA_collected_here":0,
            "phaseA_targets_are_not_independent_trait_return_events":True
        },
        "exclusions":{
            "leaf_RNA_not_floral_regulatory_competence":True,
            "2020_park_C_japonicum_transcripts_not_2026_Ryukyu_specimens":True,
            "molecular_family_ambiguity_not_biological_function_loss":True,
            "observed_W_BP_code_not_same_plant_allele_morph_link":True,
            "no_reexpression_absence_outcome_imputed":True
        },
        "per_taxon":table,
    }

def run(root=ROOT):
    return group_readiness(csv_rows("manifest",root),json_data("qc",root),
                           csv_rows("takaoense",root),json_data("family",root),
                           csv_rows("panel",root))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=ROOT)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--expected",type=Path)
    args=ap.parse_args()
    result=run(args.root)
    if args.expected:
        old=json.loads(args.expected.read_text(encoding="utf-8"))
        if old!=result:
            raise ValueError("Committed source-readiness receipt differs from recomputed source")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("AZA3_P04_SOURCE_GATE_PASS",json.dumps(result["source_summary"],sort_keys=True))

if __name__=="__main__":
    main()
