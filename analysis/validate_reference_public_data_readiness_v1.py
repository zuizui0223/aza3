#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, json

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT/"docs"/"REFERENCE_GENOME_AND_PUBLIC_DATA_READINESS_V1.md"
PREFLIGHT_DOC = ROOT/"docs"/"REFERENCE_TRANSFER_PUBLIC_PREFLIGHT_V1.md"
RUNBOOK = ROOT/"docs"/"R1B0_EXECUTION_RUNBOOK_V1.md"
CONTRACT = ROOT/"data"/"contracts"/"aza3_reference_public_data_readiness_v1.json"
PREFLIGHT = ROOT/"data"/"contracts"/"aza3_reference_transfer_public_preflight_v1.json"
HANDOFF = ROOT/"data"/"contracts"/"aza3_eazami_moreyra_locus_handoff_v1.json"
INV = ROOT/"data"/"planning"/"reference_public_resource_inventory_v1.csv"
RUNS = ROOT/"data"/"planning"/"reference_public_preflight_runs_v1.csv"
TARGET_AUDIT = ROOT/"data"/"evidence"/"comp1061_public_target_audit_v1.json"
LOCI241 = ROOT/"data"/"evidence"/"moreyra_conservative_241_no_warning_loci_v1.txt"
METRICS = ROOT/"data"/"contracts"/"r1b0_primary_metric_registry_v1.csv"
RECOVERY_SCHEMA = ROOT/"data"/"templates"/"r1b0_target_recovery_results_v1.csv"
LOCALIZATION_SCHEMA = ROOT/"data"/"templates"/"r1b0_reference_localization_results_v1.csv"
README = ROOT/"README.md"

SHA241 = "d561c6e393b1964fdd4b3acf14fda8b10f2f43923b1074cd35f86bfed07ebf73"
ARTIFACT_DIGEST = "sha256:9503762705a12bdf9bc28496cb6f5a273e388a86f3b1b255634c30012a35ad4d"

def need(text, token):
    if token not in text:
        raise AssertionError(f"missing: {token}")

def csvrows(path):
    return list(csv.DictReader(path.open(encoding="utf-8")))

def main():
    doc=DOC.read_text(encoding="utf-8")
    pdoc=PREFLIGHT_DOC.read_text(encoding="utf-8")
    runbook=RUNBOOK.read_text(encoding="utf-8")
    con=json.loads(CONTRACT.read_text(encoding="utf-8"))
    pcon=json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    hand=json.loads(HANDOFF.read_text(encoding="utf-8"))
    inv=csvrows(INV)
    runs=csvrows(RUNS)
    taudit=json.loads(TARGET_AUDIT.read_text(encoding="utf-8"))
    metrics=csvrows(METRICS)
    readme=README.read_text(encoding="utf-8")

    assert con["status"]=="R1A_RESOURCE_AUDIT_COMPLETE__R1B_TRANSFERABILITY_PILOT_PENDING"
    assert con["current_decision"]=="PROCEED_TO_R1B_TRANSFERABILITY_PILOT"
    assert con["focal_system"]=="Cirsium sieboldii"

    observed={r["accession_or_doi"] for r in inv}
    required={
        "GCA_965225835.1","GCA_965276805.1","PRJNA1127082",
        "PRJNA957074","PRJNA1158676","PRJNA1311153",
        "github:carol-siniscalchi/Comp1061-Angio353/comp1061_hybpiper_reference.fasta",
        "SRR25265669",
    }
    assert required <= observed, sorted(required-observed)

    # Public Compositae1061 compatibility reference.
    assert taudit["fasta_sequence_records"]==2597
    assert taudit["distinct_loci"]==1061
    assert set(taudit["reference_prefixes"])=={"lett","saff","sunf"}

    # EAzami -> aza3 public-author-repository handoff.
    assert hand["source_workflow_run"]==31557261287
    assert hand["source_artifact"]["id"]==9126428776
    assert hand["source_artifact"]["digest"]==ARTIFACT_DIGEST
    assert hand["historical_target_counts"]["paper_reported_mapped_target_loci"]==1064
    assert hand["historical_target_counts"]["public_named_loci_recovered_from_author_matrices"]==1061
    assert hand["historical_target_counts"]["unresolved_difference"]==3
    assert hand["exact_final_350_locus_names_recovered"] is False
    sets=hand["reproducible_locus_sets"]
    assert sets["public_1061"]["count"]==1061
    assert sets["reproducible_531"]["count"]==531
    assert sets["conservative_241"]["count"]==241
    assert sets["conservative_241"]["sha256"]==SHA241
    assert sets["manual_review_290"]["count"]==290

    vals=[x.strip() for x in LOCI241.read_text(encoding="utf-8").splitlines() if x.strip()]
    assert len(vals)==241 and len(set(vals))==241
    assert hashlib.sha256(LOCI241.read_bytes()).hexdigest()==SHA241

    # Focal public run and target provenance.
    assert pcon["status"]=="EXECUTABLE_PUBLIC_DATA_PREFLIGHT"
    assert pcon["focal_run"]["run"]=="SRR30887308"
    assert pcon["focal_run"]["biosample"]=="SAMN44017917"
    assert pcon["focal_run"]["taxon"]=="Cirsium sieboldii"
    hu=pcon["historical_target_universe"]
    assert hu["paper_reported_mapped_loci"]==1064
    assert hu["reproducible_public_named_loci"]==1061
    assert hu["unresolved_difference"]==3
    assert hu["exact_author_1064_target_recovered"] is False
    target_ids={x["id"] for x in pcon["target_reference_versions"]}
    assert target_ids=={"PUBLIC_COMP1061_1061","PUBLIC_1061_PLUS_TIOGANUM_COMPATIBILITY"}
    fixed=pcon["fixed_compatibility_layer"]
    assert fixed["path"]=="data/evidence/moreyra_conservative_241_no_warning_loci_v1.txt"
    assert fixed["count"]==241 and fixed["sha256"]==SHA241
    assert pcon["clean_locus_layers"]==[
        "PUBLIC_1061_BROAD_RECOVERY","MOREYRA_COMPATIBILITY_241","AUTO_STRICT_CLEAN"
    ]

    # Frozen comparison panel.
    assert len(runs)==10
    byrun={r["run"]:r for r in runs}
    for run,taxon in {
        "SRR30887308":"Cirsium sieboldii",
        "SRR30887271":"Cirsium japonicum",
        "SRR30887240":"Cirsium lineare",
        "SRR25265660":"Cirsium nipponicum",
        "SRR25265649":"Cirsium pendulum",
        "SRR30887259":"Cirsium dipsacolepis",
        "SRR30887291":"Cirsium tanakae",
    }.items():
        assert byrun[run]["taxon"]==taxon

    # Metrics were frozen before any empirical R1B-0 focal result.
    assert [m["metric_id"] for m in metrics]==["M01","M02","M03","M04","M05","M06","M07"]
    assert all(m["status"]=="FROZEN_BEFORE_RESULTS" for m in metrics)
    assert metrics[0]["definition"].endswith("PUBLIC_1061_BROAD_RECOVERY")
    assert "MOREYRA_COMPATIBILITY_241" in metrics[1]["definition"]
    assert metrics[6]["primary_or_sensitivity"]=="sensitivity"

    rec_header=RECOVERY_SCHEMA.read_text(encoding="utf-8").splitlines()[0].split(",")
    assert rec_header[:10]==[
        "run","taxon","target_version","hybpiper_version",
        "public_named_loci_count","public_named_loci_recovered",
        "public_named_loci_recovery_fraction","compatibility_241_count",
        "compatibility_241_recovered","compatibility_241_recovery_fraction",
    ]
    loc_header=LOCALIZATION_SCHEMA.read_text(encoding="utf-8").splitlines()[0].split(",")
    assert loc_header[:6]==["run","taxon","target_version","clean_layer","locus","genome_reference"]

    # Textual claim boundaries cannot drift back to exact-reproduction language.
    for t in (
        "PUBLIC_COMP1061_1061",
        "PUBLIC_1061_PLUS_TIOGANUM_COMPATIBILITY",
        "MOREYRA_COMPATIBILITY_241",
        "1,064",
        "1,061",
        "published final 350",
        "SRR30887308",
        "PUBLIC_TRANSFER_GREEN","PUBLIC_TRANSFER_AMBER","PUBLIC_TRANSFER_RED",
    ):
        need(pdoc,t)
    for t in (
        "MOREYRA_COMPATIBILITY_241",
        "d561c6e393b1964fdd4b3acf14fda8b10f2f43923b1074cd35f86bfed07ebf73",
        "1,064",
        "1,061",
    ):
        need(runbook,t)
    assert "Therefore a Cirsium-adapted target can be reconstructed from public materials" not in pdoc

    forbidden={
        "local ancestry","recombination blocks","haplotype age",
        "structural variant reuse","module-specific local genealogy",
        "genotype-phenotype association",
    }
    assert forbidden <= set(pcon["prohibited_inferences"])

    for t in (
        "REFERENCE_GENOME_AND_PUBLIC_DATA_READINESS_V1.md",
        "REFERENCE_TRANSFER_PUBLIC_PREFLIGHT_V1.md",
        "aza3_eazami_moreyra_locus_handoff_v1.json",
    ):
        need(readme,t)

    print(json.dumps({
        "status":"ok",
        "r1a":"resource_audit_complete",
        "r1b0":"public1061_plus_frozen241_preflight_executable",
        "historical_target_count_reported":1064,
        "reproducible_public_named_loci":1061,
        "unresolved_target_locus_difference":3,
        "primary_compatibility_loci":241,
        "public_preflight_runs":len(runs),
        "focal_run":pcon["focal_run"]["run"],
        "next":"await_or_execute_focal_hybpiper_public1061_and_241_recovery",
    },indent=2))

if __name__=="__main__":
    main()
