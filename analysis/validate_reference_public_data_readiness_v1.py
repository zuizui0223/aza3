#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "REFERENCE_GENOME_AND_PUBLIC_DATA_READINESS_V1.md"
CONTRACT = ROOT / "data" / "contracts" / "aza3_reference_public_data_readiness_v1.json"
INV = ROOT / "data" / "planning" / "reference_public_resource_inventory_v1.csv"
PREFLIGHT_DOC = ROOT / "docs" / "REFERENCE_TRANSFER_PUBLIC_PREFLIGHT_V1.md"
PREFLIGHT_CONTRACT = ROOT / "data" / "contracts" / "aza3_reference_transfer_public_preflight_v1.json"
PREFLIGHT_RUNS = ROOT / "data" / "planning" / "reference_public_preflight_runs_v1.csv"
TARGET_AUDIT = ROOT / "data" / "evidence" / "comp1061_public_target_audit_v1.json"
METRICS = ROOT / "data" / "contracts" / "r1b0_primary_metric_registry_v1.csv"
RECOVERY_SCHEMA = ROOT / "data" / "templates" / "r1b0_target_recovery_results_v1.csv"
LOCALIZATION_SCHEMA = ROOT / "data" / "templates" / "r1b0_reference_localization_results_v1.csv"
README = ROOT / "README.md"

def need(text, token):
    if token not in text:
        raise AssertionError(f"missing: {token}")

def main():
    doc = DOC.read_text(encoding="utf-8")
    con = json.loads(CONTRACT.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(INV.open(encoding="utf-8")))
    pdoc = PREFLIGHT_DOC.read_text(encoding="utf-8")
    pcon = json.loads(PREFLIGHT_CONTRACT.read_text(encoding="utf-8"))
    pruns = list(csv.DictReader(PREFLIGHT_RUNS.open(encoding="utf-8")))
    taudit = json.loads(TARGET_AUDIT.read_text(encoding="utf-8"))
    metrics = list(csv.DictReader(METRICS.open(encoding="utf-8")))
    recovery_header = RECOVERY_SCHEMA.read_text(encoding="utf-8").splitlines()[0].split(",")
    localization_header = LOCALIZATION_SCHEMA.read_text(encoding="utf-8").splitlines()[0].split(",")
    readme = README.read_text(encoding="utf-8")

    assert con["status"] == "R1A_RESOURCE_AUDIT_COMPLETE__R1B_TRANSFERABILITY_PILOT_PENDING"
    assert con["current_decision"] == "PROCEED_TO_R1B_TRANSFERABILITY_PILOT"
    assert con["focal_system"] == "Cirsium sieboldii"

    required_resources = {
        "GCA_965225835.1",
        "GCA_965276805.1",
        "PRJNA1127082",
        "PRJNA957074",
        "PRJNA1158676",
        "PRJNA1311153",
        "github:carol-siniscalchi/Comp1061-Angio353/comp1061_hybpiper_reference.fasta",
        "SRR25265669",
    }
    observed = {r["accession_or_doi"] for r in rows}
    missing = required_resources - observed
    if missing:
        raise AssertionError(f"missing required public resources: {sorted(missing)}")

    for token in (
        "R1B — C. sieboldii reference-transfer pilot",
        "GREEN_EXISTING_REFERENCE_SUFFICIENT",
        "AMBER_FOCAL_REFERENCE_REQUIRED",
        "RED_PANGENOME_OR_GRAPH_REQUIRED",
        "NOT_IDENTIFIABLE",
        "A high overall mapping rate is not enough.",
        "Moreyra LD et al. 2025",
        "PRJNA957074",
        "PRJNA1158676",
        "PRJNA1311153",
        "10.6084/m9.figshare.26927092",
    ):
        need(doc, token)

    for token in (
        "Target-capture loci from PRJNA957074 may inform broad topology and orthology but cannot by themselves establish fine local ancestry, haplotype age, recombination blocks, or structural-variant reuse.",
        "Young-leaf transcriptomes from PRJNA1158676 and PRJNA1311153 may inform coding orthology and expressed haplotypes but cannot establish intergenic regulatory architecture or genome-wide local ancestry.",
        "A high mapping rate to a congener does not by itself authorize local-genealogy inference.",
    ):
        if token not in con["claim_boundaries"]:
            raise AssertionError(f"missing claim boundary: {token}")

    focal = [r for r in rows if r["resource_id"] == "R12"]
    assert len(focal) == 1
    assert focal[0]["status"] == "MISSING_PUBLIC"
    assert "audit-bounded" in focal[0]["claim_ceiling"]

    # Audit the public Compositae1061 target universe before downstream filtering.
    assert taudit["source_blob_sha"] == "4f89e234007f367ffa8aa5e2be536bc44f31f445"
    assert taudit["fasta_sequence_records"] == 2597
    assert taudit["distinct_loci"] == 1061
    assert set(taudit["reference_prefixes"]) == {"lett", "saff", "sunf"}

    # R1B-0 public target-capture preflight must stay bounded.
    assert pcon["status"] == "EXECUTABLE_PUBLIC_DATA_PREFLIGHT"
    assert pcon["source_bioproject"] == "PRJNA957074"
    assert pcon["focal_run"]["run"] == "SRR30887308"
    assert pcon["focal_run"]["biosample"] == "SAMN44017917"
    assert pcon["focal_run"]["taxon"] == "Cirsium sieboldii"
    assert len(pruns) == 10

    target_versions = {x["id"]: x for x in pcon["target_reference_versions"]}
    assert set(target_versions) == {"ORIGINAL_COMP1061", "RECONSTRUCTED_CIRSIUM_TARGET"}
    assert target_versions["ORIGINAL_COMP1061"]["status"] == "PUBLIC"
    recon = target_versions["RECONSTRUCTED_CIRSIUM_TARGET"]
    assert recon["source_run"] == "SRR25265669"
    assert recon["source_experiment"] == "SRX21011548"
    assert recon["source_biosample"] == "SAMN34240347"
    assert recon["source_bases"] == 36939020978
    assert recon["status"] == "RECONSTRUCTABLE_PUBLICLY"

    stages = [x["stage"] for x in pcon["stages"]]
    assert stages == ["target_file_sensitivity", "genome_reference_localization"]

    orth = pcon["published_moreyra_orthology_rule"]
    assert orth["discard_if_paralog_warnings_gt"] == 10
    assert orth["intermediate_warning_range"] == "1-10"
    assert orth["final_missing_data_lt"] == 0.5
    assert orth["final_species_presence_gte"] == 0.8
    assert orth["published_final_loci"] == 350
    assert pcon["clean_locus_layers"] == [
        "ALL_1061_TARGET_LOCI",
        "PUBLISHED_RULE_COMPATIBLE",
        "AUTO_STRICT_CLEAN",
    ]

    runs = {r["run"]: r for r in pruns}
    required_runs = {
        "SRR30887308": "Cirsium sieboldii",
        "SRR30887271": "Cirsium japonicum",
        "SRR30887240": "Cirsium lineare",
        "SRR25265660": "Cirsium nipponicum",
        "SRR25265649": "Cirsium pendulum",
        "SRR30887259": "Cirsium dipsacolepis",
        "SRR30887291": "Cirsium tanakae",
    }
    for run, taxon in required_runs.items():
        assert run in runs, run
        assert runs[run]["taxon"] == taxon, (run, runs[run]["taxon"], taxon)

    for token in (
        "SRR30887308",
        "SRR25265669",
        "ORIGINAL_COMP1061",
        "RECONSTRUCTED_CIRSIUM_TARGET",
        "more than 10 HybPiper paralog warnings",
        "PUBLISHED_RULE_COMPATIBLE",
        "AUTO_STRICT_CLEAN",
        "PUBLIC_TRANSFER_GREEN",
        "PUBLIC_TRANSFER_AMBER",
        "PUBLIC_TRANSFER_RED",
        "Target-capture data interrogate a sparse, bait-defined subset of the genome.",
        "Do not expand it into a new phylogenomics project.",
    ):
        need(pdoc, token)

    assert pcon["primary_metric_registry"] == "data/contracts/r1b0_primary_metric_registry_v1.csv"
    assert pcon["result_schemas"]["target_recovery"] == "data/templates/r1b0_target_recovery_results_v1.csv"
    assert pcon["result_schemas"]["reference_localization"] == "data/templates/r1b0_reference_localization_results_v1.csv"
    assert [m["metric_id"] for m in metrics] == ["M01","M02","M03","M04","M05","M06","M07"]
    assert all(m["status"] == "FROZEN_BEFORE_RESULTS" for m in metrics)
    assert recovery_header[:4] == ["run","taxon","target_version","hybpiper_version"]
    assert localization_header[:6] == ["run","taxon","target_version","clean_layer","locus","genome_reference"]

    forbidden = {
        "local ancestry",
        "recombination blocks",
        "haplotype age",
        "structural variant reuse",
        "module-specific local genealogy",
        "genotype-phenotype association",
    }
    if not forbidden.issubset(set(pcon["prohibited_inferences"])):
        raise AssertionError("R1B-0 inference ceiling was relaxed")

    need(readme, "REFERENCE_GENOME_AND_PUBLIC_DATA_READINESS_V1.md")
    need(readme, "REFERENCE_TRANSFER_PUBLIC_PREFLIGHT_V1.md")
    need(readme, "reference_public_resource_inventory_v1.csv")
    need(readme, "reference_public_preflight_runs_v1.csv")

    print(json.dumps({
        "status": "ok",
        "r1a": "resource_audit_complete",
        "r1b0": "two_stage_public_target_capture_preflight_executable",
        "r1b": "own_wgs_transferability_pilot_pending",
        "resource_rows": len(rows),
        "public_preflight_runs": len(pruns),
        "public_target_loci": taudit["distinct_loci"],
        "public_target_records": taudit["fasta_sequence_records"],
        "public_focal_run": pcon["focal_run"]["run"],
        "cirsium_target_source_run": recon["source_run"],
        "focal_system": con["focal_system"],
        "decision": con["current_decision"],
    }, indent=2))

if __name__ == "__main__":
    main()
