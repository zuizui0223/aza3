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

    # R1B-0 public target-capture preflight must stay bounded.
    assert pcon["status"] == "EXECUTABLE_PUBLIC_DATA_PREFLIGHT"
    assert pcon["source_bioproject"] == "PRJNA957074"
    assert pcon["focal_run"]["run"] == "SRR30887308"
    assert pcon["focal_run"]["biosample"] == "SAMN44017917"
    assert pcon["focal_run"]["taxon"] == "Cirsium sieboldii"
    assert len(pruns) == 10

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
        "PUBLIC_TRANSFER_GREEN",
        "PUBLIC_TRANSFER_AMBER",
        "PUBLIC_TRANSFER_RED",
        "Target-capture data interrogate a sparse, bait-defined subset of the genome.",
        "Do not expand it into a new phylogenomics project.",
    ):
        need(pdoc, token)

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
    need(readme, "reference_public_resource_inventory_v1.csv")

    print(json.dumps({
        "status": "ok",
        "r1a": "resource_audit_complete",
        "r1b0": "public_target_capture_preflight_executable",
        "r1b": "own_wgs_transferability_pilot_pending",
        "resource_rows": len(rows),
        "public_preflight_runs": len(pruns),
        "public_focal_run": pcon["focal_run"]["run"],
        "focal_system": con["focal_system"],
        "decision": con["current_decision"],
    }, indent=2))

if __name__ == "__main__":
    main()
