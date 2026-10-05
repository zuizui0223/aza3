#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "REFERENCE_GENOME_AND_PUBLIC_DATA_READINESS_V1.md"
CONTRACT = ROOT / "data" / "contracts" / "aza3_reference_public_data_readiness_v1.json"
INV = ROOT / "data" / "planning" / "reference_public_resource_inventory_v1.csv"
README = ROOT / "README.md"

def need(text, token):
    if token not in text:
        raise AssertionError(f"missing: {token}")

def main():
    doc = DOC.read_text(encoding="utf-8")
    con = json.loads(CONTRACT.read_text(encoding="utf-8"))
    rows = list(csv.DictReader(INV.open(encoding="utf-8")))
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

    # Fail closed on the exact inferential boundaries.
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

    need(readme, "REFERENCE_GENOME_AND_PUBLIC_DATA_READINESS_V1.md")
    need(readme, "reference_public_resource_inventory_v1.csv")

    print(json.dumps({
        "status": "ok",
        "r1a": "resource_audit_complete",
        "r1b": "transferability_pilot_pending",
        "resource_rows": len(rows),
        "focal_system": con["focal_system"],
        "decision": con["current_decision"],
    }, indent=2))

if __name__ == "__main__":
    main()
