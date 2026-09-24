from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "data/planning/chapter3_p02_stickiness_radseq_contract_v1.json"
SCHEMA = ROOT / "data/planning/chapter3_p02_same_individual_schema_v1.csv"
MANIFEST = ROOT / "data/planning/chapter3_p02_recommended_sample_manifest_v1.csv"
DECISIONS = ROOT / "data/planning/chapter3_p02_decision_matrix_v1.csv"
PLAN = ROOT / "docs/CHAPTER3_P02_STICKINESS_RADSEQ_EXECUTION_V1.md"


def main() -> None:
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    assert c["contract_id"] == "chapter3_p02_stickiness_radseq_contract_v1"
    assert c["status"] == "PRE_DATA_FROZEN_MIGRATED"
    assert c["aza3_priority"] == "P02"
    assert c["aza3_issue"] == 8
    src = c["source_provenance"]
    assert src["repository"] == "zuizui0223/EAzami"
    assert src["source_issue"] == 154 and src["source_pr"] == 155
    assert src["source_pr_head"] == "b3f91d98952b4f04bd27282e33d31b80694e248a"
    assert src["original_priority_label"] == "P01"

    focal = {x["paper_concept"]: x for x in c["focal_concepts"]}
    assert set(focal) == {"JPN_06", "JPN_15"}
    for concept in ("JPN_06", "JPN_15"):
        assert focal[concept]["minimum_individuals"] == 16
        assert focal[concept]["recommended_individuals"] == 24
        assert focal[concept]["minimum_populations"] == 2
        assert focal[concept]["recommended_populations"] == 3
    assert c["total_sampling_target"]["minimum_individuals"] == 32
    assert c["total_sampling_target"]["recommended_individuals"] == 48
    assert c["current_readiness"]["collection_authorized"] is False
    assert c["current_readiness"]["minimum_current_populations_verified_for_jpn06"] is False
    assert c["current_readiness"]["minimum_current_populations_verified_for_jpn15"] is False

    with SCHEMA.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    required_fields = {r["field"] for r in rows if r["required"] == "yes"}
    assert set(c["same_individual_required_fields"]).issubset(required_fields)

    with MANIFEST.open(encoding="utf-8", newline="") as fh:
        manifest = list(csv.DictReader(fh))
    assert len(manifest) == 48
    assert Counter(r["taxon_concept"] for r in manifest) == {"JPN06": 24, "JPN15": 24}
    assert Counter(r["taxon_concept"] for r in manifest if r["minimum_core"] == "yes") == {"JPN06": 16, "JPN15": 16}
    pop_counts = Counter((r["taxon_concept"], r["population_slot"]) for r in manifest)
    for taxon in ("JPN06", "JPN15"):
        assert [pop_counts[(taxon, p)] for p in ("P1", "P2", "P3")] == [8, 8, 8]
    assert all(r["authorization_status"] == "required" for r in manifest)

    with DECISIONS.open(encoding="utf-8", newline="") as fh:
        decisions = list(csv.DictReader(fh))
    assert [r["gate_id"] for r in decisions] == [f"G{i:02d}" for i in range(1, 12)]
    by_id = {r["gate_id"]: r for r in decisions}
    assert all(by_id[g]["chapter2_effect"] == "revise" for g in ("G08", "G09", "G10"))
    assert "P02b" in by_id["G11"]["next_action"]

    text = PLAN.read_text(encoding="utf-8")
    for phrase in (
        "P02",
        "32",
        "48",
        "Current biological gate",
        "Primary falsifiers",
        "P02b",
        "cannot by itself demonstrate defence",
    ):
        assert phrase in text, phrase

    print("chapter3_p02_stickiness_radseq_contract_v1: PASS")


if __name__ == "__main__":
    main()
