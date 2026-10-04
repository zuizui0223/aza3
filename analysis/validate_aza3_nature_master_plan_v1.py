#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT=Path(__file__).resolve().parents[1]
PLAN=ROOT/"docs"/"AZA3_NATURE_SCALE_MASTER_PLAN_V1.md"
SAMP=ROOT/"data"/"planning"/"aza3_nature_sampling_priorities_v2.csv"
README=ROOT/"README.md"
BRIDGE=ROOT/"data"/"planning"/"eazami_to_aza3_hypothesis_bridge_v1.csv"

def need(t,x):
    if x not in t:
        raise AssertionError(f"missing: {x}")

def main():
    p=PLAN.read_text(encoding="utf-8")
    r=README.read_text(encoding="utf-8")
    rows=list(csv.DictReader(SAMP.open(encoding="utf-8")))
    bridge=list(csv.DictReader(BRIDGE.open(encoding="utf-8")))

    for x in (
        "36/38 sampled Japanese concepts",
        "Orientation requires 4–6 minimum changes",
        "0/3 trait pairs pass robust shared transition localization",
        "Species-tip coding loses state multiplicity in 4/4 audited polymorphic systems",
        "aza3 hypothesis H0",
        "aza3 hypothesis H1",
        "aza3 hypothesis H2",
        "aza3 hypothesis H3",
        "Gate 0",
        "G1 — repeated de novo mutation",
        "G2 — ancestral standing variation",
        "G3 — introgressive reuse",
        "G4 — polyploid, homeolog or regulatory reuse",
        "genomic combinatorial reuse",
        "Do not call the radiation exceptionally evolvable before Gate 0",
    ):
        need(p,x)

    need(r,"docs/AZA3_NATURE_SCALE_MASTER_PLAN_V1.md")
    need(r,"data/planning/aza3_nature_sampling_priorities_v2.csv")
    need(r,"data/planning/eazami_to_aza3_hypothesis_bridge_v1.csv")

    ids=[x["priority_id"] for x in rows]
    assert ids==["N00","N01","N02","N03","N04","N05","N06"], ids
    assert rows[0]["stage"]=="rate_gate"
    assert "accelerated_evolvability" in rows[0]["models_discriminated"]
    assert "genomic_combinatorial_reuse" in rows[5]["models_discriminated"]
    assert rows[6]["stage"]=="generality"

    assert len(bridge)==5
    mapping={x["eazami_result_id"]:x["aza3_hypothesis_id"] for x in bridge}
    assert mapping=={
        "E1":"H0_RAPID_INNOVATION",
        "E2":"H1_REUSABLE_VARIATION",
        "E3":"H2_COMBINATORIAL_GENOMIC_REUSE",
        "E4":"H3_SEGREGATING_SUBSTRATE",
        "E5":"H4_CAPACITY_OVER_TRIGGER",
    }

    # Explicitly preserve the logic: EAzami result -> aza3 hypothesis, not conclusion.
    forbidden=(
        "EAzami demonstrates genomic reuse",
        "EAzami proves standing variation",
        "EAzami proves introgression",
        "EAzami proves genomic modularity",
    )
    low=p.casefold()
    for x in forbidden:
        if x.casefold() in low:
            raise AssertionError(f"overclaim: {x}")

    print(json.dumps({
        "status":"ok",
        "eazami_to_aza3_bridge":"explicit_H0_H1_H2_H3_H4",
        "gate0":"required_before_exceptional_evolvability_claim",
        "genomic_models":4,
        "nature_sampling_priorities":len(rows),
        "combinatorial_reuse_test":"N05",
        "generality_test":"N06"
    },indent=2))

if __name__=="__main__":
    main()
