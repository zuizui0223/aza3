#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT=Path(__file__).resolve().parents[1]
PLAN=ROOT/"docs"/"AZA3_NATURE_SCALE_MASTER_PLAN_V1.md"
SAMP=ROOT/"data"/"planning"/"aza3_nature_sampling_priorities_v2.csv"
BRIDGE=ROOT/"data"/"planning"/"eazami_to_aza3_hypothesis_bridge_v1.csv"
README=ROOT/"README.md"
ORIGIN=ROOT/"docs"/"AZAMI_EAZAMI_TO_AZA3_NATURE_ORIGIN_HANDOFF_V1.md"

def need(t,x):
    if x not in t:
        raise AssertionError(f"missing: {x}")

def read_csv(p):
    with p.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def main():
    p=PLAN.read_text(encoding="utf-8")
    r=README.read_text(encoding="utf-8")
    o=ORIGIN.read_text(encoding="utf-8")
    srows=read_csv(SAMP)
    brows=read_csv(BRIDGE)

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
    need(r,"data/planning/eazami_to_aza3_hypothesis_bridge_v1.csv")
    need(r,"data/planning/aza3_nature_sampling_priorities_v2.csv")
    if "aza3_nature_sampling_priorities_v1.csv" in r:
        raise AssertionError("README still references redundant sampling v1")

    ids=[x["priority_id"] for x in srows]
    assert ids==["N00","N01","N02","N03","N04","N05","N06"], ids
    assert srows[0]["stage"]=="rate_gate"
    assert "accelerated_evolvability" in srows[0]["models_discriminated"]
    assert "genomic_combinatorial_reuse" in srows[5]["models_discriminated"]
    assert srows[6]["stage"]=="generality"

    assert len(brows)==5
    mapping={x["eazami_result_id"]:x["aza3_hypothesis_id"] for x in brows}
    assert mapping=={
        "E1":"H0_RAPID_INNOVATION",
        "E2":"H1_REUSABLE_VARIATION",
        "E3":"H2_COMBINATORIAL_GENOMIC_REUSE",
        "E4":"H3_SEGREGATING_SUBSTRATE",
        "E5":"H4_CAPACITY_OVER_TRIGGER",
    }

    need(o,"How could so much phenotypic change happen so quickly?")
    need(o,"Where does the capacity for rapid phenotypic innovation come from?")

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
        "eazami_to_aza3_bridge":"explicit_H0_to_H4",
        "gate0":"required_before_exceptional_evolvability_claim",
        "genomic_models":4,
        "nature_sampling_priorities":len(srows),
        "combinatorial_reuse_test":"N05",
        "generality_test":"N06",
        "sampling_v2":"canonical"
    },indent=2))

if __name__=="__main__":
    main()
