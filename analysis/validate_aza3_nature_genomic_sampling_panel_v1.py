#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"/"AZA3_NATURE_GENOMIC_SAMPLING_PANEL_V1.md"
PANEL=ROOT/"data"/"planning"/"aza3_nature_genomic_sampling_panel_v1.csv"
README=ROOT/"README.md"
MASTER=ROOT/"docs"/"AZA3_NATURE_SCALE_MASTER_PLAN_V1.md"
INTAKE=ROOT/"data"/"intake"/"aza3_nature_individual_intake_v1.csv"

def need(t,x):
    if x not in t:
        raise AssertionError(f"missing: {x}")

def main():
    d=DOC.read_text(encoding="utf-8")
    r=README.read_text(encoding="utf-8")
    m=MASTER.read_text(encoding="utf-8")
    with PANEL.open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    with INTAKE.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        intake_fields=reader.fieldnames or []

    phase_a=[x for x in rows if x["phase"]=="A"]
    counts={x["panel_id"]:int(x["recommended_n"]) for x in phase_a}
    assert counts=={
        "PEN":60,
        "SIE":40,
        "LIN":24,
        "DIP":24,
        "BREV":75,
        "IRUM":75,
    }
    assert sum(counts.values())==298

    sie=next(x for x in rows if x["panel_id"]=="SIE")
    assert "colour;orientation;phyllary;stickiness"==sie["phenotype_modules"]
    assert "genomic_combinatorial_reuse" in sie["models_discriminated"]
    assert "highest-value same-species multi-module system" in sie["genomic_role"]

    tak=next(x for x in rows if x["panel_id"]=="TAK")
    assert tak["phase"]=="B"
    assert "6 public + 40-60 new if opened"==tak["recommended_n"]
    assert "six public transcriptomes are not population replication" in tak["stop_rule"]

    for x in (
        "298 individuals",
        "highest-value combinatorial system",
        "Every Phase-A individual must link",
        "Tier 1 — all 298 Phase-A individuals",
        "Tier 2 — dense resequencing of selected systems",
        "Primary candidate: *C. sieboldii*.",
        "No RAD-only causal haplotype or allele-age claim.",
        "No field collection is authorized by this document.",
    ):
        need(d,x)

    for x in (
        "immutable individual ID",
        "voucher or voucher-linked diagnostic images",
        "standardized visible flower colour",
        "gravity-referenced capitulum orientation",
        "direct phyllary posture",
        "stickiness/gland/exudate state",
        "cytotype / relative genome size",
        "plastid companion",
    ):
        need(d,x)

    need(r,"AZA3_NATURE_GENOMIC_SAMPLING_PANEL_V1.md")
    need(r,"aza3_nature_genomic_sampling_panel_v1.csv")
    need(m,"AZA3_NATURE_GENOMIC_SAMPLING_PANEL_V1.md")

    required_intake={
        "individual_id","taxon_concept","population_id","deidentified_locality_key",
        "voucher_id","diagnostic_image_id","colour_state","orientation_deg_gravity",
        "orientation_state","phyllary_posture","stickiness_state","leaf_dna_tissue_id",
        "tier1_library_id","tier2_library_id","cytotype_x","relative_genome_size_2C",
        "flow_cytometry_record_id","plastid_sample_id","floral_rna_sample_id",
        "pigment_sample_id","access_authorization_id","collection_authorization_id",
        "conservation_review_id","phase","panel_id","admission_status","exclusion_reason"
    }
    assert required_intake.issubset(set(intake_fields))

    print(json.dumps({
        "status":"ok",
        "phaseA_total":298,
        "phaseA_systems":len(phase_a),
        "combinatorial_primary":"Cirsium sieboldii",
        "molecular_flagship_public_anchor":"takaoense_6_transcriptomes",
        "tier1_role":"ancestry_triage",
        "tier2_role":"genomic_source_discrimination",
        "field_authorization":False,
        "same_individual_intake_fields":len(intake_fields)
    },indent=2))

if __name__=="__main__":
    main()
