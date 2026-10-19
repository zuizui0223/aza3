# Source migrated from zuizui0223/anto/src/anto_sra_intake.py (anto v4-v6); empirical ownership now aza3.
"""Fail-closed acceptance of public PRJNA1311153 SRA *metadata*, not reads.

The 25 rows are ENA's publicly registered metadata for Chang et al. 2026:
- source tissue is young leaf; no flower-stage RNA can be asserted,
- taxon and code are parsed from submitted experiment titles,
- the specimen code is NOT independently verified as a herbarium voucher,
- taxon-level flower colour is NOT an individual experimental measurement,
- six takaoense phenotype labels are linked via the independently governed
  aza3 contract, one plant/locality and strongly confounded with geography.

This module deliberately uses the standard library and no network.
"""
from collections import Counter
import csv
import json
from pathlib import Path
import re

PROJECT = "PRJNA1311153"
EXPECTED_TAXON_COUNTS = {
    "Cirsium brevicaule": 3,
    "Cirsium irumtiense": 3,
    "Cirsium japonicum var. albescens": 2,
    "Cirsium japonicum var. australe": 3,
    "Cirsium japonicum var. fukienense": 4,
    "Cirsium japonicum var. japonicum": 2,
    "Cirsium japonicum var. takaoense": 6,
    "Cirsium morii": 2,
}
REQUIRED = (
    "study_accession", "run_accession", "experiment_accession",
    "biosample_accession", "experiment_title",
    "taxon_as_experiment_titled", "specimen_code_from_experiment_title",
    "tissue_from_project_description", "library_strategy",
    "library_source", "library_layout", "raw_fastq_pair_bytes",
    "raw_fastq_total_bytes", "individual_flower_colour", "source_api",
)
TITLE = re.compile(
    r"^Illumina MiniSeq sequencing: Transcriptome of (Cirsium .+)_([A-Za-z0-9]+-\d+)$"
)


def load_rows(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or any(key not in reader.fieldnames for key in REQUIRED):
            raise ValueError("Missing required manifest fields")
        rows = list(reader)
    if len(rows) != 25:
        raise ValueError(f"Expected 25 run rows, got {len(rows)}")
    for key, pattern in (
        ("run_accession",r"SRR\d+"),
        ("experiment_accession",r"SRX\d+"),
        ("biosample_accession",r"SAMN\d+"),
    ):
        seen=set()
        for row in rows:
            value=row[key]
            if re.fullmatch(pattern, value) is None or value in seen:
                raise ValueError(f"Invalid/duplicate {key}: {value}")
            seen.add(value)
    seen_codes=set()
    for row in rows:
        if row["study_accession"] != PROJECT:
            raise ValueError("Unexpected study")
        match=TITLE.fullmatch(row["experiment_title"])
        if not match or tuple(match.groups()) != (
            row["taxon_as_experiment_titled"], row["specimen_code_from_experiment_title"]
        ):
            raise ValueError(f"Taxon and specimen code must match source experiment title: {row['run_accession']}")
        if row["specimen_code_from_experiment_title"] in seen_codes:
            raise ValueError("A sample specimen code appears twice")
        seen_codes.add(row["specimen_code_from_experiment_title"])
        if row["tissue_from_project_description"] != "young leaf":
            raise ValueError("Do not classify leaf RNA as floral RNA")
        if row["individual_flower_colour"] != "NOT_MEASURED_IN_THIS_DATASET":
            raise ValueError("No individual-level floral pigmentation measurement in the primary SRA manifest")
        if (row["library_source"],row["library_strategy"],row["library_layout"]) != (
            "TRANSCRIPTOMIC","RNA-Seq","PAIRED"
        ):
            raise ValueError("Unexpected sequencing experiment design")
        parts=row["raw_fastq_pair_bytes"].split(";")
        if len(parts) != 2 or any(not p.isdigit() or int(p)<=0 for p in parts):
            raise ValueError("Required paired-end FASTQ-size declarations")
        if sum(int(p) for p in parts) != int(row["raw_fastq_total_bytes"]):
            raise ValueError("Per-run declared FASTQ byte size mismatch")
        if PROJECT not in row["source_api"] or "ebi.ac.uk/ena/" not in row["source_api"]:
            raise ValueError("Missing original ENA provenance")
    actual=dict(Counter(r["taxon_as_experiment_titled"] for r in rows))
    if actual!=EXPECTED_TAXON_COUNTS:
        raise ValueError(f"Changed taxon counts: {actual}")
    return rows


def validate_takaoense(morph_path, manifest_rows):
    with Path(morph_path).open(encoding="utf-8",newline="") as handle:
        reader=csv.DictReader(handle)
        rows=list(reader)
    if len(rows)!=6:
        raise ValueError("Expected six independently curated morphological labels")
    by_run={r["run_accession"]:r for r in manifest_rows}
    seen=set()
    counts=Counter()
    for r in rows:
        run=r["run_accession"]
        if run in seen or run not in by_run:
            raise ValueError("Unmatched or duplicate phenotype run")
        seen.add(run)
        orig=by_run[run]
        if orig["taxon_as_experiment_titled"]!="Cirsium japonicum var. takaoense":
            raise ValueError("Morph only assigned to verified takaoense run")
        if (r["specimen_code_from_experiment_title"] != orig["specimen_code_from_experiment_title"]
                or r["biosample_accession"] != orig["biosample_accession"]):
            raise ValueError("Morph/specimen BioSample crosswalk mismatch")
        if r["specimen_code_from_experiment_title"].split("-")[0] not in ("FC","WY","FB","TJ","NH","LT"):
            raise ValueError("Unknown specimen morph code")
        if r["morph_code_from_aza3"] not in ("W","BP"):
            raise ValueError("Invalid independently curated morph code")
        if not r["morphology_provenance"].startswith("aza3:"):
            raise ValueError("Missing independent morph provenance")
        for k in ("floral_rna_measured","floral_pigment_chemistry_measured","independent_population_replication"):
            if r[k] != "false":
                raise ValueError("Leaf transcriptomes do not supply functional floral assays")
        counts[r["morph_code_from_aza3"]]+=1
    if counts!={"W":3,"BP":3}:
        raise ValueError("3/3 morph labels required")
    if seen!={r["run_accession"] for r in manifest_rows if r["taxon_as_experiment_titled"]=="Cirsium japonicum var. takaoense"}:
        raise ValueError("Crosswalk must exactly cover takaoense samples")
    return rows


def validate_platform(platform_path, manifest_rows):
    p=json.loads(Path(platform_path).read_text(encoding="utf-8"))
    if p.get("n_runs") != 25 or p.get("source_project") != PROJECT:
        raise ValueError("Platform QC must match SRA project")
    runs=p.get("run_fingerprints")
    if not isinstance(runs,list) or len(runs)!=25:
        raise ValueError("Expected 25 platform run fingerprints")
    seen=set()
    for r in runs:
        if r["run_accession"] in seen:
            raise ValueError("Duplicate run fingerprint")
        seen.add(r["run_accession"])
        md5s=r["fastq_md5"].split(";")
        if len(md5s)!=2 or any(re.fullmatch(r"[a-fA-F0-9]{32}",s) is None for s in md5s):
            raise ValueError("Bad provider MD5 declaration")
        if r["instrument_model"] != "Illumina MiniSeq":
            raise ValueError("ENA platform values changed; re-audit discrepancy")
    if seen!={r["run_accession"] for r in manifest_rows}:
        raise ValueError("Fingerprint crosswalk mismatch")
    bytes_sum=sum(int(r["raw_fastq_total_bytes"]) for r in manifest_rows)
    if p["total_ena_fastq_bytes"] != bytes_sum:
        raise ValueError("Total fastq byte declarations mismatch")
    if "NovaSeq" not in p["paper_methods_instrument"] or "discrepancy" not in p["metadata_discrepancy"].lower():
        raise ValueError("Publication-vs-deposit provenance discrepancy must stay explicit")
    return p


def audit(manifest_path, morph_path, platform_path):
    rows=load_rows(manifest_path)
    morphs=validate_takaoense(morph_path,rows)
    instrument=validate_platform(platform_path,rows)
    by_taxon=dict(sorted(Counter(r["taxon_as_experiment_titled"] for r in rows).items()))
    ryukyu=[r for r in rows if r["taxon_as_experiment_titled"] in (
        "Cirsium brevicaule","Cirsium irumtiense")]
    return {
        "dataset":"PRJNA1311153_metadata_only",
        "manifest_n_runs":len(rows),
        "taxa":by_taxon,
        "ryukyu_n_runs":len(ryukyu),
        "ryukyu_samples":[
            {"taxon":r["taxon_as_experiment_titled"],
             "specimen_code":r["specimen_code_from_experiment_title"],
             "run":r["run_accession"],"biosample":r["biosample_accession"]}
            for r in ryukyu],
        "takaoense_verified_morph_crosswalk_n":len(morphs),
        "takaoense_morph_counts":dict(Counter(r["morph_code_from_aza3"] for r in morphs)),
        "declared_ena_fastq_bytes_total":instrument["total_ena_fastq_bytes"],
        "unresolved_instrument_inconsistency":True,
        "ryukyu_flower_stage_expression_available":False,
        "ryukyu_individual_pigment_assayed":False,
        "floral_reexpression_demonstrated":False,
        "evidence_tier":"SPECIMEN_METADATA_LINKED_NOT_FUNCTIONALLY_VALIDATED",
        "biological_conclusion":"Public accession-to-specimen bridge established, but trait polarity and functional reacquisition remain unmeasured.",
    }
