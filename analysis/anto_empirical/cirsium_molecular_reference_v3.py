# Source migrated from zuizui0223/anto/src/anto_molecular_screen_v5.py (anto v4-v6); empirical ownership now aza3.
"""Fail-closed molecular pre-screen for anto v5.

Reads verified Asteraceae reference CDS sequences and *source-confirmed*
Ryukyu run×locus register. The new sequences are **non-Ryukyu positive
references**, never evidence that Cirsium inherited/used the same allele.
Optionally checks user-supplied candidate CDS contigs with strict identifiers.

Crucial: read expression in young leaf is not floral developmental competence.
Presence of an apparently intact CDS is NOT evidence of protein function;
absence of a leaf transcript is NOT a deletion; intact ORFs can be paralogs;
putative disruptive sites demand genomic support/orthology and ploidy control.
"""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

REFERENCES = {
    "DFR": ("EF600682.1", "Saussurea_medusa", 61, 1089),
    "ANS": ("AY547342.1", "Saussurea_medusa", 29, 1099),
    "GMYC1": ("AJ007709.1", "Gerbera_hybrida", 310, 1911),
}
RYUKYU_TAXA = ("Cirsium brevicaule", "Cirsium irumtiense")
CODONS = (
    "TTT F,TTC F,TTA L,TTG L,TCT S,TCC S,TCA S,TCG S,TAT Y,TAC Y,TAA *,TAG *,"
    "TGT C,TGC C,TGA *,TGG W,CTT L,CTC L,CTA L,CTG L,CCT P,CCC P,CCA P,CCG P,"
    "CAT H,CAC H,CAA Q,CAG Q,CGT R,CGC R,CGA R,CGG R,ATT I,ATC I,ATA I,ATG M,"
    "ACT T,ACC T,ACA T,ACG T,AAT N,AAC N,AAA K,AAG K,AGT S,AGC S,AGA R,AGG R,"
    "GTT V,GTC V,GTA V,GTG V,GCT A,GCC A,GCA A,GCG A,GAT D,GAC D,GAA E,GAG E,"
    "GGT G,GGC G,GGA G,GGG G"
)
GENETIC_CODE = dict(part.split(" ") for part in CODONS.split(","))


def read_fasta(path):
    records = {}
    name, seq = None, []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(seq).upper()
            name = line[1:].strip()
            if not name or name in records:
                raise ValueError("Duplicate/empty FASTA identifier")
            seq = []
        else:
            if name is None:
                raise ValueError("Sequence before FASTA header")
            seq.append(line)
    if name is not None:
        if name in records:
            raise ValueError("Duplicate FASTA identifier")
        records[name] = "".join(seq).upper()
    return records


def cds_quality(seq):
    """Pure sequence-level diagnostics, never a true biological LOF call."""
    seq = seq.upper()
    if not seq or any(b not in "ACGTN" for b in seq):
        raise ValueError("Nonempty A/C/G/T/N coding sequence required")
    n_fraction = seq.count("N") / len(seq)
    frame = (len(seq) % 3 == 0)
    codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
    stop_positions = [i for i, c in enumerate(codons) if GENETIC_CODE.get(c) == "*"]
    complete = frame and bool(codons) and codons[0] == "ATG" and codons[-1] in ("TAA", "TAG", "TGA")
    internal_stop = any(i < len(codons)-1 for i in stop_positions)
    if n_fraction > .02:
        status = "AMBIGUOUS_CDS_NOT_CALLABLE"
    elif not frame or not complete:
        status = "INCOMPLETE_OR_UNCERTAIN_CDS"
    elif internal_stop:
        status = "PREMATURE_STOP_CANDIDATE_REQUIRES_GENOMIC_VALIDATION"
    else:
        status = "OPEN_READING_FRAME_CANDIDATE_ONLY"
    return {
        "nucleotide_length": len(seq),
        "ambiguous_fraction": n_fraction,
        "divisible_by_three": frame,
        "canonical_start": bool(codons) and codons[0] == "ATG",
        "canonical_terminal_stop": bool(codons) and codons[-1] in ("TAA","TAG","TGA"),
        "internal_stop_codon_count": sum(i < len(codons)-1 for i in stop_positions),
        "sequence_qc_status": status,
        "sha256_cds": hashlib.sha256(seq.encode("ascii")).hexdigest(),
        "biological_interpretation": "Coding sequence candidate; no validated orthology, phenotype, or floral function.",
    }


def validate_references(fasta, provenance):
    records = read_fasta(fasta)
    registry = json.loads(Path(provenance).read_text(encoding="utf-8"))
    features = registry.get("reference_features", [])
    if registry.get("panel_id") != "anto_asteraceae_cds_references_v5" or len(features) != 3:
        raise ValueError("Wrong reference panel")
    if len(records) != len(REFERENCES):
        raise ValueError("Reference panel must have exactly three CDS sequences")
    out = {}
    for feature in features:
        locus = feature["gene"]
        if locus not in REFERENCES:
            raise ValueError("Unknown reference locus")
        accession, taxon, start, end = REFERENCES[locus]
        if (feature["accession"], feature["taxon"], feature["start"], feature["end"]) != (accession, taxon, start, end):
            raise ValueError("Reference provenance mismatch")
        header = f"{accession}|{locus}|{taxon}|CDS_{start}_{end}"
        if header not in records:
            raise ValueError("Missing annotated reference CDS")
        seq = records[header]
        q = cds_quality(seq)
        if (len(seq) != end-start+1 or
            q["sequence_qc_status"] != "OPEN_READING_FRAME_CANDIDATE_ONLY" or
            feature["cds_length"] != len(seq)):
            raise ValueError("Reference CDS not complete, or wrong coordinates")
        out[locus] = {**q,"reference_accession":accession,
                      "reference_species":taxon,
                      "reference_only_no_focal_sequence":True}
    return out


def load_register(path, source_run_manifest):
    with Path(source_run_manifest).open(newline="", encoding="utf-8") as handle:
        source = {r["run_accession"]: r for r in csv.DictReader(handle)}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 18:
        raise ValueError("Six Ryukyu runs × three loci required")
    expected_runs = {r["run_accession"] for r in source.values()
                     if r["taxon_as_experiment_titled"] in RYUKYU_TAXA}
    if len(expected_runs) != 6:
        raise ValueError("Unexpected Ryukyu source manifest")
    seen = set()
    for r in rows:
        locus,run = r["screen_locus"],r["run_accession"]
        key=(run,locus)
        if key in seen or locus not in REFERENCES or run not in expected_runs:
            raise ValueError("Bad or duplicate declared screen")
        seen.add(key)
        orig=source[run]
        for field,canonical in (
            ("taxon_title","taxon_as_experiment_titled"),
            ("source_specimen_code","specimen_code_from_experiment_title"),
            ("biosample_accession","biosample_accession")):
            if r[field] != orig[canonical]:
                raise ValueError(f"Incorrect sample crosswalk: {field}")
        if r["reference_accession"] != REFERENCES[locus][0]:
            raise ValueError("Locus reference swap")
        for f in ("coding_sequence_reconstructed","orthology_verified",
                  "genomic_coding_variant_verified","floral_expression_measured",
                  "floral_pigment_measured","independent_event_origin_inferred"):
            if r[f] != "false":
                raise ValueError(f"Unverified biological evidence cannot be admitted: {f}")
        if r["candidate_accession"] != "NOT_RETRIEVED" or r["status"] != "NOT_TESTED":
            raise ValueError("Unassessed project cannot make candidate presence claims")
    if seen != {(run,locus) for run in expected_runs for locus in REFERENCES}:
        raise ValueError("All expected sample×locus screens required")
    return rows


def inspect_optional_candidates(candidate_fasta, registered_rows):
    """Optional *hypothetical or separately sourced* candidate CDS stage.

    The 'candidate CDS' can be inspected only at the ORF QC level; a FASTA
    header never establishes that a sample truly carries a functional gene.
    Header: >SRRnnnnnnnn|DFR|contig_id.
    """
    if candidate_fasta is None:
        return []
    known = {(r["run_accession"],r["screen_locus"]) for r in registered_rows}
    counts = Counter()
    results = []
    for head, seq in read_fasta(candidate_fasta).items():
        bits = head.split("|")
        if len(bits) != 3 or not bits[2] or (bits[0],bits[1]) not in known:
            raise ValueError("Candidate is not linked to a registered focal run×locus")
        key = (bits[0],bits[1])
        counts[key] += 1
        q = cds_quality(seq)
        results.append({
            "run_accession":bits[0], "locus":bits[1], "candidate_contig":bits[2],
            **q,"orthology_verified":False,"genomic_validation":False,
            "floral_pathway_competence_established":False
        })
    for x in results:
        if counts[(x["run_accession"],x["locus"])] > 1:
            x["multiple_paralogs_or_assembly_candidates"] = True
            x["requires_locus_tree_before_gene_assignment"] = True
        else:
            x["multiple_paralogs_or_assembly_candidates"] = False
            x["requires_locus_tree_before_gene_assignment"] = True
    return sorted(results,key=lambda x:(x["run_accession"],x["locus"],x["candidate_contig"]))


def screen(reference_fasta, provenance, register_path, run_manifest, candidate_fasta=None):
    ref=validate_references(reference_fasta,provenance)
    rows=load_register(register_path,run_manifest)
    candidates=inspect_optional_candidates(candidate_fasta,rows)
    return {
       "version":"anto_molecular_preflight_v5","focal_run_count":6,"registered_screen_count":len(rows),
       "references":ref, "n_actual_focal_sequences_reconstructed_by_this_run":0,
       "n_candidate_sequences_supplied":len(candidates),
       "candidate_orf_qc":candidates,
       "focal_gene_state":"UNMEASURED" if not candidates else "CANDIDATE_ORF_QC_ONLY",
       "ancestral_regain_demonstrated":False,
       "flower_specific_competence_demonstrated":False,
       "introgression_or_denovo_event_resolved":False,
       "stop_rule":"Do not infer Ryukyu regain, pathway function, coding loss or variant history from reference genes, missing leaf transcripts or unverified candidate CDS.",
    }
