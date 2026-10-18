"""Validate aza3 canonical Cirsium empirical ownership and anto source handoff.

Strictly source/provenance-only; NOT a test of evolutionary re-expression,
LEGO module independence, or reproductive fitness. Uses Python stdlib.
"""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import csv
import json
import re

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"data/contracts/aza3_anto_empirical_lego_bridge_v1.json"
BASE=ROOT/"data/anto_source_handoff"


def _json(name):
    return json.loads((BASE/name).read_text(encoding="utf-8"))


def _csv(name):
    with (BASE/name).open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))


def _fasta(path):
    seen={}
    name=None
    bases=[]
    for line in path.read_text(encoding="utf-8").splitlines():
        line=line.strip()
        if line.startswith(">"):
            if name is not None:
                if name in seen:raise ValueError("Duplicate transcript ID")
                seen[name]="".join(bases)
            name=line[1:].split()[0]
            bases=[]
        elif line and name:
            bases.append(line.upper())
    if name is not None:
        if name in seen:raise ValueError("Duplicate transcript ID")
        seen[name]="".join(bases)
    return seen


def run():
    contract=json.loads(CONTRACT.read_text(encoding="utf-8"))
    assert contract["version"]=="aza3_anto_empirical_lego_bridge_v1"
    assert contract["canonical_empirical_repository"]=="zuizui0223/aza3"
    assert contract["theoretical_repository"]=="zuizui0223/anto"
    assert contract["samplings"]["phase_A_collected_n_currently_admitted"]==0
    samples=contract["samplings"]["by_taxon_target"]
    assert len(samples)==6
    assert sum(x["planned_n"] for x in samples)==contract["samplings"]["phase_A_original_total"]==298
    assert sorted(x["id"] for x in contract["experimental_gates"])==["T0","T1","T2","T3","T4"]
    assert len(contract["source_snapshots"])==12
    for path in contract["source_snapshots"]:
        assert (BASE/path).is_file(),f"Required empirical handoff source absent: {path}"

    reads=_csv("public_sra/prjna1311153_run_manifest_20261017.csv")
    assert len(reads)==25
    assert len({r["run_accession"] for r in reads})==25
    assert len({r["biosample_accession"] for r in reads})==25
    assert all(r["individual_flower_colour"]=="NOT_MEASURED_IN_THIS_DATASET" for r in reads)
    ry=[r for r in reads if r["taxon_as_experiment_titled"] in ("Cirsium brevicaule","Cirsium irumtiense")]
    assert Counter(r["taxon_as_experiment_titled"] for r in ry)=={"Cirsium brevicaule":3,"Cirsium irumtiense":3}
    assert {r["specimen_code_from_experiment_title"] for r in ry}=={
        "AMM-4163","YRN-4166","OKN-4295","MMY-4078","ISG-4111","YGN-4296"}
    assert all(r["tissue_from_project_description"]=="young leaf" for r in reads)
    profile=_json("public_sra/prjna1311153_ena_platform_fingerprints_20261017.json")
    assert profile["n_runs"]==25
    assert profile["ena_instrument_models"]=={"Illumina MiniSeq":25}
    assert "NovaSeq" in profile["paper_methods_instrument"] and profile["metadata_discrepancy"]

    morph=_csv("public_sra/takaoense_aza3_morph_crosswalk_20261017.csv")
    assert len(morph)==6
    assert Counter(r["morph_code_from_aza3"] for r in morph)=={"W":3,"BP":3}
    byrun={r["run_accession"]:r for r in reads}
    for r in morph:
        src=byrun[r["run_accession"]]
        assert src["taxon_as_experiment_titled"]=="Cirsium japonicum var. takaoense"
        assert src["specimen_code_from_experiment_title"]==r["specimen_code_from_experiment_title"]
        assert src["biosample_accession"]==r["biosample_accession"]
        assert r["floral_rna_measured"]=="false"

    focus=_csv("molecular/ryukyu_gene_screen_register_v5.csv")
    assert len(focus)==18
    assert all(r["status"]=="NOT_TESTED" and r["floral_expression_measured"]=="false" for r in focus)
    assert {r["run_accession"] for r in focus}=={r["run_accession"] for r in ry}

    source=_json("molecular/cirsium_japonicum_2020_26_sequence_provenance_v6.json")
    assert source["contract"]=="anto_park2020_26_source_transcripts_v6"
    assert source["sequence_extracted_rows"]==26
    assert source["selected_annotation_rows_verified_against_original_workbook"]==26
    assert source["source_fasta_entry_count"]==104890
    assert source["source_workbook_row_count"]==54791
    assert source["focal_Ryukyu_reexpression_confirmed"] is False
    fasta=_fasta(BASE/"molecular/cirsium_japonicum_2020_26_source_sequences_v6.fasta")
    assert len(fasta)==26
    assert len({r["id"] for r in source["records"]})==26
    assert set(fasta)=={r["id"] for r in source["records"]}
    for r in source["records"]:
        seq=fasta[r["id"]]
        assert len(seq)==r["source_fasta_length"]
        assert sha256(seq.encode("ascii")).hexdigest()==r["sha256_nucleotide"]
        assert r["orthology_verified"] is False and r["Ryukyu_2026_focal_sample"] is False
        assert r["protein_function_verified"] is False
    assert sum(map(len,fasta.values()))==source["selected_total_nucleotide_chars"]==20781

    annotations=_json("molecular/cirsium_japonicum_2020_blast_candidate_hits_v6.json")
    assert len(annotations["records"])==26 and annotations["no_ryukyu_sequences_recovered"] is True
    annotated={r["unigene_id"]:r for r in annotations["records"]}
    assert set(annotated)==set(fasta)
    for r in source["records"]:
        orig=annotated[r["id"]]
        assert r["published_blast_annotation"]==orig["annotated_best_blast_hit"]
        assert r["published_blast_subject"]==orig["blast_subject_identifier"]

    organs=_json("molecular/cirsium_japonicum_organ_evidence_v6.json")
    assert len(organs["records"])==11
    assert organs["six_ryukyu_specimens_linked"] is False
    assert organs["trait_regain_event_observed"] is False
    dfr=[r for r in organs["records"] if r["molecule"].startswith("DFR_")]
    assert len(dfr)==2 and all((r["numerator_organ"],r["denominator_organ"])==("root","leaf") for r in dfr)

    ext=_json("molecular/asteraceae_reference_provenance_v5.json")
    assert {r["gene"] for r in ext["reference_features"]}=={"DFR","ANS","GMYC1"}
    assert {r["accession"] for r in ext["reference_features"]}=={"EF600682.1","AY547342.1","AJ007709.1"}
    return {
        "bridge_version":contract["version"],
        "canonical_empirical_owner":"aza3",
        "source_snapshots_verified":len(contract["source_snapshots"]),
        "rna_source_runs":len(reads),
        "ryukyu_source_runs":len(ry),
        "public_takaoense_morphs":dict(Counter(r["morph_code_from_aza3"] for r in morph)),
        "2020_annotated_transcripts":len(fasta),
        "2020_sequence_nt":sum(map(len,fasta.values())),
        "admitted_own_collected_n":0,
        "molecular_reexpression_demonstrated":False,
        "independent_genetic_LEGO_demonstrated":False,
        "fieldwork_authorized":False
    }


if __name__=="__main__":
    print(json.dumps(run(),indent=2,ensure_ascii=False))
