#!/usr/bin/env python3
"""Scientific provenance firewall: *Cirsium schantarense* != *C. shansiense*.

A nuclear ITS result in Bae 2015 describes Korean *C. shantarense*
(accepted C. schantarense). Accession MN871982 is a plastid genome
of CHINESE *C. shansiense*, voucher WNU52012, Xu et al. 2020.

A secondary 2022 paper assigns that accession to *C. schantarense*.
The original sequence/material source controls its known identity
unless independently reidentified at voucher/record level.

Input files are literature provenance declarations. This script does
not fetch NCBI nucleotides, infer phylogenetic distance, measure a
head, or assert a white flower pigment regain event.
"""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path

TARGET="Cirsium schantarense"
OTHER="Cirsium shansiense"
ITS="nuclear_ITS"
PLASTID="plastid"
ACCESSION="MN871982"

def require(ok, code):
    if not ok: raise ValueError(code)

def canonical(x):
    x=" ".join(str(x).strip().split())
    if x=="Cirsium shantarense":
        return TARGET
    if x in (TARGET,OTHER):
        return x
    raise ValueError("UNREGISTERED_CIRSIUM_SPECIES_OR_FORMAT")

def validate_source(j, external):
    require(j["target"]["accepted_species"]==TARGET,
            "TARGET_SPECIES_SOURCE_CHANGED")
    bae=j["target"]["Bae_original_publication"]
    require(bae["doi"]=="10.5352/JLS.2015.25.2.243"
            and canonical(j["target"]["observed_spelling_in_Bae_2015"])==TARGET,
            "BAE_2015_PRIMARY_TARGET_MISMATCH")
    require("ITS" in bae["nuclear_marker"]
            and bae["sample_n_of_target"]==1
            and bae["accession"] is None
            and bae["accession_status"]=="NOT_INDEPENDENTLY_VERIFIED",
            "UNVERIFIED_BAE_ITS_ACCESSION_FALSELY_PROMOTED")
    require(bae["collection_site"]=="Dooryoon Mountain, Korea",
            "BAE_PRIMARY_SPECIMEN_SITE_DRIFT")
    a=j["accession"]
    require(a["accession"]==ACCESSION
            and canonical(a["original_species"])==OTHER
            and a["original_genome_compartment"]==PLASTID
            and a["reported_length_bp"]==152396
            and a["voucher"]=="WNU52012"
            and a["collecting_region"]=="Qingyang, Gansu, China"
            and a["primary_source_doi"]=="10.1080/23802359.2020.1721032",
            "MN871982_NOT_THE_ORIGINAL_C_SHANSIENSE_PLASTOME")
    require(not a["genbank_sequence_bytes_downloaded_in_this_turn"],
            "SEQ_DOWNLOAD_CLAIM_WITHOUT_PRIMARY_FASTA")
    second=a["secondary_source_discordance"]
    require(second["printed_label"]==TARGET
            and second["printed_accession"]==ACCESSION
            and not second["independent_voucher_confirmation"]
            and second["assessment"]==
            "CONFLICTS_WITH_2020_PRIMARY_VOUCHER_SPECIES_DO_NOT_REASSIGN",
            "SECONDARY_ERROR_WAS_FALSELY_PROMOTED")
    tax=j["taxonomy"]
    require({x["name"] for x in tax["distinct_current_accepted_species"]}==
            {TARGET,OTHER},"DISTINCT_ACCEPTED_SPECIES_FALSELY_MERGED")
    require(tax["white_form"]["name"]==
            "Cirsium schantarense f. albiflorum"
            and tax["white_form"]["accepted_taxon"]==TARGET
            and tax["white_form"]["status"]==
            "SYNONYM_NO_INDEPENDENT_COLOUR_REGAIN_EVENT_VERIFIED"
            and not tax["white_form"]["genetic_cause_known"]
            and not tax["white_form"]["pigment_reinduction_assayed"]
            and not tax["white_form"]["heritable_reappearance_demonstrated"],
            "WHITE_FORM_NAME_WAS_FALSELY_PROMOTED_TO_REEVOLUTION")
    for key,val in j["evidence_gates"].items():
        require(val is False,"GENOMIC_OR_COEVOLUTION_GATE_FALSELY_PASSED:"+key)

    literature=external["independent_lower_resolution_nuclear_literature"]
    require(literature["doi"]==bae["doi"]
            and literature["source_publication_has_molecular_observation"]
            and not literature["is_moreyra_2025_1061_locus_nuclear_sample"]
            and not literature["verified_GenBank_accession_or_sequence_recovered_here"]
            and not literature["shared_comparable_branch_length_ensemble_with_2025_Moreyra"]
            and not literature["trait_evolution_direction_or_genetic_resource_regain_inferred"],
            "PREVIOUS_THREE_QUESTION_CROSSWALK_SCOPE_CONTRADICTION")
    require(external["n_downward_candidate_species_with_both"]==0
            and external["n_species_with_both_safe_nuclear_sample_and_primary_occurrence_support"]==1,
            "PRIOR_PINNED_NUCLEAR_GBIF_INTERSECTION_DRIFT")
    return {
      "status":"PASS_PRIMARY_VOUCHER_AND_GENOME_COMPARTMENT_SOURCE_GATES",
      "BAE_2015_Korean_schantarense_ITS_exists_in_publication":True,
      "BAE_2015_Korean_schantarense_ITS_accession_recovered":False,
      "MN871982_original_species":OTHER,
      "MN871982_original_genome_compartment":PLASTID,
      "MN871982_voucher":"WNU52012",
      "2022_secondary_species_label_disagrees_with_2020_original":True,
      "Kew_two_species_are_distinct_accepted_names":True,
      "white_form_taxon_name_is_not_heritable_regain":True,
      "Q1_phylogenetic_memory_estimate":None,
      "Q3_reexpression_or_genetic_reuse_evidence":None,
      "genuine_new_nucleotide_sequences_recovered":0,
    }

def check_claim(species,compartment,accession=None,source="primary"):
    species=canonical(species)
    require(compartment in {ITS,PLASTID},"GENOME_COMPARTMENT_UNCLASSIFIED")
    require(source in {"primary","secondary"},"EVIDENCE_SOURCE_NOT_DEFINED")
    if accession==ACCESSION:
        require(source=="primary",
                "SECONDARY_PLASTOME_LABEL_NOT_INDEPENDENT_GENBANK_VALIDATION")
        require(species==OTHER,
                "MN871982_WRONG_SPECIES_C_SCHANTARENSE_NOT_SHANSIENSE")
        require(compartment==PLASTID,
                "MN871982_IS_PLASTID_NOT_NUCLEAR_ITS")
        return "ORIGINAL_XU_2020_SHANSIENSE_PLASTOME_ONLY"
    require(accession is None,
            "UNVALIDATED_ACCESSION_MUST_NOT_BE_INFERRED")
    if species==TARGET and compartment==ITS:
        return "BAE_2015_NUCLEAR_ITS_LITERATURE_NO_ACCESSION_RECOVERED"
    return "NO_SPECIMEN_LINK_PROVEN_BY_THIS_LITERATURE_MANIFEST"

def claim_regained_white_form(name,genotype_linked,independent_origin):
    if name=="Cirsium schantarense f. albiflorum":
        require(genotype_linked is True and independent_origin is True,
                "WHITE_FORM_SYNONYM_NOT_INDEPENDENT_GENETIC_REGAIN")
        return "FUTURE_INDEPENDENT_GENETIC_EVIDENCE_REQUIRED_SOURCE_REGISTRY_CANNOT_VERIFY"
    raise ValueError("UNREGISTERED_COLOUR_FORM_NO_REEVOLUTION_TEST")

def tests(j,external):
    good=validate_source(j,external)
    require(check_claim("Cirsium shantarense",ITS) ==
            "BAE_2015_NUCLEAR_ITS_LITERATURE_NO_ACCESSION_RECOVERED",
            "2015_ORIGINAL_SPELLING_FALSIFIED")
    require(check_claim(OTHER,PLASTID,ACCESSION)==
            "ORIGINAL_XU_2020_SHANSIENSE_PLASTOME_ONLY",
            "TRUE_PRIMARY_ACCESSION_REJECTED")
    bad_claims=[
      (TARGET,PLASTID,ACCESSION,"primary","MN871982_WRONG_SPECIES"),
      (TARGET,ITS,ACCESSION,"primary","MN871982_WRONG_SPECIES"),
      (OTHER,ITS,ACCESSION,"primary","MN871982_IS_PLASTID"),
      (OTHER,PLASTID,ACCESSION,"secondary","SECONDARY_PLASTOME_LABEL"),
      (TARGET,ITS,"AB123456","primary","UNVALIDATED_ACCESSION"),
    ]
    for species,comp,acc,source,needle in bad_claims:
        try:check_claim(species,comp,acc,source)
        except ValueError as e:require(needle in str(e),
                                      "INCORRECT_REJECTION_REASON_"+str(e))
        else:raise AssertionError("INVALID_GENBANK_CIRSIUM_CLAIM_ACCEPTED")
    try:claim_regained_white_form("Cirsium schantarense f. albiflorum",
                                  False,False)
    except ValueError as e:
        require("WHITE_FORM_SYNONYM_NOT_INDEPENDENT_GENETIC_REGAIN" in str(e),
                "WHITE_FORM_GATE_INCORRECT")
    else:raise AssertionError("A_WHITE_TAXON_NAME_WAS_TREATED_AS_REGATTRAINED")

    corruptions=[
      (lambda x:x["accession"].update(original_species=TARGET),
       "MN871982_NOT_THE_ORIGINAL_C_SHANSIENSE_PLASTOME"),
      (lambda x:x["accession"].update(original_genome_compartment=ITS),
       "MN871982_NOT_THE_ORIGINAL_C_SHANSIENSE_PLASTOME"),
      (lambda x:x["target"]["Bae_original_publication"].update(accession="MN871982"),
       "UNVERIFIED_BAE_ITS_ACCESSION_FALSELY_PROMOTED"),
      (lambda x:x["taxonomy"]["white_form"].update(genetic_cause_known=True),
       "WHITE_FORM_NAME_WAS_FALSELY_PROMOTED_TO_REEVOLUTION"),
      (lambda x:x["evidence_gates"].update(Q1_memory_slope_estimated=True),
       "GENOMIC_OR_COEVOLUTION_GATE_FALSELY_PASSED"),
    ]
    for modify,expected in corruptions:
        bad=copy.deepcopy(j);modify(bad)
        try:validate_source(bad,external)
        except ValueError as e:
            require(expected in str(e),
                    "SOURCE_CORRUPTION_GATE_WRONG_REASON_"+str(e))
        else:raise AssertionError("IMPROPER_SOURCE_CLAIM_WAS_ACCEPTED")
    return {
      "status":"PASS_SYNTHETIC_BAD_ACCESSION_AND_FALSE_REGAIN_CONTROLS",
      "rejected_invalid_accession_claims":len(bad_claims),
      "rejected_source_manifest_fabrications":len(corruptions),
      "rejected_taxonomic_form_as_regain":True,
      "genuine_new_GenBank_sequences_retrieved":0,
      "phylogenomic_memory_estimates_created":0,
      "independent_Cirsium_regain_events_verified":0,
      "original_source_guard":good["status"],
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-manifest",type=Path,default=Path(
      "data/evidence/cirsium_schantarense_vs_shansiense_primary_accession_audit_v1.json"))
    p.add_argument("--previous-crosswalk",type=Path,default=Path(
      "data/evidence/cirsium_q1_nuclear_occurrence_source_overlap_v1.json"))
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    m=json.loads(a.source_manifest.read_text(encoding="utf-8"))
    earlier=json.loads(a.previous_crosswalk.read_text(encoding="utf-8"))
    receipt={
       "version":"cirsium_schantarense_primary_ITS_plastome_source_integrity_v1",
       "as_of":"2026-10-19",
       "original_public_literature_and_taxonomy_audit":validate_source(m,earlier),
       "synthetic_negative_controls":tests(m,earlier),
       "actual_sequence_data_reanalyzed":None,
       "conclusion_ceiling":"SOURCE_MISATTRIBUTION_REMOVED_BUT_AUTHENTIC_TARGET_ITS_FASTA_AND_MULTILOCUS_TREE_NOT_RECOVERED",
    }
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",
                         encoding="utf-8")
    print(json.dumps({
        "status":receipt["original_public_literature_and_taxonomy_audit"]["status"],
        "misattributed_accession":ACCESSION,
        "true_species":OTHER,
        "true_compartment":PLASTID,
        "accession_target_ITS_recovered":False,
        "Q1_memory_estimated":False,
        "negative_tests":receipt["synthetic_negative_controls"]["status"],
    },indent=2))

if __name__=="__main__":
    main()
