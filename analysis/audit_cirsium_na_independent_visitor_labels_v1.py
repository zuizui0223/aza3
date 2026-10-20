#!/usr/bin/env python3
"""Verify Cirsium author-coded pollination syndrome vs independent visitor reports.

Strict source claim scope:
  (1) published 2023 Cirsium genus poll.csv calls Wrightii 'Insect';
  (2) Flora of North America reports direct hummingbird visits to Wrightii;
  (3) source bird-coded taxa also receive insects.
Independent taxon-flora evidence was MANUALLY CURATED and source-linked; this
script validates source label/trait identity and scientific claim boundaries.
It does NOT scrape source pages, obtain individual visit events, independently
assay pollen transfer, or estimate effect of orientation/spine on fitness.

Sources: Siniscalchi et al. (2023), DOI 10.1086/724310, GitHub fixed main
commit c925d4af86160064448f794a9715f70b917a2786; source FNA taxon
treatments cited in JSON; C. wrightii USFWS species status evidence.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

SOURCE_REPO="carol-siniscalchi/cirsium"
SOURCE_COMMIT="c925d4af86160064448f794a9715f70b917a2786"
BLOBS={
   "traits.csv":"fa127f0587b676e6e004ea121424d26069586bcb",
   "poll.csv":"641ad617b92857a4db0744f799c7ac66ba483ce7",
}
KNOWN_SERVICE_CEILINGS={
    "UNKNOWN",
    "UNKNOWN_QUANTITATIVE",
    "POLLEN_VECTOR_REPORTED_BUT_NOT_SOURCE_QUANTIFIED",
}
MORPH_FIELDS=("CorollaColor","GlutinousRidge","HeadPosition",
              "InnerPhyllaryColor","OuterPhyllaryPosition")
EXPECTED_TAXA=("Cirsium_wrightii","Cirsium_vinaceum","Cirsium_andersonii",
               "Cirsium_occidentale_candidissimum")


def require(ok,error):
    if not ok:raise ValueError(error)


def blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()


def parse_pinned(path,source_name,idcol):
    src=Path(path).read_bytes()
    require(blob(src)==BLOBS[source_name],"ORIGINAL_GIT_BLOB_SHA_DRIFT")
    rows=list(csv.DictReader(io.StringIO(src.decode("utf-8-sig"))))
    ids=[z[idcol] for z in rows]
    require(len(rows)==63 and len(ids)==len(set(ids)) and all(ids),
            "SOURCE_63_TAXON_COUNT_OR_UNIQUE_IDS_CHANGED")
    return {z[idcol]:z for z in rows}


def claim_limits():
    return {
        "bird_visit_means_bird_effective_pollination":False,
        "actual_same_head_visit_pollen_and_fitness":False,
        "head_architecture_as_mechanical_bird_access_barrier_identified":False,
        "evolutionary_gain_or_reexpression_identified":False,
        "author_bird_insect_code_is_exclusive_observed_visitor_list":False,
        "manual_flora_review_is_complete_63_taxon_external_survey":False,
        "external_observation_dates_and_raw_effort_recovered":False,
    }


def audit_four(catalog,traits,poll,expect_full=True):
    cases=catalog["cases"]
    names=[z["taxon"] for z in cases]
    require(len(names)==len(set(names)),"DUPLICATE_TAXON")
    require(set(names)==set(EXPECTED_TAXA),
            "UNEXPECTED_PURPOSEFUL_SOURCE_CASES")
    require(catalog["broader_complex_warning"][
            "may_not_assign_to_both_original_subspecies_as_independent_field_visits"]
            is True,"ARIZONICUM_AGGREGATE_IS_NOT_SUBSPECIES_REPLICATION")
    count_b=count_i=bc_other=ic_bird=0
    taxa_summaries=[]
    for z in cases:
        name=z["taxon"]
        require(name in poll and name in traits,"ORIGINAL_SOURCES_DISAGREE")
        require(poll[name] in ("Hummingbird","Insect"),
                "ORIGINAL_SOURCES_DISAGREE")
        require(z["author_syndrome"]==poll[name],
                "ORIGINAL_SOURCES_DISAGREE")
        for k in MORPH_FIELDS:
            require(str(z["author_trait_codes"].get(k))==str(traits[name].get(k)),
                    "ORIGINAL_SOURCES_DISAGREE")
        require(z["visitor_bird_observed"] is True and
                z["visitor_insect_observed"] is True,
                "EXTERNAL_VISITOR_PRESENCE_SOURCE_CHANGED")
        require(z["reproductive_service_by_bird"] in KNOWN_SERVICE_CEILINGS,
                "UNSUPPORTED_POLLINATION_EFFECTIVENESS")
        urls=z.get("external_urls",[])
        require(isinstance(urls,list) and any(
            x.startswith("https://www.efloras.org/florataxon.aspx?") for x in urls),
            "MISSING_PRIMARY_FLORA_REFERENCE")
        if poll[name]=="Hummingbird":
            count_b+=1
            bc_other+=int(z["visitor_insect_observed"])
        else:
            count_i+=1
            ic_bird+=int(z["visitor_bird_observed"])
        taxa_summaries.append({
            "taxon":name,
            "author_syndrome":poll[name],
            "author_head_position":traits[name]["HeadPosition"],
            "author_outer_phyllary_position":traits[name]["OuterPhyllaryPosition"],
            "author_inner_phyllary_colour":traits[name]["InnerPhyllaryColor"],
            "bird_visitor_in_flora":True,
            "insect_visitor_in_flora":True,
            "bird_pollination_service":z["reproductive_service_by_bird"],
            "source_urls":urls,
        })
    s=catalog["summary"]
    require(s["independent_bird_pollination_efficiency_experiments_in_curated_records"]==0,
            "NO_INDEPENDENT_EFFECTIVENESS_EXPERIMENT")
    require(s["n_taxa_targeted"]==len(cases) and
            s["n_original_bird_coded"]==count_b and
            s["n_original_insect_coded"]==count_i and
            s["n_original_insect_coded_with_external_bird_visits"]==ic_bird and
            s["n_original_bird_coded_with_external_insect_visits"]==bc_other,
            "FROZEN_CURATED_SPOTCHECK_COUNT_CHANGED")
    require(catalog["required_interpretation"][
                "do_not_automatically_recode_original_syndrome_using_visit_record"]
            is True,"MUST_PRESERVE_ORIGINAL_SYNDROME_CODE")
    require(catalog["required_interpretation"][
                "all_63_species_observationally_verified"] is False,
            "FOUR_SPOT_CHECKS_CANNOT_REPRESENT_ALL_63")
    if expect_full:
        require(len(traits)==len(poll)==63 and set(traits)==set(poll),
                "SOURCE_ALL_63_TIP_LABELS_NOT_MATCHED")
        require(s["n_original_source_taxa_total"]==63,
                "SOURCE_63_TAXA_CHANGED")
        require(count_b==3 and count_i==1 and ic_bird==1 and bc_other==3,
                "SOURCE_SPOTCHECK_REPRODUCTION_DRIFT")
    return {
        "n_taxa_spotchecked":len(cases),
        "n_original_bird_coded":count_b,
        "n_original_insect_coded":count_i,
        "n_insect_coded_with_external_hummingbird_visits":ic_bird,
        "n_bird_coded_with_external_insect_visits":bc_other,
        "no_overlap_in_visitor_records_asserted":False,
        "direct_effective_pollination_confirmed":False,
        "spot_checked_taxa":taxa_summaries,
        "independent_source_scope":"MANUALLY_CURATED_TAXON_FLORA_SUMMARY_NO_RAW_VISIT_IDENTIFIERS_OR_REPLICATE_EFFORT",
        "conclusion":"Bird/insect syndrome codes are non-exclusive with respect to visiting animal guilds. Bird service or mechanical access fitness is not identified.",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--poll",type=Path,required=True)
    p.add_argument("--traits",type=Path,required=True)
    p.add_argument("--catalog",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    q=p.parse_args()
    poll_rows=parse_pinned(q.poll,"poll.csv","species")
    trait_rows=parse_pinned(q.traits,"traits.csv","Species")
    labels={k:z["pollinator"] for k,z in poll_rows.items()}
    record=json.loads(q.catalog.read_text(encoding="utf-8"))
    assert record["classification_source"]["commit"]==SOURCE_COMMIT
    assert record["classification_source"]["poll_csv_git_blob"]==BLOBS["poll.csv"]
    assert record["classification_source"]["traits_csv_git_blob"]==BLOBS["traits.csv"]
    case=audit_four(record,trait_rows,labels,expect_full=True)
    result={
        "version":"cirsium_na_independent_visitor_syndrome_code_audit_v1",
        "original_author_repo":SOURCE_REPO,
        "original_author_commit":SOURCE_COMMIT,
        "original_author_blob_hashes":BLOBS,
        "n_original_author_classified_taxa":len(labels),
        "n_original_author_hummingbird_syndrome":sum(
            x=="Hummingbird" for x in labels.values()),
        "case_evidence":case,
        "claim_limits":claim_limits(),
        "not_independent_effectiveness_measurement":True,
    }
    q.out.parent.mkdir(parents=True,exist_ok=True)
    q.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",
                     encoding="utf-8")
    print(json.dumps({
        "n_source_taxa":len(labels),
        "n_checked":case["n_taxa_spotchecked"],
        "source_insect_coded_hummingbird_visited":
            case["n_insect_coded_with_external_hummingbird_visits"],
        "source_bird_coded_insect_visited":
            case["n_bird_coded_with_external_insect_visits"],
        "functional_bird_pollen_assay":False,
    },indent=2))

if __name__=="__main__":
    main()
