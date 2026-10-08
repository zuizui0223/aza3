#!/usr/bin/env python3
"""Auditable GloBI Cirsium interaction discovery. Descriptive, not ecological efficacy."""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import io
import json
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

API = "https://api.globalbioticinteractions.org/interaction.csv"
SPECIES = re.compile(r"^[A-Z][a-zA-Z-]+ [a-z][a-z.-]+(?: [a-z]+)?$")
CATS = {
    "pollinates": "pollination_claim", "pollinatedby": "pollination_claim",
    "visitsflowersof": "flower_visit", "flowersvisitedby": "flower_visit",
    "visits": "general_visit", "visitedby": "general_visit",
    "eats": "feeding", "eatenby": "feeding",
    "preyson": "predation", "preyeduponby": "predation",
    "kills": "mortality_relation", "killedby": "mortality_relation",
    "parasiteof": "parasitism", "hasparasite": "parasitism",
    "endoparasiteof": "parasitism", "hasendoparasite": "parasitism",
    "ectoparasiteof": "parasitism", "hasectoparasite": "parasitism",
    "parasitoidof": "parasitoid_claim", "hasparasitoid": "parasitoid_claim",
    "hostof": "host_association", "hashost": "host_association",
    "pathogenof": "pathogen_claim", "haspathogen": "pathogen_claim",
    "providesnutrientsfor": "nutritional_association",
    "acquiresnutrientsfrom": "nutritional_association",
    "symbiontof": "symbiotic_association", "mutualistof": "mutualist_claim",
    "commensalistof": "commensal_claim",
    "cooccurswith": "cooccurrence_only",
    "interactswith": "nonspecific_interaction",
    "ecologicallyrelatedto": "nonspecific_interaction",
    "dispersalvectorof": "dispersal_claim",
    "hasdispersalvector": "dispersal_claim",
}
OUT_COLS = [
    "record_key", "focal_taxon", "partner_taxon", "partner_species_rank_candidate",
    "focal_position", "interaction_type", "interaction_category",
    "source_taxon", "target_taxon", "source_taxon_id", "target_taxon_id",
    "study_title", "study_url", "study_doi", "study_citation",
    "study_source_citation", "study_source_doi", "study_source_id", "source_namespace",
    "source_taxon_path", "target_taxon_path",
    "source_specimen_body_part", "target_specimen_body_part",
    "source_specimen_life_stage", "target_specimen_life_stage",
    "source_specimen_occurrence_id", "target_specimen_occurrence_id",
    "event_date", "locality", "latitude", "longitude", "source_queries",
    "head_specific_evidence", "ecological_effectiveness",
]

def field(row, *keys):
    for k in keys:
        v = str(row.get(k, "") or "").strip()
        if v:
            return " ".join(v.split())
    return ""

def normtype(value):
    return "".join(x for x in str(value).lower() if x.isalnum())

def is_cirsium(name, genus, taxon_path):
    parts = [x.strip().casefold() for x in taxon_path.split("|")]
    return ("cirsium" in parts or genus.casefold()=="cirsium"
            or name.casefold()=="cirsium" or name.casefold().startswith("cirsium "))

def parse_row(row, side):
    source = field(row, "source_taxon_name", "sourceTaxonName")
    target = field(row, "target_taxon_name", "targetTaxonName")
    sgenus = field(row, "source_taxon_genus_name", "sourceTaxonGenusName")
    tgenus = field(row, "target_taxon_genus_name", "targetTaxonGenusName")
    spath = field(row, "source_taxon_path", "sourceTaxonPath")
    tpath = field(row, "target_taxon_path", "targetTaxonPath")
    s = is_cirsium(source, sgenus, spath)
    t = is_cirsium(target, tgenus, tpath)
    if not (s or t):
        return None
    if s and t:
        return None   # exclude Cirsium–Cirsium relationships from partner inventory
    source_id = field(row, "source_taxon_external_id", "sourceTaxonExternalId", "source_taxon_id", "sourceTaxonId")
    target_id = field(row, "target_taxon_external_id", "targetTaxonExternalId", "target_taxon_id", "targetTaxonId")
    interaction = field(row, "interaction_type", "interactionTypeName", "interaction_type_name")
    typ = normtype(interaction)
    study_id = field(row, "study_source_id", "studySourceId", "study_external_id", "studyExternalId")
    study_url = field(row, "study_url", "studyUrl")
    citation = field(row, "study_citation", "studyCitation")
    src_citation = field(row, "study_source_citation", "studySourceCitation")
    doi = field(row, "study_doi", "studyDoi")
    dataset_doi = field(row, "study_source_doi", "studySourceDoi")
    source_specimen_id = field(row, "source_specimen_occurrence_id", "sourceSpecimenOccurrenceId")
    target_specimen_id = field(row, "target_specimen_occurrence_id", "targetSpecimenOccurrenceId")
    s_part = field(row, "source_specimen_body_part", "sourceSpecimenBodyPart")
    t_part = field(row, "target_specimen_body_part", "targetSpecimenBodyPart")
    s_stage = field(row, "source_specimen_life_stage", "sourceSpecimenLifeStage")
    t_stage = field(row, "target_specimen_life_stage", "targetSpecimenLifeStage")
    event = field(row, "event_date", "eventDate")
    locality = field(row, "locality")
    lat = field(row, "latitude", "decimalLatitude")
    lon = field(row, "longitude", "decimalLongitude")
    focal, partner = (source, target) if s else (target, source)
    partner = partner.strip()
    if not partner:
        return None
    pos = "source" if s else "target"
    # Identity is occurrence/provenance oriented, not collapsed to one species edge.
    identity = (source_id or source, target_id or target, typ,
                study_id, study_url, doi, dataset_doi, citation, src_citation, 
                source_specimen_id, target_specimen_id, event, locality, lat, lon)
    key = hashlib.sha256("\x1f".join(identity).encode()).hexdigest()[:24]
    namespace = study_id or field(row, "source_namespace", "sourceNamespace")
    return dict(
        record_key=key, focal_taxon=focal, partner_taxon=partner,
        partner_species_rank_candidate=str(bool(SPECIES.match(partner))).lower(),
        focal_position=pos,
        interaction_type=interaction, interaction_category=CATS.get(typ, "other_relation"),
        source_taxon=source, target_taxon=target,
        source_taxon_id=source_id, target_taxon_id=target_id,
        study_title=field(row, "study_title", "studyTitle"),
        study_url=study_url, study_doi=doi, study_citation=citation,
        study_source_citation=src_citation, study_source_doi=dataset_doi,
        study_source_id=study_id, source_namespace=namespace,
        source_taxon_path=spath, target_taxon_path=tpath,
        source_specimen_body_part=s_part, target_specimen_body_part=t_part,
        source_specimen_life_stage=s_stage, target_specimen_life_stage=t_stage,
        source_specimen_occurrence_id=source_specimen_id,
        target_specimen_occurrence_id=target_specimen_id,
        event_date=event, locality=locality,
        latitude=lat, longitude=lon, source_queries=side,
        head_specific_evidence="not_assessed", ecological_effectiveness="not_assessed",
    )

def fetch_csv(side, page_size, offset, timeout, attempts=3):
    param = "sourceTaxon" if side == "source" else "targetTaxon"
    url = API + "?" + urlencode({param:"Cirsium","includeObservations":"true",
                                  "limit":page_size, "offset":offset})
    for attempt in range(attempts):
        try:
            req = Request(url, headers={"User-Agent":"aza3-Cirsium-GloBI-audit/1.0 (research; read-only)",
                                        "Accept":"text/csv"})
            with urlopen(req, timeout=timeout) as response:
                raw = response.read()
                content_type = response.headers.get("Content-Type", "")
            decoded = raw.decode("utf-8-sig", "replace")
            reader = csv.DictReader(io.StringIO(decoded, newline=""))
            header = set(reader.fieldnames or [])
            if not ({"source_taxon_name", "target_taxon_name", "interaction_type"} <= header
                    or {"sourceTaxonName", "targetTaxonName", "interactionTypeName"} <= header):
                raise ValueError("GloBI response lacks expected source/target/relation CSV header")
            rows = list(reader)
            if rows and any(not isinstance(x, dict) for x in rows):
                raise ValueError("CSV parser failed")
            if not rows and not decoded.strip():
                raise ValueError("GloBI returned empty response, not verified no-data CSV")
            return rows, url, hashlib.sha256(raw).hexdigest(), content_type
        except (HTTPError, URLError, OSError, ValueError) as exc:
            if attempt == attempts - 1:
                raise RuntimeError(f"{side} offset={offset}: {type(exc).__name__}: {exc}") from exc
            time.sleep(2 ** (attempt + 1))
    raise AssertionError("unreachable")

def run(args):
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    fetched_at = datetime.now(timezone.utc).isoformat()
    records = {}
    pages = []
    failures = []
    raw_rows = 0
    unmatched = 0
    terminal = {}
    for side in ("source", "target"):
        reached_end = False
        for page in range(args.max_pages):
            offset = page * args.page_size
            try:
                rows, url, digest, content_type = fetch_csv(side, args.page_size,
                                                            offset, args.timeout)
            except RuntimeError as exc:
                failures.append(str(exc))
                break
            raw_rows += len(rows)
            kept = 0
            for row in rows:
                parsed = parse_row(row, side)
                if parsed is None:
                    unmatched += 1
                    continue
                old = records.get(parsed["record_key"])
                if old:
                    old["source_queries"] = ";".join(sorted(set(
                        old["source_queries"].split(";") + [side])))
                else:
                    records[parsed["record_key"]] = parsed
                kept += 1
            pages.append({"side":side,"offset":offset,"limit":args.page_size,
                          "rows":len(rows),"cirsium_partner_rows":kept,
                          "url":url,"sha256":digest,"content_type":content_type})
            print(f"FETCH {side} offset={offset} rows={len(rows)} eligible={kept}", flush=True)
            if len(rows) < args.page_size:
                reached_end = True
                break
            time.sleep(args.pause)
        terminal[side] = "EXHAUSTED" if reached_end else (
            "ERROR" if any(f.startswith(side + " ") for f in failures) else "PAGE_CAP")
    fields = OUT_COLS
    datafile = output / "cirsium_globi_observations.csv"
    with datafile.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(sorted(records.values(),
                                key=lambda r:(r["interaction_category"],r["partner_taxon"],
                                              r["focal_taxon"],r["record_key"])))
    (output / "pages.json").write_text(json.dumps(pages,indent=2)+"\n")
    per_type = collections.Counter(r["interaction_type"] for r in records.values())
    per_cat = collections.Counter(r["interaction_category"] for r in records.values())
    names = collections.defaultdict(set)
    species_names = collections.defaultdict(set)
    sources = collections.defaultdict(set)
    species_frequency = collections.defaultdict(collections.Counter)
    for r in records.values():
        cat = r["interaction_category"]
        names[cat].add(r["partner_taxon"])
        if r["partner_species_rank_candidate"] == "true":
            species_names[cat].add(r["partner_taxon"])
            species_frequency[cat][r["partner_taxon"]] += 1
        marker = r["source_namespace"] or r["study_source_citation"] or r["study_url"] or r["study_doi"]
        if marker:
            sources[cat].add(marker)
    groups = []
    for cat in sorted(per_cat):
        groups.append({
            "category":cat,
            "deduplicated_observation_rows":per_cat[cat],
            "distinct_partner_taxon_strings":len(names[cat]),
            "species_rank_name_candidates":len(species_names[cat]),
            "source_identifiers":len(sources[cat]),
            "leading_partner_species_by_rows":species_frequency[cat].most_common(20),
        })
    # Distinct relation + focal name + partner name + dataset provenance = cautious edge,
    # not replicated studies or biological effect sizes.
    edges = {}
    for r in records.values():
        signature = (r["focal_taxon"],r["partner_taxon"],r["interaction_type"],
                     r["source_namespace"],r["study_source_citation"] or r["study_url"])
        if signature not in edges:
            edges[signature] = {k:r[k] for k in fields}
    edgefile = output / "cirsium_partner_relation_sources.csv"
    with edgefile.open("w",encoding="utf-8",newline="") as f:
        writer = csv.DictWriter(f,fieldnames=fields)
        writer.writeheader()
        writer.writerows(sorted(edges.values(), key=lambda x:(x["interaction_category"],
                               x["partner_taxon"],x["interaction_type"],x["focal_taxon"])))
    status = "COMPLETE_DISCOVERY_QUERY" if set(terminal.values())=={"EXHAUSTED"} else "PARTIAL_NOT_EXHAUSTIVE"
    if not records or failures:
        status = "RETRIEVAL_INCOMPLETE_DO_NOT_INFER_ABSENCE"
    summary = {
        "version":"cirsium_globi_bidirectional_observation_audit_v1",
        "retrieved_at_utc":fetched_at,"source":"GloBI unstable live API",
        "source_doi_recommended":"10.5281/zenodo.3950589",
        "query_taxon":"Cirsium", "both_directions_queried":True,
        "terminal_by_query":terminal, "status":status,
        "query_errors":failures,
        "query_pages":len(pages), "raw_rows_retrieved":raw_rows,
        "nonfocal_or_self_or_unnamed_rows":unmatched,
        "deduplicated_observation_rows":len(records),
        "distinct_partner_taxon_strings":len({r["partner_taxon"] for r in records.values()}),
        "partner_species_rank_name_candidates":len({
            r["partner_taxon"] for r in records.values()
            if r["partner_species_rank_candidate"]=="true"}),
        "distinct_focal_taxon_strings":len({r["focal_taxon"] for r in records.values()}),
        "unique_focal_partner_relation_provenance_edges":len(edges),
        "relation_type_counts":dict(per_type.most_common()),
        "categories":groups,
        "csv_sha256":hashlib.sha256(datafile.read_bytes()).hexdigest(),
        "edge_csv_sha256":hashlib.sha256(edgefile.read_bytes()).hexdigest(),
        "hard_limits":[
            "The two directions can expose reciprocal exports of the same claim; deduplication is a heuristic, not guaranteed specimen-level independence.",
            "GloBI live API is mutable, not the recommended frozen publication snapshot.",
            "Distinct partner species strings are not verified accepted taxa or global species richness.",
            "A flower visitor is not necessarily a pollinator; a host association is not head florivory.",
            "An interacting species does not imply a head-specific interaction, selection, defence, or adaptation.",
            "Absence of records is not evidence of ecological absence.",
            "Only records explicitly marked parasitoidOf/hasParasitoid support a parasitoid-type claim, not the impact on Cirsium seed success.",
            "Species interaction counts are strongly confounded by study effort and provenance.",
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    print("BEGIN_CIRSIUM_GLOBI_SUMMARY")
    print(json.dumps({k:summary[k] for k in (
        "status","terminal_by_query","query_errors","raw_rows_retrieved",
        "deduplicated_observation_rows","distinct_partner_taxon_strings",
        "partner_species_rank_name_candidates","unique_focal_partner_relation_provenance_edges",
        "relation_type_counts","categories")},ensure_ascii=False,indent=2))
    print("END_CIRSIUM_GLOBI_SUMMARY",flush=True)
    if failures or not records:
        raise SystemExit(2)
    if status != "COMPLETE_DISCOVERY_QUERY":
        print("WARNING: query truncated; all counts are lower bounds of retrieved claims.",file=sys.stderr)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output-dir",default="outputs/globi_cirsium_v1")
    p.add_argument("--page-size",type=int,default=1000)
    p.add_argument("--max-pages",type=int,default=32)
    p.add_argument("--pause",type=float,default=0.25)
    p.add_argument("--timeout",type=int,default=60)
    a=p.parse_args()
    if not 1<=a.page_size<=1024 or not 1<=a.max_pages<=500 or a.pause<0:
        p.error("invalid page size, page count, or pause")
    run(a)
if __name__=="__main__":
    main()
