#!/usr/bin/env python3
"""Can actual Cirsium capitulum traits and independently observed visitors be linked?

Two ORIGINAL data systems, rather than syndrome label recycling:
(1) original 63 phylogenetic-tip Cirsium capitulum trait panel,
    Siniscalchi et al (2023) DOI 10.1086/724310, author traits.csv,
    commit c925d4af86160064448f794a9715f70b917a2786;
(2) original aza3 2026-10-08 archived GloBI/iNaturalist adult explicit
    visitsFlowersOf source events, raw 104162 and admitted 3243 UNIQUE
    observation IDs, using existing aza3 v1 exact admission.

This reports exact original binomial overlap, independent observational
visitor order records, per-species sampling counts, and *whether* a
between-species original-head-morphology contrast has adequate comparable
observation units in same coarse region. It NEVER infers legitimate pollen
delivery, plant-specific fitness, physical spine contact, true head-level
traits, adaptations, heritability or retained ancestral toolkit.

Missing native head-visitor links are missing, not zero. Source 63 tips
include named infraspecific lineages; truncate NO subspecies to species
because this would fabricate a trait observation. GloBI names are original
names, not independently accepted/vouchered taxon resolutions.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import hashlib
import io
import itertools
import json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

from cirsium_flower_partner_regional_turnover_v1 import (
    ARCHIVE_SHA256, OBSERVATION_MEMBER, prepare_events,
)

ORIGINAL_AUTHOR_REPO="carol-siniscalchi/cirsium"
ORIGINAL_AUTHOR_COMMIT="c925d4af86160064448f794a9715f70b917a2786"
ORIGINAL_TRAITS_BLOB="fa127f0587b676e6e004ea121424d26069586bcb"
ORIGINAL_EVENT_CSV_SHA256="6d2732565c4f158a0f531b4671832738dafe0b06af7f25b635a5d26b7804deff"
ORDERS=("Lepidoptera","Hymenoptera","Coleoptera","Diptera",
        "Hemiptera","Orthoptera","Thysanoptera")
TRAIT_AXES={
    "head_pending":"HeadPosition",
    "outer_phyllaries_outward":"OuterPhyllaryPosition",
    "sticky_ridge":"GlutinousRidge",
}
MIN_GROUPS_FOR_TAXON_CONTRAST=3
LAT_RANGE=(18.,72.)
LON_RANGE=(-170.,-50.)

def require(ok,msg):
    if not ok:raise ValueError(msg)

def blob(b):
    return hashlib.sha1(b"blob "+str(len(b)).encode("ascii")+b"\0"+b).hexdigest()

def trait_panel(raw):
    require(blob(raw)==ORIGINAL_TRAITS_BLOB,
            "AUTHOR_NA_63_TRAIT_SOURCE_GIT_BLOB_CHANGED")
    panel=pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False)
    require(len(panel)==63 and panel.Species.is_unique,
            "AUTHOR_63_PHYLOGENETIC_TIPS_UNIQUE_DRIFT")
    require({"Species",*TRAIT_AXES.values()}.issubset(panel),
            "AUTHOR_HEAD_TRAIT_COLUMNS_UNAVAILABLE")
    for field,vals in (("HeadPosition",{"0","1"}),
                       ("OuterPhyllaryPosition",{"0","1","2"}),
                       ("GlutinousRidge",{"0","1"})):
        require(set(panel[field]).issubset(vals),
                "UNSUPPORTED_OR_UNKNOWN_SOURCE_MORPHOLOGY_CODE_"+field)
    # Strict source-name identity. 'Cirsium_x_variety' must NOT be truncated
    # to 'Cirsium x' without an independent voucher-linked crosswalk.
    strict=panel[panel.Species.str.fullmatch(r"Cirsium_[a-z]+")].copy()
    strict["source_binomial"]=strict.Species.str.replace("_"," ",regex=False)
    require(strict.source_binomial.is_unique,
            "SOURCE_SPECIES_RENAMING_COLLISION")
    return panel,strict

def original_events(archive):
    require(hashlib.sha256(archive).hexdigest()==ARCHIVE_SHA256,
            "GLOBI_ORIGINAL_GIT_ACTIONS_SOURCE_SHA256_CHANGED")
    with ZipFile(io.BytesIO(archive)) as z:
        meta=json.loads(z.read("summary.json"))
        source=z.read(OBSERVATION_MEMBER)
    require(hashlib.sha256(source).hexdigest()==ORIGINAL_EVENT_CSV_SHA256,
            "GLOBI_UNMODIFIED_OBSERVATION_SOURCE_SHA256_CHANGED")
    require(meta["csv_sha256"]==ORIGINAL_EVENT_CSV_SHA256,
            "GLOBI_SOURCE_MANIFEST_CSV_DRIFT")
    raw=pd.read_csv(io.BytesIO(source),dtype=str,keep_default_na=False)
    events,audit=prepare_events(raw)
    require(len(raw)==104162 and len(events)==3243 and
            audit["retained_unique_observations"]==3243 and
            audit["retained_rows_without_event_date"]==3243 and
            audit["retained_host_name_strings"]==42,
            "GLOBI_ORIGINAL_20261008_ADMISSION_DRIFT")
    return events,audit

def insect_order(path):
    text=str(path)
    found=[o for o in ORDERS
           if ("|"+o+"|" in "|"+text.replace(" | ","|")
               .replace("| ","|").replace(" |","|")+"|")]
    # Observations without EXACTLY ONE order have no comparable order label.
    return found[0] if len(found)==1 else None

def na_region(x):
    return x.lat.between(*LAT_RANGE) & x.lon.between(*LON_RANGE)

def source_prepared(tax,event):
    events=event.copy()
    events["original_insect_order"]=[
        insect_order(x) for x in events.source_taxon_path
    ]
    events["source_binomial"]=events.focal_taxon
    joined=events[events.source_binomial.isin(
        set(tax.source_binomial))].copy()
    allowed=joined[na_region(joined)].copy()
    return events,joined,allowed

def species_summary(panel,eligible):
    meta={r.source_binomial:r for r in panel.itertuples()}
    rows=[]
    for taxon in panel.source_binomial:
        subset=eligible[eligible.source_binomial.eq(taxon)]
        all_order=Counter(subset.original_insect_order.dropna())
        total=int(len(subset))
        n_classified=int(sum(all_order.values()))
        head=meta[taxon]
        # Distinct 2° cells: sampling opportunity proxy only, NOT effort.
        cells={
            (int(np.floor((float(z.lat)+90)/2)),
             int(np.floor((float(z.lon)+180)/2)))
            for z in subset.itertuples()
        }
        rows.append({
            "source_taxon":taxon,
            "n_source_admitted_adult_flower_visit_events_NA":total,
            "n_source_order_classified_visits":n_classified,
            "n_coarse_2deg_cells":len(cells),
            "source_head_pending":int(head.HeadPosition),
            "source_outer_phyllaries_outward":
                int(head.OuterPhyllaryPosition) in (1,2),
            "source_glutinous_ridge":int(head.GlutinousRidge),
            "source_order_record_counts":{
                x:int(all_order.get(x,0)) for x in ORDERS},
            "recorded_hymenoptera_fraction":
                all_order["Hymenoptera"]/n_classified
                if n_classified else None,
            "recorded_lepidoptera_fraction":
                all_order["Lepidoptera"]/n_classified
                if n_classified else None,
        })
    return rows

def coverage(rows,minimum):
    return [r for r in rows if
            r["n_source_admitted_adult_flower_visit_events_NA"]>=minimum
            and r["n_coarse_2deg_cells"]>=2
            and r["n_source_order_classified_visits"]>=minimum]

def contrasts(species):
    out={}
    for label in (
        "source_head_pending","source_outer_phyllaries_outward",
        "source_glutinous_ridge"
    ):
        negative=[r for r in species if not r[label]]
        positive=[r for r in species if r[label]]
        report={"n_species_state0":len(negative),
                "n_species_state1":len(positive),
                "eligibility_minimum_taxa_per_state":
                    MIN_GROUPS_FOR_TAXON_CONTRAST}
        if min(len(negative),len(positive))<MIN_GROUPS_FOR_TAXON_CONTRAST:
            report.update({
                "status":"NOT_IDENTIFIABLE_INDEPENDENT_MORPHOLOGY_GROUPS",
                "mean_recorded_hymenoptera_fraction_state1_minus_state0":None,
                "mean_recorded_lepidoptera_fraction_state1_minus_state0":None,
            })
        else:
            def avg(group,key):
                return sum(r[key] for r in group)/len(group)
            report.update({
                "status":"DESCRIPTIVE_SPECIES_LEVEL_ONLY",
                "mean_recorded_hymenoptera_fraction_state1_minus_state0":
                    avg(positive,"recorded_hymenoptera_fraction")-
                    avg(negative,"recorded_hymenoptera_fraction"),
                "mean_recorded_lepidoptera_fraction_state1_minus_state0":
                    avg(positive,"recorded_lepidoptera_fraction")-
                    avg(negative,"recorded_lepidoptera_fraction"),
            })
        out[label]=report
    return out

def matched_cells(eligible,source_traits,minimum=5):
    """True original host pairs WITH contrasting morphology in same coarse NA cells.

    No species pair is equated with an independent evolution event. Shared
    region must have >=5 adult visit records for EACH head-state host.
    """
    if eligible.empty:return {"n_shared_region_host_pairs":0,
                              "n_species_pairs_any_region":0,
                              "n_species_pairs_at_least_3_regions":0,
                              "by_trait":{}}
    x=eligible.copy()
    x["cell"]=[
        str(int(np.floor((float(r.lat)+90)/2)))+"_"+
        str(int(np.floor((float(r.lon)+180)/2)))
        for r in x.itertuples()
    ]
    counts=x.groupby(["source_binomial","cell"],sort=True).size()
    good={(s,c) for (s,c),n in counts.items() if n>=minimum}
    by_host=defaultdict(set)
    for s,c in good:by_host[s].add(c)
    states={r.source_binomial:r for r in source_traits.itertuples()}
    groups={}
    for trait in TRAIT_AXES:
        pair_list=[]
        for a,b in itertools.combinations(sorted(by_host),2):
            sa,sb=states[a],states[b]
            def state(row):
                if trait=="head_pending":
                    return int(row.HeadPosition)
                if trait=="sticky_ridge":
                    return int(row.GlutinousRidge)
                return int(row.OuterPhyllaryPosition) in (1,2)
            if state(sa)==state(sb):continue
            shared=by_host[a]&by_host[b]
            if shared:pair_list.append({
                "host_pair":[a,b],
                "n_original_coarse_regions_same_both_hosts":
                    len(shared),
            })
        groups[trait]={
            "n_species_pairs_at_least_one_matched_region":len(pair_list),
            "n_species_pairs_at_least_three_matched_regions":
                sum(p["n_original_coarse_regions_same_both_hosts"]>=3
                    for p in pair_list),
            "candidate_pairs":pair_list,
        }
    union=set()
    for item in groups.values():
        for p in item["candidate_pairs"]:
            union.add(tuple(p["host_pair"]))
    return {
        "minimum_records_per_taxon_coarse_region":minimum,
        "eligible_taxon_coarse_region_units":len(good),
        "n_species_pairs_any_trait_state_contrast":len(union),
        "by_trait":groups,
    }

def analyze(trait_bytes,archive_bytes):
    full,strict=trait_panel(trait_bytes)
    event,audit=original_events(archive_bytes)
    all_events,joined,eligible=source_prepared(strict,event)
    native=species_summary(strict,eligible)
    admit=coverage(native,10)
    counts={threshold:len(coverage(native,threshold))
            for threshold in (2,5,10,20)}
    report={
        "version":"cirsium_native_head_x_independent_recorded_visitor_audit_v1",
        "computed_date":"2026-10-24",
        "source":{
            "original_NA_author_repo":ORIGINAL_AUTHOR_REPO,
            "original_NA_author_commit":ORIGINAL_AUTHOR_COMMIT,
            "original_NA_traits_git_blob":ORIGINAL_TRAITS_BLOB,
            "original_GloBI_archive_sha256":ARCHIVE_SHA256,
            "original_GloBI_observations_csv_sha256":
                ORIGINAL_EVENT_CSV_SHA256,
            "original_GloBI_source_artifact_id":11534639952,
            "original_GloBI_observation_snapshot":"2026-10-08",
        },
        "source_admission":audit,
        "source_panel":{
            "n_original_phylogenetic_tip_labels":len(full),
            "n_exact_unambiguous_binomial_labels":len(strict),
            "n_nonbinomial_subspecies_or_multiword_tip_labels":
                len(full)-len(strict),
            "do_not_collapse_infraspecies_labels_to_a_binomial":True,
            "n_event_host_binomial_strings":int(event.focal_taxon.nunique()),
            "n_original_GloBI_eligible_events":len(event),
            "n_events_joined_exact_binomial_anywhere":len(joined),
            "n_events_joined_exact_binomial_NA_bounds":len(eligible),
            "n_joined_taxa_anywhere":int(joined.source_binomial.nunique()),
            "n_joined_taxa_in_NA_bounds":int(eligible.source_binomial.nunique()),
            "NA_window_degrees":{
                "lat_min":LAT_RANGE[0],"lat_max":LAT_RANGE[1],
                "lon_min":LON_RANGE[0],"lon_max":LON_RANGE[1],
            },
            "n_original_admitted_visits_not_morphologically_linkable":
                len(event)-len(joined),
            "n_source_visits_without_identifiable_insect_order_NA":
                int(eligible.original_insect_order.isna().sum()),
        },
        "species_all_original_native_binomial_tips":native,
        "taxa_original_records_at_least_N_in_2_or_more_regions":
            {str(k):n for k,n in counts.items()},
        "primary_ten_records_per_species_2_regions":{
            "n_species":len(admit),
            "source_taxa":[r["source_taxon"] for r in admit],
            "total_original_admitted_visits_NA":sum(
                r["n_source_admitted_adult_flower_visit_events_NA"]
                for r in admit),
            "trait_distribution":contrasts(admit),
        },
        "actual_coarse_region_matching":matched_cells(
            eligible,strict,minimum=5),
        "limits":{
            "visitor_record_is_not_pollination_success_or_oviposition":True,
            "original_recording_effort_not_measured":True,
            "original_all_retained_gloBI_event_dates_missing":True,
            "species_trait_2013_taxon_mean_not_observed_on_iNat_individual":True,
            "no_taxonomic_synonym_resolution_without_vouchers":True,
            "NA_geographic_bounds_are_heuristic_not_species_native_range":True,
            "source_taxa_are_phylogenetically_dependent":True,
            "no_guild_effect_if_true_morphological_groups_missing":True,
            "no_population_head_angle_or_spine_size_measured":True,
            "no_new_organismal_field_data":0,
        },
    }
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--traits",type=Path,required=True)
    parser.add_argument("--archive",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    result=analyze(a.traits.read_bytes(),a.archive.read_bytes())
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                                allow_nan=False)+"\n")
    # Intentionally emit a second single-line AGGREGATED receipt in CI so
    # reviewers can source-freeze its 56 taxon coverage outcomes without
    # copying sensitive original event IDs, observer IDs or coordinates.
    print("NATIVE_CIRSIUM_EVIDENCE_JSON="+json.dumps(
        result,separators=(",",":"),ensure_ascii=False,allow_nan=False))
    print(json.dumps({
        "panel":result["source_panel"],
        "coverage":result["taxa_original_records_at_least_N_in_2_or_more_regions"],
        "primary":result["primary_ten_records_per_species_2_regions"],
        "matched":result["actual_coarse_region_matching"],
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
