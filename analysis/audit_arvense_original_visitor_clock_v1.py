#!/usr/bin/env python3
"""Original 2022 Cirsium arvense visitor-bout clock audit for SLK × ANTO.

This reads only two provenance-pinned author CSV files. It does NOT estimate
multi-generation evolution, environmental Markov persistence, pollen delivery,
head architectural fitness, or genomic re-expression competence.

The insect Order field is a taxonomic order, not a measured pollination role.
A bout can contain visits to multiple heads on a plant; two orders in the same
bout are not evidence of simultaneous access to the SAME flower head.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import hashlib
import io
import json
from pathlib import Path

SOURCE_REPO = "Yogita25/Multitrophic_Interactions"
SOURCE_COMMIT = "8a151c62ed069f775a16b9bc5fbdce209cbb8691"
SOURCE_BLOBS = {
    "Pollinators.csv": "1723e96b8fc8ea755832daa786c387c5f3c188a8",
    "Pollinator_Id.csv": "875f2046e94b5deca62220e605ed60c4efe232bb",
}
BOUT_HEADER = (
    "Start", "End", "Site_Code", "Block_Code", "Treatment_Code",
    "Number_openFlowerhead", "Remarks", "Instance_ID", "Instance_Name", "Key",
)
VISITOR_HEADER = (
    "Common_Name", "Visited_Flowers", "Key", "Feeding_Guild",
    "Pollination", "Order", "Family", "Genus", "Species", "Species_Code",
)
ALLOWED_ORDERS = {"Hymenoptera", "Diptera", "Coleoptera", "Colleoptera", "Lepidoptera"}


def require(condition: bool, reason: str):
    if not condition:
        raise ValueError(reason)


def git_blob_sha1(source: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(source)).encode("ascii") +
                        b"\\0" + source).hexdigest()


def parse_source_csv(source: bytes, columns: tuple,
                     expected_git_sha1: str | None = None) -> list[dict]:
    if expected_git_sha1 is not None:
        require(git_blob_sha1(source) == expected_git_sha1,
                "SOURCE_SHA_MISMATCH")
    raw = source.decode("utf-8-sig")
    parser = csv.DictReader(io.StringIO(raw, newline=""))
    require(tuple(parser.fieldnames or ()) == tuple(columns),
            "SOURCE_HEADER_MISMATCH")
    rows = []
    for row in parser:
        require(None not in row and row[None] is None if None in row else True,
                "SOURCE_ROW_WIDTH_MISMATCH")
        require(all(v is not None for v in row.values()),
                "SOURCE_ROW_MISSING_FIELD")
        rows.append({k: v.strip() for k, v in row.items()})
    return rows


def parse_bout_time(text: str) -> dt.datetime:
    for fmt in ("%d/%m/%Y %H:%M", "%d/%m/%y %H:%M"):
        try:
            return dt.datetime.strptime(text.strip(), fmt)
        except (ValueError, TypeError):
            continue
    raise ValueError("SOURCE_DATE_PARSE_FAILURE")


def audit(bout_rows: list[dict], visitor_rows: list[dict],
          require_original_counts: bool = True) -> dict:
    if require_original_counts:
        require(len(bout_rows) == 120 and len(visitor_rows) == 649,
                "SOURCE_ROW_COUNT_DRIFT")
    require(len(bout_rows) > 0, "SOURCE_NO_BOUTS")
    by_key: dict[str, list[dict]] = collections.defaultdict(list)
    bouts = []
    site_days = set()
    site_ids = set()
    sites_seen_dates: dict[tuple, set] = collections.defaultdict(set)
    all_dates = set()
    duplicate_duration = 0.0
    zero_duration_count = 0

    for raw in bout_rows:
        require(all(k in raw for k in BOUT_HEADER), "SOURCE_BOUT_FIELD_MISSING")
        start, end = parse_bout_time(raw["Start"]), parse_bout_time(raw["End"])
        delta = (end - start).total_seconds() / 60
        require(delta >= 0, "SOURCE_NEGATIVE_BOUT_DURATION")
        zero_duration_count += int(delta == 0)
        duplicate_duration += delta
        site = raw["Site_Code"].strip().lower()
        block = raw["Block_Code"].strip().lower()
        treatment = raw["Treatment_Code"].strip().lower()
        require(bool(site and block and treatment and raw["Key"]),
                "SOURCE_REQUIRED_BOUT_ID_MISSING")
        x = {
            "date": start.date().isoformat(), "site": site,
            "block": block, "treatment": treatment, "duration_min": delta,
            "key": raw["Key"],
        }
        bouts.append(x)
        by_key[x["key"]].append(x)
        all_dates.add(x["date"])
        site_days.add((site, x["date"]))
        site_ids.add(site)
        sites_seen_dates[(site, block, treatment)].add(x["date"])

    events: dict[str, dict] = {}
    ambiguous = unmatched = blank_order = 0
    for raw in visitor_rows:
        require(all(k in raw for k in VISITOR_HEADER),
                "SOURCE_VISITOR_FIELD_MISSING")
        order = raw["Order"].strip()
        if not order:
            blank_order += 1
            continue
        require(order in ALLOWED_ORDERS, "SOURCE_UNRECOGNIZED_INSECT_ORDER")
        located = by_key.get(raw["Key"], [])
        if len(located) == 0:
            unmatched += 1
            continue
        if len(located) != 1:
            ambiguous += 1
            continue
        key = raw["Key"]
        event = events.setdefault(key, {
            "site": located[0]["site"], "date": located[0]["date"],
            "counts": collections.Counter(),
        })
        event["counts"][order] += 1

    strata: dict[tuple, dict] = {}
    both = 0
    for event in events.values():
        count = event["counts"]
        hym = count["Hymenoptera"]
        dip = count["Diptera"]
        both += int(hym > 0 and dip > 0)
        key = (event["site"], event["date"])
        item = strata.setdefault(key, {
            "site": event["site"], "date": event["date"],
            "order_record_hym": 0, "order_record_dip": 0,
            "named_bouts": 0, "bouts_with_hym_or_dip": 0,
            "bout_hym_fraction_sum": 0.0,
        })
        item["named_bouts"] += 1
        item["order_record_hym"] += hym
        item["order_record_dip"] += dip
        if hym + dip:
            item["bouts_with_hym_or_dip"] += 1
            item["bout_hym_fraction_sum"] += hym/(hym+dip)

    by_site: dict[str, list[dict]] = collections.defaultdict(list)
    strata_public = []
    for (_, _), item in sorted(strata.items()):
        h, d = item["order_record_hym"], item["order_record_dip"]
        record_share = h/(h+d) if h+d else None
        n = item["bouts_with_hym_or_dip"]
        bout_share = item["bout_hym_fraction_sum"]/n if n else None
        entry = {
            "site": item["site"], "date": item["date"],
            "named_order_bouts": item["named_bouts"],
            "order_records_Hymenoptera": h,
            "order_records_Diptera": d,
            "record_weighted_Hymenoptera_share": record_share,
            "equal_bout_Hymenoptera_share": bout_share,
        }
        strata_public.append(entry)
        by_site[item["site"]].append(entry)

    comparisons = []
    for site, collection in sorted(by_site.items()):
        ordered = sorted(collection, key=lambda x: x["date"])
        if len(ordered) < 2:
            continue
        first, last = ordered[0], ordered[-1]
        if any(first[k] is None or last[k] is None
               for k in ("record_weighted_Hymenoptera_share",
                         "equal_bout_Hymenoptera_share")):
            continue
        rdelta = (last["record_weighted_Hymenoptera_share"] -
                  first["record_weighted_Hymenoptera_share"])
        bdelta = (last["equal_bout_Hymenoptera_share"] -
                  first["equal_bout_Hymenoptera_share"])
        comparisons.append({
            "site": site, "first_date": first["date"], "last_date": last["date"],
            "observed_site_dates": len(ordered),
            "record_weighted_final_minus_first": rdelta,
            "equal_bout_final_minus_first": bdelta,
        })

    return {
        "schema": "aza3_slk_anto_cirsium_original_visitor_clock_v1",
        "source": {
            "study": "Cirsium arvense 2022 original field observer CSV",
            "source_repo": SOURCE_REPO,
            "source_commit": SOURCE_COMMIT,
            "source_study_status": "Karpate_et_al_2025_research_square_preprint",
            "raw_bout_git_blob_sha1": SOURCE_BLOBS["Pollinators.csv"],
            "raw_visitor_git_blob_sha1": SOURCE_BLOBS["Pollinator_Id.csv"],
            "new_plant_samples_collected": 0,
        },
        "source_rows": {
            "bouts": len(bout_rows),
            "visitor_lines": len(visitor_rows),
            "distinct_calendar_dates": len(all_dates),
            "survey_site_count": len(site_ids),
            "site_date_strata": len(site_days),
            "distinct_bout_keys": len(by_key),
            "duplicated_bout_keys": sum(len(v)>1 for v in by_key.values()),
            "nominal_bout_minutes_not_unique_effort": duplicate_duration,
            "zero_duration_bout_rows": zero_duration_count,
            "repeated_site_block_treatment_across_dates":
                sum(len(x)>1 for x in sites_seen_dates.values()),
        },
        "matched_taxonomic_observations": {
            "named_order_visitor_lines": sum(sum(v["counts"].values())
                                             for v in events.values()),
            "named_order_bouts": len(events),
            "mixed_Hymenoptera_Diptera_bouts": both,
            "ambiguous_bout_key_visitor_lines": ambiguous,
            "unmatched_bout_key_visitor_lines": unmatched,
            "blank_order_visitor_lines": blank_order,
            "mixed_bout_fraction": both/len(events) if events else None,
        },
        "site_date_record_composition": strata_public,
        "site_time_comparisons": comparisons,
        "first_last_comparisons": {
            "site_count": len(comparisons),
            "positive_record_weighted_sites":
                sum(x["record_weighted_final_minus_first"] > 0 for x in comparisons),
            "positive_equal_bout_sites":
                sum(x["equal_bout_final_minus_first"] > 0 for x in comparisons),
        },
        "scientific_boundaries": {
            "Markov_environment_transition_probability_identified": False,
            "same_plant_heritable_functional_modules_linked": False,
            "same_head_reproductive_fitness_from_visits_identified": False,
            "genealogy_regain_or_fixation_inferred": False,
        },
        "interpretation": {
            "ecologically_observed": "Hymenoptera and Diptera share short observational bouts in the same plant survey frame",
            "ecological_effect": "Different orders and bout compositions, not measured pollination effectiveness",
            "site_first_last_aggregate_not_inferential_test": True,
            "botanical_population_transitions_not_inferred": True,
            "original_bout_key_ambiguity_excluded_not_imputed": True,
            "within_bout_same_head_coaccess_not_proven": True,
            "dates_are_within_one_2022_flowering_season_not_generations": True,
        },
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--bouts",type=Path,required=True)
    parser.add_argument("--visitors",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    bouts=parse_source_csv(args.bouts.read_bytes(),BOUT_HEADER,
                           SOURCE_BLOBS["Pollinators.csv"])
    visitors=parse_source_csv(args.visitors.read_bytes(),VISITOR_HEADER,
                             SOURCE_BLOBS["Pollinator_Id.csv"])
    outcome=audit(bouts,visitors)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(outcome,ensure_ascii=False,
                                     sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("SOURCE_VERIFIED_ARVENSE_VISITOR_CLOCK",outcome["source_rows"]["bouts"],
          outcome["source_rows"]["visitor_lines"],
          outcome["matched_taxonomic_observations"]["named_order_bouts"])
    return outcome


if __name__=="__main__":
    main()
