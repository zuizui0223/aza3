#!/usr/bin/env python3
"""Real Cirsium arvense visitor-count × maternal seed-outcome LABEL audit.

Two different original author tables share site x block x treatment labels,
NOT source-verified unique head IDs. Therefore all contrasts are exploratory
SOURCE-LABEL associations. They cannot measure actual structural SLK fitness,
individual insect pollination effectiveness, ANTO genetic-resource retention,
nor environmental Markov transition rates.

Source data are downloaded separately by CI, byte-verified at pinned Git
blob SHA values, and are deliberately not duplicated in this repository.
"""
from __future__ import annotations
import argparse
import collections
import json
import math
from pathlib import Path

from audit_arvense_original_visitor_clock_v1 import (
    BOUT_HEADER, VISITOR_HEADER, SOURCE_BLOBS, parse_source_csv,
    parse_bout_time,
)
from reanalyse_cirsium_arvense_2022_real_head_seed_bookkeeping_v1 import (
    GIT_BLOBS, parse_head_rows, source_csv,
)

SITES = {"campus", "construction", "ditch", "rapseed"}
ARMS = {"control", "predator exclusion", "predator model"}
ALL_SEED_ARMS = ARMS | {"pollinator exclusion"}


def need(ok, msg):
    if not ok:
        raise ValueError(msg)


def normalize_arm(value):
    a = str(value).strip().lower()
    return "predator exclusion" if a == "exclude predators" else a


def label_key(site, block, arm):
    return (str(site).strip().lower(), str(block).strip().lower(),
            normalize_arm(arm))


def pearson(xs, ys):
    if len(xs) != len(ys):
        raise ValueError("CORRELATION_LENGTH_MISMATCH")
    if len(xs) < 4:
        return None
    if any(not math.isfinite(x) for x in xs+ys):
        raise ValueError("CORRELATION_NONFINITE")
    mx, my = sum(xs)/len(xs), sum(ys)/len(ys)
    a = sum((x-mx)**2 for x in xs)
    b = sum((y-my)**2 for y in ys)
    if a < 1e-24 or b < 1e-24:
        return None
    return sum((x-mx)*(y-my) for x,y in zip(xs,ys))/math.sqrt(a*b)


def association(rows, exposure, outcome="viable_fraction"):
    candidates = [r for r in rows if r.get(exposure) is not None
                  and r.get(outcome) is not None]
    return {
        "n": len(candidates),
        "r": pearson([float(r[exposure]) for r in candidates],
                     [float(r[outcome]) for r in candidates]),
        "analysis": "exploratory_unadjusted_pearson_not_causal",
    }


def controlled_correlation(rows, exposure, outcome="viable_fraction"):
    """Demeaned descriptive Pearson by original site × treatment groups.

    This removes additive group means; it does not identify pollinator causal
    benefits, solve selection into visited plants, or control head phenology.
    """
    subsets = collections.defaultdict(list)
    for r in rows:
        if r.get(exposure) is not None and r.get(outcome) is not None:
            subsets[(r["site"],r["treatment"])].append(r)
    xx, yy = [], []
    kept = 0
    for vals in subsets.values():
        if len(vals) < 2:
            continue
        kept += 1
        mx = sum(float(r[exposure]) for r in vals)/len(vals)
        my = sum(float(r[outcome]) for r in vals)/len(vals)
        xx.extend(float(r[exposure])-mx for r in vals)
        yy.extend(float(r[outcome])-my for r in vals)
    return {
        "n": len(xx), "n_strata": kept, "r": pearson(xx,yy),
        "analysis": "exploratory_site_and_treatment_demeaned_pearson_not_causal",
    }


def join_source_labels(heads, bouts, visits, strict=True):
    """Join ORIGINAL experiment source labels, not exact biological heads.

    Accept parsed head rows from the verified original seed reader.
    The site/block/treatment triple is at best an INFERRED maternal group.
    Visitor UUIDs duplicated across bout rows are rejected, not matched to
    the first arbitrary row (the original author's match() convention).

    Only positive-duration bouts with >=1 open head contribute nominal
    exposure effort. Output internal group rows for summary use; never
    distribute their original labels or unique event IDs in public files.
    """
    if strict:
        need((len(heads),len(bouts),len(visits)) == (615,120,649),
             "SOURCE_ROW_COUNT_DRIFT")
    groups = {}
    valid_heads = 0
    for h in heads:
        site = str(h["site"]).strip().lower()
        arm = normalize_arm(h["treatment"])
        if site not in SITES or arm not in ALL_SEED_ARMS:
            continue
        valid_heads += 1
        v = [int(h[n]) for n in ("Viable_Seeds", "NViable_Seeds",
                                  "Damaged_Seeds", "Missing_Seeds", "Flowers")]
        need(v[0]+v[1]+v[2]+v[3]==v[4], "SEED_MASS_BALANCE_MISMATCH")
        need(v[4]>0 and min(v[:3])>=0, "INVALID_OBSERVED_SEED_TOTAL")
        k = label_key(site,h["block"],arm)
        g = groups.setdefault(k,{"site":site,"treatment":arm,
                 "n_heads":0,"flowers":0,"viable":0,"direct_damaged":0,
                 "emerged_insects":0})
        g["n_heads"]+=1
        g["flowers"]+=v[4]
        g["viable"]+=v[0]
        g["direct_damaged"]+=v[2]
        g["emerged_insects"]+=int(h["Insects"])

    if strict:
        need(valid_heads==612 and len(groups)==164,
             "SOURCE_SEED_GROUP_UNIVERSE_CHANGED")

    visitor_groups=collections.defaultdict(list)
    for r in visits:
        visitor_groups[str(r.get("Key","")).strip()].append(r)
    bout_key_rows=collections.defaultdict(list)
    for b in bouts:
        bout_key_rows[str(b["Key"]).strip()].append(b)

    candidates = collections.defaultdict(list)
    for b in bouts:
        k=label_key(b["Site_Code"],b["Block_Code"],b["Treatment_Code"])
        if k in groups and k[2] in ARMS:
            candidates[k].append(b)
    results=[]
    valid=ambiguous=zero_duration=zero_open=0
    for k,source_bouts in candidates.items():
        g=groups[k]
        hcount=dcount=ocount=hlines=dlines=0
        n_bouts=0
        open_head_hours=0.0
        bout_shares=[]
        for b in source_bouts:
            key=str(b["Key"]).strip()
            if len(bout_key_rows[key])!=1:
                ambiguous+=1
                continue
            start,end=parse_bout_time(b["Start"]),parse_bout_time(b["End"])
            duration=(end-start).total_seconds()/60
            need(duration>=0, "NEGATIVE_OBSERVATION_DURATION")
            try:
                open_heads=float(b["Number_openFlowerhead"])
            except (ValueError,TypeError):
                raise ValueError("NONNUMERIC_OPEN_HEAD_COUNT") from None
            need(math.isfinite(open_heads) and open_heads>=0,
                 "INVALID_OPEN_HEAD_COUNT")
            if duration==0:
                zero_duration+=1
                continue
            if open_heads==0:
                zero_open+=1
                continue
            h=d=o=nh=nd=0
            for v in visitor_groups.get(key,()):
                order=str(v.get("Order","")).strip()
                raw_visit=str(v.get("Visited_Flowers","")).strip()
                if not raw_visit:
                    continue
                try:
                    value=float(raw_visit)
                except ValueError:
                    raise ValueError("VISITOR_COUNT_NONNUMERIC") from None
                need(math.isfinite(value) and value>=0,
                     "VISITOR_NEGATIVE_COUNT")
                if not order:
                    continue
                if order=="Hymenoptera":
                    h+=value;nh+=1
                elif order=="Diptera":
                    d+=value;nd+=1
                else:
                    o+=value
            valid+=1
            n_bouts+=1
            open_head_hours+=open_heads*duration/60
            hcount+=h;dcount+=d;ocount+=o;hlines+=nh;dlines+=nd
            if h+d>0:
                bout_shares.append(h/(h+d))
        if n_bouts==0:
            continue
        need(open_head_hours>0, "INVALID_NOMINAL_EXPOSURE_SUM")
        combined=hcount+dcount
        results.append({
            "site":g["site"], "treatment":g["treatment"],
            "n_heads":g["n_heads"], "n_bouts":n_bouts,
            "viable_fraction":g["viable"]/g["flowers"],
            "direct_damaged_fraction":g["direct_damaged"]/g["flowers"],
            "emerged_insects_per_head":g["emerged_insects"]/g["n_heads"],
            "hym_count":hcount, "dip_count":dcount, "other_order_count":ocount,
            "hym_visit_share":hcount/combined if combined else None,
            "hym_line_share":hlines/(hlines+dlines) if hlines+dlines else None,
            "equal_bout_hym_share":
                sum(bout_shares)/len(bout_shares) if bout_shares else None,
            "rate_visit_per_nominal_head_hour":combined/open_head_hours,
            "rate_hym_per_nominal_head_hour":hcount/open_head_hours,
        })
    return {
        "group_records":results,   # internal only: no source identity keys
        "source_join":{
            "original_bout_rows":len(bouts),
            "original_visitor_rows":len(visits),
            "original_seed_head_rows":len(heads),
            "original_seed_label_groups_all_four_arms":len(groups),
            "matched_label_groups_before_effort_filter":len(candidates),
            "matched_label_groups_after_effort_filter":len(results),
            "valid_nonzero_head_bouts":valid,
            "excluded_ambiguous_key_bouts":ambiguous,
            "excluded_zero_duration_bouts":zero_duration,
            "excluded_zero_open_head_bouts":zero_open,
            "group_label_join_not_verified_same_head":True,
            "source_bout_key_duplicate_groups":
                sum(len(r)>1 for r in bout_key_rows.values()),
        },
    }


def summarize_source_biological_join(heads,bouts,visitors,strict=True):
    linked=join_source_labels(heads,bouts,visitors,strict=strict)
    vals=linked["group_records"]
    sites=sorted({v["site"] for v in vals})
    arms=sorted({v["treatment"] for v in vals})
    share="hym_visit_share"
    within=controlled_correlation(vals,share)
    main=association(vals,share)
    site_only={site:association([r for r in vals if r["site"]==site],share)
               for site in sites}
    leave_site={
        site:controlled_correlation([r for r in vals if r["site"]!=site],share)
        for site in sites
    }
    groups_one=[r for r in vals if r["n_bouts"]==1]
    groups_repeated=[r for r in vals if r["n_bouts"]>=2]
    proxy="rate_visit_per_nominal_head_hour"
    sorted_rates=sorted(r[proxy] for r in vals)
    cutoff=sorted_rates[min(len(sorted_rates)-1,int(.95*len(sorted_rates)))]
    no_extreme=[r for r in vals if r[proxy]<=cutoff]
    return {
        "schema":"aza3_slk_anto_source_bout_seed_group_20261020_v1",
        "source":{
            "study":"Cirsium arvense source-authenticated 2022 field observation",
            "author_repo":"Yogita25/Multitrophic_Interactions",
            "source_commit":"8a151c62ed069f775a16b9bc5fbdce209cbb8691",
            "seed_sha1":GIT_BLOBS["seeds.csv"],
            "bout_sha1":SOURCE_BLOBS["Pollinators.csv"],
            "visitor_sha1":SOURCE_BLOBS["Pollinator_Id.csv"],
            "raw_seed_rows":len(heads),"raw_bouts":len(bouts),
            "raw_visitor_rows":len(visitors),
            "source_document":"Karpate et al 2025 Research Square preprint",
        },
        "source_link":linked["source_join"],
        "results":{
            "total_groups_with_valid_nominal_rate":len(vals),
            "groups_with_Hym_Dip_share":main["n"],
            "groups_by_site":{s:sum(r["site"]==s for r in vals) for s in sites},
            "groups_by_treatment":{a:sum(r["treatment"]==a for r in vals)
                                   for a in arms},
            "groups_one_bout":len(groups_one),
            "groups_repeated_bouts":len(groups_repeated),
            "visitor_share_vs_viable_fraction":{
                "unadjusted":main,
                "within_site_and_treatment":within,
                "site_leave_one_out_controlled":leave_site,
                "one_bout_only":association(groups_one,share),
                "repeated_bouts_only":association(groups_repeated,share),
                "alternate_line_count":association(vals,"hym_line_share"),
                "alternate_equal_bout":association(vals,"equal_bout_hym_share"),
                "direct_damaged_outcome":association(vals,share,
                                                     "direct_damaged_fraction"),
            },
            "site_specific_visitor_share_vs_viable_fraction":site_only,
            "nominal_visit_count_rate_vs_viable_fraction":{
                "all":association(vals,proxy),
                "within_site_and_treatment":controlled_correlation(vals,proxy),
                "remove_top_5_percent_rate":association(no_extreme,proxy),
                "remove_top_5_percent_controlled":
                    controlled_correlation(no_extreme,proxy),
                "rate_cutoff":cutoff,
                "rate_is_not_independent_pollinator_service":True,
            },
            "actual_plant_genomic_data_records":0,
        },
        "scientific_boundaries":{
            "experimentally_verified_pollination_effectiveness":False,
            "same_biological_head_contact_and_seed_fate":False,
            "anatomical_SLK_division_of_labour_profit_estimated":False,
            "functional_ANTO_reexpression_competence_measured":False,
            "intergenerational_environment_markov_rate_inferred":False,
            "heritable_genetic_recovery_or_population_fixation_inferred":False,
        },
        "interpretation":{
            "analysis_status":"EXPLORATORY_POST_HOC_SOURCE_LABEL_ASSOCIATION",
            "comparison_unit":"inferred site_x_block_x_treatment mother group, not an identified head",
            "source_visits":"reported visited-flower counts, not pollen effectiveness",
            "survey_denominator":"nominal sum of open heads x bout hours; duplicated and unquantified survey effort may remain",
            "site_specific_sign_heterogeneity_not_causal_adaptive_reversal":True,
            "no_authentic_head_morphology_measured":True,
            "raw_individual_keys_not_exported":True,
            "source_overlap_with_existing_aza3_615_head_analysis":True,
            "source_link_not_independent_replicate_dataset":True,
        },
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--seeds",type=Path,required=True)
    p.add_argument("--bouts",type=Path,required=True)
    p.add_argument("--visitors",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    _,raw_heads=source_csv(a.seeds,"seeds.csv")
    heads=parse_head_rows(raw_heads)
    bout_rows=parse_source_csv(a.bouts.read_bytes(),BOUT_HEADER,
                               SOURCE_BLOBS["Pollinators.csv"])
    visitor_rows=parse_source_csv(a.visitors.read_bytes(),VISITOR_HEADER,
                                  SOURCE_BLOBS["Pollinator_Id.csv"])
    report=summarize_source_biological_join(heads,bout_rows,visitor_rows)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True,
                                    ensure_ascii=False)+"\n",encoding="utf-8")
    print("REAL_ORIGINAL_ARVENSE_VISITOR_SEED_LABEL_LINK",
          report["source_link"]["matched_label_groups_after_effort_filter"],
          report["results"]["groups_with_Hym_Dip_share"],
          report["results"]["visitor_share_vs_viable_fraction"]
                ["within_site_and_treatment"]["r"])
    return report


if __name__=="__main__":
    main()
