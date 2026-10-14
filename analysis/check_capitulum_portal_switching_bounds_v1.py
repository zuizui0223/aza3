#!/usr/bin/env python3
"""Cirsium capitulum portal switching: restricted-model tests + empty intake audit.

All numbers in self-tests are deliberately synthetic, not field observations.
This file does NOT identify the causal effect of morphology or adaptation.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

FIELDS = [
    "record_id", "individual_id", "population_id", "capitulum_id",
    "observation_bout_id", "approach_episode_id", "video_clip_id",
    "attempt_index", "attempt_time_s", "head_stage", "contact_zone",
    "contact_world", "obstruction_observed", "access_observed",
    "continuous_coverage", "evidence_uri", "scorer_blinded",
]
STAGES = {
    "bud", "pre_anthesis", "early_anthesis", "full_anthesis",
    "late_anthesis", "postanthesis", "fruit_development", "mature_head",
}
ZONES = {
    "involucre_outer", "involucre_gap", "floret_disc",
    "receptacle_base", "stem_entry", "undetermined",
}
WORLD = {"above", "side", "below", "stem_crawl", "undetermined"}
BINARY_MISSING = {"0", "1", "NA"}

def entry_probability(p_outer, outer_access, floret_access, reroute):
    """Two-portal success per identifiable arrival, not seed or pollen fitness."""
    values = (p_outer, outer_access, floret_access, reroute)
    if not all(math.isfinite(v) and 0 <= v <= 1 for v in values):
        raise ValueError("Probabilities must be finite in [0,1]")
    return (p_outer * (outer_access +
                       (1 - outer_access) * reroute * floret_access)
            + (1 - p_outer) * floret_access)

def synthetic_tests():
    # Invariant p,b,s: lowering outer access cannot improve total access.
    combinations = 0
    for p in (0., .2, .8, 1.):
        for b in (0., .4, .95, 1.):
            for s in (0., .6, .8, 1.):
                low = entry_probability(p, .3, b, s)
                high = entry_probability(p, .7, b, s)
                derivative = p * (1 - s * b)
                assert abs(high - low - .4 * derivative) < 1e-12
                assert low <= high + 1e-12
                combinations += 1

    bud_control = entry_probability(1, .7, 0, .8)
    bud_defence = entry_probability(1, .3, 0, .8)
    anth_control = entry_probability(1, .7, .95, .8)
    anth_defence = entry_probability(1, .3, .95, .8)
    assert abs(bud_control - .7) < 1e-12
    assert abs(bud_defence - .3) < 1e-12
    assert abs(anth_control - .928) < 1e-12
    assert abs(anth_defence - .832) < 1e-12
    # Inversion requires a treatment change in another process parameter.
    cue_control = entry_probability(.8, .5, .95, 0.)
    cue_defence = entry_probability(.2, .1, .95, 0.)
    assert abs(cue_control - .59) < 1e-12
    assert abs(cue_defence - .78) < 1e-12
    assert cue_defence > cue_control
    # Complete rescue can entirely hide outer-wall blockage.
    assert entry_probability(1, 0, 1, 1) == 1
    try:
        entry_probability(-.1, .2, .2, .2)
    except ValueError:
        pass
    else:
        raise AssertionError("Illegal probability accepted")
    return {
        "status": "PASS_SYNTHETIC_BOUNDARY_TESTS_NOT_FIELD_DATA",
        "grid_combinations": combinations,
        "fixed_params_monotone": True,
        "bud_p1": {"control": bud_control, "barrier": bud_defence},
        "anthesis_p1": {"control": anth_control, "barrier": anth_defence},
        "initial_contact_reallocation_counterexample":
            {"control": cue_control, "barrier": cue_defence,
             "changed_parameter": "p_outer"},
    }

def audit(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None or reader.fieldnames != FIELDS:
            raise ValueError("Attempt CSV must have exactly specified ordered header")
        rows = list(reader)
    if not rows:
        return {
            "status": "NO_REAL_ATTEMPT_SEQUENCES",
            "n_rows": 0, "n_episodes": 0, "qualified_descriptive_switches": None,
            "ecological_claim": "NOT_IDENTIFIABLE",
        }
    record_ids = set()
    groups = defaultdict(list)
    for row in rows:
        record_id = row["record_id"]
        if not record_id or record_id in record_ids:
            raise ValueError("Duplicate or empty record_id")
        record_ids.add(record_id)
        for k in FIELDS:
            if not row[k]:
                raise ValueError(f"Missing required value {k}")
        if row["head_stage"] not in STAGES or row["contact_zone"] not in ZONES:
            raise ValueError("Unrecognized head stage or contact interface")
        if row["contact_world"] not in WORLD:
            raise ValueError("Invalid world direction")
        if row["obstruction_observed"] not in BINARY_MISSING:
            raise ValueError("Invalid obstruction state")
        if row["access_observed"] not in BINARY_MISSING:
            raise ValueError("Invalid access state")
        if row["continuous_coverage"] not in {"0", "1"} or row["scorer_blinded"] not in {"0", "1"}:
            raise ValueError("Invalid coverage/blinding flag")
        if row["obstruction_observed"] == "1" and row["access_observed"] != "0":
            raise ValueError("Confirmed obstruction is incompatible with successful same-attempt entry")
        idx = int(row["attempt_index"])
        time = float(row["attempt_time_s"])
        if idx < 1 or not math.isfinite(time) or time < 0:
            raise ValueError("Invalid attempt index/time")
        row["_idx"], row["_time"] = idx, time
        episode = (
            row["population_id"], row["individual_id"], row["capitulum_id"],
            row["observation_bout_id"], row["approach_episode_id"],
        )
        groups[episode].append(row)
    qualified = 0
    unknown = 0
    for _, attempts in groups.items():
        attempts.sort(key=lambda r: r["_idx"])
        if [r["_idx"] for r in attempts] != list(range(1, len(attempts)+1)):
            raise ValueError("Nonconsecutive or repeated attempt indices")
        if any(b["_time"] <= a["_time"] for a,b in zip(attempts, attempts[1:])):
            raise ValueError("Attempt order is inconsistent with elapsed time")
        for field in ("head_stage", "video_clip_id"):
            if len({r[field] for r in attempts}) != 1:
                raise ValueError("One episode must retain its stage and video ID")
        first = attempts[0]
        complete = all(r["continuous_coverage"] == "1" for r in attempts)
        any_changed = any(r["contact_zone"] != first["contact_zone"] and
                          r["contact_zone"] != "undetermined"
                          for r in attempts[1:])
        if len(attempts) >= 2 and first["obstruction_observed"] == "1" and any_changed:
            if complete and first["contact_zone"] != "undetermined":
                qualified += 1
            else:
                unknown += 1
    return {
        "status": "DESCRIPTIVE_CONTACT_SEQUENCES_NOT_CAUSAL",
        "n_rows": len(rows), "n_episodes": len(groups),
        "qualified_descriptive_switches": qualified,
        "unknown_due_to_incomplete_coverage": unknown,
        "ecological_claim": "NEEDS_PARENT_LEDGER_JOIN_AND_HEAD_LEVEL_REPLICATION",
    }

def intake_negative_tests():
    """Fail closed on missing evidence, impossible blocking and incomplete routes."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "test.csv"
        def row(index, zone, blocked, accessed, coverage="1"):
            return dict(zip(FIELDS, [
                f"r{index}", "p1", "pop1", "h1", "b1", "e1", "clip1",
                str(index), str(index), "full_anthesis", zone, "side",
                blocked, accessed, coverage, "synthetic://fixture", "1",
            ]))
        def write(rows):
            with path.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=FIELDS)
                writer.writeheader()
                writer.writerows(rows)
        first = row(1, "involucre_outer", "1", "0")
        second = row(2, "floret_disc", "0", "1")
        write([first, second])
        assert audit(path)["qualified_descriptive_switches"] == 1
        write([first, {**second, "continuous_coverage": "0"}])
        r = audit(path)
        assert r["qualified_descriptive_switches"] == 0
        assert r["unknown_due_to_incomplete_coverage"] == 1
        for corrupt in [
            [{**first, "access_observed": "1"}, second],
            [{**first, "evidence_uri": ""}, second],
            [first, {**second, "attempt_index": "3"}],
        ]:
            write(corrupt)
            try:
                audit(path)
            except ValueError:
                pass
            else:
                raise AssertionError("Invalid filmed contact accepted")
        write([])
        assert audit(path)["status"] == "NO_REAL_ATTEMPT_SEQUENCES"
        assert audit(path)["qualified_descriptive_switches"] is None
    return "PASS_EVIDENCE_REVIEW_ROUTE_GAPS_AND_NO_STRUCTURAL_ZERO"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--intake", type=Path,
                    default=Path("data/intake/capitulum_contact_attempt_sequence_v1.csv"))
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    payload = {
        "version": "capitulum_portal_switching_bounds_v1",
        "date": "2026-10-08",
        "model": "S=p[a+(1-a)sb]+(1-p)b; dS/da=p(1-sb)>=0 if p,b,s fixed",
        "synthetic_only": synthetic_tests(),
        "semantics_test": intake_negative_tests(),
        "intake": audit(args.intake),
        "not_estimated": [
            "arrival rates", "causal barrier effects", "seed damage",
            "pollen deposition", "fitness", "natural selection",
            "phylogenetic transition counts",
        ],
    }
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=2) + "\n",
                            encoding="utf-8")
    print(json.dumps(payload, indent=2))

if __name__ == "__main__":
    main()
