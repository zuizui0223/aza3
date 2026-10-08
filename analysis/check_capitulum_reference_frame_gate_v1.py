#!/usr/bin/env python3
"""Orientation reference-frame identifiability for Cirsium head arrivals.

OBSERVATION DESIGN ONLY. Simulations are not field results. Head axis and
insect's *pre-contact* approach ray must be independently calibrated in the
same 3D world coordinates; no single-camera image-vertical extrapolation.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import tempfile
from collections import defaultdict
from pathlib import Path

FIELDS = [
    "record_id", "population_id", "taxon_label", "individual_id",
    "capitulum_id", "observation_bout_id", "approach_episode_id",
    "video_clip_id", "phenological_stage", "head_rank", "pre_entry_guild",
    "head_axis_world_x", "head_axis_world_y", "head_axis_world_z",
    "insect_ray_world_x", "insect_ray_world_y", "insect_ray_world_z",
    "approach_shell_diameters", "geometry_method",
    "calibration_evidence_uri", "episode_evidence_uri",
    "scorer_blinded_to_head_traits",
]
STAGES = {
    "bud", "pre_anthesis", "early_anthesis", "full_anthesis",
    "late_anthesis", "postanthesis", "fruit_development", "mature_head",
}
RANKS = {"terminal", "secondary", "tertiary", "other", "not_determined"}
GUILDS = {
    "legitimate_pollinator_candidate", "seed_feeder_candidate",
    "seed_feeder_parasitoid_candidate", "florivore_candidate",
    "other_known", "unresolved",
}
METHODS = {"synchronized_multiview", "calibrated_stereo_video"}
ELIGIBLE = GUILDS - {"unresolved"}

def _finite_float(value, name):
    try:
        f = float(value)
    except (ValueError, TypeError):
        raise ValueError(f"Not numeric: {name}") from None
    if not math.isfinite(f):
        raise ValueError(f"Not finite: {name}")
    return f

def _unit3(row, prefix):
    v = tuple(_finite_float(row[f"{prefix}_{d}"], f"{prefix}_{d}")
              for d in ("x", "y", "z"))
    norm = math.sqrt(sum(x*x for x in v))
    if not 0.97 <= norm <= 1.03:
        raise ValueError(f"{prefix}: not a calibrated 3D unit vector")
    return tuple(x/norm for x in v)

def _angle_bin(z):
    if z >= 0.5:
        return "upright"
    if z <= -0.5:
        return "pendant"
    return "lateral"

def _average(values):
    return sum(values) / len(values)

def evaluate_rows(rows, min_plants_per_bin=3):
    """Descriptive, pre-matched frame contrasts. No p values, no causal fit."""
    if not rows:
        return {
            "status": "NO_REAL_STEREO_GEOMETRY",
            "n_rows": 0, "n_eligible_arrival_rows": 0,
            "n_matched_strata": 0, "contrasts": None,
            "claim_ceiling": "NOT_IDENTIFIABLE",
        }
    if min_plants_per_bin < 2:
        raise ValueError("Need multiple independent plants per orientation bin")
    record_ids, episode_ids = set(), set()
    groups = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    unblinded = 0
    for row in rows:
        if set(row) != set(FIELDS):
            raise ValueError("Row column set differs from contract")
        if any(not row[k] for k in FIELDS):
            raise ValueError("Blank geometry/candidate evidence is not a zero")
        if row["record_id"] in record_ids:
            raise ValueError("Duplicate record_id")
        record_ids.add(row["record_id"])
        episode_key = (
            row["population_id"], row["individual_id"],
            row["capitulum_id"], row["observation_bout_id"],
            row["approach_episode_id"],
        )
        if episode_key in episode_ids:
            raise ValueError("More than one pre-contact ray per approach episode")
        episode_ids.add(episode_key)
        if row["phenological_stage"] not in STAGES:
            raise ValueError("Unknown reproductive stage")
        if row["head_rank"] not in RANKS:
            raise ValueError("Unknown head rank")
        if row["pre_entry_guild"] not in GUILDS:
            raise ValueError("Unknown independently assigned guild")
        if row["geometry_method"] not in METHODS:
            raise ValueError("3D inference from one photograph/uncalibrated camera invalid")
        if row["scorer_blinded_to_head_traits"] not in {"0", "1"}:
            raise ValueError("Blinding must be explicitly 0 or 1")
        unblinded += int(row["scorer_blinded_to_head_traits"] == "0")
        shell = _finite_float(row["approach_shell_diameters"], "shell")
        if not 1 <= shell <= 3:
            raise ValueError("Shell distance must be 1–3 head diameters")
        # This version fixes the approach shell to avoid a post-hoc radius
        # linked to insect guild, head stage or success.
        if abs(shell - 1.5) > 1e-6:
            raise ValueError("v1 requires pre-contact ray at exactly 1.5 head diameters")
        n = _unit3(row, "head_axis_world")
        ray = _unit3(row, "insect_ray_world")
        world_up = ray[2]  # positive world z is gravity-opposing up
        head_front = sum(a*b for a,b in zip(n, ray))
        head_orientation = _angle_bin(n[2])
        if row["pre_entry_guild"] == "unresolved":
            continue
        stratum = (
            row["population_id"], row["taxon_label"],
            row["phenological_stage"], row["head_rank"], row["pre_entry_guild"],
        )
        # Repeated approaches of a plant are explicitly pooled; they do not
        # become independent replicates.
        groups[stratum][head_orientation][row["individual_id"]].append(
            (world_up, head_front)
        )
    contrasts = []
    n_eligible_rows = sum(
        len(observations)
        for bins in groups.values()
        for plants in bins.values()
        for observations in plants.values()
    )
    for stratum, bins in groups.items():
        valid_bins = {
            k: plants for k, plants in bins.items()
            if len(plants) >= min_plants_per_bin
        }
        if len(valid_bins) < 2:
            continue
        for orient_a, orient_b in [
            ("upright", "lateral"), ("upright", "pendant"),
            ("lateral", "pendant"),
        ]:
            if orient_a not in valid_bins or orient_b not in valid_bins:
                continue
            def plant_means(orientation):
                return [
                    (_average([x[0] for x in records]),
                     _average([x[1] for x in records]))
                    for records in valid_bins[orientation].values()
                ]
            a, b = plant_means(orient_a), plant_means(orient_b)
            world_a, head_a = _average([x[0] for x in a]), _average([x[1] for x in a])
            world_b, head_b = _average([x[0] for x in b]), _average([x[1] for x in b])
            contrasts.append({
                "stratum": list(stratum), "orientation_bins": [orient_a, orient_b],
                "n_independent_plants_by_bin": [len(a), len(b)],
                "n_unique_plants_across_bins": len(
                    set(valid_bins[orient_a]) | set(valid_bins[orient_b])
                ),
                "world_vertical_mean_difference_a_minus_b": world_a-world_b,
                "head_front_mean_difference_a_minus_b": head_a-head_b,
                "interpretation": "DESCRIPTIVE_ONLY; requires head/camera/effort join, exposure and phenotype controls",
            })
    return {
        "status": "MATCHED_DESCRIPTIVE_FRAME_CONTRAST_ONLY"
                 if contrasts else "GEOMETRY_PRESENT_BUT_FRAME_NOT_COMPARABLE",
        "n_rows": len(rows),
        "n_eligible_arrival_rows": n_eligible_rows,
        "n_matched_strata": len(set(tuple(c["stratum"]) for c in contrasts)),
        "unblinded_rows": unblinded,
        "contrasts": contrasts if contrasts else None,
        "claim_ceiling": "NO_CAUSAL_OR_ADAPTIVE_INFERENCE_FROM_OBSERVATIONAL_FRAMES",
    }

def evaluate(path):
    with path.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != FIELDS:
            raise ValueError("Exact ordered CSV header required")
        return evaluate_rows(list(reader))

def synthetic_tests():
    def make(i, orientation, insect_ray, guild="seed_feeder_candidate",
             stage="full_anthesis"):
        n = (0., 0., 1.) if orientation == "upright" else (1., 0., 0.)
        return dict(zip(FIELDS, [
            str(i), "p1", "Cirsium synthetic", f"plant{i}",
            f"head{i}", f"bout{i}", f"e{i}", f"clip{i}",
            stage, "terminal", guild, *map(str,n),
            *map(str,insect_ray), "1.5", "calibrated_stereo_video",
            "synthetic://calibration", "synthetic://episode", "1",
        ]))
    a = [
        make(i, "upright" if i<3 else "lateral", (0,0,1))
        for i in range(6)
    ]
    b = [
        make(i, "upright" if i<3 else "lateral",
             (0,0,1) if i<3 else (1,0,0))
        for i in range(6)
    ]
    x = evaluate_rows(a)
    y = evaluate_rows(b)
    assert x["n_matched_strata"] == y["n_matched_strata"] == 1
    assert abs(x["contrasts"][0]["world_vertical_mean_difference_a_minus_b"]) < 1e-12
    assert abs(x["contrasts"][0]["head_front_mean_difference_a_minus_b"]-1) < 1e-12
    assert abs(y["contrasts"][0]["world_vertical_mean_difference_a_minus_b"]-1) < 1e-12
    assert abs(y["contrasts"][0]["head_front_mean_difference_a_minus_b"]) < 1e-12
    assert evaluate_rows(a[:3])["contrasts"] is None
    assert evaluate_rows([make(i, "upright", (0,0,1)) for i in range(6)])["contrasts"] is None
    assert evaluate_rows([make(i, "upright" if i<3 else "lateral",
                          (0,0,1), guild="unresolved") for i in range(6)])["contrasts"] is None
    assert evaluate_rows([])["contrasts"] is None
    bad = [
        {**a[0], "head_axis_world_z": "4"},
        {**a[0], "geometry_method": "single_photo_guess"},
        {**a[0], "calibration_evidence_uri": ""},
        {**a[0], "approach_shell_diameters": "2.0"},
    ]
    for invalid in bad:
        try:
            evaluate_rows([invalid])
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid geometry admitted")
    with tempfile.TemporaryDirectory() as d:
        fp = Path(d)/"empty.csv"
        fp.write_text(",".join(FIELDS)+"\n", encoding="utf-8")
        assert evaluate(fp)["status"] == "NO_REAL_STEREO_GEOMETRY"
    return {
        "status": "PASS_SYNTHETIC_REFERENCE_FRAME_GATES",
        "world_anchored_example": x["contrasts"][0],
        "head_anchored_example": y["contrasts"][0],
        "real_observations": 0,
        "not_phylogenetic_or_experimental": True,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--intake", type=Path,
        default=Path("data/intake/capitulum_world_head_3d_approach_v1.csv"),
    )
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    result = {
        "version": "capitulum_world_head_frame_gate_v1",
        "date": "2026-10-08",
        "synthetic_tests": synthetic_tests(),
        "real_intake": evaluate(args.intake),
        "requirements": [
            "Approach ray precedes first portal contact and is at fixed 1.5-head-diameter shell",
            "Same 3D world calibration for insect ray and floral disc-normal",
            "Within population, taxon, head stage, rank, independently identified guild",
            "At least three independent plants per contrasted natural orientation bin",
            "Parent event/effort joins and genuine effective fitness results required separately",
            "Stage/rain and floral-display confounding do not vanish with a 3D camera",
        ],
    }
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n",
                            encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()
