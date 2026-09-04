#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "data" / "planning" / "chapter3_p03_orientation_neighbourhood_registry_v3.csv"


def rows() -> list[dict[str, str]]:
    with REGISTRY.open(encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h))


def validate_registry() -> list[dict[str, str]]:
    r = rows()
    if [x["neighbourhood_id"] for x in r] != [
        "O1_PENDULUM_BRANCH",
        "O2_NORIKURENSE_KAMTSCHATICUM_CLUSTER",
        "O3_YEZOENSE_REGION",
        "O4_OPEN_NONOVERLAPPING_SLOT",
    ]:
        raise AssertionError("P03 provisional neighbourhood registry drift")
    status = [x["status"] for x in r]
    if status != ["BACKBONE_CONTRAST_READY", "BACKBONE_CONTRAST_READY", "TRAIT_GAP_FIRST", "OPEN_DO_NOT_FORCE"]:
        raise AssertionError("P03 readiness was over-promoted")
    if r[3]["field_rad_status"] != "CLOSED":
        raise AssertionError("open fourth neighbourhood was prematurely opened")
    if "JPN36" not in r[3]["why_not_double_count"]:
        raise AssertionError("JPN36 anti-double-count rule missing")
    if "JPN_22" not in r[1]["focal_downward_concepts"] or "JPN_37" not in r[1]["focal_downward_concepts"]:
        raise AssertionError("O2 local downward cluster drift")
    if "JPN_04" not in r[0]["current_upward_bracket"] or "JPN_23" not in r[0]["current_upward_bracket"]:
        raise AssertionError("O1 current upward bracket drift")
    if "JPN_10" not in r[2]["current_gap"] or "JPN_24" not in r[2]["current_gap"]:
        raise AssertionError("O3 trait-gap gate drift")
    return r


def main() -> int:
    if not REGISTRY.exists() or REGISTRY.stat().st_size == 0:
        raise AssertionError("missing P03 v3 neighbourhood registry")
    r = validate_registry()
    print("p03_orientation_neighbourhood_registry_v3_valid=true")
    print(f"provisional_slots={len(r)}")
    print(f"backbone_contrast_ready={sum(x['status']=='BACKBONE_CONTRAST_READY' for x in r)}")
    print("fourth_slot_open=true")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
