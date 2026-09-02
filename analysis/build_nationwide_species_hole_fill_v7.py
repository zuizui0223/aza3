#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
import pandas as pd

REGIONS = [
    ("OGASAWARA", ["小笠原"]),
    ("IZU_ISLANDS", ["伊豆七島", "八丈", "三宅", "御蔵", "利島", "新島", "神津"]),
    ("SADO", ["佐渡"]),
    ("RISHIRI", ["利尻"]),
    ("RYUKYU", ["沖縄", "奄美", "宮古", "八重山", "西表", "石垣", "与那国", "琉球"]),
    ("OSUMI_ISLANDS", ["屋久島", "種子島"]),
    ("TSUSHIMA", ["対馬"]),
    ("HOKKAIDO", ["北海道"]),
    ("TOHOKU", ["東北", "青森", "岩手", "宮城", "秋田", "山形", "福島"]),
    ("KANTO", ["関東", "東京", "神奈川", "千葉", "埼玉", "茨城", "栃木", "群馬"]),
    ("CHUBU", ["中部", "長野", "新潟", "富山", "石川", "福井", "岐阜", "山梨", "静岡", "愛知"]),
    ("KINKI", ["近畿", "滋賀", "京都", "奈良", "大阪", "兵庫", "和歌山", "三重"]),
    ("CHUGOKU", ["中国", "山口", "岡山", "広島", "鳥取", "島根"]),
    ("SHIKOKU", ["四国", "徳島", "香川", "愛媛", "高知"]),
    ("KYUSHU", ["九州", "福岡", "佐賀", "長崎", "熊本", "大分", "宮崎", "鹿児島"]),
]

PRIORITY_RANK = {
    "H1A_EMPTY_TAXONOMIC_BLOCK": 1,
    "H2A_BROKEN_PUBLIC_TIP": 2,
    "H1B_SPARSE_BLOCK_LE25PCT": 3,
    "H1C_PARTIAL_BLOCK_LE50PCT": 4,
    "H1D_LOCAL_GAP_GT50PCT": 5,
    "H3_PUBLIC_SINGLE_ADD_ONE_OWN": 6,
    "H0_PUBLIC_TWO_PLUS_TREE_FILLED": 7,
}


def bucket(text: str) -> str:
    labels = []
    for label, keys in REGIONS:
        if any(key in text for key in keys):
            labels.append(label)
    return "|".join(labels) if labels else "REVIEW_DISTRIBUTION"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--universe", required=True)
    p.add_argument("--crosswalk", required=True)
    p.add_argument("--repair-ledger", required=True)
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    u = pd.read_csv(args.universe, dtype=str, keep_default_na=False)
    c = pd.read_csv(args.crosswalk, dtype=str, keep_default_na=False)
    r = pd.read_csv(args.repair_ledger, dtype=str, keep_default_na=False)

    required_u = {
        "authority_record_id", "japanese_name", "source_species_string", "species_binomial",
        "name_status", "taxonomic_block", "distribution_summary",
    }
    if required_u.difference(u.columns):
        raise ValueError(f"universe schema missing {sorted(required_u.difference(u.columns))}")

    # Aggregate Moreyra biological samples at current NMNS species level. Multiple
    # paper concepts may collapse to one current species (notably JPN01/JPN20).
    current = c[c["current_nmns_species_binomial"].ne("")].copy()
    current["moreyra_public_biological_samples"] = current["moreyra_public_biological_samples"].astype(int)
    public_n = current.groupby("current_nmns_species_binomial")["moreyra_public_biological_samples"].sum().to_dict()
    jpn_map = current.groupby("current_nmns_species_binomial")["member_id"].apply(lambda x: "|".join(x)).to_dict()

    repair_members = set(r["member_id"])
    repair_species = set(current.loc[current["member_id"].isin(repair_members), "current_nmns_species_binomial"])

    # Collapse authority rows to species-binomial screen units without claiming that
    # this equals a final nomenclaturally accepted species list.
    species_rows = []
    for species, g in u.groupby("species_binomial", sort=True):
        statuses = set(g["name_status"])
        if "PUBLISHED_OR_NO_NEW_LABEL" in statuses:
            name_status = "PUBLISHED_CORE_SCREEN"
        elif "PROVISIONAL_NAME" in statuses:
            name_status = "PROVISIONAL_NAME"
        else:
            name_status = "UNPUBLISHED_NEW_NAME"
        species_rows.append({
            "species_binomial": species,
            "source_species_strings": " | ".join(sorted(set(g["source_species_string"]))),
            "japanese_names": " | ".join(sorted(set(g["japanese_name"]))),
            "name_status": name_status,
            "nmns_record_ids": " | ".join(sorted(set(g["authority_record_id"]))),
            "nmns_taxonomic_block": " | ".join(sorted(set(g["taxonomic_block"]))),
            "region_bucket": bucket(" | ".join(sorted(set(g["distribution_summary"])))),
            "n_nmns_records": int(len(g)),
            "moreyra_jpns": jpn_map.get(species, ""),
            "moreyra_public_biological_n": int(public_n.get(species, 0)),
        })
    species = pd.DataFrame(species_rows)

    core = species[species["name_status"].eq("PUBLISHED_CORE_SCREEN")].copy()
    core["moreyra_covered"] = core["moreyra_public_biological_n"].gt(0)

    block_stats = {}
    block_anchors = {}
    for block, g in core.groupby("nmns_taxonomic_block"):
        total = int(len(g))
        covered = int(g["moreyra_covered"].sum())
        block_stats[block] = {
            "total": total,
            "covered": covered,
            "missing": total - covered,
            "ratio": covered / total,
        }
        anchors = []
        for _, row in g[g["moreyra_covered"]].sort_values("species_binomial").iterrows():
            anchors.append(f"{row['moreyra_jpns']}:{row['species_binomial']}")
        block_anchors[block] = " | ".join(anchors)

    def slot(row: pd.Series) -> str:
        sp = row["species_binomial"]
        if sp in repair_species:
            return "H2A_BROKEN_PUBLIC_TIP"
        if not row["moreyra_covered"]:
            stat = block_stats[row["nmns_taxonomic_block"]]
            if stat["covered"] == 0:
                return "H1A_EMPTY_TAXONOMIC_BLOCK"
            if stat["ratio"] <= 0.25:
                return "H1B_SPARSE_BLOCK_LE25PCT"
            if stat["ratio"] <= 0.50:
                return "H1C_PARTIAL_BLOCK_LE50PCT"
            return "H1D_LOCAL_GAP_GT50PCT"
        if int(row["moreyra_public_biological_n"]) >= 2:
            return "H0_PUBLIC_TWO_PLUS_TREE_FILLED"
        return "H3_PUBLIC_SINGLE_ADD_ONE_OWN"

    core["slot_priority"] = core.apply(slot, axis=1)
    core["priority_rank"] = core["slot_priority"].map(PRIORITY_RANK)
    core["new_tree_samples_needed"] = core["slot_priority"].map(
        lambda x: 2 if x.startswith("H1") or x.startswith("H2") else (1 if x.startswith("H3") else 0)
    )
    core["own_trait_link_sample_recommended"] = core["slot_priority"].eq("H0_PUBLIC_TWO_PLUS_TREE_FILLED").astype(int)
    core["block_total_species"] = core["nmns_taxonomic_block"].map(lambda x: block_stats[x]["total"])
    core["block_moreyra_covered_species"] = core["nmns_taxonomic_block"].map(lambda x: block_stats[x]["covered"])
    core["block_missing_species"] = core["nmns_taxonomic_block"].map(lambda x: block_stats[x]["missing"])
    core["same_block_moreyra_anchors"] = core["nmns_taxonomic_block"].map(block_anchors)

    core_cols = [
        "priority_rank", "slot_priority", "species_binomial", "source_species_strings", "japanese_names",
        "nmns_record_ids", "region_bucket", "nmns_taxonomic_block", "block_total_species",
        "block_moreyra_covered_species", "block_missing_species", "same_block_moreyra_anchors",
        "moreyra_jpns", "moreyra_public_biological_n", "new_tree_samples_needed",
        "own_trait_link_sample_recommended", "n_nmns_records",
    ]
    core_out = core[core_cols].sort_values(["priority_rank", "nmns_taxonomic_block", "species_binomial"])

    block_rows = []
    for block, stat in block_stats.items():
        g = core[core["nmns_taxonomic_block"].eq(block)]
        if stat["covered"] == 0:
            bp = "B1_EMPTY_BLOCK"
        elif stat["ratio"] <= 0.25:
            bp = "B2_SPARSE_LE25PCT"
        elif stat["ratio"] <= 0.50:
            bp = "B3_PARTIAL_LE50PCT"
        elif stat["missing"]:
            bp = "B4_LOCAL_GAPS"
        else:
            bp = "B5_COVERED"
        block_rows.append({
            "block_priority": bp,
            "nmns_taxonomic_block": block,
            "published_core_species": stat["total"],
            "moreyra_covered_species": stat["covered"],
            "missing_species": stat["missing"],
            "coverage_fraction": round(stat["ratio"], 3),
            "same_block_moreyra_anchors": block_anchors[block],
            "missing_species_binomials": " | ".join(sorted(g.loc[~g["moreyra_covered"], "species_binomial"])),
            "broken_anchor_species": " | ".join(sorted(g.loc[g["species_binomial"].isin(repair_species), "species_binomial"])),
        })
    blocks = pd.DataFrame(block_rows).sort_values(["missing_species", "published_core_species"], ascending=False)

    extension = species[~species["name_status"].eq("PUBLISHED_CORE_SCREEN")].copy()
    extension["slot_priority"] = extension["name_status"].map({
        "PROVISIONAL_NAME": "P1_PROVISIONAL_NAME",
        "UNPUBLISHED_NEW_NAME": "P2_UNPUBLISHED_NEW_NAME",
    })
    extension["priority_rank"] = extension["slot_priority"].map({"P1_PROVISIONAL_NAME": 1, "P2_UNPUBLISHED_NEW_NAME": 2})
    extension["new_samples_if_extension_admitted"] = extension["moreyra_public_biological_n"].map(lambda n: 1 if int(n) > 0 else 2)
    extension_cols = [
        "priority_rank", "slot_priority", "species_binomial", "source_species_strings", "japanese_names",
        "nmns_record_ids", "region_bucket", "nmns_taxonomic_block", "moreyra_jpns",
        "moreyra_public_biological_n", "new_samples_if_extension_admitted",
    ]
    extension_out = extension[extension_cols].sort_values(["priority_rank", "species_binomial"])

    conflicts = c[c["current_nmns_species_binomial"].eq("")].copy()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    core_out.to_csv(out / "chapter3_nmns_published_core_hole_fill_v7.csv", index=False, encoding="utf-8")
    blocks.to_csv(out / "chapter3_nmns_block_gap_summary_v7.csv", index=False, encoding="utf-8")
    extension_out.to_csv(out / "chapter3_nmns_provisional_extension_v7.csv", index=False, encoding="utf-8")
    conflicts.to_csv(out / "chapter3_moreyra_taxonomic_conflicts_v7.csv", index=False, encoding="utf-8")

    counts = core_out["slot_priority"].value_counts().to_dict()
    summary = {
        "contract_version": "chapter3_nationwide_species_hole_fill_v7",
        "nmns_authority_records": int(len(u)),
        "nmns_unique_species_binomials": int(species["species_binomial"].nunique()),
        "published_core_screen_species": int(len(core_out)),
        "unpublished_new_name_species": int(extension_out["slot_priority"].eq("P2_UNPUBLISHED_NEW_NAME").sum()),
        "provisional_name_species": int(extension_out["slot_priority"].eq("P1_PROVISIONAL_NAME").sum()),
        "moreyra_current_published_species_coverage": int(core["moreyra_covered"].sum()),
        "h1_absent_published_species": int((core_out["slot_priority"].str.startswith("H1")).sum()),
        "h2_broken_current_species": int(core_out["slot_priority"].str.startswith("H2").sum()),
        "h3_public_single_species": int(core_out["slot_priority"].str.startswith("H3").sum()),
        "h0_public_two_plus_species": int(core_out["slot_priority"].str.startswith("H0").sum()),
        "new_tree_samples_for_published_core": int(core_out["new_tree_samples_needed"].sum()),
        "additional_own_trait_link_samples_recommended": int(core_out["own_trait_link_sample_recommended"].sum()),
        "new_samples_if_all_provisional_newname_extension_admitted": int(extension_out["new_samples_if_extension_admitted"].sum()),
        "unresolved_moreyra_paper_concepts_not_in_current_nmns_screen": int(len(conflicts)),
        "slot_counts": counts,
        "claim_boundary": "Current NMNS species-binomial hole-fill screen. Published-core means no 新称/仮称 label in this snapshot, not a final nomenclatural acceptance decision. Collection remains blocked until the operational species census and permits are frozen.",
    }
    (out / "chapter3_nmns_hole_fill_summary_v7.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
