#!/usr/bin/env python3
"""Build the authority-wide Japanese Cirsium candidate universe for v3.

This script is intentionally not a live CI dependency. It is run when a fresh
NMNS snapshot is requested. It stores derived categorical screening states and
hashes rather than redistributing authority catchphrase prose.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
from pathlib import Path

import pandas as pd
import requests

DEFAULT_URL = "https://www.kahaku.go.jp/research/db/botany/azami/list.html?word=all"
INFRA_MARKERS = {"var.", "subsp.", "ssp.", "f."}


def clean(value: object) -> str:
    text = "" if value is None or pd.isna(value) else str(value)
    text = text.replace("\u3000", " ").replace("（", "(").replace("）", ")")
    return re.sub(r"\s+", " ", text).strip()


def normalize_taxon(species: object, infra: object = "") -> str:
    s = clean(species)
    i = clean(infra)
    raw = f"{s} {i}".strip()
    raw = re.sub(r"^C\.\s+", "Cirsium ", raw)
    tokens = raw.split()
    if len(tokens) < 2:
        return raw
    out = [tokens[0], tokens[1].strip(",;()")]
    for idx, token in enumerate(tokens[2:], start=2):
        if token in INFRA_MARKERS and idx + 1 < len(tokens):
            out += [token, tokens[idx + 1].strip(",;()")]
            break
    return " ".join(out)


def classify_orientation(text: str) -> str:
    t = clean(text)
    states: list[str] = []
    if any(k in t for k in ("下向き", "点頭", "懸垂")):
        states.append("downward_or_nodding")
    if "横向き" in t:
        states.append("lateral")
    if any(k in t for k in ("上向き", "直立")):
        states.append("upward_or_erect")
    if any(k in t for k in ("斜め上向き", "斜上")):
        states.append("ascending")
    return "|".join(dict.fromkeys(states)) if states else "unknown_from_index"


def classify_phyllary(text: str) -> str:
    t = clean(text)
    states: list[str] = []
    if "圧着" in t:
        states.append("appressed")
    if "斜上" in t:
        states.append("ascending")
    if any(k in t for k in ("開出", "張り出")):
        states.append("spreading")
    if any(k in t for k in ("反曲", "反り返")):
        states.append("recurved")
    return "|".join(dict.fromkeys(states)) if states else "unknown_from_index"


def classify_stickiness(text: str) -> str:
    t = clean(text)
    if any(k in t for k in ("粘らない", "粘らず", "ほとんど粘らない")):
        return "nonsticky_or_nearly_nonsticky"
    if any(k in t for k in ("著しく粘る", "良く粘る", "よく粘る", "粘る", "粘着")):
        return "sticky"
    return "unknown_from_index"


def find_table(html: str) -> pd.DataFrame:
    for table in pd.read_html(io.StringIO(html)):
        cols = [clean(c) for c in table.columns]
        if "種名" in cols and "キャッチフレーズ" in cols:
            table = table.copy()
            table.columns = cols
            return table
    raise ValueError("NMNS thistle table not found")


def load_moreyra_map(path: str | None) -> set[str]:
    if not path:
        return set()
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    if "paper_taxon_concept" not in frame.columns:
        raise ValueError("Moreyra map lacks paper_taxon_concept")
    concepts = set()
    for value in frame["paper_taxon_concept"]:
        tokens = clean(value).split()
        if len(tokens) >= 2:
            concepts.add(" ".join(tokens[:2]))
    return concepts


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--url", default=DEFAULT_URL)
    p.add_argument("--moreyra-map")
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    response = requests.get(args.url, timeout=60)
    response.raise_for_status()
    table = find_table(response.text)
    moreyra = load_moreyra_map(args.moreyra_map)

    rows = []
    for idx, row in table.iterrows():
        species = clean(row.get("種名", ""))
        infra = clean(row.get("変種名", ""))
        phrase = clean(row.get("キャッチフレーズ", ""))
        concept = normalize_taxon(species, infra)
        binomial = " ".join(concept.split()[:2])
        rows.append({
            "authority_record_id": f"NMNS_{idx+1:03d}",
            "japanese_name": clean(row.get("和名", "")),
            "authority_taxon_concept": concept,
            "species_binomial": binomial,
            "infraspecific_record": any(m in concept.split() for m in INFRA_MARKERS),
            "distribution_summary": clean(row.get("分布", "")),
            "orientation_screen": classify_orientation(phrase),
            "phyllary_screen": classify_phyllary(phrase),
            "stickiness_screen": classify_stickiness(phrase),
            "source_catchphrase_sha256": hashlib.sha256(phrase.encode("utf-8")).hexdigest(),
            "represented_in_moreyra_binomial_screen": binomial in moreyra if moreyra else "UNKNOWN_MAP_NOT_SUPPLIED",
            "screening_claim_boundary": "authority categorical screening only; requires taxonomic reconciliation and individual validation before Chapter 3 inference",
        })

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "nmns_transition_candidate_universe_v3.csv", index=False, encoding="utf-8")
    summary = {
        "contract_version": "nmns_transition_candidate_universe_v3",
        "source_url": args.url,
        "source_html_sha256": hashlib.sha256(response.content).hexdigest(),
        "n_authority_records": int(len(frame)),
        "record_count_is_species_count": False,
        "n_unique_normalized_concepts": int(frame["authority_taxon_concept"].nunique()),
        "n_unique_species_binomials": int(frame["species_binomial"].nunique()),
        "moreyra_map_supplied": bool(args.moreyra_map),
        "claim_boundary": "candidate universe for transition-first screening, not a phylogeny or a set of independent evolutionary origins",
    }
    (out / "nmns_transition_candidate_universe_v3.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
