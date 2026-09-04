#!/usr/bin/env python3
"""Build the authority-wide Japanese Cirsium candidate universe from NMNS data.js.

The current NMNS list page is populated client-side from data/data.js, so this
builder reads that source directly. It stores only short derived fields and hashes;
source prose is not redistributed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
import requests

DEFAULT_DATA_JS_URL = "https://www.kahaku.go.jp/research/db/botany/azami/data/data.js"
INFRA_MARKERS = {"var.", "subsp.", "ssp.", "f."}


def clean(value: object) -> str:
    text = "" if value is None or pd.isna(value) else str(value)
    text = text.replace("\u3000", " ").replace("（", "(").replace("）", ")")
    return re.sub(r"\s+", " ", text).strip()


def normalize_species_seed(value: object) -> str:
    text = clean(value)
    if text.startswith("Cirsiumi "):
        text = "Cirsium " + text.split(" ", 1)[1]
    return text


def normalize_taxon(species: object, infra: object = "") -> str:
    s = normalize_species_seed(species)
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


def parse_data_js(content: bytes) -> list[dict]:
    text = content.decode("utf-8")
    match = re.search(r"var\s+data\s*=\s*(\[.*\])\s*;?\s*$", text, re.S)
    if not match:
        raise ValueError("NMNS data.js does not contain expected `var data = [...]` payload")
    records = json.loads(match.group(1))
    if not isinstance(records, list):
        raise ValueError("NMNS data.js payload is not a list")
    return records


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
    p.add_argument("--data-js-url", default=DEFAULT_DATA_JS_URL)
    p.add_argument("--moreyra-map")
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    response = requests.get(
        args.data_js_url,
        timeout=60,
        headers={"User-Agent": "Mozilla/5.0 aza3-research-snapshot/1.0"},
    )
    response.raise_for_status()
    records = parse_data_js(response.content)
    moreyra = load_moreyra_map(args.moreyra_map)

    rows = []
    for rec in records:
        seed_raw = clean(rec.get("seed", ""))
        species_seed = normalize_species_seed(seed_raw)
        infra = clean(rec.get("var", ""))
        phrase = clean(rec.get("catch", ""))
        concept = normalize_taxon(species_seed, infra)
        tokens = species_seed.split()
        binomial = " ".join(tokens[:2]) if len(tokens) >= 2 else species_seed
        name = clean(rec.get("name", ""))
        if "新称" in name:
            name_status = "UNPUBLISHED_NEW_NAME"
        elif "仮称" in name:
            name_status = "PROVISIONAL_NAME"
        else:
            name_status = "PUBLISHED_OR_NO_NEW_LABEL"
        rows.append({
            "authority_record_id": f"NMNS_{clean(rec.get('no',''))}",
            "japanese_name": name,
            "source_species_string": seed_raw,
            "authority_taxon_concept": concept,
            "species_binomial": binomial,
            "name_status": name_status,
            "infraspecific_record": bool(infra),
            "taxonomic_block": clean(rec.get("class", "")),
            "distribution_summary": clean(rec.get("dist", "")),
            "orientation_screen": classify_orientation(phrase),
            "phyllary_screen": classify_phyllary(phrase),
            "stickiness_screen": classify_stickiness(phrase),
            "source_catchphrase_sha256": hashlib.sha256(phrase.encode("utf-8")).hexdigest(),
            "represented_in_moreyra_binomial_screen": binomial in moreyra if moreyra else "UNKNOWN_MAP_NOT_SUPPLIED",
            "screening_claim_boundary": "authority categorical screening only; taxonomic reconciliation and individual validation required before Chapter 3 inference",
        })

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "nmns_transition_candidate_universe_v3.csv", index=False, encoding="utf-8")
    summary = {
        "contract_version": "nmns_transition_candidate_universe_v3_data_js",
        "source_url": args.data_js_url,
        "source_data_js_sha256": hashlib.sha256(response.content).hexdigest(),
        "n_authority_records": int(len(frame)),
        "record_count_is_species_count": False,
        "n_unique_species_binomials": int(frame["species_binomial"].nunique()),
        "n_published_or_no_new_label_species_binomials": int(frame.loc[frame["name_status"].eq("PUBLISHED_OR_NO_NEW_LABEL"), "species_binomial"].nunique()),
        "n_unpublished_new_name_species_binomials": int(frame.loc[frame["name_status"].eq("UNPUBLISHED_NEW_NAME"), "species_binomial"].nunique()),
        "n_provisional_name_species_binomials": int(frame.loc[frame["name_status"].eq("PROVISIONAL_NAME"), "species_binomial"].nunique()),
        "moreyra_map_supplied": bool(args.moreyra_map),
        "claim_boundary": "current NMNS species-binomial screening universe, not a final accepted-species census or phylogeny",
    }
    (out / "nmns_transition_candidate_universe_v3.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
