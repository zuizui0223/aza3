# Chapter 3 master sample list v11

Status: **AUTHORITATIVE SAMPLING ROUTER; NO COLLECTION AUTHORIZATION IMPLIED**

## Purpose

This document replaces the scattered interpretation of older `38 x 3`, `core190`, `222`, `298`, nominal `120 species`, and calendar-first sampling plans with one current sampling router.

The governing hierarchy is:

> exact 128-species nationwide tree hole fill -> nested P01-P05 focal history tests -> nested population RAD where discriminating -> embedded M01 molecular depth case -> conditional function/fitness

The master list is generated from current authoritative sources only:

- `chapter3_exact_site_freeze_queue_v8.csv` generated from the v7 128-species hole-fill plus v10 conservation screen;
- `chapter3_p02_recommended_sample_manifest_v1.csv`;
- `m01_operational_population_ledger_v1.csv`;
- the conditional P01b source protocol migrated from EAzami PR #114.

Older EAzami sampling designs remain under `legacy/eazami_ch3_sources/` and are provenance, not current sample-count authority.

## The number to use

### Active minimum: **314 unique physical representation sources/plants**

This is not `228 + 32 + 60 = 320`, because focal samples are deliberately nested into the nationwide tree.

| Layer | Raw requirement | Already counted in tree | Additional unique sources | Running total |
|---|---:|---:|---:|---:|
| nationwide tree | 228 | — | 228 | **228** |
| P02 minimum RAD | 32 | 2 | 30 | **258** |
| M01 four discovery populations | 60 | 4 | 56 | **314** |

The six overlaps are:

- one `C. dipsacolepis` P02 individual = the required P6 own tree complement;
- one `C. lineare` P02 individual = the required P6 own tree complement;
- one `C. brevicaule` Okinawa-Honto discovery individual = one P3 tree slot;
- one `C. brevicaule` Amami-Oshima discovery individual = the second P3 tree slot;
- one `C. irumtiense` Miyako discovery individual = one P3 tree slot;
- one `C. irumtiense` Ishigaki discovery individual = the second P3 tree slot.

The NMNS tree ledger currently spells the last species `Cirsium irimtiense`; the M01 ledger uses `Cirsium irumtiense`. The master builder carries an explicit equivalence only for cross-ledger matching and does not silently rewrite either source.

## Recommended active depth

- active minimum unique sources: **314**;
- complete the P02 recommended third-population depth: **330**;
- also collect the 12 M01 discovery reserves: **342**;
- also obtain the optional P7 own phenotype-linked `C. incomptum` complement: **343**.

These are staged counts, not one mandatory field campaign.

## Conditional layers — do not collect by default

The full prospective manifest contains **451 rows** only so every possible currently declared source has one identity slot. It does **not** mean 451 plants should be collected now.

After the active core, conditional rows are:

- M01 E3 expansion: **60**;
- M01 E1 bracket material if suitable existing/collaborator material is unavailable: **18**;
- optional M01 positive controls: **6**;
- P01b JPN36 function pilot: **24**.

These open only after their own preregistered history/material/authorization gates.

## Nationwide tree component

The 228 required tree representation slots remain:

| Priority | Meaning | Species | Required representation slots |
|---|---|---:|---:|
| P1 | zero-Moreyra taxonomic block | 33 | 66 |
| P2 | broken/weak public tip | 6 | 12 |
| P3 | <=25% block coverage | 44 | 88 |
| P4 | <=50% block coverage | 16 | 32 |
| P5 | local gap | 2 | 4 |
| P6 | public singleton + one own complement | 26 | 26 |
| P7 | already >=2 public representatives | 1 | 0 required |

A tree representation slot can be filled by authorized fresh material, minimal authorized tissue, authenticated herbarium/historical material, or remain unresolved depending on conservation status. Therefore 228 is not a destructive-wild-collection target.

## What every fresh focal plant must carry

At minimum, one immutable plant ID must link:

- taxon determination and population/deidentified-locality key;
- voucher or diagnostic image record;
- target-capture and/or RAD tissue as applicable;
- orientation;
- phyllary posture/direct dimensions where assessable;
- stickiness plus gland/exudate documentation where focal;
- flower colour and developmental stage;
- cytotype/genome-size status;
- authorization record key.

For M01 flowering individuals, pigment and floral RNA samples remain material IDs nested under the same biological individual.

## Public-locality rule

The master manifest is intentionally public-safe:

- exact latitude/longitude: **prohibited**;
- rare-species micro-locality directions: **prohibited**;
- public fields: broad range sector, deidentified locality ID, acquisition mode and gate state only;
- exact sites and permits remain in a private field ledger.

## Generated artifacts

Run:

```bash
python analysis/build_core_conservation_screen_v10.py
python analysis/build_exact_site_freeze_queue_v8.py
python analysis/build_chapter3_master_sample_manifest_v11.py
```

Outputs:

- `data/planning/chapter3_master_sample_manifest_v11.csv` — one row per unique prospective physical source/plant;
- `data/planning/chapter3_master_sample_summary_v11.json` — stage counts and de-duplication totals;
- `data/planning/chapter3_master_sample_blocks_v11.csv` — human-auditable block-level source plan.

The builder fails if the full manifest is not 451 unique rows, if the active minimum is not 314, if the six tree/focal overlaps are not found exactly, or if exact-coordinate publication is enabled.

## Operational reading

If planning the next field/material campaign, filter the generated manifest in this order:

1. `ACTIVE_TREE_REQUIRED`;
2. `ACTIVE_P02_MINIMUM`;
3. `ACTIVE_M01_DISCOVERY`;
4. `RECOMMENDED_P02_EXTENSION` if P02 third populations are feasible;
5. `RECOMMENDED_M01_RESERVE`;
6. everything marked `CONDITIONAL_*`, `MATERIAL_GAP_DEPENDENT_*`, or `OPTIONAL_*` only when its gate explicitly opens.

Within `ACTIVE_TREE_REQUIRED`, P1/P2/P3/P4/P5/P6 priority and conservation mode determine which holes should be resolved first. Route efficiency may reorder samples within a priority class but must not redefine the hole.

## Claim boundary

A completed sample bank is not a biological result. Target capture, RAD, molecular assays and functional experiments each retain separate QC and inferential gates. No sample count establishes an independent origin, adaptation, defence, pollinator causation, regain or convergence.
