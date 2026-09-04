# Chapter 3 exact backbone hole-fill ledger v7

## Purpose

This document converts the Japan-wide species-tree sampling problem into an explicit **backbone hole-fill problem** against the existing Moreyra Japanese target-capture skeleton.

The unit is the current NMNS species-binomial screen, not the old paper-tip count and not the raw NMNS row count. Exact phylogenetic insertion is **not** assigned from morphology, geography or NMNS taxonomy. New samples are target-captured and their placement is inferred from the nuclear data.

## Current NMNS / Moreyra reconciliation snapshot

The current NMNS data source contains:

- 161 authority records;
- 154 unique species-binomial strings;
- 128 species in the `PUBLISHED_OR_NO_NEW_LABEL` core screen;
- 25 `新称` species in a separate unpublished-new-name extension lane;
- 1 `仮称` species in a separate provisional lane.

At the current species-level crosswalk, the Moreyra Japanese material covers 33 of the 128 core species. The previous arithmetic `120 - 38 = 82 missing species` is therefore superseded.

## Exact published-core hole-fill counts

| Priority | Hole type | Species | New tree samples |
| --- | --- | ---: | ---: |
| P1 | taxonomic block has zero Moreyra species | 33 | 66 |
| P2 | current species has a broken/weak public tip | 6 | 12 |
| P3 | block coverage <=25% | 44 | 88 |
| P4 | block coverage <=50% | 16 | 32 |
| P5 | local gap in an otherwise better-covered block | 2 | 4 |
| P6 | adequate public singleton; add one own trait-linked representative | 26 | 26 |
| P7 | public tree already has >=2 representatives | 1 | 0 for tree |
| **Total** |  | **128** | **228** |

For P7 (`Cirsium incomptum`), one additional own wild individual is still recommended for the same-individual phenotype scaffold, but it is not required to fill the nuclear tree slot.

## Sampling order

The biological priority is:

> **P1 empty blocks -> P2 broken tips -> P3 sparse blocks -> P4 partial blocks -> P5 local gaps -> P6 own trait-linked complements -> P7 optional trait link.**

Field season and route efficiency can reorder taxa **within** one priority class, but they do not promote a lower-priority tree slot above an unfilled higher-priority hole unless conservation/permission constraints make the higher slot unavailable.

## What “two samples” means

For P1-P5 and P2 repair taxa, the default is **two independent wild Japanese biological individuals**.

- If the species spans multiple range sectors, Sample A and Sample B should come from different verified populations/range sectors and maximize geographic separation where feasible.
- If the species is narrow endemic or only one authorized population is available, use two spatially separated independent plants and flag the tree slot as `SINGLE_POPULATION_LIMITED` rather than inventing a second locality.
- Exact sites are not frozen from NMNS database prose; current occurrence, conservation status and permissions must be checked before collection.

For P6, retain the usable Moreyra public accession as representative A and add one own wild, voucher-linked, phenotype-linked representative B.

## Priority 1: block-level holes with zero Moreyra coverage

These 33 species occupy NMNS taxonomic blocks for which the current Moreyra species-level skeleton contributes no species. They are the first hole-fill targets because entire classification regions are missing from the skeleton.

- `C. hanamakiense`, `C. tamastoloniferum`
- `C. ashinokurense`
- `C. apoense`, `C. charkeviczii`, `C. pectinellum`, `C. yezoalpinum`
- `C. borealinipponense`, `C. chokaiense`, `C. hachimantaiense`, `C. maruyamanum`, `C. occidentalinipponense`, `C. okamotoi`, `C. shimae`
- `C. hyugamontanum`, `C. kirishimense`, `C. kujuense`, `C. tenuisquamatum`, `C. unzenense`
- `C. homolepis`, `C. inundatum`, `C. ugoense`
- `C. babanum`
- `C. ganjuense`
- `C. boreale`, `C. hidakamontanum`
- `C. chikushiense`
- `C. nambuense`
- `C. hachijoense`, `C. toyoshimae`
- `C. shinanense`
- `C. tashiroi`
- `C. purpratum`

The source spellings above are preserved from the current NMNS snapshot; nomenclatural corrections are not silently introduced in this ledger.

## Priority 2: six broken/weak Moreyra tips

These species already have a Moreyra species-level anchor, but the current public material is inadequate as the primary Japanese wild representative and is replaced by two new wild Japanese samples:

- `C. pendulum` — foreign reference for a Japanese-distributed taxon;
- `C. sieboldii` — cultivated reference; draw both replacement representatives from the P01 population bank;
- `C. nipponicum` — cultivated reference;
- `C. kamtschaticum` — foreign reference for a Japanese-distributed taxon;
- `C. buergeri` — cultivated reference;
- `C. microspicatum` — cultivated reference.

The old accessions remain provenance/geographic sensitivity samples when interpretable; they do not count as primary Japanese wild representatives.

## Priority 3: sparse blocks

Five blocks contain Moreyra anchors but <=25% of their current core species are represented. The exact missing species are stored in block-specific ledgers:

1. Kaganoazami block: 17 missing species around JPN08/JPN11/JPN13/JPN14/JPN28.
2. Sawaazami block: 8 missing species around JPN24/JPN30.
3. Norikura-series block: 7 missing species around JPN18.
4. Hamaazami block: 4 missing species around JPN17; this contains `C. boninense`, `C. brevicaule`, `C. irimtiense`, and `C. spinosum`.
5. Yamaazami block: 8 missing species around JPN34/JPN26.

These anchor labels identify the **coverage hole** only. They do not assert that the new species will attach as sisters to those named tips.

## Priority 4 and Priority 5

Priority 4 contains 16 species in partially covered blocks (<=50% coverage). Priority 5 contains two local gaps, `C. katoanum` and `C. takahashii`, in the currently better-covered Himeazami block.

## Priority 6 and Priority 7

Priority 6 contains 26 species with one adequate public Moreyra representative. Add one own wild representative linked to voucher, orientation, phyllary, stickiness, colour and cytotype/genome-size status.

Priority 7 contains only `C. incomptum`, which already has >=2 public biological representatives at the current species-level crosswalk. No new tree sample is required; one own phenotype-linked plant remains desirable for Chapter 3 trait linkage.

## Taxonomy extension lane

The 25 `新称` species and 1 `仮称` species are **not silently merged into the 128-species core**. If all 26 are later admitted to the operational species census, they require 51 additional new samples: two each for the 25 unpublished-new-name species, and one additional representative for provisional `C. ishidatense` because JPN10 already supplies one public biological sample.

## Taxonomy gates, not collection slots

Three Moreyra paper concepts are not currently matched to a species in the NMNS species screen:

- JPN29 `C. verutum`;
- JPN31 `C. yuki-uenoanum`;
- JPN33 `C. effusum`.

Do not assign two-sample collection quotas to these names until their taxonomic concepts are reconciled. A convenience substitute is not permitted.

## Focal population studies are nested into the hole fill

- `C. sieboldii`: P2 replacement representatives come from the P01 RAD population bank.
- `C. dipsacolepis` and `C. lineare`: their P6 own representatives come from their RAD banks.
- `C. brevicaule`: the two P3 tree samples are one Okinawa-Honto and one Amami-Oshima discovery individual.
- `C. irimtiense`: the two P3 tree samples are one Miyako and one Ishigaki discovery individual.

No separate tree-only plants are collected for these focal systems when the population bank can supply the required representative.

## Authoritative machine-readable files

- `data/planning/chapter3_exact_hole_fill_summary_v7.json`
- `data/planning/chapter3_holefill_priority1_empty_blocks_v7.csv`
- `data/planning/chapter3_holefill_priority2_broken_tips_v7.csv`
- `data/planning/chapter3_holefill_priority3_kaganoazami_block_v7.csv`
- `data/planning/chapter3_holefill_priority3_sawaazami_block_v7.csv`
- `data/planning/chapter3_holefill_priority3_norikura_series_v7.csv`
- `data/planning/chapter3_holefill_priority3_hamaazami_block_v7.csv`
- `data/planning/chapter3_holefill_priority3_yamaazami_block_v7.csv`
- `data/planning/chapter3_holefill_priority4_partial_blocks_v7.csv`
- `data/planning/chapter3_holefill_priority5_local_gaps_v7.csv`
- `data/planning/chapter3_holefill_priority6_add_own_trait_link_v7.csv`
- `data/planning/chapter3_holefill_priority7_tree_filled_v7.csv`
- `data/planning/chapter3_holefill_taxonomy_extension_v7.csv`
- `data/planning/chapter3_holefill_taxonomy_conflicts_v7.csv`

Collection and sequencing remain fail-closed until the operational species census, current occurrence and permission gates are frozen.
