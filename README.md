# aza3 — Japan-wide Cirsium tree to population and mechanism

This repository starts Chapter 3 from the uncertainty left by the completed public-data Chapter 2. Future target-capture, RAD-seq, transcriptomic or field data are not treated as retroactive confirmation.

## Locked starting point

EAzami Chapter 2 remains frozen at merge `62fa8c5c913c2b236e710f6bad366e80676aa78f`.

Within its admitted public topology ensemble:

- orientation requires **at least four** state changes;
- phyllary posture requires at least three;
- stickiness requires at least five;
- minimum counts are better resolved than individual event placements;
- **species-tip compression** hides within-species state multiplicity.

## one question, breadth and depth

> What is the own-data evolutionary history of the Japanese Cirsium radiation, where did capitulum traits change on that history, and how far can selected transitions be decomposed from population history to molecular mechanism?

### P01-P05 breadth

P01-P05 ask where apparent transitions in orientation, phyllary posture, stickiness and related morphology survive when the analysis is rebuilt on own individual-linked phenotype, cytotype and nuclear ancestry.

### M01 depth — embedded worked case, not P06

The white versus bluish-purple floral-pigmentation contrast in `Cirsium brevicaule` and `C. irumtiense` remains an **embedded worked case**, **not P06**. It follows one retained contrast through history, pigment chemistry, corolla expression and selection-consistent genomic evidence without presupposing loss, regain or a selective agent.

## Primary genomic product

The first Chapter 3 genomic product is:

> **a Japan-wide nuclear Cirsium species-tree/network built with Comp1061-compatible target capture.**

The existing Moreyra Japanese target-capture material is reused as a public skeleton. New acquisition is therefore a **hole-fill problem**, not a blanket `all species x2 new plants` recollection.

## Exact hole-fill v7 — authoritative sampling ledger

`docs/CHAPTER3_EXACT_BACKBONE_HOLE_FILL_LEDGER_V7.md` and `data/planning/chapter3_exact_hole_fill_summary_v7.json` are authoritative for Level-1 acquisition priority.

The current NMNS snapshot contains:

- **161 authority records**;
- **154 unique species-binomial strings**;
- **128 species** in the current `PUBLISHED_OR_NO_NEW_LABEL` core screen;
- **25 新称** species in a separate unpublished-new-name extension lane;
- **1 仮称** species in a separate provisional lane.

At the current species-level crosswalk, the Moreyra Japanese material covers **33 of the 128 core species**. Thus the older `120 - 38 = 82` and `211 new individuals` arithmetic is superseded.

### Exact published-core acquisition priorities

| Priority | Tree hole | Species | New tree representation slots |
| --- | --- | ---: | ---: |
| **P1** | taxonomic block has zero Moreyra species | 33 | **66** |
| **P2** | broken/weak public species tip | 6 | **12** |
| **P3** | block coverage <=25% | 44 | **88** |
| **P4** | block coverage <=50% | 16 | **32** |
| **P5** | local gap in a better-covered block | 2 | **4** |
| **P6** | adequate public singleton; add one own trait-linked representative | 26 | **26** |
| **P7** | public tree already has >=2 representatives | 1 | **0** |
| **Total** |  | **128** | **228 required representation slots** |

These 228 slots are **not automatically 228 wild collections**. `data/planning/chapter3_conservation_override_v9.csv` can replace wild collection with authenticated historical/herbarium material or minimal authorized tissue sampling while preserving the tree-representation goal.

For P7 `C. incomptum`, one additional own phenotype-linked individual is recommended for Chapter 3 trait linkage but is not required to fill the nuclear tree slot.

The biological priority is:

> **P1 empty blocks -> P2 broken tips -> P3 sparse blocks -> P4 partial blocks -> P5 local gaps -> P6 own trait-linked complements -> P7 optional trait link.**

Field season and route efficiency may reorder taxa within one priority class, but they do not redefine the tree holes. Conservation status overrides travel convenience and default wild-collection logic.

## Priority 1 — entire classification blocks missing

P1 contains 33 core species in NMNS taxonomic blocks with zero current Moreyra species-level representation. The exact list is frozen in `data/planning/chapter3_holefill_priority1_empty_blocks_v7.csv`.

Examples include the entire current blocks containing `C. apoense / C. pectinellum / C. yezoalpinum`, the seven-species Oni-azami block, the five-species Kirishima-azami block, `C. hachijoense / C. toyoshimae`, and the Fuji-azami section entry in the current source snapshot.

These block names define coverage gaps only. They do **not** predetermine the nuclear placement of the new samples.

Conservation overrides are mandatory where needed. For example, `C. toyoshimae` is treated as a historical-material-only tree slot rather than a wild collection target, and highly threatened taxa such as `C. apoense` enter a conservation gate before any wild tissue acquisition is considered.

## Priority 2 — six broken public tips

Two new Japanese representatives are required for:

- `C. pendulum`;
- `C. sieboldii`;
- `C. nipponicum`;
- `C. kamtschaticum`;
- `C. buergeri`;
- `C. microspicatum`.

Wild material is preferred where safe, legal and biologically justified. Old cultivated/foreign accessions remain provenance or geographic sensitivity samples when interpretable; they do not automatically occupy the primary Japanese slots.

## Priority 3 — five sparse blocks

P3 contains 44 missing species in five blocks with <=25% current Moreyra coverage:

- Kaganoazami block: **17** missing species;
- Sawaazami block: **8**;
- Norikura-series block: **7**;
- Hamaazami block: **4** (`C. boninense`, `C. brevicaule`, `C. irimtiense`, `C. spinosum` in the current source snapshot);
- Yamaazami block: **8**.

Each block has its own exact CSV in `data/planning/`.

## What two representation slots mean

For P1-P5 and P2 repair taxa, the default is **two independent Japanese biological representatives**, but acquisition mode is conservation-aware.

- For ordinary broad-range taxa, use different verified wild populations/range sectors and maximize separation where feasible.
- For narrow endemic or conservation-sensitive taxa, prefer authenticated herbarium/ex-situ material where it can answer the tree question; minimal wild tissue is considered only after explicit conservation review and permission.
- For extinct taxa, wild collection is impossible and the slot must be filled from authenticated historical material or remain unresolved.
- Exact sites are not inferred from NMNS prose; current occurrence, conservation and permission checks must be frozen before field collection.

For P6, retain the usable Moreyra public accession as representative A and add one own phenotype-linked representative B only after auditing the public accession locality. If that locality remains too coarse, do not guess the complementary sector.

## Taxonomy is a separate gate

The 25 `新称` species and one `仮称` species are not silently promoted into the 128-species core. If all 26 extension units are later admitted, they require **51 additional representation slots** under the current arithmetic.

Three Moreyra paper concepts are not assigned new collection slots until concept reconciliation:

- JPN29 `C. verutum`;
- JPN31 `C. yuki-uenoanum`;
- JPN33 `C. effusum`.

No convenience replacement is permitted.

## Every own fresh tree plant is also a trait plant

Each own fresh target-capture individual links one immutable ID to taxon/population, voucher/diagnostic images, orientation, phyllary posture/calibrated image, stickiness/gland state, flower colour, developmental stage, DNA source and cytotype/genome-size evidence status.

Historical/herbarium acquisitions cannot provide the same live phenotype bundle and are marked accordingly rather than treated as equivalent trait plants.

Thus the expanded tree does not recreate species-tip compression silently.

## Focal RAD is nested into the hole fill

Do not collect separate tree-only plants when the focal population bank already supplies the needed representative.

- `C. sieboldii`: its two P2 replacement representatives come from the P01 RAD bank if conservation/permission gates permit;
- `C. dipsacolepis`: its P6 own representative comes from its RAD bank;
- `C. lineare`: same;
- `C. brevicaule`: its two P3 tree samples are one Okinawa-Honto and one Amami-Oshima discovery individual;
- `C. irimtiense`: its two P3 tree samples are one Miyako and one Ishigaki discovery individual.

RAD remains a population-ancestry/reticulation instrument below the nationwide target-capture scaffold. It is not the method used to define the nationwide species tree.

## Analysis hierarchy

1. **Level 0 — taxonomic census:** reconcile the NMNS authority universe into operational species concepts.
2. **Level 1 — exact backbone hole fill + target capture:** combine compatible Moreyra public anchors with P1-P7 own/historical acquisitions to build the Japanese nuclear tree/network.
3. **Level 2 — trait history:** remap orientation, phyllary, stickiness and cross-module histories on the own-data nationwide ensemble.
4. **Level 3 — population RAD:** deepen only focal lineages/transition neighbourhoods where shallow ancestry, introgression or morph history discriminates competing histories.
5. **Level 4 — M01:** add pigment/RNA and later candidate-region confirmation to the embedded floral-pigmentation case.

## Start here

1. `docs/CHAPTER3_EXACT_BACKBONE_HOLE_FILL_LEDGER_V7.md` — authoritative exact species-level hole-fill plan.
2. `data/planning/chapter3_exact_hole_fill_summary_v7.json` — exact 128-species / 228-slot contract.
3. `data/planning/chapter3_conservation_override_v9.csv` — conservation-aware acquisition overrides.
4. `docs/CHAPTER3_EXACT_SITE_FREEZE_V8.md` — field/historical acquisition freeze rules.
5. `data/planning/chapter3_holefill_priority1_empty_blocks_v7.csv` — 33 species in fully empty blocks.
6. `data/planning/chapter3_holefill_priority2_broken_tips_v7.csv` — six broken public tips.
7. `data/planning/chapter3_holefill_priority3_*_v7.csv` — five sparse-block ledgers, 44 species total.
8. `data/planning/chapter3_holefill_priority4_partial_blocks_v7.csv` — 16 partial-block holes.
9. `data/planning/chapter3_holefill_priority5_local_gaps_v7.csv` — two local gaps.
10. `data/planning/chapter3_holefill_priority6_add_own_trait_link_v7.csv` — 26 own complements.
11. `data/planning/chapter3_holefill_priority7_tree_filled_v7.csv` — current tree-filled slot.
12. `data/planning/chapter3_holefill_taxonomy_extension_v7.csv` — 25 new-name + one provisional extension units.
13. `data/planning/chapter3_holefill_taxonomy_conflicts_v7.csv` — JPN29/JPN31/JPN33 concept gates.
14. `data/planning/chapter3_sampling_priorities_v1.csv` — frozen P01-P05 inherited from Chapter 2.
15. `data/planning/chapter3_mechanistic_demonstration_v1.json` — M01 evidence ladder.
16. `docs/CHAPTER3_RADSEQ_END_TO_END_V1.md` — population-RAD safeguards.

v1-v6 inverse-sampling, 38x3, zero-baseline, transition-first, calendar-first and nominal `120/82/211` files remain design history only and do not govern new collection.

## Current state

- own biological data admitted: **0**;
- physical samples: **0**;
- operational nationwide species census frozen: **false**;
- nationwide collection authorized: **false**;
- target-capture sequencing authorized: **false**;
- population RAD authorized: **false**;
- M01 expansion authorized: **false**;
- regain claim authorized: **false**;
- selection claim authorized: **false**;
- pollinator-agent claim authorized: **false**.

Run the current Chapter 3 validators before accepting changes to these boundaries.
