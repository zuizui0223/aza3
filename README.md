# aza3 — Japan-wide Cirsium tree to population and mechanism

This repository starts Chapter 3 from the uncertainty left by the completed public-data Chapter 2. It does not treat future target-capture, RAD-seq, transcriptomic or field data as retroactive confirmation.

## Locked starting point

EAzami Chapter 2 remains frozen at merge `62fa8c5c913c2b236e710f6bad366e80676aa78f`.

Within its admitted public topology ensemble:

- orientation requires **at least four** state changes;
- phyllary posture requires at least three;
- stickiness requires at least five;
- minimum counts are better resolved than individual event placements;
- **species-tip compression** hides within-species state multiplicity.

## One question, breadth and depth

> What is the own-data evolutionary history of the Japanese Cirsium radiation, where did capitulum traits change on that history, and how far can selected transitions be decomposed from population history to molecular mechanism?

### P01-P05 breadth

P01-P05 ask where apparent transitions in orientation, phyllary posture, stickiness and related morphology survive when the analysis is rebuilt on own individual-linked phenotype, cytotype and nuclear ancestry.

### M01 depth — embedded worked case, not P06

The white versus bluish-purple floral-pigmentation contrast in `Cirsium brevicaule` and `C. irumtiense` remains an **embedded worked case**, **not P06**. It follows one retained contrast through history, pigment chemistry, corolla expression and selection-consistent genomic evidence without presupposing loss, regain or a selective agent.

## Nationwide species tree first

The primary genomic product is:

> **a Japan-wide nuclear Cirsium species-tree/network built with Comp1061-compatible target capture.**

The current planning baseline is approximately **120 primary species concepts**, but the final operational list must be frozen from the NMNS authority universe before destructive collection.

The nationwide analysis still targets two primary representatives/species, but **new collection is now hole-filled against the existing Moreyra 38-species public skeleton rather than collecting 240 new plants blindly**.

## Backbone hole-fill v6 — authoritative Level-1 acquisition priority

`docs/CHAPTER3_BACKBONE_HOLE_FILL_V6.md` defines which new samples are actually needed.

At the nominal 120-species baseline:

### H1 — species absent from the Moreyra public skeleton

- approximately 82 species are absent from the 38-tip Japanese public scaffold;
- collect **2 new wild Japanese individuals/species**;
- nominal new material = **164 individuals**.

These are true empty species slots. Their exact branch position is inferred from target capture; it is never assigned from geography or morphology in advance.

### H2 — weak or broken Moreyra tips

Nine current Moreyra concepts need high-priority repair:

- identity blocked: JPN29 `C. verutum`, JPN31 `C. yuki-uenoanum`, JPN33 `C. effusum`;
- cultivated/foreign provenance requiring wild-Japan replacement: JPN32 `C. buergeri`, JPN34 `C. microspicatum`, JPN35 `C. nipponicum`, JPN36 `C. sieboldii`, JPN37 `C. kamtschaticum`, JPN38 `C. pendulum`.

Collect **2 new wild Japanese individuals/concept** = **18 individuals**. Old accessions remain labelled sensitivity samples where appropriate; they do not occupy a primary wild-Japan slot.

### H3 — adequate public tip but no own trait-linked individual

The remaining nominal 29 adequate public tips keep one Moreyra representative and receive **one new voucher-linked wild Japanese individual** each.

Nominal H3 = **29 individuals**.

### Level-1 arithmetic

- Phase A: H1 + H2 = **182 new individuals**;
- Phase B: H3 = **+29**;
- total new material = **211 individuals**;
- retain 29 adequate Moreyra public representatives;
- primary nationwide matrix = **240 representatives = two/species at 120 species**.

General formula:

`new = 2*(N-38) + (38-R) + 2R = 2N - 38 + R`,

where `N` is the frozen Japanese species count and `R` is the number of public concepts requiring two-new-sample repair. At `N=120, R=9`, new = **211**.

This is the governing acquisition rule: **fill an empty or broken species slot before adding redundant depth to a tip that is already represented.**

## Geographic hole order

Within H1, prioritize lineage-space that is poorly represented in the 38-tip scaffold:

1. **isolated island lineages** — e.g. `C. boninense`, `C. hachijoense`, `C. sadoense`, `C. brevicaule`, `C. irumtiense`, `C. umezawanum`;
2. **Hokkaido/northern radiation gaps** — e.g. `C. apoense`, `C. yezoalpinum`, `C. pectinellum`, `C. albrechtii`, `C. boreale`, `C. iito-kojianum`;
3. **alpine/narrow mainland gaps** — e.g. `C. zawoense`, `C. ugoense`, `C. fauriei`, `C. babanum`, `C. furusei`;
4. **western regional radiations** — e.g. `C. calcicola`, `C. ashiuense`, `C. taishakuense`, `C. kirishimense`, `C. nishimeraense`, `C. austrokiusianum`.

These are verified examples, not the final operational census.

## Every new tree sample is also a trait sample

Every new nationwide-tree individual links one immutable ID to taxon/population, voucher/images, orientation, phyllary posture/calibrated image, stickiness/gland state, flower colour, developmental stage, DNA source and cytotype/genome-size evidence status.

Thus the expanded tree does not recreate species-tip compression.

## Focal RAD is nested inside the hole fill

Do not collect separate tree-only plants when the population bank already supplies the needed Level-1 representative.

- `C. sieboldii`: its two H2 repair individuals come from the P01 RAD bank;
- `C. dipsacolepis`: its one H3 own tree representative comes from its RAD bank;
- `C. lineare`: same;
- `C. brevicaule`: its two H1 tree representatives are one from Okinawa Honto and one from Amami Oshima;
- `C. irumtiense`: its two H1 tree representatives are one from Miyako and one from Ishigaki.

Full focal RAD target remains 168 primary individuals, but it is downstream of Level-1 coverage and is not a substitute for species-tree hole filling.

## Target-capture staging

The first technical qualification batch is drawn from actual H1/H2 samples and remains in the nationwide matrix. No separate pilot taxa are collected.

After H1/H2 coverage is established, add H3 own trait-linked representatives. Third/fourth target-capture individuals are H4 enrichment and open only after the first nationwide tree/network identifies non-monophyly, strong discordance, ploidy/hybrid problems or major geographic structure.

## Analysis hierarchy

1. **Level 0 — taxonomic census:** reconcile the NMNS authority universe into operational species concepts.
2. **Level 1 — backbone hole fill + nationwide target capture:** combine compatible Moreyra public anchors with H1/H2/H3 own samples to build the Japanese nuclear tree/network.
3. **Level 2 — trait history:** remap orientation, phyllary, stickiness and cross-module histories on the expanded nationwide ensemble.
4. **Level 3 — population RAD:** deepen only focal lineages/transition neighbourhoods where shallow ancestry, introgression or morph history can discriminate competing histories.
5. **Level 4 — M01:** add pigment/RNA and later candidate-region confirmation to the embedded floral-pigmentation case.

## Start here

1. `docs/CHAPTER3_BACKBONE_HOLE_FILL_V6.md` — authoritative new-sample hole-fill priority.
2. `data/planning/chapter3_backbone_hole_fill_v6.json` — machine-readable hole-fill contract.
3. `data/planning/chapter3_moreyra_anchor_repair_v6.csv` — nine broken/weak public tips requiring repair.
4. `docs/CHAPTER3_NATIONWIDE_SPECIES_TREE_SAMPLING_V4.md` — nationwide species-tree architecture.
5. `docs/CHAPTER3_FIELD_PRIORITY_CALENDAR_V5.md` — seasonal logistics only; subordinate to the hole-fill priority.
6. `data/planning/chapter3_sampling_priorities_v1.csv` — frozen P01-P05 inherited from Chapter 2.
7. `data/planning/chapter3_mechanistic_demonstration_v1.json` — M01 evidence ladder.
8. `docs/CHAPTER3_RADSEQ_END_TO_END_V1.md` — RAD safeguards.

Older v1-v3 inverse-sampling, 38×3, zero-baseline and transition-first plans remain design history only.

## Current state

- own biological data admitted: **0**;
- physical samples: **0**;
- operational nationwide species census frozen: **false**;
- nationwide collection authorized: **false**;
- H1/H2/H3/H4 acquisition open: **false**;
- target-capture sequencing authorized: **false**;
- population RAD authorized: **false**;
- M01 expansion authorized: **false**;
- regain claim authorized: **false**;
- selection claim authorized: **false**;
- pollinator-agent claim authorized: **false**.

Run the Chapter 3 validators before accepting changes to these boundaries.
