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

## Authoritative sampling v4: nationwide species tree first

New sampling is governed by `docs/CHAPTER3_NATIONWIDE_SPECIES_TREE_SAMPLING_V4.md` and `data/planning/chapter3_nationwide_species_tree_sampling_v4.json`.

The primary genomic product is now:

> **a Japan-wide nuclear Cirsium species-tree/network built with target capture.**

RAD is nested underneath that tree for population ancestry and reticulation; it is not the method used to define the nationwide cross-species scaffold.

## How many species and individuals?

Current NMNS public resources use different scopes: a 2025 article describes approximately **120 species**, while the thistle database currently returns **161 authority records**, including infraspecific and unpublished/new-name entries. Therefore 161 database rows are not treated as 161 independent species.

Planning baseline:

- **120 primary species concepts** as the working nationwide count;
- final operational count must be frozen by taxonomic reconciliation before destructive collection;
- if the final primary-species count exceeds 125, reopen the individual budget rather than reducing replication.

### Nationwide target-capture bank

Every admitted species receives **two independent biological individuals** where possible.

At 120 species:

- base: 120 ×2 = **240 individuals**;
- up to 30 widespread/variable taxa receive +1 = **+30**;
- up to 10 complex/polyploid/hybrid-suspected taxa receive +2 = **+20**;
- nominal full panel = **290 target-capture candidates**.

At 125 species with the same enrichment envelope, the planning ceiling is **300 individuals**.

The two base individuals should come from different populations or clearly separated geographic occurrences whenever feasible. A conservation-limited species may remain one-sample with an explicit `SINGLE_SAMPLE_LIMITED` flag; it is not replaced by a convenience relative.

### Target-capture sequencing waves

**Wave 1:** sequence the two base individuals/species first — nominally **240 individuals** at 120 species.

**Wave 2:** sequence up to **50 preregistered enrichment individuals** for widespread/variable or complex taxa.

The intended full nationwide nuclear panel is therefore **290–300 individuals**, not 38×3 RAD samples.

Use a **Comp1061-compatible nuclear target-capture assay**, or an explicitly crosswalked equivalent, so the own-data panel remains compatible with the existing Moreyra/EAzami nuclear scaffold.

## Every tree sample is also a trait sample

The nationwide tree must not recreate species-tip compression. Every newly collected sequenced individual links one immutable ID to:

- taxon and population;
- voucher/diagnostic images;
- capitulum orientation;
- phyllary posture plus calibrated image;
- stickiness/gland state;
- flower colour;
- developmental stage;
- silica/fresh leaf DNA;
- cytotype or genome-size evidence status;
- deidentified authorization/conservation IDs.

Thus Level 1 provides both the own nuclear tree and the individual-linked phenotype scaffold for P01-P05.

## Analysis hierarchy

1. **Level 0 — taxonomic census:** reconcile the NMNS authority universe into operational species concepts.
2. **Level 1 — nationwide target capture:** build the Japanese nuclear species-tree/network ensemble.
3. **Level 2 — trait history:** remap orientation, phyllary, stickiness and cross-module histories on the own-data nationwide ensemble.
4. **Level 3 — population RAD:** expand only focal lineages/transition neighbourhoods where shallow ancestry, introgression or morph history can discriminate competing histories.
5. **Level 4 — M01:** add pigment/RNA and later candidate-region confirmation to the embedded floral-pigmentation case.

## Nested focal RAD programme

Current target population-RAD design:

- `C. sieboldii`: 4 populations ×12 = **48**;
- `C. dipsacolepis`: 3 ×12 = **36**;
- `C. lineare`: 3 ×12 = **36**;
- `C. brevicaule`: 2 ×12 = **24**;
- `C. irumtiense`: 2 ×12 = **24**.

Total initial focal RAD = **168 primary individuals**.

These are not added as a separate tree-only collection. Up to four RAD individuals per focal species also serve as nationwide target-capture representatives. Under the nominal full design:

- nationwide bank = 290;
- focal RAD = 168;
- same-individual overlap = 20;
- additional focal plants beyond the nationwide bank ≈148;
- full national-tree + initial focal-RAD programme ≈ **438 unique physical plants**.

The species tree can be completed before the full 438-plant programme: the first publishable nationwide genomic wave is the 240-individual two-per-species target-capture panel.

## RAD assay boundary

The focal RAD enzyme/complexity pilot remains nested within the five mandatory systems:

- Stage A: 5 systems ×1 DNA ×3 candidates = **15 shallow libraries**;
- Stage B: 5 systems ×2 independent individuals ×top-2 protocols =20 primary +5 cross-batch repeats = **25 libraries**;
- genotype concordance gate ≥0.95;
- core-locus recovery gate ≥0.90.

If later Level 2 selects high-ploidy taxa for population RAD, run a separate ploidy/complexity-stratum assay check. Mixed-ploidy RAD is never allowed to define the nationwide species tree.

## M01 remains nested

`C. brevicaule` and `C. irumtiense` are already included in the nationwide target-capture panel and initial focal RAD programme. M01 adds floral material in Okinawa Honto, Amami Oshima, Miyako and Ishigaki: six developmentally matched floral collections/population, five primary RNA replicates and one reserve.

Further islands remain conditional on E1/E2 promotion. RAD-only outliers are not final selection evidence, and genomics alone cannot identify the selective agent.

## Start here

1. `docs/CHAPTER3_NATIONWIDE_SPECIES_TREE_SAMPLING_V4.md` — authoritative species/individual sampling architecture.
2. `data/planning/chapter3_nationwide_species_tree_sampling_v4.json` — machine-readable v4 contract.
3. `data/planning/chapter3_nationwide_species_tree_budget_v4.csv` — 120-species / 290-individual planning budget and nested RAD totals.
4. `docs/CHAPTER3_SCOPE_AND_HANDOFF_V1.md` — frozen Chapter 3 claim boundaries.
5. `data/planning/chapter3_sampling_priorities_v1.csv` — frozen P01-P05 inherited from Chapter 2.
6. `data/planning/chapter3_mechanistic_demonstration_v1.json` — M01 evidence ladder.
7. `docs/CHAPTER3_RADSEQ_END_TO_END_V1.md` — RAD laboratory/analysis safeguards.
8. `data/intake/chapter3_individual_intake_v1.csv` — empty biological intake ledger.
9. `data/intake/chapter3_radseq_library_intake_v1.csv` — empty RAD library/QC ledger.

v1-v3 inverse-sampling, 38×3, zero-baseline, transition-first, eight-anchor and inventory-first files remain design history only. Transition-neighbourhood logic from v3 remains useful at Level 2/3 after the nationwide tree is built.

## Current state

- own biological data admitted: **0**;
- physical samples: **0**;
- operational nationwide species census frozen: **false**;
- nationwide collection authorized: **false**;
- target-capture Wave 1 authorized: **false**;
- target-capture Wave 2 authorized: **false**;
- population RAD authorized: **false**;
- M01 expansion authorized: **false**;
- regain claim authorized: **false**;
- selection claim authorized: **false**;
- pollinator-agent claim authorized: **false**.

Run the Chapter 3 validators before accepting changes to these boundaries.
