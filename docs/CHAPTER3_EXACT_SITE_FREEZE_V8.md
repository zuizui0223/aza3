# Chapter 3 exact-site freeze workflow v8

The v7 hole-fill plan answers **which species and which range sector**. v8 answers **when a specific acquisition locality or non-field material source is sufficiently verified to represent that tree slot**.

## Queue size

`analysis/build_exact_site_freeze_queue_v8.py` expands the v7 species/sector ledgers into one execution row per planned sample:

- P1-P5: 101 two-representative species = **202 required slots**;
- P6: 26 public-singleton complements = **26 required slots**;
- P7: `Cirsium incomptum` = **1 optional own trait-linked slot**;
- total = **229 execution slots**, of which **228 are required for the published-core tree**.

The 228 slots are **tree-representation slots, not automatically 228 wild collections**.

## Conservation override

`data/planning/chapter3_conservation_override_v9.csv` can replace wild collection with safer acquisition modes without deleting the corresponding tree slot.

Examples:

- extinct taxa: authenticated historical/herbarium DNA only; wild collection is blocked;
- CR/EN or otherwise highly conservation-sensitive taxa: prefer authenticated herbarium/ex-situ material and allow minimal wild tissue only after explicit conservation review and permission;
- ordinary taxa: wild tissue may proceed only after the normal occurrence and permission gates.

The queue therefore distinguishes `NOT_FROZEN`, `WILD_COLLECTION_BLOCKED`, and `CONSERVATION_GATE_REQUIRED`. None of these means that collection has been authorized.

## Range sector is not an exact site

A v7 statement such as `Shiretoko + Daisetsu`, `Aomori + Iwate`, or `Amami Oshima + Okinawa Honto` is an acquisition constraint, not proof that a currently collectable population exists at a particular point.

A wild-field slot can become `FROZEN_FOR_FIELD_COLLECTION` only after all of the following are resolved:

1. current occurrence is supported by dated evidence;
2. conservation restrictions are reviewed;
3. the relevant land manager is identified or shown not to apply;
4. collection permission is approved or documented as not required;
5. tissue collection permission is approved or documented as not required;
6. the flowering/field window is checked;
7. a deidentified `target_locality_id` is assigned;
8. the exact site is stored in a private field record outside the public repository.

Historical/herbarium slots instead require authenticated specimen identity, institution/material availability, destructive-sampling permission, and sufficient DNA feasibility for the chosen target-capture workflow.

## No sensitive coordinates in the public repository

The public repository must never contain precise coordinates or rare-species microhabitat directions. The public queue records only range-sector evidence, permission/review status, acquisition mode, and `target_locality_id`. Coordinates remain in a private field ledger linked by that ID.

This is especially important for narrow endemic, island, alpine and conservation-sensitive taxa.

## P6 is deliberately delayed

For P6 the Moreyra tree already contains one usable public individual. The purpose of the new sample is geographic and phenotype-linked complementarity. Therefore the public accession locality must be audited first. Only then is the new own range sector chosen to maximize nonredundant coverage within the NMNS-documented distribution.

If the existing public locality is only `Japan` or otherwise too coarse, an opposite sector must not be guessed. The taxon remains unresolved until voucher provenance is recovered or a two-sector own sensitivity design is justified.

`C. dipsacolepis` and `C. lineare` remain operational exceptions in that their own representative is drawn from the focal RAD population banks; the selected RAD population must still pass current-occurrence and permission gates.

## Freeze order

Freeze acquisition slots in biological hole-fill priority, not travel convenience:

1. P1 whole-block gaps;
2. P2 broken public tips;
3. P3 sparse blocks;
4. P4 partial blocks;
5. P5 local gaps;
6. P6 singleton complements after public-locality audit;
7. P7 optional trait-linked complement.

Within every priority, conservation overrides supersede convenience. An extinct or critically endangered taxon is not promoted to wild collection merely because its range sector is logistically convenient.

## Generated queue fields

The queue records:

- immutable slot ID;
- priority/species/sample slot;
- required versus optional status;
- acquisition mode;
- v7 range sector and source distribution;
- whether public-locality audit is required;
- current-occurrence evidence status and dated evidence ID;
- conservation and land-manager review;
- collection/tissue permission status;
- deidentified target-locality ID;
- status of the private exact-site record;
- field-window status;
- collector/collaborator;
- final freeze/block status.

The exact queue is reproducibly generated in CI and audited to contain 229 rows, 228 required slots, no public coordinate fields, and explicit conservation blocks/gates where required.
