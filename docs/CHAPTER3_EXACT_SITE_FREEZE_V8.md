# Chapter 3 exact-site freeze workflow v8

The v7 hole-fill plan answers **which species and which range sector**. v8 answers **when a specific field locality is sufficiently verified to be used for collection**.

## Queue size

`analysis/build_exact_site_freeze_queue_v8.py` expands the v7 species/sector ledgers into one execution row per planned sample:

- P1-P5: 101 two-sample species = **202 required slots**;
- P6: 26 public-singleton complements = **26 required slots**;
- P7: `Cirsium incomptum` = **1 optional own trait-linked slot**;
- total = **229 execution slots**, of which **228 are required for the published-core tree**.

Every slot begins `NOT_FROZEN`.

## Range sector is not an exact site

A v7 statement such as `Shiretoko + Daisetsu`, `Aomori + Iwate`, or `Amami Oshima + Okinawa Honto` is an acquisition constraint, not proof that a currently collectable population exists at a particular point.

A slot can become `FROZEN_FOR_FIELD_COLLECTION` only after all of the following are resolved:

1. current occurrence is supported by dated evidence;
2. conservation restrictions are reviewed;
3. the relevant land manager is identified or shown not to apply;
4. collection permission is approved or documented as not required;
5. tissue collection permission is approved or documented as not required;
6. the flowering/field window is checked;
7. a deidentified `target_locality_id` is assigned;
8. the exact site is stored in a private field record outside the public repository.

## No sensitive coordinates in the public repository

The public repository must never contain precise coordinates or rare-species microhabitat directions. The public queue records only range-sector evidence, permission/review status, and `target_locality_id`. Coordinates remain in a private field ledger linked by that ID.

This is especially important for narrow endemic, island, alpine and conservation-sensitive taxa.

## P6 is deliberately delayed

For P6 the Moreyra tree already contains one usable public individual. The purpose of the new sample is geographic and phenotype-linked complementarity. Therefore the public accession locality must be audited first. Only then is the new own range sector chosen to maximize nonredundant coverage within the NMNS-documented distribution.

`C. dipsacolepis` and `C. lineare` remain operational exceptions in that their own representative is drawn from the focal RAD population banks; the selected RAD population must still pass current-occurrence and permission gates.

## Freeze order

Freeze sites in biological hole-fill priority, not travel convenience:

1. P1 whole-block gaps;
2. P2 broken public tips;
3. P3 sparse blocks;
4. P4 partial blocks;
5. P5 local gaps;
6. P6 singleton complements after public-locality audit;
7. P7 optional trait-linked complement.

Logistics may combine nearby frozen sites into one trip, but travel efficiency cannot promote a lower-priority sample ahead of an unresolved higher-value hole.

## Generated queue fields

The queue records:

- immutable slot ID;
- priority/species/sample slot;
- required versus optional status;
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

The exact queue is reproducibly generated in CI and audited to contain 229 rows, 228 required slots, no pre-frozen records, and no coordinate fields.
