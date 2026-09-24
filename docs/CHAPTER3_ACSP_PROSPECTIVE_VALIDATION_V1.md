# Chapter 3 ACSP prospective validation cohort v1

Status: **frozen before candidate-patch construction and before any 2026+ Cirsium field outcome**

This file mirrors the prospectively selected 13-unit Cirsium validation cohort from `zuizui0223/acsp` into the biological execution repository. The mirror is intentionally public-safe: it contains range sectors and deidentified validation-unit IDs, but no exact candidate or population coordinates.

## Why this layer exists

Chapter 3 and ACSP answer different questions.

- **aza3** decides which biological tree/RAD/M01 representation is needed, whether collection is legal and appropriate, and which assay a verified plant enters.
- **ACSP** decides whether occurrence evidence justifies a bounded search patch, which ecological structural family is allowed, and which comparator is evaluated.

A successful ACSP search does not automatically authorize tissue collection. A failed access attempt is not a biological non-detection.

## Frozen cohort

The cohort contains 13 distinct species and 13 distinct aza3 slots:

- 9 `LOCAL_CONTINUATION` units;
- 4 `SENTINEL` units;
- anchor replication: 4 zero, 4 single, 5 multiple;
- method arms: 8 structural-local, 4 structural-sentinel, 1 spatial-baseline-only.

The structural families were assigned before field outcomes from source-backed habitat statements:

- wetland/moisture structure;
- alpine topographic structure;
- open grassland structure;
- coastal/island structure;
- forest-edge structure;
- one explicit general-spatial negative-control lane.

No taxon, structural family or method arm may be replaced after seeing field success.

## Three-way prospective comparison

For local-continuation units the frozen comparison is:

`source-backed structural selector` vs `annular nearest-known` vs `deterministic spatial balance`.

For sentinel units there is no eligible local anchor, so the frozen comparison is:

`source-backed structural sentinel selector` vs `validated broad robust support` vs `deterministic spatial balance`.

`Cirsium tamastoloniferum` remains the negative-control lane with no added structural ecological family; only the annular and deterministic spatial controls are retained. An ecological feature cannot be added to that unit after its result is known.

## Anchor adequacy is not raw GBIF abundance

The frozen cohort explicitly distinguishes:

- `ZERO_PRIMARY_ANCHOR`;
- `SINGLE_PRIMARY_ANCHOR`;
- `MULTIPLE_PRIMARY_ANCHORS`.

This is a key methodological result of the pre-field audit. Large retained-record counts do not imply a usable local anchor set: `C. dipsacolepis` and `C. lineare`, for example, enter as sentinel problems because no record survives the frozen local-anchor eligibility rules for their declared P02 OWN search sector. Conversely, one precise eligible coordinate is kept distinct from genuinely replicated local-anchor evidence.

## Files

- `data/planning/chapter3_acsp_prospective_validation_cohort_v1.csv` — upstream cohort mirror.
- `data/planning/chapter3_acsp_prospective_validation_contract_v1.json` — upstream provenance and freeze rules.
- `data/planning/chapter3_acsp_field_handoff_prefill_v1.csv` — field-log prefill with patch/outcome/site fields intentionally unopened.
- `analysis/validate_acsp_prospective_validation_v1.py` — fail-closed validator.
- `tests/test_validate_acsp_prospective_validation_v1.py` — regression check.

## Next execution boundary

Candidate patches may now be generated in ACSP using only the frozen pre-2026 occurrence evidence and the frozen structural-family contract. Before any field outcome is opened, the resulting method/comparator patch identifiers must be written back to the handoff table. Exact sensitive coordinates remain private.

Only after that patch freeze may a field search record detection/non-detection, search effort, access/permission failures and identity verification. P02/M01 tissue reuse is downstream of verified current occurrence and aza3 authorization.

## Claim boundary

This cohort is a prospective method test, not evidence that the 13 species occupy the proposed patches, that those patches are accessible, that ACSP outperforms its comparators, or that any Chapter 3 evolutionary hypothesis is supported.
