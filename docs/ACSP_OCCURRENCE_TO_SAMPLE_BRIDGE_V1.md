# ACSP occurrence-to-sample bridge v1

Status: **prospective bridge; no exact site frozen**

This document defines how ACSP candidate search patches enter the aza3 sampling workflow without allowing a survey-planning algorithm to redefine the biological priorities of Chapter 3.

## Direction of information

```text
aza3 master sample / tree slot
        ↓
species + hole priority + range sector + conservation mode
        ↓
ACSP occurrence audit and candidate-patch generation
        ↓
patch / abstain + evidence/feasibility diagnostics
        ↓
field current-occurrence verification
        ↓
aza3 conservation + permission + private exact-site freeze
        ↓
target capture / RAD / M01 acquisition
```

ACSP never promotes a patch directly into `target_locality_id`. A patch is only a bounded place where a search is scientifically justified.

## Frozen first-pass occurrence window

The first Cirsium ACSP bridge uses public occurrence evidence with event date <= **2025-12-31**. The purpose is to keep the initial candidate construction independent of the prospective 2026+ field era.

The ACSP contract types occurrence evidence by taxon match, spatial precision, event age and provenance. Only recent, sufficiently precise reconciled records can define the primary local-continuation anchor set. Older, obscured, uncertain and region-only records remain context/sensitivity evidence rather than pseudo-exact current populations.

## What aza3 supplies

For each required representation slot:

- immutable `aza3_slot_id`;
- operational species concept;
- P1-P6 tree-hole priority;
- A/B/OWN range sector;
- conservation/acquisition mode;
- no public sensitive exact coordinates.

The current authoritative nationwide tree requirement remains **228 required representation slots**. P02 and M01 physical plants remain nested as defined by master sample manifest v11.

## What ACSP returns

Public-safe handoff fields include:

- `acsp_patch_id`;
- `discovery_regime` (`LOCAL_CONTINUATION`, `DETACHED_COMPONENT`, `SENTINEL`, `ABSTAIN_LOCAL_PATCH`);
- range sector;
- occurrence-anchor evidence counts/classes;
- ecological support status;
- survey feasibility status;
- comparator assignment for prospective method evaluation;
- abstention/block reason.

Exact candidate coordinates for sensitive species are not committed to either public repository.

## Ecological support versus operational feasibility

ACSP v2 separates:

- `G_E`: ecological support/continuity;
- `G_F`: survey feasibility/access/permission/conservation.

A patch may be ecologically supported but inaccessible. That is not a biological non-detection and does not justify moving support toward the nearest trail.

## Field outcome handoff

The primary ACSP validation endpoint is **verified current occurrence**, not tissue collection.

A field visit must distinguish:

- completed search + verified detection;
- completed search + non-detection;
- detected but taxonomically unresolved;
- access failure;
- permission block;
- phenology not evaluable;
- incomplete search.

Search minutes and observer count are mandatory; route length/searched area are added when defensibly measured. Exact coordinates of verified populations stay in the private aza3 field ledger.

A verified plant found where tissue removal is prohibited is:

- ACSP current-occurrence success;
- aza3 collection blocked/pending.

Those two outcomes must never be collapsed.

## Concrete sampling plan after ACSP

For each aza3 slot:

1. audit pre-2026 occurrences;
2. classify the slot as local-continuation-input available or sentinel/abstain candidate;
3. build ecological candidate patches under the frozen structural family;
4. apply feasibility/conservation mask separately;
5. freeze method/comparator patch before field outcome;
6. conduct standardized current-occurrence search;
7. after verified detection, complete local conservation/permission review;
8. create private exact-site record and only then promote the slot toward collection.

For P1-P5 A/B slots, the two samples should remain independent range-sector representations whenever the species distribution allows it. For P6, the ACSP search targets the own complement after public-accession provenance auditing. P02/M01 focal population depth does not override missing nationwide-tree holes.

## Algorithm-development role

Campanula showed that generic NDVI filtering did not provide stable added recovery beyond deterministic spatial balance. The Cirsium system therefore tests whether **predeclared structural ecological continuity** — such as wetland/moisture, alpine terrain, coastal/island structure, or source-supported geology/substrate — adds prospective discovery value beyond annular distance and spatial balance.

This is an external empirical development/validation of the occurrence-to-patch line. It is not a retrospective addition to the frozen ACSP Ecological Informatics result.
