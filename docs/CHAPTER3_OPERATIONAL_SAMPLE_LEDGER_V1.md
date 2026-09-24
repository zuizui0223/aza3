# Chapter 3 operational sample acquisition ledger v1

Status date: 2026-09-02

## Purpose

This document turns the claim-backward sample counts into an acquisition ledger without inventing exact localities that are not supported by the frozen sources. Exact coordinates, permit documents and sensitive locality details remain outside the public repository.

The machine-readable ledgers are:

- `data/planning/chapter3_core_operational_sample_ledger_v1.csv` — one row for every frozen Japan38 paper concept;
- `data/planning/m01_operational_population_ledger_v1.csv` — focal Ryukyu populations, ancestral-state brackets and optional molecular controls.

## Rule 1 — target count is not admission

The core target remains 167 minimum primary individuals and 193 recommended primary individuals, but those counts are not automatically admissible. Every individual must pass identity, current-occurrence, rights, conservation and tissue-quality gates.

Three Japan38 concepts have explicit pre-collection identity blocks in the current ledger:

- `JPN_29 Cirsium verutum` — an explicitly Japanese voucher exists, but the taxon-concept identity remains contradictory;
- `JPN_31 Cirsium yuki-uenoanum` — the paper-Japan concept conflicts with a Ukraine BioSample locality and cannot be used as a clean Japanese tip;
- `JPN_33 Cirsium effusum` — cultivated Japanese provenance exists, but the name/voucher concept remains taxonomically problematic.

If any of these remain unresolved, no convenience substitute is added after the topology consequences are known. The breadth analysis reports the missing concept and runs the corresponding sensitivity rather than silently restoring n=38.

## Rule 2 — not all Japan38 concepts have the same acquisition priority

### Tier A — focal population replication

These require deliberate multi-population sampling beyond the 3-individual breadth floor:

- `JPN_36 C. sieboldii`: minimum 30, recommended 40, across 3–4 wild Japanese populations; P01 phyllary history;
- `JPN_06 C. dipsacolepis`: minimum 16, recommended 24, across 2–3 verified populations; P02 stickiness history;
- `JPN_15 C. lineare`: minimum 16, recommended 24, across 2–3 verified populations; P02 stickiness history.

These samples also carry P04 same-individual orientation, phyllary and stickiness measurements.

### Tier B — wild-Japan provenance repair

The published reference exists, but it does not satisfy the new wild-Japan breadth floor:

- `JPN_32 C. buergeri` — current reference is cultivated; source audit documents Japanese distribution in Honshu (Shiga/Yamaguchi) and Kyushu;
- `JPN_34 C. microspicatum` — current reference is cultivated; wild Honshu sample required;
- `JPN_35 C. nipponicum` — current reference is cultivated; wild northern-Honshu sample required and cytotype must be explicit;
- `JPN_36 C. sieboldii` — current reference is cultivated outside Japan; already Tier A for P01;
- `JPN_37 C. kamtschaticum` — current reference is Russian; a wild Hokkaido sample is required for Japanese population history;
- `JPN_38 C. pendulum` — current reference is Russian; a wild Japanese sample is required for the breadth floor.

### Tier C — identity-aware routine breadth sampling

Most JPN01–JPN30 concepts already have a direct-Japan Moreyra sample or public locality. They still need three newly linked wild individuals for the own-data sensitivity panel, but exact populations are not guessed from the current repository. For these taxa, the operational order is:

1. authority-backed paper-concept reconciliation;
2. current-occurrence verification;
3. land/tissue/conservation authorization;
4. freeze one or more deidentified population IDs;
5. collect the three linked individuals.

Minor name or tree-code discrepancies for JPN02, JPN10, JPN12, JPN14, JPN19, JPN20 and JPN21 are recorded in the ledger and require voucher-linked confirmation rather than silent normalization.

## Rule 3 — one individual package for the breadth panel

Every admitted breadth individual uses the same immutable ID for:

- taxon concept and population ID;
- voucher or diagnostic images;
- standardized head orientation;
- phyllary posture and calibrated phyllary image;
- stickiness state, with gland/exudate documentation for JPN06/JPN15;
- standardized flower-colour image;
- cytotype status, including explicit missingness;
- leaf DNA / genomic library ID;
- deidentified authorization and conservation IDs.

Technical library duplicates are not biological replicates.

## M01 discovery — the four populations to prioritize personally

The discovery layer is deliberately smaller than the final selection panel.

### Cirsium brevicaule

- `M01_BREV_OKI` — Okinawa Honto, 15 primary + 3 predeclared reserves;
- `M01_BREV_AMAMI` — Amami Oshima, 15 primary + 3 predeclared reserves.

### Cirsium irumtiense

- `M01_IRUM_MIYAKO` — Miyako, 15 primary + 3 predeclared reserves;
- `M01_IRUM_ISHIGAKI` — Ishigaki, 15 primary + 3 predeclared reserves.

Ishigaki is frozen as the discovery Yaeyama population to reduce the temptation to choose between Ishigaki and Iriomote after seeing an outcome. Iriomote remains an E3 expansion population. A change is allowed only before focal data inspection and only for access/authorization or current-occurrence failure, with the reason recorded prospectively.

Each discovery population banks 15 same-individual phenotype/photo/DNA samples. Nested within that bank:

- pigment chemistry: 3 primary individuals/population = 12 focal chemistry samples;
- floral RNA: collect 6/population and sequence 5 primary/population = 20 focal RNA-seq libraries at one frozen developmental stage.

No bagging, artificial-flower experiment, single-visit pollination assay or long camera deployment is required.

## M01 E1 ancestral-state bracket — exact taxa

Chang et al. (2026; DOI `10.1186/s12870-026-08097-6`) resolves subsect. Arenicola as sister to subsect. Nipponocirsium and both as sister to Sinocirsium. The minimum six-concept nonfocal bracket is therefore frozen as follows:

### Immediate sister-clade coverage — Nipponocirsium

- `C. morii`;
- `C. pengii`;
- `C. tatakaense`;
- `C. kawakamii`.

Each requires at least three phenotype-confirmed floral individuals for the E1 chemistry/colour layer. Existing transcriptomic topology can be reused; new field sequencing is not required merely to rebuild the backbone.

### Deeper state-balanced bracket — Sinocirsium

- `C. japonicum var. albescens` — published floral diagnosis is white;
- `C. japonicum var. fukienense` — published floral diagnosis is bluish-purple / lighter purple.

These two are deliberately state-balanced bracketing taxa, not a random sample used to estimate colour prevalence or transition rate.

Optional pipeline controls are the white and bluish-purple morphs of `C. japonicum var. takaoense`, three phenotype-confirmed individuals per morph. They are not counted toward the six-concept E1 minimum. Their role is to check whether the pigment/expression workflow can recover a published within-lineage colour contrast before over-interpreting the Ryukyu species comparison.

The Taiwan material should preferentially be obtained through a collaborator or documented fresh material source; it does not require the investigator to add a separate Taiwan field campaign.

## M01 E3 expansion — only after the E2 gate

The final E3 primary panel is 120 focal individuals:

### C. brevicaule — 4 populations × 15

1. Amami Oshima — already in discovery;
2. one verified intermediate Amami-group population;
3. one verified southern Amami-group population;
4. Okinawa Honto — already in discovery.

### C. irumtiense — 4 populations × 15

1. Miyako — already in discovery;
2. Ishigaki — already in discovery;
3. Iriomote;
4. Yonaguni.

Thus the discovery bank contributes 60 of the final 120 primary focal individuals. If E3 is promoted, the default incremental burden is approximately 60 additional primary leaf-DNA/photo individuals, not a second independent 120-individual collection.

E3 sequencing does not open until history is interpretable, the E2 expression mechanism replicates across two populations per lineage, and orthology/reference/mapping-bias qualification is credible. RAD-only two-species FST is not sufficient evidence of selection.

## Investigator versus collaborator burden

The investigator's personal effort is prioritized to samples where handling consistency matters most:

1. four M01 discovery populations for same-stage floral RNA and pigment material;
2. difficult P01/P02 focal populations where phyllary/stickiness phenotyping must be standardized;
3. any identity-gate population that requires expert voucher documentation.

Routine breadth-floor leaf DNA and image packages, and E3 terminal-island leaf collections, may be supplied by standardized collaborators if individual identity, authorization and chain-of-custody requirements are met.

## What is deliberately not frozen here

For most Japan38 concepts the current repository does not provide a defensible population-level locality. This document therefore does not invent one. Exact sites are admitted only through a later current-occurrence/permission ledger and remain deidentified in the public repository.
