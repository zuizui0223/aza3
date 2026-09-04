# Chapter 3 nationwide species-tree sampling v4

Status date: 2026-09-02

This is the authoritative sampling architecture for new Chapter 3 acquisition. It starts from zero physical samples and treats the **Japan-wide Cirsium species tree/network as the first genomic product**, not as a late sensitivity analysis.

## 1. Why the design changes again

The previous v3 correctly rejected `Japan38 = all Japanese Cirsium`, but it over-corrected by making transition-first population genomics the first sampling product.

For Chapter 3, the broader scientific object is the Japanese radiation itself. The correct hierarchy is therefore:

1. taxonomic census;
2. nationwide nuclear species tree/network;
3. capitulum-trait history mapping on that own-data tree ensemble;
4. population RAD in selected focal lineages/transition neighbourhoods;
5. M01 mechanistic decomposition as an embedded depth case.

Moreyra's 38 Japanese species remain the Chapter 2 hypothesis-origin subset and an existing nuclear scaffold, but they are not the Chapter 3 sampling universe.

## 2. How many Japanese species?

Two current NMNS public descriptions use different scopes:

- the 2025 NMNS article describes approximately **120 species** of Japanese Cirsium;
- the NMNS thistle database currently returns **161 authority records** and its top page describes more than 150 thistles, while explicitly noting unpublished/new-name entries.

Therefore `161 records = 161 independent species` is not admitted.

### Operational rule

Use **120 species as the planning baseline**, not as a frozen taxonomic truth.

Before destructive collection, reconcile the current authority records into one operational census with explicit statuses for:

- accepted primary species concept;
- infraspecific concept;
- synonym / duplicate treatment;
- unpublished or `new name` concept;
- taxonomic conflict;
- non-Japanese/shared continental taxon;
- conservation-limited concept.

If the frozen operational primary-species count exceeds **125**, reopen the individual budget rather than reducing replication below the rules below.

## 3. Nationwide sampling target: 290–300 target-capture individuals

The simplest defensible baseline is **two independent individuals per admitted species**.

At a nominal 120 species:

> 120 species × 2 individuals = **240 base individuals**

The two individuals should come from different populations or clearly separated geographic occurrences whenever feasible. Two adjacent plants from one patch are not geographic replication.

### Enrichment tiers

Not every species needs the same additional replication.

#### T1 — standard species

Default: **2 individuals/species**.

#### T2 — widespread or variable species

Nominally up to **30 species** receive one additional individual:

> 30 × +1 = **+30**

Triggers include broad range, documented within-species capitulum polymorphism, suspected multiple cytotypes/genome-size classes, or first-wave gene-tree discordance.

#### T3 — taxonomically/genomically complex species

Nominally up to **10 species** receive two additional individuals:

> 10 × +2 = **+20**

Triggers include hybrid/introgression suspicion, polyploid/allopolyploid uncertainty, unresolved taxonomic complexes, or unusually high leverage for P01–P05/M01.

Therefore the nominal nationwide bank is:

> **240 + 30 + 20 = 290 individuals**

If the frozen operational census is 125 species, the same enrichment envelope gives:

> 125 × 2 + 50 = **300 individuals**

Hence the planning target is **about 120–125 species and 290–300 physical target-capture candidates**.

This is a capacity envelope. It is not permission to force the census to 120 or 125 species.

### Rare-species exception

If conservation or legal constraints allow only one individual, retain it as `SINGLE_SAMPLE_LIMITED`. Do not replace a rare species with a convenient relative simply to keep the matrix rectangular.

## 4. Target capture, not RAD, builds the nationwide species tree

Primary cross-species genomic method:

> **Comp1061-compatible nuclear target capture**, or an explicitly crosswalked equivalent.

This preserves compatibility with the existing Moreyra/EAzami nuclear scaffold.

### Wave 1

Sequence the two base individuals for every admitted species:

- nominal 120-species case: **240 target-capture individuals**.

Keep individual sequences separate through QC and gene-tree inference. Do not collapse both individuals to one tip before checking whether they agree.

### Wave 2

After technical QC and before inspecting focal trait-history outcomes, sequence the preregistered enrichment set:

- up to **50 additional individuals**;
- nominal full target-capture panel = **290 individuals**.

### Species-tree/network inference boundaries

The nationwide product must:

- use a coalescent-aware species-tree ensemble;
- retain gene-tree discordance;
- not assume species monophyly merely because two individuals share a name;
- retain reticulation/network sensitivities for replicate non-monophyly or strong genomic conflict;
- keep plastid history separate from the nuclear organismal scaffold;
- treat mixed ploidy explicitly rather than allowing a pooled diploid RAD matrix to define the nationwide tree.

## 5. Every tree individual is also a phenotype individual

The species-tree panel must not recreate species-tip compression.

Every newly collected sequenced individual links one immutable ID to:

- taxon determination;
- population/locality key;
- voucher or diagnostic images;
- capitulum orientation;
- phyllary posture plus calibrated image;
- stickiness / gland state;
- flower colour;
- phenological/developmental stage;
- silica/fresh leaf DNA source;
- cytotype or genome-size evidence status;
- deidentified authorization/conservation IDs.

Thus Level 1 creates the own-data tree and the individual-linked trait scaffold simultaneously.

## 6. Level 2: remap P01–P05 on the own nationwide tree

Only after the nationwide nuclear tree/network ensemble is admitted do we re-estimate:

- orientation transition count and placement;
- phyllary transition placement;
- stickiness transition placement;
- shared-history / cross-module overlap;
- continuous phyllary dimensions.

The old Japan38 placements remain priors/design history, not constraints on the own-data answer.

## 7. Level 3: population RAD is nested inside the nationwide tree

RAD answers shallower questions:

- population ancestry;
- admixture / reticulation;
- morph-history within species;
- whether a species-level transition actually decomposes into population structure or introgression.

The current target focal RAD design is:

| system | populations | primary RAD/pop | total primary RAD |
|---|---:|---:|---:|
| `C. sieboldii` | 4 | 12 | 48 |
| `C. dipsacolepis` | 3 | 12 | 36 |
| `C. lineare` | 3 | 12 | 36 |
| `C. brevicaule` | 2 | 12 | 24 |
| `C. irumtiense` | 2 | 12 | 24 |
| **total** | **14** |  | **168** |

For each focal species, choose up to four of these RAD individuals as the nationwide target-capture representatives. Do not collect a second set of `tree-only` plants if a RAD individual already has the full voucher/phenotype package.

Under the nominal full design:

- nationwide bank = **290** physical plants;
- target focal RAD = **168** plants;
- nominal overlap = **20** plants (4 per five focal species);
- additional focal plants beyond the nationwide bank ≈ **148**;
- total unique physical plants for the full national-tree + initial focal-RAD programme ≈ **438**.

This 438 is a planning total, not a requirement to complete every downstream population module before publishing the species tree.

## 8. RAD assay pilot

The RAD pilot remains nested inside the five focal systems rather than collecting unrelated assay-only taxa.

- Stage A: 5 systems ×1 DNA ×3 enzyme candidates = **15 shallow libraries**;
- Stage B: 5 systems ×2 independent individuals ×top-2 protocols =20 primary libraries +5 cross-batch technical repeats = **25 libraries**;
- callable-overlap genotype concordance gate ≥0.95;
- core-locus recovery gate ≥0.90.

If later transition mapping selects high-ploidy taxa for population RAD, run a second ploidy/complexity-stratum assay check before sequencing them.

## 9. Level 4: M01 remains embedded, not a separate dissertation

`C. brevicaule` and `C. irumtiense` are already part of the nationwide species-tree sample and the initial focal RAD panel.

M01 uses the same individuals and adds floral material in four discovery populations:

- Okinawa Honto;
- Amami Oshima;
- Miyako;
- Ishigaki.

Collect six developmentally matched floral samples/population, sequence five RNA replicates and retain one reserve.

Further islands/populations remain conditional on E1/E2 success and are not required to complete the Japanese species tree.

## 10. What is frozen now

Planning baseline:

- **~120 species**;
- **2 individuals/species minimum**;
- **240 first-wave target-capture individuals**;
- **+50 enrichment individuals**;
- **290 nominal / 300 upper-envelope target-capture individuals**;
- **168 initial focal population-RAD individuals** nested underneath the tree;
- **~438 unique physical plants** only if the full national tree and initial focal RAD programme are both completed.

What is not frozen:

- the final accepted operational species count;
- which 30 taxa receive T2 enrichment;
- which 10 taxa receive T3 enrichment;
- exact collection localities;
- final species-tree topology;
- trait-transition placements;
- any regain, selection or selective-agent claim.

## 11. Supersession

For new field acquisition, this v4 supersedes the transition-first v3, zero-baseline v2 and earlier 38×3/inventory-first designs.

The useful parts of v3 remain downstream: after the nationwide tree is built, transition-neighbourhood logic is used to decide which additional lineages deserve population RAD. It is no longer the rule for deciding which Japanese species enter the species tree.
