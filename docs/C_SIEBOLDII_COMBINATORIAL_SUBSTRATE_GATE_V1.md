# aza3 Gate CS — Cirsium sieboldii combinatorial-substrate verification v1

**Status:** REQUIRED BEFORE OWN-WGS SAMPLE SELECTION — OBSERVATION-ONLY DESIGN  
**Date:** 2026-10-05  
**Focal taxon:** *Cirsium sieboldii*

## Why this gate exists

The active Nature question is:

> **Must an integrated complex organ evolve as an integrated unit?**

For *C. sieboldii* to provide the decisive within-species test, more is required than species-level authority records showing that several states occur somewhere in the taxon.

The Nature mechanism requires a **segregating/reassemblable substrate**:

- at least two capitulum components must vary among directly observed wild individuals;
- their variation must not be explained by flower developmental stage alone;
- the components must not be perfectly confounded with distant geography, taxonomic misidentification or cytotype;
- preferably, different component combinations must coexist within one population or within a tightly matched local population set.

At present, public evidence establishes potential, not this condition.

## Public evidence ceiling

Current authoritative/public evidence supports:

1. *C. sieboldii* normally flowers with nodding capitula, but upright-flowering plants have been reported from multiple places;
2. normal plants become upright after flowering, making phenology a critical confound;
3. a white-flowered form, *C. sieboldii* f. *leucanthum*, is formally recognized;
4. herbarium comparison shows intraspecific variation in capitulum size and phyllary-row number;
5. no source recovered in the current audit demonstrates that white/coloured, anthesis-upright/nodding, phyllary or stickiness states independently co-segregate within the same natural population.

Therefore:

> **species-wide state multiplicity is not evidence of a combinatorial population substrate.**

## Public-locality prescreen and CS0 priority

A bounded public-source prescreen now ranks localities for **observation only**. It does not classify Gate CS.

### Priority 1 — Kurumayama Moor / Kirigamine, Nagano

Reason:
- type locality of *C. sieboldii* f. *leucanthum*;
- current public locality information still lists Shirobana-kiseruazami from the Kurumayama area.

Primary CS0 task:
- verify whether white and non-white flowering individuals coexist;
- measure A1–A2 orientation continuously on all accessible flowering individuals;
- do **not** infer orientation polymorphism from the white-form record itself.

### Priority 2 — Tsukude Plateau wetland local set, Aichi

Reason:
- multiple historical *C. sieboldii* herbarium records from former Tsukude-mura;
- local field accounts describe abundant wetland populations and within-region morphological variation.

Primary CS0 task:
- measure A1–A2 orientation distribution;
- search prospectively for colour variation;
- treat hairiness variation only as evidence that local individual variation exists, not as a capitulum-module result.

Public-source ceiling:
**neither locality currently demonstrates same-population colour × anthesis-orientation co-segregation.**

The machine-readable prescreen is:
`data/evidence/c_sieboldii_public_locality_prescreen_v1.json`.

## Gate question

> Does one natural *C. sieboldii* population, or a tightly matched local population set, contain directly observed individuals that independently vary in at least two capitulum components at the same flowering stage?

This is an observation/verification gate, not a fitness experiment and not a genomic association study.

## Stage standardization — mandatory

Orientation is only admissible when the head is in a frozen flowering-stage class.

### A0 — bud
No open florets. Not admissible for the primary orientation state.

### A1 — early anthesis
Fresh florets open; no visible mature pappus/fruiting transition.

### A2 — full anthesis
Substantial fresh floret presentation; no obvious post-flowering head reorientation.

### A3 — late anthesis
Some senescence but still active flowering; no mature dispersal structures.

### P — post-anthesis / fruiting
Flowering largely finished, head reorientation toward upright may occur.

**Primary orientation inference uses A1–A2 only.**

A3 is sensitivity-only.

P is never coded as evidence for an upright-flowering morph.

Where practical, mark and re-photograph the same head through time. A head that is nodding at A1–A2 and upright only at P is phenological, not an upright-flowering phenotype.

## Direct component measurements

Every observed individual receives an immutable observational ID and same-individual records.

### C — colour
Primary:
- standardized visible photograph with calibration target where feasible;
- floral CIELAB/chroma/lightness from fresh florets;
- binary white/non-white state retained only as a descriptive derivative.

### O — orientation
Primary:
- gravity-referenced capitulum axis angle at A1–A2;
- photograph must show gravity reference or a standardized vertical frame.

Frozen angular convention:

- **0° = vertically upward**;
- **90° = horizontal**;
- **180° = vertically downward**.

The continuous angle is primary.

For the transparent 2×2 gate summary only:

- descriptive `upright` = angle <= **60°**;
- descriptive `nodding` = angle >= **120°**;
- 60–120° = intermediate and excluded from binary-cell occupancy while retained in continuous analyses.

Do not replace angle with observer labels such as “upright” or “nodding” in the primary analysis.

### P — phyllary / involucre
Record, without requiring it to vary:
- phyllary row number;
- outer-phyllary posture/angle;
- involucre diameter and length;
- standardized diagnostic image.

### S — stickiness
Record:
- predeclared touch/adhesion score or standardized surface assay;
- visible gland/exudate status;
- do not infer stickiness from species authority descriptions.

## Predeclared component-pair hierarchy

The gate does not search all possible trait pairs after seeing the data.

### Primary pair
**C × O: floral colour × anthesis orientation**

Reason:
both state dimensions have prior public evidence of species-wide variation and can be scored on the same flowering individual.

### Secondary pair 1
**O × P: anthesis orientation × quantitative phyllary/involucre state**

Opened if P shows sufficient directly measured individual variation under the frozen measurement protocol.

### Secondary pair 2
**C × P: floral colour × quantitative phyllary/involucre state**

Opened under the same condition.

### Exploratory only
Stickiness combinations.

Stickiness cannot rescue a failed primary/secondary gate unless an independent future protocol upgrades it prospectively.

## Population verification design

### Phase CS0 — locality verification

Before DNA collection, visit candidate populations/local population sets and record:

- taxonomic identity;
- population extent;
- flowering individual count;
- flowering-stage distribution;
- conservation/access status;
- whether white/non-white individuals are present;
- whether A1–A2 orientation spans meaningful variation;
- whether phyllary quantitative variation is observable.

No tissue collection is implied by this gate.

### Phase CS1 — combinatorial census

For a candidate population:

- target at least **30 flowering genets/individuals** if population size permits;
- primary gate requires at least **20 A1–A2 individuals** with complete C and O records;
- avoid counting multiple ramets from one obvious genet as independent when clonality can be recognized;
- spatially spread observations across the population rather than choosing striking morphs.

If no single population is sufficiently large, use a predeclared local population set only when sites are geographically/ecologically close enough to avoid turning the gate into a broad geographic comparison.

## Quantitative gate logic

### Step 1 — does each component actually vary?

For C and O separately, the gate requires nontrivial individual-level spread.

Orientation:
- A1–A2 angle range must exceed **30 degrees** after excluding obvious measurement failures;
- at least 3 A1–A2 individuals must satisfy the frozen descriptive upright threshold (<=60°) and at least 3 must satisfy the nodding threshold (>=120°) for orientation to contribute to the binary combination gate.

Colour:
- continuous colour spread must exceed the frozen technical-error envelope;
- a white/non-white contrast is only counted if each state has at least **3 individuals** in the same candidate population/local set.

These are substrate-verification thresholds, not evolutionary effect sizes.

### Step 2 — are multiple combinations present?

For C × O, a decisive within-population substrate requires:

- both C and O pass Step 1;
- at least **3 of the 4 descriptive C × O cells** occur;
- every counted occupied cell contains at least **3 individuals**;
- no occupied-cell pattern is created solely by A3/P phenology.

The continuous C and O measurements remain the primary stored data; the 2×2 display is only a transparent gate summary. Intermediate-orientation individuals are retained in the continuous correlation and range calculations but do not create binary cells.

### Step 3 — is coupling incomplete?

If C × O variation exists but the two dimensions are effectively deterministic functions of one another, the system is not yet a strong reassembly substrate.

Primary criterion:
- absolute Spearman correlation between continuous C and O < **0.8**;
- and neither colour state is confined to a single orientation class with complete separation.

This is an engineering gate, not a null-hypothesis test of independence.

### Step 4 — obvious confounds

Automatic downgrade if combinations track:

- flowering stage;
- taxonomic mixture;
- gross habitat discontinuity within the supposed population;
- cytotype/genome-size class once measured;
- sampling date/batch in a way that cannot be separated.

## Decision classes

### CS_GREEN — within-population combinatorial substrate verified

Required:

- one natural population passes primary C × O Step 1–3;
- or a predeclared secondary pair passes equivalently while C × O remains informative but rare;
- phenotype combinations are not explained by phenology/taxonomic mixture.

**Action:** use that population as the first-choice source for the 8-individual R1B own-WGS pilot and later natural-reassembly tests.

Nature implication:
a same-population segregating substrate exists, so recombinant/local-ancestry tests can directly address reconfigurable integration.

### CS_GREEN_LOCAL_SET — tightly matched local substrate

No one population reaches CS_GREEN, but two or more nearby/ecologically matched populations jointly contain the component variation and combinations without complete trait–population confounding.

**Action:** own-WGS pilot may proceed, but the paper may not call this within-population natural reassembly until genomics demonstrates exchange/segregation.

### CS_AMBER — species-level variation only

Multiple component states are verified, but they are segregated among distant populations or strongly confounded with geography.

**Action:**
- *C. sieboldii* can still test repeated/local genomic histories;
- it is weak for the Nature Fig.4 natural-reassembly kill shot;
- verify an alternative focal system before scaling the full panel.

### CS_RED_PHENOLOGY

Apparent orientation polymorphism is explained by anthesis-to-fruiting reorientation.

**Action:** orientation is removed as a segregating component in *C. sieboldii* for the Nature mechanism unless a separately verified anthesis-upright population is found.

### CS_RED_NO_COMBINATORIAL_SUBSTRATE

After standardized observation, only one component varies, or component states remain inseparable from taxon/geography.

**Action:** do not center Nature Fig.3–4 on *C. sieboldii*. Retain it as a historical/phylogenetic system and move the within-species reassembly test to another verified system.

### NOT_IDENTIFIABLE

Population sizes, access or flowering overlap are insufficient.

**Action:** do not infer absence. Verify another locality/system.

## Relationship to the 8-individual WGS pilot

Reference transferability and biological substrate are separate gates.

The WGS reference pilot is technically feasible after R1B-0, but **sample selection for the biological Nature test is not opened until Gate CS is classified**.

If WGS must be performed before CS for logistical reasons, it remains a reference-transferability dataset only and cannot be promoted retrospectively to a combinatorial-substrate test without satisfying the phenotype gate.

## Why this gate is cheap and decisive

The critical failure mode is biological, not computational:

> the focal species may contain several states somewhere, yet no population may actually contain independently recombinable component variation.

A short, standardized observation census can detect that failure before expensive WGS, transcriptomics, fitness manipulation or focal-reference construction.

## Stop rules

- Do not code post-flowering upright heads as upright-flowering morphs.
- Do not use species-level authority states as individual phenotypes.
- Do not balance phenotype classes by sampling distant populations and then call the result segregation.
- Do not let stickiness rescue the gate post hoc.
- Do not sequence 40 *C. sieboldii* merely because the species is morphologically interesting.
- Do not claim natural reassembly until directly observed combinations and genomic segregation are both demonstrated.
