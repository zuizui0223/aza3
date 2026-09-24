# Chapter 3 RAD-seq pilot execution v1

Status date: 2026-09-02

## Decision

The production RAD programme does not start by choosing one enzyme pair from convention and applying it to the full biological panel. The first executable genomic experiment is a technical pilot designed to answer one question:

> Can one restriction protocol recover reproducible, sufficiently shared nuclear RAD loci across the focal Japanese Cirsium systems, including the mixed-ploidy and large genome-size contrasts needed by Chapter 3?

This is an assay-development question, not a biological hypothesis test. Pilot templates do not increase P01-P05 biological sample size and no trait-history result may be inspected to choose the protocol.

## Stage A — three-protocol screen

Three enzyme combinations are preregistered as **candidates, not winners**:

1. `AvaII + MspI`;
2. `EcoRI + MspI`;
3. `PstI + MspI`.

Plant ddRAD comparisons report all three as plausible plant combinations, with `AvaII + MspI` providing high tag recovery across a broad simulated angiosperm genome-size range. The initial 400–700 bp window follows that plant-method screen only as a **pilot starting window**. It is not yet the production size-selection contract.

The eight assay anchors already frozen in `chapter3_radseq_pilot_anchor_ledger_v1.csv` are:

- JPN03 `C. alpicola` — high-ploidy extreme;
- JPN05 `C. aomorense` — diploid benchmark;
- JPN06 `C. dipsacolepis` — P02 focal;
- JPN15 `C. lineare` — P02 focal counterpart;
- JPN35 `C. nipponicum` — tetraploid benchmark;
- JPN36 `C. sieboldii` — P01 focal;
- `C. brevicaule` — larger-genome M01 lineage;
- `C. irumtiense` — smaller-genome M01 lineage.

Stage A uses **one qualified DNA template from each anchor under each of the three candidate protocols = 24 shallow pilot libraries**.

All 24 are sequenced together or otherwise balanced so protocol cannot be confused with lane or run. Protocols are compared only after equal-read downsampling.

A candidate is hard-failed if it repeatedly fails library formation in a mandatory anchor, shows persistent adapter/off-window dominance, or selectively loses one of the two M01 genome-size anchors under otherwise adequate DNA QC.

The surviving protocols are ranked using technical metrics only:

- retained reads;
- usable loci after equal-depth downsampling;
- dispersion in locus recovery among the eight anchors;
- loci shared by 6/8, 7/8 and 8/8 anchor concepts;
- fragment-size consistency;
- depth evenness and the high-depth repeat/paralog tail.

The best two technically viable protocols advance. If fewer than two are viable, Stage B does not open: the candidate set or size window is redesigned first.

## Stage B — reproducibility and one-protocol viability

Stage B uses **two independent biological DNA templates per anchor × two advanced protocols = 32 primary pilot libraries**.

Add **8 preregistered technical repeat libraries**, one sentinel repeat for each anchor system. Four repeats are assigned to the Stage-A TOP1 protocol and four to TOP2, and every repeat is prepared in an independent library batch from the original. Total planned Stage-B libraries = **40**.

The 16 biological templates remain assay-development material. They are not a population sample and cannot be counted toward the 167/193 biological design.

### Technical reproducibility gates

For each technical repeat pair, after the frozen pilot genotype/depth filters:

- overlap genotype concordance must be **>= 0.95**;
- recovery of the original core loci must be **>= 0.90**.

These thresholds are assay-QC gates rather than biological effect thresholds. Published ddRAD validation shows that call-rate/depth filtering can raise genotype concordance to approximately 95%, and repeated ddRAD runs using the same laboratory protocol have recovered >90% of core loci across runs.

A pair that fails because of a documented rebuildable laboratory incident may be rebuilt once. Repeated unexplained failure is a protocol/batch failure, not a reason to lower the gate.

## Choosing one protocol versus stratifying

A single Chapter 3 production protocol is admitted only if:

1. all eight anchor systems produce interpretable libraries;
2. all preregistered technical repeats pass or have a single documented, successful technical rebuild;
3. the larger-genome, smaller-genome and high-ploidy classes do not show a protocol-specific collapse after equal-depth normalization;
4. the 6/8, 7/8 and 8/8 anchor shared-locus curves justify attempting a strict all-Japan matrix.

If exactly one protocol passes, freeze it for production.

If both pass, choose the one with stronger shared-locus recovery and lower between-anchor dispersion. Trait states, P01/P02 topology, M01 colour and any candidate biological result are prohibited decision inputs.

If neither protocol works across all anchor classes but coherent subsets work, use **protocol stratification**. For example, a Japan38 core stratum and an M01 Ryukyu stratum may use different protocols, but loci from different restriction protocols are not concatenated into one SNP matrix. The target-capture scaffold then provides cross-stratum organismal history.

## Sequencing depth is learned, not guessed

Stage-B reads are downsampled to 20%, 40%, 60%, 80% and 100% of available reads. Track:

- usable loci;
- r80 loci inside focal/population-like subsets;
- anchor shared-locus counts;
- median locus depth;
- technical concordance.

A production read target is frozen only when a locus-recovery plateau is visible. Operationally, choose a depth at or above the first point where further sequencing adds <10% usable loci for the relevant analysis stratum and technical concordance is already above the gate.

If no plateau is approached, the correct response is to top up the pilot or reduce assay complexity, not to sequence hundreds of production individuals at an unvalidated depth.

## Stacks optimization during the pilot

The Paris et al. r80 logic is used **within focal/population-like subsets**, where broad locus replication has a sensible biological meaning.

Evaluate:

- `m = 2–6`;
- `M = 1–6`;
- after choosing M, `n = M-1, M, M+1` where nonnegative.

For the deeply divergent Japan38 anchors, do not maximize a single r80 number across all concepts. Instead report homologous catalog-locus recovery by **anchor-concept occupancy**: 6/8, 7/8 and 8/8. This separates population-level Stacks optimization from the separate question of whether an all-Japan shared-locus sensitivity is even feasible.

The primary Stacks settings, paralog/homeolog filters and analysis-specific occupancy matrices are frozen before any P01-P05 trait-history result is inspected.

## Batch design

Every production stratum inherits the pilot anti-confounding rules:

- population, taxon, colour/trait state and cytotype are interleaved across extraction and library batches;
- focal contrasts cannot occupy separate lanes by design;
- 5–10% of production libraries are outcome-blind technical duplicates in different batches;
- low-yield resequencing is allowed only for a documented technical yield failure;
- technical replacement never changes biological sample identity.

`data/intake/chapter3_radseq_library_intake_v1.csv` records the full chain from individual/tissue and DNA QC to enzyme protocol, batch, read yield, locus recovery and technical concordance.

## Production opening checklist

Production RAD remains closed until all of the following are frozen:

- Stage A completed;
- Stage B completed;
- technical-repeat concordance passed;
- enzyme pair and production size window;
- production read target;
- Stacks m/M/n and primary filters;
- focal and all-Japan occupancy-matrix rules;
- ploidy handling for every admitted concept or explicit exclusion of unresolved cytotypes.

## What this pilot cannot prove

A successful pilot proves that a technical RAD protocol is reproducible enough for the specified ancestry analyses. It does not prove:

- a definitive Japanese species tree;
- absence of restriction-site dropout;
- adaptation;
- selection on an M01 candidate gene;
- any pollinator or climate selective agent.

Those remain bounded by the parent RAD and M01 contracts.

## Current state

- Stage A authorized: **false**;
- Stage A completed: **false**;
- Stage B authorized: **false**;
- Stage B completed: **false**;
- winner protocol frozen: **false**;
- production open: **false**.
