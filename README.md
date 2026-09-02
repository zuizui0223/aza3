# aza3 — own-data discrimination after EAzami Chapter 2

This repository starts Chapter 3 from the uncertainty left by the completed public-data Chapter 2. It is not a second copy of the Chapter 2 paper and does not treat future RAD-seq, genomic, transcriptomic or field data as retroactive confirmation.

## Locked starting point

EAzami Chapter 2 and its complete meta/simulation disposition are frozen at merge `62fa8c5c913c2b236e710f6bad366e80676aa78f`.

Within its admitted public topology ensemble:

- orientation requires **at least four** state changes;
- phyllary posture requires **at least three**;
- stickiness requires **at least five**;
- minimum counts are better resolved than individual event placements;
- **species-tip compression** hides state multiplicity in 4/4 audited polymorphic systems, while direct morph-genotype linkage exists in only 1/4.

These are topology-conditioned lower bounds, not counts of independent origins, convergence events, rates or adaptations.

## Chapter 3 architecture — one question, breadth and depth

> Which histories and causal paths that remain compatible with Chapter 2 can be discriminated by own ancestry data linked to phenotype, cytotype and separately authorized molecular or functional evidence?

### Breadth — P01-P05

The frozen P01-P05 core asks whether apparent species-tip transitions in orientation, phyllary posture, stickiness and related morphology survive when analysis is returned to linked individuals, populations and nuclear ancestry.

The claim-backward biological panel is:

- 38 concepts × 3 wild individuals = **114** all-Japan floor;
- JPN36 = 30 minimum / 40 recommended;
- JPN06 = 16 minimum / 24 recommended;
- JPN15 = 16 minimum / 24 recommended;
- total = **167 minimum** or **193 recommended** primary individuals;
- P04 reuses those individuals;
- P05 initially reuses calibrated images rather than opening another fresh-collection quota.

### Depth — M01 embedded worked case, not P06

M01 is an **embedded worked case** of the same identifiability problem, **not P06** and not a separate Ryukyu-thistle dissertation theme.

The white versus bluish-purple floral-pigmentation contrast in `Cirsium brevicaule` and `C. irumtiense` is decomposed through:

1. ancestral-state / transition history;
2. pigment chemistry;
3. same-individual corolla expression;
4. population-genomic evidence consistent with selection after background controls.

Loss, regain/re-evolution and selective-agent attribution are output-dependent claims, not assumptions.

M01 discovery uses Okinawa Honto + Amami Oshima for `C. brevicaule` and Miyako + Ishigaki for `C. irumtiense`, 15 primary individuals/population. RNA-seq is nested at 5 primary individuals/population = **20 libraries**. E3 expands only after E1/E2 gates pass, to 8 populations ×15 = **120 primary focal individuals**.

## RAD-seq is an instrument, not the conclusion

`docs/CHAPTER3_RADSEQ_END_TO_END_V1.md` freezes the RAD workflow from DNA to inference.

The principal correction is that **one same-library all-Japan RAD matrix is conditional, not guaranteed**.

Before production, an empirical pilot uses **8 anchor systems / 16 DNA templates** spanning P01/P02, diploid and polyploid concepts, and the `C. brevicaule`–`C. irumtiense` genome-size contrast. At least three candidate enzyme combinations are screened; the best two are taken through a reproducibility pilot. Enzyme pair, size window, read target, Stacks parameters and QC gates are frozen before trait outcomes are inspected.

If one protocol does not work across the anchor panel, production is stratified and different RAD protocols are **not concatenated into one SNP matrix**. The Moreyra/Comp1061 target-capture framework remains the primary cross-species scaffold.

RAD roles are bounded:

- P01/P02/M01: population ancestry, structure and admixture/network sensitivity;
- P03/P04: a strict all-Japan RAD topology/network **secondary sensitivity** only if shared-locus gates pass;
- mixed ploidy: within-cytotype or ploidy-aware probabilistic inference, never one pooled diploid caller;
- M01 E3: RAD supplies population structure/demographic background, but **RAD-only outlier scans are not final selection evidence**. Replicated E2 candidates should move to a qualified non-RAD confirmation route such as target capture, amplicon resequencing or low-pass/whole-genome resequencing.

Primary de novo assembly is Stacks 2. `m/M/n`, missing-data filters and locus filters are optimized on outcome-blind technical metrics and then frozen. Trait states cannot be used to tune the RAD pipeline.

The blank `data/intake/chapter3_radseq_library_intake_v1.csv` schema records DNA QC, extraction/library batch, protocol/index, sequencing run, read retention, locus depth/missingness and technical-replicate concordance. It must remain empty until production RAD is separately opened.

## Operational acquisition

`docs/CHAPTER3_OPERATIONAL_SAMPLE_LEDGER_V1.md` contains one row for every Japan38 concept plus the M01 focal/outgroup populations.

`JPN_29`, `JPN_31` and `JPN_33` remain identity-blocked and cannot be replaced by convenience taxa. `JPN_32`, `JPN_34`, `JPN_35`, `JPN_36`, `JPN_37` and `JPN_38` require wild-Japan provenance repair for the Chapter 3 biological panel.

## Start here

1. `docs/CHAPTER3_SCOPE_AND_HANDOFF_V1.md` — Chapter 3 question and claim boundaries.
2. `docs/CHAPTER3_INVERSE_SAMPLING_DESIGN_V1.md` — claim-backward biological sample counts and assay gates.
3. `docs/CHAPTER3_OPERATIONAL_SAMPLE_LEDGER_V1.md` — concept/population acquisition ledger.
4. `docs/CHAPTER3_RADSEQ_END_TO_END_V1.md` — RAD pilot, library design, QC, Stacks, missingness, ploidy and inference boundaries.
5. `data/planning/chapter3_radseq_end_to_end_v1.json` — machine-readable RAD contract.
6. `data/planning/chapter3_radseq_pilot_anchor_ledger_v1.csv` — eight-system RAD assay-pilot anchors.
7. `data/intake/chapter3_radseq_library_intake_v1.csv` — empty production-library QC/batch/read/locus ledger.
8. `data/planning/chapter3_inverse_sampling_design_v1.json` — machine-readable inverse sampling design.
9. `data/planning/chapter3_core_operational_sample_ledger_v1.csv` — Japan38 operational ledger.
10. `data/planning/m01_operational_population_ledger_v1.csv` — M01 focal/outgroup ledger.
11. `data/planning/chapter3_sampling_priorities_v1.csv` — frozen P01-P05.
12. `data/planning/chapter3_mechanistic_demonstration_v1.json` — M01 evidence ladder and claim ceilings.
13. `docs/M01_FLORAL_PIGMENTATION_MECHANISTIC_DEMONSTRATION_V1.md` — readable M01 design.
14. `data/intake/chapter3_individual_intake_v1.csv` — currently empty same-individual biological intake ledger.

## Current authorization state

Core:
- own biological data admitted: **0**;
- field execution authorized: **false**;
- tissue collection authorized: **false**;
- sensitive coordinates permitted in this repository: **false**;
- definitive Japan-wide species-tree claim authorized: **false**.

M01:
- own M01 biological data admitted: **0**;
- tissue collection authorized: **false**;
- field manipulation authorized: **false**;
- regain claim authorized: **false**;
- selection claim authorized: **false**;
- pollinator-agent claim authorized: **false**.

RAD production:
- RAD pilot completed: **false**;
- enzyme pair frozen: **false**;
- production RAD authorized: **false**;
- RAD library records admitted: **0**;
- M01 RAD-only selection claim authorized: **false**.

Run the Chapter 3 validators before accepting any change to these boundaries.
