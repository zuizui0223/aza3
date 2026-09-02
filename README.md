# aza3 — own-data discrimination after EAzami Chapter 2

This repository starts Chapter 3 from the uncertainty left by the completed public-data Chapter 2. It is not a second copy of the Chapter 2 paper and does not treat future RAD-seq, genomic, transcriptomic or field data as retroactive confirmation.

## Locked starting point

EAzami Chapter 2 and its complete meta/simulation disposition are frozen at merge `62fa8c5c913c2b236e710f6bad366e80676aa78f` ([core PR #129](https://github.com/zuizui0223/EAzami/pull/129); [completeness PR #130](https://github.com/zuizui0223/EAzami/pull/130)). Within its admitted public topology ensemble:

- orientation requires **at least four** state changes;
- phyllary posture requires **at least three**;
- stickiness requires **at least five**;
- minimum counts are better resolved than individual event placements;
- species-tip coding hides state multiplicity in 4/4 audited polymorphic systems, while direct morph-genotype linkage exists in only 1/4.

These are topology-conditioned lower bounds, not counts of independent origins, convergence events, rates or adaptations.

## Chapter 3 architecture — one question, breadth and depth

The primary Chapter 3 question remains:

> Which histories and causal paths that remain compatible with Chapter 2 can be discriminated by own Japan-wide ancestry data linked to phenotype, cytotype and separately authorized functional or molecular evidence?

Chapter 3 now uses the same identifiability problem at two scales rather than creating a separate Ryukyu-thistle project.

### Breadth: core history-discrimination programme P01-P05

The first core product is an all-Japan same-library RAD-seq topology/network sensitivity, not an unconditional species tree. Every genomic record must link to an immutable individual, a voucher or diagnostic image, phenotype states, cytotype status and deidentified authorization records.

The five ranked history discriminators in `chapter3_sampling_priorities_v1.csv` remain frozen as P01-P05. They ask whether apparent species-tip transitions in orientation, phyllary posture, stickiness and related morphology survive when the analysis is returned to linked individuals, populations and nuclear ancestry. They are not reordered or replaced by later mechanistic examples.

### Depth: M01 embedded mechanistic case — floral pigmentation

M01 is an **embedded worked case of the same Chapter 3 problem**, **not P06** and not a sixth Chapter 2-derived history priority. The core asks **where a trait-state change remains identifiable** after species-tip compression is reduced; M01 asks **what one such trait-state contrast consists of biologically** when followed below the topology level.

M01 decomposes one conspicuous contrast into:

1. ancestral-state and transition history;
2. homologous pigment chemistry;
3. same-individual corolla gene expression and regulatory mechanism;
4. population-genomic evidence that is consistent with selection after demographic and mapping-bias controls.

The worked example is the white versus bluish-purple corolla contrast between `Cirsium brevicaule` and `C. irumtiense`. The taxa are therefore an empirical example, not the primary taxonomic subject of Chapter 3. Published 2026 phylotranscriptomic work is treated only as a literature premise for their lineage separation and reported genome-size contrast; M01 does not make basic species delimitation its primary estimand.

M01 freezes three competing histories rather than assuming re-evolution: pigmented ancestor followed by loss/reduction, unpigmented ancestor followed by secondary gain/regain, or a more complex/unresolved history. `regain` is therefore an output-dependent claim, not the starting hypothesis.

M01 does **not** require manipulative pollination experiments or pollinator observation for its E0-E3 evidence ladder. Genomics alone cannot identify butterflies, flies, bees, climate or another ecological factor as the selective agent. Failure, ambiguity or cancellation of M01 limits only the depth demonstration; it does not invalidate or reorder P01-P05.

## Claim-backward sampling

Sampling is frozen by inference rather than by a single maximum field count. `docs/CHAPTER3_INVERSE_SAMPLING_DESIGN_V1.md` works backward from P01-P05 and M01 claim ceilings:

- all-Japan breadth floor: **38 concepts × 3 = 114** wild primary individuals;
- focal P01/P02 top-ups produce a **167-individual minimum core** or **193 recommended** core;
- P04 reuses the same individuals and P05 initially reuses calibrated images instead of opening another fresh-collection quota;
- M01 E1/E2 starts with **two populations per lineage**, 15 primary individuals per population, a 12-individual focal pigment subset and **20 primary corolla RNA-seq libraries**;
- only after E1/E2 promotion gates pass does E3 expand to **4 populations × 15 individuals × 2 lineages = 120 primary population-genomic individuals**.

Cheap leaf/photo material may be banked early, but expensive population-genomic sequencing remains closed until the previous inference layer survives. This is a distributed collection design; the investigator is not required to visit every population personally.

## Start here

1. `docs/CHAPTER3_SCOPE_AND_HANDOFF_V1.md` — the breadth-to-depth Chapter 3 question, work packages, M01 relationship and claim boundary.
2. `docs/CHAPTER3_INVERSE_SAMPLING_DESIGN_V1.md` — claim-backward core and M01 sampling, assay gates and field-burden rules.
3. `data/planning/chapter3_inverse_sampling_design_v1.json` — machine-readable inverse sampling contract.
4. `data/contracts/chapter3_eazami_handoff_contract_v1.json` — fail-closed source and authorization contract.
5. `data/planning/chapter3_sampling_priorities_v1.csv` — five ranked core history discriminators P01-P05.
6. `data/planning/chapter3_bounded_prior_registry_v1.csv` — 14 meta-analysis, programme-routing and simulation boundaries.
7. `data/planning/chapter3_protocol_registry_v1.csv` — experiment readiness without field authorization.
8. `data/planning/chapter3_mechanistic_demonstration_v1.json` — embedded M01 competing histories, evidence ladder, selection gates and stop rules.
9. `docs/M01_FLORAL_PIGMENTATION_MECHANISTIC_DEMONSTRATION_V1.md` — readable M01 research design.
10. `data/intake/chapter3_individual_intake_v1.csv` — currently empty same-individual intake ledger.

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

Run `python analysis/validate_chapter3_handoff_v1.py`, `python analysis/validate_mechanistic_demonstration_v1.py` and `python analysis/validate_inverse_sampling_design_v1.py` before accepting any change to this starting state.
