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

## Chapter 3 architecture

### Core history-discrimination programme: P01-P05

The primary Chapter 3 question remains:

> Which histories and causal paths that remain compatible with Chapter 2 can be discriminated by own Japan-wide ancestry data linked to phenotype, cytotype and separately authorized functional experiments?

The first core product is an all-Japan same-library RAD-seq topology/network sensitivity, not an unconditional species tree. Every genomic record must link to an immutable individual, a voucher or diagnostic image, phenotype states, cytotype status and deidentified authorization records.

The five ranked history discriminators in `chapter3_sampling_priorities_v1.csv` remain frozen as P01-P05. They are not reordered or replaced by later mechanistic examples.

### M01 mechanistic demonstration: floral pigmentation

M01 is an independent **mechanistic demonstration**, **not P06** and not a sixth Chapter 2-derived history priority. It asks whether a clearly observed trait-state contrast can be decomposed into:

1. ancestral-state and transition history;
2. homologous pigment chemistry;
3. same-individual corolla gene expression and regulatory mechanism;
4. population-genomic evidence that is consistent with selection after demographic and mapping-bias controls.

The worked example is the white versus bluish-purple corolla contrast between `Cirsium brevicaule` and `C. irumtiense`. Published 2026 phylotranscriptomic work is treated only as a literature premise for their lineage separation and reported genome-size contrast; M01 does not make basic species delimitation its primary estimand.

M01 freezes three competing histories rather than assuming re-evolution: pigmented ancestor followed by loss/reduction, unpigmented ancestor followed by secondary gain/regain, or a more complex/unresolved history. `regain` is therefore an output-dependent claim, not the starting hypothesis.

M01 does **not** require manipulative pollination experiments or pollinator observation for its E0-E3 evidence ladder. Genomics alone cannot identify butterflies, flies, bees, climate or another ecological factor as the selective agent. Failure, ambiguity or cancellation of M01 does not invalidate P01-P05.

## Start here

1. `docs/CHAPTER3_SCOPE_AND_HANDOFF_V1.md` — core question, work packages, M01 relationship and claim boundary.
2. `data/contracts/chapter3_eazami_handoff_contract_v1.json` — fail-closed source and authorization contract.
3. `data/planning/chapter3_sampling_priorities_v1.csv` — five ranked core history discriminators P01-P05.
4. `data/planning/chapter3_bounded_prior_registry_v1.csv` — 14 meta-analysis, programme-routing and simulation boundaries.
5. `data/planning/chapter3_protocol_registry_v1.csv` — experiment readiness without field authorization.
6. `data/planning/chapter3_mechanistic_demonstration_v1.json` — independent M01 competing histories, evidence ladder, selection gates and stop rules.
7. `docs/M01_FLORAL_PIGMENTATION_MECHANISTIC_DEMONSTRATION_V1.md` — readable M01 research design.
8. `data/intake/chapter3_individual_intake_v1.csv` — currently empty same-individual intake ledger.

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

Run `python analysis/validate_chapter3_handoff_v1.py` and `python analysis/validate_mechanistic_demonstration_v1.py` before accepting any change to this starting state.
