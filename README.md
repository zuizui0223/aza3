# aza3 — own-data discrimination after EAzami Chapter 2

This repository starts Chapter 3 from the uncertainty left by the completed public-data Chapter 2. It is not a second copy of the Chapter 2 paper and does not treat future RAD-seq or field data as retroactive confirmation.

## Locked starting point

EAzami Chapter 2 is frozen at merge `4fc03f128a7ec05ce9e16e1daedef23b61104b89` ([PR #129](https://github.com/zuizui0223/EAzami/pull/129)). Within its admitted public topology ensemble:

- orientation requires **at least four** state changes;
- phyllary posture requires **at least three**;
- stickiness requires **at least five**;
- minimum counts are better resolved than individual event placements;
- species-tip coding hides state multiplicity in 4/4 audited polymorphic systems, while direct morph-genotype linkage exists in only 1/4.

These are topology-conditioned lower bounds, not counts of independent origins, convergence events, rates or adaptations.

## Chapter 3 question

> Which histories and causal paths that remain compatible with Chapter 2 can be discriminated by own Japan-wide ancestry data linked to phenotype, cytotype and separately authorized functional experiments?

The first product is an all-Japan same-library RAD-seq topology/network sensitivity, not an unconditional species tree. Every genomic record must link to an immutable individual, a voucher or diagnostic image, phenotype states, cytotype status and deidentified authorization records.

## Start here

1. `docs/CHAPTER3_SCOPE_AND_HANDOFF_V1.md` — question, work packages and claim boundary.
2. `data/contracts/chapter3_eazami_handoff_contract_v1.json` — fail-closed source and authorization contract.
3. `data/planning/chapter3_sampling_priorities_v1.csv` — five ranked history discriminators.
4. `data/planning/chapter3_bounded_prior_registry_v1.csv` — meta-analysis priors and simulation boundaries.
5. `data/planning/chapter3_protocol_registry_v1.csv` — experiment readiness without field authorization.
6. `data/intake/chapter3_individual_intake_v1.csv` — currently empty same-individual intake ledger.

## Current authorization state

- own biological data admitted: **0**;
- field execution authorized: **false**;
- tissue collection authorized: **false**;
- sensitive coordinates permitted in this repository: **false**;
- definitive Japan-wide species-tree claim authorized: **false**.

Run `python analysis/validate_chapter3_handoff_v1.py` before accepting any change to this starting state.
