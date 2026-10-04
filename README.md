# aza3 — source of phenotypic evolvability in a young thistle radiation

## Nature-scale programme

aza3 now preserves two nested roles:

1. the original Chapter-3 role: use own same-individual ancestry, phenotype and cytotype data to discriminate histories left unresolved by EAzami;
2. the higher-order Nature-scale question that emerged from Azami + EAzami:

> **Why could a young thistle radiation repeatedly generate so much reproductive-phenotype diversity in so little evolutionary time?**

The current high-risk working hypothesis is that rapid capitulum diversification was enabled by **reuse, sorting, exchange and regulatory redeployment of pre-existing genomic variation**, rather than by repeated de novo invention alone.

Competing genomic explanations are retained explicitly:

- repeated de novo mutation;
- ancestral standing variation;
- introgressive reuse;
- polyploid/homeolog/regulatory reuse.

This is a hypothesis programme, not a result. EAzami establishes repeated phenotypic reassembly, not its genomic cause.

### Current Nature-scale execution plan

1. `docs/AZA3_NATURE_SCALE_MASTER_PLAN_V1.md` — active Nature plan; Gate 0 is split into **0A dated-tree recovery** and **0B phenotype-rate testing** before the four genomic source models.
2. `data/planning/eazami_to_aza3_hypothesis_bridge_v1.csv` — explicit EAzami-result → aza3-hypothesis map, including falsifiers and claim ceilings.
3. `docs/GATE0A_DATED_TREE_READINESS_AND_RECOVERY_V1.md` + `data/planning/aza3_gate0a_recovery_routes_v1.csv` — current first executable gate; author-source recovery before independent rebuild.
4. `data/planning/aza3_nature_sampling_priorities_v2.csv` — canonical sampling order ranked by hypothesis discrimination rather than taxonomic coverage.
5. The original `data/planning/chapter3_sampling_priorities_v1.csv` remains preserved as the narrower history-resolution plan and is not the current strategic ranking.

### Read the origin of this question first

1. `docs/AZAMI_EAZAMI_TO_AZA3_NATURE_ORIGIN_HANDOFF_V1.md` — narrative research-origin handoff: how the question changed from spatial phenotype ecology to evolutionary depth and finally to the source of rapid phenotypic innovation.
2. `data/contracts/aza3_nature_origin_handoff_v1.json` — machine-readable source commits, turning points, competing models and claim boundaries.
3. `docs/CHAPTER3_SCOPE_AND_HANDOFF_V1.md` — original own-data Chapter-3 scope, retained as the operational ancestry/phenotype/cytotype layer.

The programme-level evidence chain is:

`phenotype space (Azami) → repeated reassembly in a young radiation (EAzami) → genomic source of evolvability (aza3) → representative fitness consequences`

---

# aza3 — own-data discrimination after EAzami Chapter 2

This repository starts Chapter 3 from the uncertainty left by the completed public-data Chapter 2. It is not a second copy of the Chapter 2 paper and does not treat future RAD-seq or field data as retroactive confirmation.

## Locked starting point

EAzami Chapter 2 and its complete meta/simulation disposition are frozen at merge `62fa8c5c913c2b236e710f6bad366e80676aa78f` ([core PR #129](https://github.com/zuizui0223/EAzami/pull/129); [completeness PR #130](https://github.com/zuizui0223/EAzami/pull/130)). Within its admitted public topology ensemble:

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
4. `data/planning/chapter3_bounded_prior_registry_v1.csv` — 14 meta-analysis, programme-routing and simulation boundaries.
5. `data/planning/chapter3_protocol_registry_v1.csv` — experiment readiness without field authorization.
6. `data/intake/chapter3_individual_intake_v1.csv` — currently empty same-individual intake ledger.

## Current authorization state

- own biological data admitted: **0**;
- field execution authorized: **false**;
- tissue collection authorized: **false**;
- sensitive coordinates permitted in this repository: **false**;
- definitive Japan-wide species-tree claim authorized: **false**.

Run `python analysis/validate_chapter3_handoff_v1.py` before accepting any change to this starting state.
