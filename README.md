# aza3 — source of phenotypic evolvability in a young thistle radiation

## Q2 new biological stage finding: capitulum good seed loss is primarily a per-Flowers *rate* contrast (23 October 2026)

**Original 2022 Cirsium arvense 615-head source, four original sites and the same 36 complete four-arm blocks used in PR #82.** The author-defined good-seed production / sampled head is the EXACT product of (source `Flowers` / sampled head) × (author `Viable_Seeds` / `Flowers`). The symmetric within-block two-factor decomposition of organza-bagged minus control **−41.59 viable seeds/head** gives **−0.33 from differing reported Flowers/head** and **−41.25 from differing author-defined good-seed fraction** (99.2% of the descriptive observed mean difference). Four-site-stratified original block bootstrap 95% for the two components: Flowers count [−2.16,+1.36], good-seed fraction [−46.56,−35.77]; the rate component is negative in all four sites.

**Second exact accounting check:** bagged-minus-control good-seed fraction −0.523 is offset by +0.389 author `NViable_Seeds` fraction, −0.036 visually `Damaged_Seeds` fraction, and +0.170 *signed accounting residual* `Missing_Seeds` fraction. This is NOT evidence of more real larval feeding, identified pollen deprivation, anatomical defence, or independent causal path shares: source Flowers and seed categories are observed/constructed after compound physical bagging, and no genuine floret/head morphology or effective pollination was assayed. Predator exclusion shows a *different descriptive composition* (−3.82 viable seeds/head ≈ −3.77 from count component and −0.05 from rate component), but treatment arm differences may reflect plant/sample composition, not changing head structure.

[Biology-first source report](docs/CIRSIUM_ARVENSE_ORIGINAL_HEAD_FLORET_VS_SEED_CONVERSION_20261023.md), [original Git-blob-pinned Python analysis](analysis/audit_cirsium_arvense_head_floret_conversion_v1.py), [independent JavaScript 36-block numerical receipt](data/evidence/cirsium_arvense_head_floret_conversion_independent_js_v1.json), [replay CI](.github/workflows/audit-cirsium-arvense-head-stage-decomposition-v1.yml). **The C. arvense functional source and 114-individual head morphology source are different taxa and unlinked individuals: no fabricated same-head trait×partner fitness demonstration.**

## Q2 new original-data joint response: fecundity and antagonists can *both* decline (23 October 2026)

**Original 2022 field data — 615 Cirsium arvense head records / 36 complete site-matched four-arm blocks (2025 preprint).** Bagged pollinator-exclusion heads relative to controls show **viable seed fraction −0.523** (descriptive four-site bootstrap 95% [−0.587,−0.456]), **emerged insects/head −1.450** ([−2.539,−0.607]); **25 of 36 blocks show both seed output and insect emergence declining**. In the same 36 blocks, predator exclusion versus control shows **viable seed fraction −0.006** ([−0.064,+0.052]) despite decreased emerged insects/head −0.928 ([−1.901,−0.150]). Within identical matched blocks, bagged minus predator-excluded viable seed fraction is **−0.517** ([−0.579,−0.454]), negative **36/36**. Thus *fewer insects* is not a sufficient proxy for *better plant reproduction*, and multitrophic responses need distinct measured endpoints.

**Crucial mechanism boundary:** The organza bag directly changes *both pollinator and other insect physical access*, plus potentially microclimate; predator exclusion uses adhesive barriers and manual removals; plastic spider is only a cue. Counts of emerged insects do not establish true larval seed-feeder identity. This source supports a joint-access ecological phenomenon, **not a demonstration of Cirsium spine/head-angle adaptive LEGO, pollinator-driven florivore reproduction, or causal partner-specific selection**.

See [original source and biological analysis](docs/CIRSIUM_ARVENSE_2022_JOINT_SEED_INSECT_ACCESS_RESPONSE_20261023.md), [reproducible code](analysis/audit_cirsium_arvense_joint_access_response_v1.py), [independently calculated original-JS 36-block receipt](data/evidence/cirsium_arvense_four_arm_joint_response_independent_js_v1.json) and [CI workflow](.github/workflows/audit-cirsium-arvense-joint-access-v1.yml). No aza3 own field sampling was added.

## New Q2 source interpretation gate: 14/106 selected mothers and year-changed protocols (23 October 2026)

**The biological denominator changes the conclusion.** Although 20 genuine x00/x10 branch-code pairs show unequal seed production, they come from only **14 of 106 original USDA *C. undulatum* mothers**: 7/39 in 2007 and 7/67 in 2008. Those selectable mothers already have higher observed total good seed (102.57 vs 54.31 in 2007; 147.71 vs 81.75 in 2008) and larger recorded rosettes, so the paired subset cannot be interpreted as the randomized population. Taking only the lowest-coded pair per mother changes the 2008 terminal-only vs intermediate-only comparison to **4 vs 0 across 7 mothers, exploratory exact two-sided p=.125**; a highest-code rule yields 5 vs 0, p=.0625. Neither establishes an adaptive position effect.

**Newly verified original-method confound:** West's [2012 University of Nebraska dissertation](https://digitalcommons.unl.edu/bioscidiss/41/), chapter 2, records **acephate every 7–10 days in 2007 versus bifenthrin every 14 days in 2008**, additional mechanically damaged apex treatments in 2008, and flowerhead observations every 2 versus 4 weeks. Year differences cannot be attributed uniquely to environment, insect ecology or architecture. The same original source identifies subordinate head positions as below a particular branch terminal, strengthening the position-code meaning without proving an individual photographed branch map. These are human-checked 2012 documentary facts, not automatically derived from the 2007–2008 source CSV.

[New original-data selection and methodological-confounding audit](docs/CIRSIUM_UNDULATUM_2012_SOURCE_YEAR_PROTOCOL_SELECTION_20261023.md) · [replay code](analysis/audit_cirsium_undulatum_year_protocol_selection_gate_v3.py) · [independent 106-mother/arm numerical receipt](data/evidence/cirsium_undulatum_year_protocol_selection_independent_js_v3.json). **Biological endpoint remains observational good seed, not demonstrated visitor-mediated protection, heritable LEGO or regained ancestral colour.**

## New Q2 result: exact branch-code pairs split reproductive success from seed yield (23 October 2026)

The original USDA 2007–2008 *Cirsium undulatum* data contain **20 recorded branch-code terminal/intermediate pairs across only 14 real mother plants** (11 pairs in 2007, 9 in 2008). In 2008, **5/9** pairs had original author-scored good seed only at the branch terminal, **3/9** at both positions, **1/9** at neither; individual good-seed counts were lower at intermediates in **7/9**, higher in **1/9** and equal in **1/9**. Thus the earlier 7/7 negative *mother-mean* effect is **not** a universal within-branch effect. The positive 2008 WR56 branch-code pair had 68 vs 85 good seeds at equal recorded injury (0.63). Injury contrasts are not uniformly worse for intermediate heads. Code-block identity is not independently verified physical branch connectivity, flowering age or manipulation. This is exploratory and does not demonstrate adaptive architecture, partner causation or genetic modularity.

See the [source-faithful stage/hurdle ecological report](docs/CIRSIUM_UNDULATUM_NESTED_BRANCH_SEED_HURDLE_20261023.md), [original-data Python audit](analysis/audit_cirsium_undulatum_nested_branch_hurdle_v2.py), [independently reconstructed 20-pair JS receipt](data/evidence/cirsium_undulatum_nested_branch_hurdle_independent_js_v2.json), and [source-matched Actions reproducibility](.github/workflows/audit-cirsium-undulatum-nested-branch-hurdle-v2.yml). The previous 2008 six-exploratory-test **BH q=.09375** still bounds inferential claims.

## New source-backed 2026 Ryukyu thistle RNA observations (19 October 2026)

**Six original public young-leaf RNA runs analysed (3 *C. brevicaule*, 3 *C. irumtiense*).** A capped 8 MiB-first-R1-per-source reanalysis of **747,234 real reads** against the **entire 104,890-transcript *C. japonicum* 2020 publisher reference** (plus composition-shuffled decoys) produced **215 strict selected-candidate-mapped fragments** and zero strict shuffled-decoy matches. The previously identified **DFR-like TR39766** source transcript had **66** compatible fragments, and **ANS-annotated TR41106** had **106**; both were detected in **all six** taxon-labelled source samples. This observation survives competition against the previously omitted ~104k author transcripts, not merely a 26-candidate index.

Read this as **shared young-leaf transcript-fragment homology between taxon-labelled white and coloured Cirsium lineages**, NOT intact catalytic DFR/ANS genes, known flower pigment concentration, genetic ancestral retention, heritable colour reactivation or adaptive genetic LEGO modules. The reference transcript IDs belong to *C. japonicum* 2020, not sequence-confirmed Ryukyu allele IDs. Original per-plant flower colour is not assayed. Raw 8/16 MiB prefixes are nonrandom, and read counts are not expression contrasts.

- [Actual six-run whole-author-transcriptome source result and experimental meaning](docs/RYUKYU_FULL_AUTHOR_TRANSCRIPTOME_RESULTS_V7.md)
- [Machine-frozen public evidence with accession joins](data/anto_source_handoff/public_sra/ryukyu_full_104890_reference_competition_8MiB_v7.json)
- [Independent source/biology claim validator](analysis/validate_ryukyu_104890_full_reference_v7.py)
- [Underlying capped read/whole-reference competitive analysis](analysis/anto_empirical/ryukyu_full_unigene_competition_v7.py)

The next biological gate is to verify real DFR/ANR/ANS/FLS locus orthology and catalytic activity in focal plants, and compare stage-matched floret/involucre/nonfloral metabolism and developmental regulatory function. Evolutionary-regain history remains unpolarized. **anto remains general theoretical work; aza3 exclusively owns these biological observations and follow-up experiments.**



## anto × aza3: empirical ownership and LEGO hypothesis (18 October 2026)

**aza3 is now the exclusive active home for Cirsium empirical study design** underlying the cross-angiosperm [anto theoretical programme](https://github.com/zuizui0223/anto): Ryukyu *brevicaule*–*irumtiense* source data, var. *takaoense* molecular contrasts, organ-scoped colour and pigment tests, population ancestry, cytotype, all same-individual capitulum modules, functional manipulations and fitness.

- **anto upstream:** What modules can plant lineages *carry* across dispersal/colonization, and how long can lost floral phenotypes remain evolutionarily reachable?
- **aza3 downstream LEGO:** Which retained/shared/exchanged modules can be jointly assembled into newly expressed capitulum configurations, conditional on ancestry, developmental compatibility and environment?
- The two claims are not equivalent. Repeat appearance of a coloured flower can be immigration/introgression, and recurrent head shapes need not use the same genetic building blocks.

**New v2 genomic LEGO test (19 October 2026):** [Higher-order module measurement and power plan](docs/ANTO_HIGHER_ORDER_LEGO_MEASUREMENT_PLAN_V2.md) operationalizes anto's exact result that identical marginal and pairwise module frequencies **do not identify 3-way genetic resource availability**. The first within-species focal system is the **prospective 40-individual C. sieboldii target**, with same-individual colour, orientation and involucre phenotypes measured separately from **functional genetic resources**. Phenotypic correlation is not yet genetic LEGO. [Machine contract](data/contracts/aza3_capitulum_higher_order_lego_v2.json) and [outcome-blind feasibility checks](analysis/validate_aza3_capitulum_higher_order_lego_v2.py) require population/ancestry controls, held-out replication, mechanistic compatibility, and separate viable-achene fitness evidence. All relevant individual datasets will live here, not in anto.

**Authoritative biological execution and decision gates:** [ANTO_ACCESSIBILITY_TO_AZA3_HEAD_LEGO_EMPIRICAL_MASTER_V1.md](docs/ANTO_ACCESSIBILITY_TO_AZA3_HEAD_LEGO_EMPIRICAL_MASTER_V1.md), [machine-readable contract](data/contracts/aza3_anto_empirical_lego_bridge_v1.json), and [fail-closed validator](analysis/validate_anto_lego_empirical_bridge_v1.py).

**New actual Ryukyu molecular observation (19 October 2026):** For the FIRST time, the [six source-verified 2026 Cirsium young-leaf RNA R1 prefixes](docs/RYUKYU_REAL_LEAF_RNA_PREFIX_PILOT_V6.md) have been analyzed, not just registered. Among six original 2026 *C. brevicaule* / *C. irumtiense* source libraries, a larger [104,890-source-transcript competition](docs/RYUKYU_FULL_2020_UNIGENE_COMPETITIVE_MAPPING_V7.md) still finds 215 stringently mapped R1 fragments to 26 2020 *C. japonicum* BLAST pigment-pathway gene-family candidates, including DFR-like TR39766 (66) and ANS-BLAST-like TR41106 (106), detected in all six; 0 shuffled decoy matches. [Actual source-frozen data and accession-by-accession counts](data/anto_source_handoff/public_sra/ryukyu_six_full_104890_unigene_competitive_mapping_v7.json). **This is young-leaf cross-study sequence homology after full published 2020 transcriptome competition, NOT a verified gene's function, whole locus, flower anthocyanin, white-ancestral re-expression or colour regain.** No comparative expression estimation was performed from these biased first-in-file R1 prefixes. The [offline verifier](analysis/validate_ryukyu_full_104890_unigene_mapping_v7.py) requires exact biological source accession, read-count conservation and stringent claim limits.

**v5 source gene-family trees (19 October 2026):** [Real MAFFT/NJ uncertainty audit](docs/CIRSIUM_FLAVONOID_GENE_FAMILY_UNCERTAINTY_V5.md) retains seven tested DFR/ANR or ANS/FLS candidate placements (only five distinct *2020 Cirsium japonicum* transcripts). **All seven identity-difference 95% intervals include zero; no candidate is confirmed as a functional anthocyanin enzyme or Ryukyu re-expression locus.** A second approximate-ML FastTree comparison is implemented with independent regression tests and must be interpreted as methodology sensitivity, not biochemical assay evidence.

**New biochemical-family specificity screen:** The [seven-reference DFR/ANS-versus-ANR/FLS panel](docs/ANTO_CIRSIUM_PARALOG_SPECIFICITY_AND_HEAD_CHEMISTRY_V4.md) directly tests whether the previous five *C. japonicum* transcript candidates preferentially match pigment-pathway enzymes rather than related reductase/oxygenase paralogs. This is still source-protein homology, **not** a documented Ryukyu floral reactivation. Same-individual floret and phyllary metabolite measurements in aza3 will be needed to distinguish pigment-pathway suppression from alternative pathway flux.

**The empirical analysis itself also runs in aza3 (19 October 2026):** [single Cirsium public-source pipeline](analysis/run_anto_cirsium_empirical_sources_v3.py) and [source-homology/metadata/organ modules](analysis/anto_empirical/) now test the 25 public RNA accession records, source 2020 Cirsium transcripts and GenBank reference CDSs here. Earlier anto code is archived provenance, not a second active Cirsium pipeline. The [full code migration and reproducibility note](docs/ANTO_CIRSIUM_EMPIRICAL_CODE_MIGRATION_V3.md) preserves the rule that sequence similarity or a leaf-only transcript is NOT a flower-pigment re-expression assay.

**Real publicly accessible external evidence copied under aza3:** [Cirsium and Asteraceae source handoff](data/anto_source_handoff/) includes 25 source-registered young-leaf runs (6 focal Ryukyu and 6 var. *takaoense*), 2020 within-*Cirsium* multi-organ phenylpropanoid comparisons, and 26 actual source transcript sequences annotated as DFR/ANS-like by published BLAST. None of these represents a confirmed genotype, regulatory regain or causal colour allele for the 2026 Ryukyu sample pair. Original aza3 own-sample authorization, prospective 298 Phase-A targets and dated-tree Gate 0 remain unchanged. No new biological field collection is authorized.



## Direct measured capitulum component integration (20 October 2026)

**Original-source same-plant head morphometrics (new aza3 reanalysis):** In 114 published *Cirsium* individual plants across three parental taxa and two sex strata, corolla length and involucre length show within-taxon/sex log residual **r=+0.670** (plant-bootstrap 95% CI +0.597 to +0.741), while corolla length and relative phyllary-tip spread **r=−0.053** (−0.185 to +0.079). Their correlation difference **Δr=+0.723** (+0.580 to +0.871); after descriptive width/size adjustment **+0.450** (+0.285 to +0.624), positive when leaving out any one species. This is **phenotypic integration of two reproductive dimensions contrasted with a relatively uncoupled external width proxy**, NOT calibrated spine length, per-capitulum verified physical matching, heritable genetic LEGO, flower visitor access, selection, or restored ancestral pigment function. [Real source-head result and inference ceiling](docs/CIRSIUM_HEAD_CORE_VS_BRACT_TIP_SPREAD_EMPIRICAL_20261020.md), [numerical receipt](data/evidence/cirsium_head_core_external_morphology_result_v1.json), and [replay workflow](.github/workflows/audit-cirsium-head-core-external-v1.yml). No own field/head manipulations have yet been authorized.

## Three linked Cirsium evolutionary questions (19 October 2026)

**New actual Q2 visitor guild resolution (20 October 2026):** [Same-record observer × insect-order Cirsium visitor identity decomposition](docs/CIRSIUM_Q2_WITHIN_ORDER_VISITOR_IDENTITY_V4.md) holds the previous five region/month units and 191 publicly recorded flower visits constant. The observed within-host visitor-name identity excess above the observer-only reference is +0.02466; conditioning on both original observer and insect order still leaves +0.01398. All reference comparisons are observational/photo sampling diagnostics, not independently confirmed insect preference, pollination effectiveness or adaptive head morphology. Same-observer/day support remains only 15 exchangeable records. [Source-pinned independent numeric receipt](data/evidence/cirsium_within_order_visitor_identity_v4.json) and [replay workflow](.github/workflows/cirsium-within-order-visitor-identity-v4.yml) keep the source-checksum and original ecological inference limits.

The new [three-question empirical integration](docs/THREE_QUESTIONS_CIRSIUM_MEMORY_PAYOFF_REACHABILITY_MASTER_20261019.md) connects `turnover + chun + EAzami` (lineage portability of trait memory), `bita + sch + balance + slk` (partner-dependent functional payoff and conflict), and `anto + egwe + izu-core + aza3` (re-expression accessibility after visible phenotype loss). It shares real Cirsium capitulum measurements across three separately falsifiable endpoints. The [first evidence-constrained JPN36/Sieboldii execution checkpoint](docs/THREE_QUESTIONS_JPN36_MINIMUM_EMPIRICAL_MILESTONE_V1.md) cross-checks the EAzami 30-plant observation-only reference, 24-head/12-pair pilot, the 2% conservation quota (at least 1,200 counted flowering plants would be necessary for a full 24-head cohort), and the **separate, not-yet-collected 40-plant** aza3 genomic plan. No field intervention is authorized. **This is a prospective programme, not proof of historical colour regain, partner-dependent adaptation or generality across angiosperm families.** All actual Cirsium empirical execution remains owned by aza3; frozen external repository results are not changed.


**New real-data Q2 check (20 October 2026):** The [four-site/39-matched-block Cirsium arvense sensitivity reanalysis](docs/CIRSIUM_ARVENSE_2022_FOUR_SITE_BLOCK_ROBUSTNESS_V2.md) uses 615 authentic author seed-head records, independently computes source-backed site and leave-out contrasts, and distinguishes observed viable seeds, visually damaged seeds, signed bookkeeping residuals and emerged insects. Bagged versus open heads have lower scored viable fraction in all 39 matched blocks, while visible damage changes in mixed directions across individual blocks. The sign of author-defined estimated Eaten differs from direct damage because Eaten includes a signed accounting residual; this is **not** evidence that a true capitulum barrier reverses its fitness benefit with insect partner. [Python code](analysis/audit_cirsium_arvense_2022_site_block_robustness_v2.py), [independent JS numeric receipt](data/evidence/cirsium_arvense_2022_site_block_robustness_v2.json) and [source-pinned CI](.github/workflows/audit-cirsium-arvense-2022-site-block-robustness-v2.yml) are in this branch. Original collection/functional head trait manipulations remain unperformed.

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
4. `docs/GATE0B_PHENOTYPE_RATE_READINESS_V1.md` + `data/contracts/aza3_gate0b_rate_estimand_v1.json` — frozen phenotype numerator and total-rate ratio; currently blocked by Japanese trait coverage.
5. `data/planning/aza3_nature_sampling_priorities_v2.csv` — canonical sampling order ranked by hypothesis discrimination rather than taxonomic coverage.
6. `docs/AZA3_NATURE_GENOMIC_SAMPLING_PANEL_V1.md` + `data/planning/aza3_nature_genomic_sampling_panel_v1.csv` — active 298-individual Phase-A genomic hypothesis-discrimination panel; *C. sieboldii* is the primary same-species combinatorial-reuse system.
7. `docs/ADAPTIVE_VARIANT_IDENTIFICATION_PROGRAM_V1.md` + `data/planning/adaptive_variant_evidence_ladder_v1.csv` — colour-led variant-to-fitness-to-reuse ladder; *takaoense* is the Phase-B molecular anchor.
8. `docs/TAKAOENSE_PHASEB_PREFIELD_CANDIDATE_REDUCTION_V1.md` + `data/contracts/takaoense_phaseb_prefield_contract_v1.json` — one-pass public 6-sample coding/haplotype candidate reduction before new W/BP population sampling.
9. The original `data/planning/chapter3_sampling_priorities_v1.csv` remains preserved as the narrower history-resolution plan and is not the current strategic ranking.

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
