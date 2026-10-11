# 頭花構造は「虫の訪問」ではなく「繁殖成功」まで選別するか：主解析を一つに絞る

2026-10-11 | `aza3` PR #44, go/no-go synthesis. **This is a prospective decision contract, not a new biological finding or field authorization.** The frozen Azami 1,734-image and EAzami analysis branches remain unchanged.

## Primary biological claim we could actually test

> In naturally comparable Cirsium capitula exposed to the same locally present insect species and head stage, does authentic gravity-referenced head orientation × measured involucral armature change the transition from **arrival** to confirmed seed-feeder oviposition without proportionally impeding independent effective pollen deposition, and does the difference predict filled/viable achene output?

This is a specific form of a functional interaction filter, not a generic pollination–defence trade-off or a theory of module evolution. It must survive direct alternatives: post-pollination floral scent, clock-time/diel activity, head rank/size, conspecific and alternative Cirsium display, microclimate, host race, camera coverage, and the age/sex of each head.

### Three distinct denominators; do not collapse them

1. **Head × verified observation effort (primary sampling frame)**: enumerate head and plant IDs before any insect is seen, including intervals with **verified** zero approaches. Observe head morphology, true gravity orientation, botanical spine dimensions and final seed fate independently from whether any larva or parasitoid emerged. All preselected heads contribute to head-level summaries; unscorable coverage and missing seed fate remain **missing, never zero**.
2. **Independently identified insect approach episode (behavioural process)**: count every approach from the source event ledger, including individuals that hover and never have an ordered contact attempt. Preserve unresolved guild identity. A candidate seed-feeder identification must be independent of whether an egg is eventually laid.
3. **Continuously filmed attempted contact (mechanical process)**: every eligible parent attempted contact requires a botanical annotation, including successful and unsuccessful ones, floret-disc access, no contact and unresolved surfaces. `verified physical blockage` is a **per-attempt** state; the later entrance of the same insect is a separate `bypass` event, not a second animal.

The internally audited `validate_capitulum_true_armature_contact_v1.py` now reports `n_independently_identified_candidate_approach_episodes` separately from its attempt counts, `n_candidate_episodes_without_recorded_attempt_sequence`, `n_episodes_with_verified_block_then_later_access`, and verified `oviposition_confirmed`/independent `pollen_deposition_assay` categories. A real positive egg requires species/role-specific oviposition evidence; pollen delivery needs real stigmatic contact and an independent assay. Neither is inferred from temporary or final entry. These are counts of **recorded process evidence**, not natural selection coefficients.

### Prespecified decision ladder

| Evidence achieved | Maximum defensible claim | Not justified |
|---|---|---|
| Validated 2-D photos alone | Nonuniform visible phenotype integration among sampled taxon labels | Real spine measurement, gravity angle, mechanistic access or coadaptation |
| Pre-entry identified insect + complete head video | Rates/route distributions of observed approaches, failed and alternative-portal attempts | Repulsion or causal mechanical trait effect from success-selected clips |
| Actual spine/phyllary contact with caliper-supported geometry | Head/guild/stage-specific *descriptive* access-filter mechanism; bypass versus obstruction | Adaptive module, fitness advantage, genuine exclusion merely from first blockage |
| Independent egg or pollen assay | Confirmed antagonistic/reproductive service linked to the same head exposure | Seed saving or female fitness if mature achenes are absent |
| All preselected heads with comparable viable/filled achenes and permitted trait variation | Association of head structural modules with maternal seed fitness under clearly stated sampling and ecological rivals | Evolutionary adaptive origin without randomization or independent heritable replication |
| Approved randomized manipulation or strong quasi-experimental identification, replicated fitness and ancestry | Candidate causal structure→interaction→fitness and eventual regime-dependent selection | Universality or recurrent evolution unless independently supported |

**Primary go/no-go:** Before testing an orientation × true armature interaction, require (a) within-population true armature/angle variation and measured morphology, (b) independent candidate guild exposure with valid continuous coverage, (c) successful and failed attempts plus final entry, and (d) filled or germinable achene fate for the same *prespecified head cohort* (or explicitly matched authorized intact heads). If (a)–(c) are available but (d) is not, publish **interaction ecology**, not plant adaptation. If the guild never touches the armature or all successful entrances bypass it, report the obstacle as unused/compensated in that observed ecological context rather than inventing a trade-off.

### Critical causal cautions

- Conditioning on heads visited by both guilds, insects that touched a spine, heads where larvae established, or those that produced an emerged adult is conditioning on post-morphology/post-exposure events. It may create a collider or change the estimand. Keep conditional process results separate from **all registered head** fitness comparisons.
- The head is a compound floral organ with repeated attempts by possibly the same insect. Ten video frames or ten contacts on one head are **not** ten biological replicates. Other insect species and distinct *Cirsium* populations are not phylogenetically independent evolutionary transitions by declaration.
- Parent `oviposition_confirmed=0` is at most not-confirmed in a reviewed encounter, not necessarily a demonstrated true lack of eggs in a matured head. `NA` is unevaluated. Likewise `pollen_deposition_assay=negative` must come from an actual measured assay.
- Earlier source-backed examples: *C. heterophyllum* parasitoid-free structural refuge (Romstöck-Völkl 1990, DOI 10.1111/j.1365-2311.1990.tb00814.x); *C. pitcheri* infestation and seed rank (Gijsman et al. 2020 DOI 10.1016/j.gecco.2020.e00945); and koinobiont stimulation of host feeding in **experimental Saussurea** with observational *C. setosum* component (Xi et al. 2015 DOI 10.1111/1365-2656.12361). These are important prior art, **not** direct evidence of the targeted orientation × spine effect.

## Status and near-term decision

- **Source-backed facts:** within-head insect/parasitoid exploitation is real and can affect seed production, but not necessarily in a uniformly beneficial sign.
- **Your numerical results:** 1,734 image observations over 42 Cirsium taxon labels show nonuniform visible phenotype covariance; 0/12 predeclared climatic/photo moderation comparisons pass BH. These image properties are not true botanical thorn anatomy.
- **Direct tests currently absent:** the current `aza3` real true-spine contact, independently observed insect approach and mature viable-achene intake files contain **zero rows**. The new tests that show bypass/egg/pollen evidence are **synthetic invariants only**.
- **Actionable next milestone:** one independently identifiable insect approaches one continuously filmed, botanically measured head; preserve source footage, time/stage, all contacts (including no touch), final access, confirmed egg/pollen observation and the intact head's authorized reproductive fate. If none is available, the scientific result remains **NOT IDENTIFIABLE**; no new mathematical counterfactual should be promoted as a discovery.

Until these observations exist, new derivations of a selection sign are subordinate to validating real field evidence. Existing JPN36 24-head feasibility cap and site/collection permissions remain unchanged.
