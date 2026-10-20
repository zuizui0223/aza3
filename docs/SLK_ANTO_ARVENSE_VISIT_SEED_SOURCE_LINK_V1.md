# Real Cirsium visitor composition and measured seed outcomes: the empirical SLK × ANTO boundary

**20 October 2026 | original author data reanalysis, exploratory after source inspection | biological empirical owner: aza3 | theory owner: SLK and ANTO separately.**

## Biological question

The SLK × ANTO v1–v6 theory asks whether ecological partner environments alter the benefit of a multifunctional floral architecture and whether a lineage retains the genetic/developmental ability to express it later. This has produced exact and simulated conditional predictions, but so far no **one-individual linkage** of floral functional structure, pollinator/enemy activity, mature achenes and developmental competence.

Here we do not add another synthetic threshold. We ask a small, real ecological question: **in the available Cirsium field data, can recorded visitation order composition be linked to final reproductive seed outcomes at all, and is the direction consistent across places and sampling definitions?**

### Original biological source and provenance

The dataset belongs to the original *Cirsium arvense* 2022 four-site field experiment of Karpate et al., 2025 **Research Square preprint** DOI [10.21203/rs.3.rs-6634172/v1](https://doi.org/10.21203/rs.3.rs-6634172/v1) and associated author-deposited [GitHub source commit](https://github.com/Yogita25/Multitrophic_Interactions/tree/8a151c62ed069f775a16b9bc5fbdce209cbb8691). This is a **secondary computation from genuine experimental field data**, not newly sampled plants or a randomization of visitor species.

We re-download and independently Git-blob-SHA-verify all three original CSVs:

- `seeds.csv` (615 measured head rows; SHA1 `ca99c762162d15f006f0c73b4603ca71fd320a29`)
- `Pollinators.csv` (120 timed observational bouts; SHA1 `1723e96b8fc8ea755832daa786c387c5f3c188a8`)
- `Pollinator_Id.csv` (649 visitor identification lines; SHA1 `875f2046e94b5deca62220e605ed60c4efe232bb`)

The earlier original-author seed-only analysis and 39 matched organza-bag versus control blocks already exist in [aza3's real 615-head study](CIRSIUM_ARVENSE_REAL_HEAD_SEED_REANALYSIS_AND_RESIDUAL_RISK_V1.md). These are **the same source plants and observations**; this result does not create an independent second biological sample or field replication. The original publisher source continues to be identified as a preprint, not a confirmed final peer-reviewed publication.

### Biological correspondence and hard join ceiling

Original `Site_Code × Block_Code × Treatment_Code` labels are shared by the two original source tables. The author's published R code also defines plant-like `Individual` identifiers from block and treatment. Here we match all THREE labels after case standardization and the documented visit treatment-name alias (`Exclude predators` ↔ `Predator Exclusion`). This is an **inferred maternal treatment-plant group**, not a source-authenticated single unique head with matching visit time, genotype, or voucher. The seed table has multiple head rows per such source label.

Of 120 bouts, 92 distinct site×block×treatment groups appear; 84 of those labels also occur in the original seed-source maternal groups. We deliberately exclude ambiguous duplicated bout UUIDs, zero-duration bouts and zero-open-head bouts from the nominal count-rate denominator. This leaves **98 valid visitor bouts on 81 inferred biological groups**; the Hymenoptera share of the **reported number of visited flowers** relative to Hymenoptera plus Diptera can be computed for **68 of the 81 groups**.

No eligible visit or eligible head is converted to a genotype, effective pollination event or independent site replicate. A row with `Pollination=Yes` is an original author's classification, **not measured stigmatic deposition or paternity**. Likewise `Feeding_Guild=Herbivore` labels in `Pollinator_Id.csv` do not imply that every bee is a seed predator. Source keys with repeated UUIDs are rejected for bout assignment, not arbitrarily assigned to the first appearance.

### Original numeric outcome: modest pooled association, opposite local signs

The primary measure is the proportion of reported visited-flower counts attributed to Hymenoptera among the two named orders. The original `Viable_Seeds / Flowers` ratio is pooled over each inferred maternal group. These are both observational ratios and neither identifies an architecture-specific SLK payoff.

| Noncausal descriptive estimand | Eligible group count | Pearson r |
|---|---:|---:|
| Hymenoptera reported visit share vs viable seed fraction | 68 | **+0.21035** |
| Same variables after removing site×treatment group means | 68, across 12 strata | **+0.28494** |
| Visitor-row-count rather than visited-flower-count share | 68 | +0.15564 |
| Visitor-bout equal rather than aggregate count share | 68 | +0.15763 |
| One valid observation bout only | 54 | +0.10418 |
| More than one valid bout | 14 | +0.64913 |

These are **post-hoc exploratory source-label correlations**. In particular, 14 multiply observed groups cannot be treated as randomized repeats of 14 new evolutionary histories; within-site×treatment demeaning removes only additive means, not phenology, flower size, head condition, actual time-of-day, missing visits, source-selection or causal treatment mediation.

The pooled sign hides site-specific heterogeneity:

| Field site | Group count with defined order share | Pearson r (order share vs seed fraction) |
|---|---:|---:|
| campus | 21 | **−0.01404** |
| construction | 17 | **−0.17515** |
| ditch | 16 | **+0.50277** |
| rapseed | 14 | **+0.44601** |

The signs vary from negative/near-zero to positive. **This is not evidence that pollinators selectively harm the plant at one site and help it at another.** The same insect order can contain species with unequal pollination function; site and observation selection can generate correlation-sign differences even if causal head-level functional selection never changes.

A separate descriptive intensity proxy—the sum of reported visited-flower counts for both orders per nominal open-head-hour—has a pooled correlation with viable seed fraction of **+0.11272** across 81 groups. Excluding the highest ~5% source-label values of that rate changes it to **−0.02692**. Its site×treatment centered correlation changes from **+0.26364** to **+0.11573** under the same exclusion. Thus a positive intensity–seed association should not be promoted to a generalized ecological mechanism; rare very high count-rate labels affect its magnitude and pooled sign.

### The two real outcomes are not interchangeable

The source head archive explicitly distinguishes `Viable_Seeds` from directly recorded `Damaged_Seeds`; the original paper's inferred `Eaten_Seeds` is an accounting residual `Flowers−Viable_Seeds−NViable_Seeds`, **not a third independently scored biological predation measure**. Negative residual values occur, as shown in aza3's source-only seed-reanalysis. Do not substitute that residual for direct insect feeding or use it to invent a pollination→florivore mediation coefficient.

The observational visitor/seed data can therefore motivate future **same-head pollen deposition, herbivore entry, direct physical damage and final filled achene** measurements, but they do not themselves establish a trade-off across those routes. Even an actually higher viable seed fraction need not establish demographic recruitment or genetic fixation, as independently illustrated by *Cirsium canescens* field experiments in earlier aza3 audits (Russell et al. 2025, DOI 10.1002/ecs2.70310).

## Explicit SLK × ANTO scientific inference

**What is source-backed:** In real *Cirsium arvense*, multiple insect orders are recorded during the same timed floral survey bouts (48 of 92 attributable visit-bearing bouts in the earlier timing analysis), and the recorded hymenopteran share shows modest, place-heterogeneous, weighting-sensitive association with later measured plant-group seed fractions. These are observations on a real multi-trophic plant system, not generated traits from the theoretical model.

**What SLK still does not have:** an experimentally identified common physical structural coordinate (e.g. head inclination, bract gap, spine, protective exudate), opposing *verified functions* on the **same** structure, an alternative integrated/differentiated candidate's net reproductive benefit, local variant access and rare phenotype performance. Visit-order correlations cannot be assigned to any of those.

**What ANTO still does not have:** functional same-plant joint genomic/developmental module inventories, genuinely repressed versus nonfunctional states, regulated re-expression by biological cue, source-resolved ancestral/introgressed alleles, or transmission across offspring. The real Ryukyu *C. brevicaule–irumtiense* young-leaf RNA fragments analyzed by aza3 are on OTHER taxa and individuals and cannot fill any of these measurements for the *C. arvense* source.

**What stochastic v6 still does not have:** longitudinal whole-environment Markov persistence `P(E_t=E_(t+1))` **over plant generations**. Ten intermittently sampled days across four sites in one summer are not ten generations, known continuous residence times, nor a stationary exogenous Markov chain. This data set is an ecological *exposure anchor* only, not a validated 14-generation dynamic process.

### Scientific decision

This evidence does **not** confirm the joint SLK–ANTO hypothesis and is not a standalone publishable new mechanism at this stage. It does show which crucial empirical bridge is missing: **same-lineage co-measurement of functional-structure performance, the ecological partner's time series, and inherited ability to re-express the identical architecture**.

The least conflated next biological comparison requires full factorial record keeping before analysis: independently surveyed pollinators and seed enemies (including failed entries), repeated head identity and real pollen/seed outcomes under known and authorized anatomical manipulations, environmental history sampled over generations, genotype/functional component assays, and controlled lineage parent-offspring links. Head morphology in the present dataset is missing. In aza3's *C. sieboldii* programme, 40 planned molecular individuals remain PROSPECTIVE, with no assured overlap to EAzami's separately planned functional-manipulation plants. Do not change authorization or claim such a join has already been achieved.

## Machine reproducibility and protection from claim inflation

- Analysis: `analysis/audit_slk_anto_arvense_bout_seed_join_v1.py` — uses the existing source-validated `parse_head_rows` and dated visitor parser instead of duplicating or altering v1 ecological data.
- Dedicated unit tests: `tests/test_slk_anto_arvense_bout_seed_join_v1.py` — nine tests for real source label meaning, deliberately ambiguous bout keys, valid effort, alternative weighting, exact seed mass balance, site/treatment-controlled association and explicit inability to infer genetic capacity.
- Original-source workflow: `.github/workflows/slk-anto-arvense-visit-seed-link-v1.yml` — re-download and Git-blob-SHA verify all three public author CSV files; recompute the ecological aggregates; compare against the source-locked numeric receipt `data/evidence/slk_anto_arvense_bout_seed_link_v1.json`.
- Review status: [aza3 draft PR #62](https://github.com/zuizui0223/aza3/pull/62), based on the original-clock [PR #58](https://github.com/zuizui0223/aza3/pull/58). [Successful original-source verification](https://github.com/zuizui0223/aza3/actions/runs/39415739190).
- Source-data evidence and source code remain in aza3 only. SLK and ANTO continue to own distinct questions; these correlations cannot be counted as **two** independently confirming studies if both repositories cite them.

**Conclusion ceiling:** a modest overall order-composition/viable-seed association exists in the source-matched label groups, but it is not uniformly signed between sites and is not robustly interpretable as a direct biological pollination benefit, structural functional conflict, capacity regain or an evolutionary adaptive landscape.
