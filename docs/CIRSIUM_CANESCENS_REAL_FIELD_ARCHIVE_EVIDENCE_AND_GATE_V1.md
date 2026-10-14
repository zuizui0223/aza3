# A real Cirsium head-level external archive: what can actually be tested next

Date: 2026-10-14. Evidence priority in PR #44. Do not confuse public archive *existence* with a completed CSV reanalysis.

## Archival positive control in the focal genus

Russell, Taylor & Louda (2025), *Ecosphere* DOI [10.1002/ecs2.70310](https://doi.org/10.1002/ecs2.70310), reported insect-exclusion experiments on *Cirsium canescens* in Colorado (five sites, 2017–2019) with individual plant and individual capitulum outcomes. Public CC0 Dryad dataset DOI [10.5061/dryad.kh18932dc](https://doi.org/10.5061/dryad.kh18932dc), version 3 dated 2025-05-06, includes:
- `Plattethistleseedproduction.csv`, 348,096 bytes, SHA256 `d8c584b0c67a80f7e637e95df2262c7cb16788b0bd09d7e6e4ae702e00508b50`, file ID 4056201;
- `Plattethistleseedlingestablishment.csv`, 178,451 bytes, SHA256 `d5d13c621ada942ad343a5bf337fd69f2f92ec8f79496fa095d36e60677f6163`, file ID 4056202;
- `Plattethistlefirstflowering.csv`, 10,685 bytes, SHA256 `f784f1c56d3f555e4ba3bfbd0ef9b6e0cbdae112f75284cead50679df5b53013`, file ID 4056203;
- `README.md`, 34,555 bytes, SHA256 `393c3b9aa2e00bee878ab74ff77fcf0e977fca25c507bb27c37bd84ccbd7b71b`, file ID 4056213.

These metadata and variables were verified from the publisher's archive landing page and the public version-3 Dryad REST metadata. **The raw files were NOT downloaded or parsed this turn**: web-facing `/downloads/file_stream/4056201` returned 403; unauthenticated `/api/v2/files/4056201/download` returned 401. Dryad's documented file-download API requires bearer authentication even for public datasets; unauthenticated metadata access succeeds. No API credential is held or requested by this repository. An endpoint returning HTTP 401/403 is not evidence of absent organisms or absent data. When the CC0 CSV is made available through a legitimate source, check the exact SHA256 before analyzing it.

## Why this archive changes the research strategy

The original 2025 study already established that reducing floral/seed herbivory increased *Cirsium* viable seed output, while recruitment into seedlings did **not** increase proportionally; elevation influenced flowering time but not the herbivore effect. Thus even **real intervention → higher maternal seed output** does not mean a proportional **plant population recruitment response**. A head defence effect, a maternal fitness effect and recruitment limitation are distinct ecological endpoints. Avoid making any claim that seed-level defence necessarily explains population growth or adaptive evolutionary origin.

Crucially, the archived **per-head** data directly record stages of an antagonistic head interaction that the local Azami photographs lack:

| Original source field | Biological meaning / inferential limit |
|---|---|
| `Project`, `Year`, `Site`, `Plant`, `HdN` | Actual grain is **individual head** nested in reproductive plant; the plant, not head, was randomized to insecticide within site |
| `FnlTrt` | Plant-randomized insecticide `I` vs pooled water/untreated control `C`; stored on FIRST row of a plant, not necessarily repeated on every head |
| `HdN` | 100 = main stem apical capitulum; 200/300 = lateral branch terminal heads; 210/220 = axillary heads on a lateral branch (position, NOT gravitational inclination) |
| `Stg` | Stages 1–6 assessed at **last mid-August census**, hence post-treatment/developmental state, not automatically a baseline covariate |
| `HdDam` | Externally visible insect damage category: 0 none, 1 ambiguous, 2–4 increasing phyllary damage, 5 stem miner, 6 chewed hole into the head. **It is NOT a continuous ordinal severity 0–6**: categories 5/6 are different attack routes, not degrees 5/6 |
| `RcEg` | Number of ovipositions by exotic *Rhinocyllus conicus* on the focal head, not all herbivore guilds, effective entry or proven seed loss |
| `FinPRGD` | Estimated **proportion** of viable seeds in head, not necessarily independent germination assay |
| `EstSdProd` vs `EstSdProdTang` | Raw estimate vs estimate corrected for unintended `Tanglefoot` damage in 2018. Published README says corrections affect **68 heads**. Do not silently combine both columns or treat Tanglefoot-killed heads as naturally insect-damaged |
| `TPltSDs` | Plant-level estimated viable seed total; may only appear on apical/first row. Do **not** sum or treat it as a head-level response |

Original README itself lists some bud-size unit phrasing oddly (`>10 cm`/ `<10 cm` for capitula), so do not assume the specified 10-cm value is anatomically correct without inspecting methodology. No actual spine length, rigidity, head inclination, authenticated insect entry attempts, or confirmed pollinator service are in this archive.

## Predeclared analysis, not an already observed result

The next archive analysis, when legitimate source bytes are available, is **one design-correct positive control**:
- Preserve the unique head key `Project × Year × Site × Plant × HdN`.
- Carry the insecticide treatment down within that **same plant only**; never use a simple entire-CSV forward fill, which would carry treatment across plants.
- Report each of the three projects (when applicable), each study year and site explicitly; Project 1 and Project 2 use different head-viability estimation protocols. Treatment assigned to **plants**; do not treat all heads as independent randomized replicates.
- Independent process readouts: fraction with unambiguous external damage (`HdDam` 2–6), actual `R. conicus` oviposition counts (`RcEg`) among eligible assessed heads, and **seed output** `EstSdProdTang` separately by intervention and head rank. Neither damage nor eggs may be used as baseline confounders in the total insecticide-effect analysis because they are outcomes/mediators.
- Main response: **plant-summed corrected viable seeds** using every enrolled head, including complete zeros but **not** converting missing values to zero. Use plant as the replicate with site/year/project blocks.
- Stage/rank descriptive readouts must not be elevated to direct causal spine effects. `Stg` is a final-stage observation, not a pretreatment covariate, and rank/HdN proxy head position but not orientation relative to gravity.
- Evaluate whether the published Tanglefoot correction changes plant/seed totals. If corrected and uncorrected columns differ, report discrepancy and the 2018 intervention artifact; do not attribute it to naturally occurring defence.

**Target interpretive boundary:** This is direct evidence about Cirsium head enemy biology and its reproductive consequences, but it does **not** test the proposed orientation × authentic-spine mechanical-access hypothesis. If analysis reproduces known herbivore effects and measures rank/stage differences, that is an independently audited positive control, not new head trait selection.

## Source access status and honest result

- Public archive: **VERIFIED** (DOI, version, file lists, field definitions, SHA256 metadata).
- Original individual-head CSV content: **NOT FETCHED / AUTH REQUIRED**.
- Recomputed Cirsium head-level contrasts: **NOT AVAILABLE**, not zero.
- New `aza3` field video/contact/viable-seed ledgers: **zero real rows**, so the focal orientation × spine fitness analysis remains not identified.
- Nature-level novelty: **not claimed** for broad herbivory seed benefits or seedling limitation (already published); remains specific **external head morphology × guild-specific physical contact/entry versus final seed fate** across independent lineages.

Original published source: Russell, F. L., Taylor, M. R. & Louda, S. M. (2025) *Ecosphere* 16 e70310, DOI 10.1002/ecs2.70310; Dryad DOI 10.5061/dryad.kh18932dc. No external field sampling permission granted by this audit.
