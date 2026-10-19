# Q2の重要な既存実験：送粉者を排除すると種子食者の羽化数も減る

2026-10-20. **Source review and independent reanalysis readiness**, not a newly established flower-head armature effect.

## 原著に基づく発見（同一Cirsium種の野外操作）

Karpate, Y., Blažek, P., Kollross, J., Rech, A. R., Sam, K. & Jorge, L. R. (2025), *Pollinators, not predators, affect florivore survival in Cirsium arvense (Asteraceae) reveals an experiment based on Structural Causal Modelling*. Research Square **preprint**, DOI [10.21203/rs.3.rs-6634172/v1](https://doi.org/10.21203/rs.3.rs-6634172/v1), publicly posted June 30, 2025. Peer-reviewed journal publication was **not independently established** in this audit.

The work randomized/manipulated access/exposure in naturally occurring female plants of **functionally dioecious** *Cirsium arvense* in **four sites around České Budějovice, Czechia**, beginning June **2022**. Its four treatments were:
1. unmanipulated control;
2. predator exclusion;
3. pollinator exclusion;
4. predator-model treatment.

A generative structural causal simulation was used to assess the proposed experimental analysis and design sensitivity; **those synthetic results are not field measurements**. The paper reports around **three heads per biological mother plant** as the field-scale sampling unit for seed counting, nested in treatment/block/site.

**Direct field observations reported by authors**:
- predator-exclusion/model arms did not exhibit detectable differences from control in pollinator visitation, florivore abundance, the fraction of fertilized seeds, or the fraction of attacked seeds, under this study's detection/design;
- pollinator exclusion **substantially reduced fertilized-seed fraction**, consistent with functional female/male separation and insect pollination dependence;
- pollinator exclusion **reduced the number of florivores emerging from heads**, consistent across four field sites in their reported study;
- paradoxically the **proportion of predated seeds was higher** when pollinators were excluded, even though fewer florivorous adults emerged.

This last result is a genuinely **nontrivial multitrophic sign divergence in the same thistle species**. It is an association within a pollinator-exclusion intervention; it does not by itself demonstrate that pollinators directly facilitate florivore survival through edible fertilized achenes rather than changed plant condition, floral access, life stages or sampling denominators. The paper's phrase "florivores depend on pollinators to complete their life cycle" is the authors' biological interpretation. **Adult emergences, damaged-seed fraction, absolute damaged achenes and absolute viable achenes are not interchangeable response variables**. Sex of mother plant, possible mechanical insect exclusion side effects and four treatment levels must be kept distinct.

## What this means for bita × sch × balance × slk

**Critical correction to Q2 novelty:** The hypothesis "changing pollinator access can change seed-predator abundance and plant fitness in Cirsium" has **already been experimentally examined** (2025 preprint, and earlier work on pollination/florivory). Therefore a future paper cannot claim to discover this generic cross-guild interdependence just by re-presenting it in aza3.

A genuinely unresolved *capitulum architecture* hypothesis remains:

> Within the SAME naturally exposed head and independently identified insect taxa, does an experimentally or quasi-experimentally distinguishable **real phyllary/spine/stickiness/orientation interface** control the pollinator's access and the seed-feeder's access **separately**; and does the resulting change in fertilization affect the *feeding quality or establishment success of enemies* as well as the plant's absolute filled achenes?

A simple monotone "barrier protects plant from seed predators" model becomes too narrow if seed-feeder recruitment or survival **also depends on pollinator-created viable seed supply**. Head phenotype may change enemy pressure by affecting **pollen transfer first**, independent of direct mechanical exclusion. In this mechanism, a trait with positive pollen benefit can increase a potential resource base for enemies. The resource-mediated feedback could oppose or amplify any direct head armature defence, but **neither the feedback magnitude nor trait-dependence has been measured for the aza3 morphological modules**.

## Specific real-data causal and denominator precautions

- The four arms are **not** a fully crossed factorial of independent actual pollinator presence × predator presence × florivore identity × head morphology. Thus they cannot estimate `W00(z), W10(z), W01(z), W11(z)` as a factorial morphological fitness response surface or prove selected trait optima.
- A reduction in total emerged seed-feeder adults and an increase in **fraction** of attacked seeds can coexist without a direct per-insect feeding-rate increase. To distinguish mechanisms, the required absolute head-level counts are `N_florets`, `N_fertilized`, `N_damaged`, `N_unfertilized`, `N_emerged by taxon`, and `N_filled viable achenes`, plus their **exact named denominators**.
- If a pollinator exclusion treatment physically covers the head, it might ALSO affect ovipositor access, weather, microclimate and predator accessibility. The paper's actual treatment details and sham/control definitions must be retained before attributing florivore emergence differences uniquely to seed provisioning.
- The study's experimental origin `2022 / four Czech sites` and female-plant sex cannot be pooled as observations from *C. purpuratum* 1993/94 or *C. discolor* 1980/82.
- Flower-stage-specific access and independent guild IDs must come from observed organisms; broad insect order is insufficient.
- Head armature true spine length, orientation to gravity, and actual contact route were **not documented** in the preprint as randomized traits. Even if this dataset becomes re-analyzable, a physical-spine adaptation coefficient remains **NOT_IDENTIFIABLE**.

## Public data archive: verify before claiming any recomputation

The manuscript explicitly names a publicly deposited data/code archive [Zenodo DOI 10.5281/zenodo.15364152](https://doi.org/10.5281/zenodo.15364152). During this session, regular web direct access to that URL was unavailable, and no raw source file has yet been downloaded or inspected. A new **read-only bounded probe** at `analysis/probe_cirsium_arvense_2025_zenodo_v1.py` uses a GitHub Actions runner to retrieve metadata, verify DOI/record ID, list checksums, and report permitted small openly accessible field-file metadata. Absence/access denial is `ARCHIVE_NOT_FETCHED` / `METADATA_UNAVAILABLE`, not a zero experimental effect.

No raw specimen head file, re-fitted effect, new Cirsium field accession, genome, pollen count or spine-contact outcome is asserted on the basis of this record review.

### Novelty decision and immediate question

Q2 has **an actual head-level multi-guild intervention with a reported sign-divergent outcome in the same Cirsium species**. The key missing causal link is not "do pollinators and seed predators affect one another?" but **which head interface/stage produces that coupling, and whether a physical module shifts it toward a net gain or loss in reproductive output**.

Before a broad new field factorial, test whether the open Zenodo release contains the absolute per-head fertilized/unfertilized/predated-seed counts and emerged insect numbers needed for a **descriptive** denominator audit and field effect replication. If the data are available, preserve original experimental randomization and source licenses. If not, retain only the authors' published outcome directions and note access failure.
