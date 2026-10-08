# Azami as an empirical baseline for the leaf–capitulum contrast

**Date:** 2026-10-08. **Status:** observational, exploratory bridge; does not alter frozen `zuizui0223/azami` GEB mainline.

## Direct answer

The frozen **Azami** data are usable *now* to evaluate **within-capitulum differential lability and integration**; they are **not** already a leaf-versus-head comparison. The current `azami` numerical archive has 22/27 measured image-derived capitulum endpoints, with leaf-only photos excluded from positive-head measurement tables. Its numerical reproduction begins from frozen head measurements rather than raw photographs. It is essential not to imply leaf segmentation already exists or that the GEB submission should be reopened.

### Frozen biological surface

- Broad environmental/trait atlas: 46,276 observations, 259 taxon labels (not necessarily all suitable for any one morphological contrast).
- Strict nine-construct covariance comparison: 1,734 observations, 42 taxa, 36 pair relations. Relations among taxa typically stronger than within taxa (33/36).
- Frozen result: scale-alignment Spearman rho=0.4391; module cohesion P=0.0013 within and P=0.0365 among (cohesion describes visible measured modules, **not developmental or genetic modules**).
- Direct pairwise RV (`azami/reproducibility/current_reference/upgrade/complete18_construct_pairwise.csv`):
  - image presentation angle × outline elongation: within 0.002343, among 0.339617.
  - image presentation angle × projection prominence: within 0.002798, among 0.003174.
  - projection prominence × projection pattern: within 0.350738, among 0.099930.
  - corolla lightness × hue: within 0.184336, among 0.465806.
This is **not** evidence that all head architecture is fixed. Correlations vary strongly by pairing and biological scale.

### New exploratory reinspection of frozen Azami numerical input (outside Chapter 1)

**Source:** `azami` Actions artifact 9612943217 (`universe/continuous_trait_universe_observation_long.csv`, 631,353,407 uncompressed bytes; exact input file SHA256 `d775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344`). Only `analysis_eligible=True` finite values; for each endpoint select species labels with at least 5 observations. Between fraction is `Var(unweighted_taxon_means) / [Var(unweighted_taxon_means) + mean(within_taxon_sample_variance)]`. This is **descriptive scale-dependent image phenotypic variance**, not genetic heritability, evolutionary rate, causal phenotypic plasticity, adaptation, or leaf phenotype. The cohorts differ by endpoint.

| Endpoint | n observations in admitted taxa | admitted taxa | between fraction |
|---|---:|---:|---:|
| corolla L* | 40,290 | 143 | 0.3406 |
| corolla C* | 40,290 | 143 | 0.3276 |
| head-outline aspect ratio | 38,289 | 141 | 0.1128 |
| head-outline circularity | 38,289 | 141 | 0.1191 |
| head-outline solidity | 38,289 | 141 | 0.0971 |
| head-outline width profile CV | 38,289 | 141 | 0.0879 |
| image-vertical presentation angle | 37,042 | 142 | 0.1223 |
| involucre length/width | 2,182 | 44 | 0.1416 |
| bract projection roughness | 2,182 | 44 | 0.1240 |
| bract spread fraction | 2,182 | 44 | 0.1369 |

**Meaning:** Some colour axes show greater between-taxon fraction than outline/projection axes in these differing observational cohorts. This is consistent with a conserved recognizable capitulum outline with evolvable ecological interfaces, but **does not test the structural conservation of floret arrangement, meristem identity, or leaf lability**. Unbalanced public photographs, lighting, image reference, taxon labels and sampling can shift these numbers.

A defensible direct comparison should FIRST reconstruct a matched individual/species cohort with leaf measurements, then test whether the same qualitative ranking holds. Current Azami data lack that contrast.

## How to reuse Azami rather than restart collection

**Path A: existing Azami numerical data, immediately admitted.** Treat flower-head outline, involucral form, projection geometry, visible colour, and image presentation angle as *separately changing capitulum dimensions*. Compare within/among taxa and within-head functional units using the frozen 18-endpoint common cohort. Do not rename RV / phenotypic integration as trait constraint, genetic modularity, or adaptive lability.

**Path B: upstream full-photo discovery, requires renewed image rights/availability.** The historical `legacy/ch1_global/v2/73_merge_exhaustive_continuous_shards.py` retains detector screening rows for head-absent images (which can include leaf-only photographs), while output measurements contain only head-positive images. The original photo IDs and observation IDs are often present. They permit **a candidate join**, but:
- photo/observation membership `!=` same individual plant unless independently validated visually/voucher;
- one full photo may lack enough leaf area for segmentation;
- photographed basal/middle/upper stem leaves have different ontogenetic morphology;
- a leaf-only image in an observation with a head-only image need not depict the same genotype or plant;
- image copyrights and user photo availability remain separate from MIT code release.
The original photographs are **not** redistributed in Azami's current v3 numerical archive. Reprocessing leaves from already **head crops** is not a meaningful leaf assay. Qualify one photograph as *same plant/leaf/head visible* and then annotate leaf length/width, lobe count/depth, serration, leaf marginal spine length/density, life stage and image pose; register an independent leaf-image detection failure/QC table.

**Path C: direct independent leaf/head material already published.** Michálková et al. (2023), DOI `10.1007/s00606-023-01854-2`, report 129 plants from three *Cirsium* parental species and hybrid crosses, including same-individual middle-leaf shape/lobing/spinosity and involucre width/shape/spinosity from flowering heads, plus plant sex and genomic admixture. Their open supplementary file (Online Resource 6) permits quantitative paired-organ variation comparisons. This is an independent **cross-organ validity/calibration dataset**; don't merge raw raw traits or treat hybrid groups as independent species. Account for sex, hybrid ancestry, geography, size and developmental stage; absent genetic causal manipulation, distinguish hybrid mosaic morphology from genetic architecture.

### Stronger evolutionary–developmental hypothesis, with real prior-art boundary

The conserved unit is potentially **a developmental patterning architecture** of repeated florets encircled by involucre, not every numeric leaf/head shape. In Asteraceae, auxin, LFY/UFO and meristem-patterning networks specify floret/phyllary identity and relative placement (Zhang et al. 2024 `10.1111/nph.19590`; 2019 auxin study PMID:30459264). Involucral phyllaries are modified bracts/leaves, a natural boundary for asking why the *same developmental currency* displays diverse leaf lobing but different constraints once recruited to a reproductive unit.

**Competing predictions:**
- `DEVELOPMENTAL_CONSTRAINT`: leaf-to-phyllary switch reduces accessible shape variation **even when** pollinator/enemy regime differs; stable meristem/floret arrangement and strong developmental covariance.
- `REPRODUCTIVE_STABILIZING_SELECTION`: matched genetic/genotypic variation exists but extreme floral-head morphology lowers effective pollen/seed production; trait fitness optimum and loss steepness vary with pollinator/enemy regimes.
- `FUNCTIONAL_MODULE_REDEPLOYMENT`: outer orientation/spines/stickiness change while internal floret/receptacle organization is retained. Across independently evolving lineages, only the outer modules change in a way predicted by ecological access route and viable-seed fitness.
- `OBSERVATION/ONTOGENY`: the apparent distinction disappears after leaf level (basal vs stem), capitulum phenophase, photo angle, individual identity and sampling selection are matched.

**Falsification:** no reliable same-individual cross-organ leaf/head measurements -> `NOT_IDENTIFIABLE_LEAF_GREATER_LABILITY`. If a matched comparison shows similar variation in leaves and phyllaries, reject the central differential-constraint prediction. If shape divergence persists but its link to seed fitness is absent, adaptation claim fails.

### Repository ownership / preservation

- `azami` stays frozen as the submitted GEB numerical paper; no import of leaf analyses into the primary manuscript.
- This cross-organ test belongs in `aza3` as a separate *test of how complex reproductive organs can preserve architecture while exposing ecological interfaces*.
- `EAzami` supplies ancestry and repeatedly changing trait histories, not direct leaf/head ratio.
