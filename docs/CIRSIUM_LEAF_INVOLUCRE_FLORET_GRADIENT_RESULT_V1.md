# Cirsium leaf–involucre–floret hierarchy: independent quantitative result

**As of 2026-10-08.** This is a new, bounded **empirical reanalysis**, separate from Azami's frozen GEB Chapter 1. It **does not** establish why the hierarchy exists, natural selection, genetic canalization, developmental constraint or an adaptive evolutionary origin.

## Direct question

Is the widespread impression that thistle leaves vary much more than floral heads supported when **leaf, involucre and individual floret traits are measured on the same plants**? Or is a more specific organ hierarchy hidden behind the coarse notion of "stable capitulum"?

## Source 1: published European Cirsium individual-level table, newly reanalysed

Michálková et al. (2023), *Plant Systematics and Evolution*, DOI [10.1007/s00606-023-01854-2](https://doi.org/10.1007/s00606-023-01854-2), **Online Resource 6**, raw XLSX https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs00606-023-01854-2/MediaObjects/606_2023_1854_MOESM6_ESM.xlsx

- Original dataset: **129** individually identified plants (3 parental species plus 2 hybrid categories); the primary test uses **114 non-hybrid plants**: *C. acaulon* n=40, *C. bertolonii* n=34, *C. erisithales* n=40.
- Six separately retained sex × species strata: 10 female and 30 hermaphrodite *C. acaulon*, 10 female and 24 hermaphrodite *C. bertolonii*, 10 female and 30 hermaphrodite *C. erisithales*. All compared endpoints observed on the same individual plant, without missing measurements.
- Raw archive obtained from Springer Online Resource 6 through aza3 GitHub Actions **37751715445** / artifact **11537953289**, SHA256 `c06423cc6f9800ba4775a4428afa072ba1b3dc6fb48ccd05fdb7e42a2f03ee79`. The third-party workbook is retained as a transient data artifact, not silently included in GitHub source files.
- Read using `artifact_tool`, the source's own character IDs: middle leaf length `Lf_L`, leaf width `Lf_W`, leaf length/width ratio `Lf_Sh`; involucre length `Inv_L`, width `Inv_W`, length/width ratio `Inv_Sh`; corolla length `Cor_L`.
- Predeclared estimand: `CV = sample SD / arithmetic mean`, then leaf/other-organ **CV ratio**, geometric-mean-aggregated over the six strata with equal weights. Paired individual resampling inside fixed stratum, **3,000 replicates**, random seed `20261008`. Confidence intervals are conditional individual-sampling intervals for **these three taxa and these sampling localities**, not cross-taxon/phylogenetic replicated uncertainty.
- Raw count or log size is not directly compared across centimetres and millimetres; CV is unit invariant. Floret length is the **corolla length**, not whole capitulum length.

### New main results

| Compared morphological features | Leaf/other relative CV ratio | Conditional 95% bootstrap interval | Positive in all six species×sex strata? |
|---|---:|---:|---|
| Leaf length / involucre length | **2.362** | **1.882–3.136** | YES |
| Leaf width / involucre width | **1.685** | **1.303–2.042** | YES |
| Leaf aspect ratio / involucre aspect ratio | **1.492** | **1.226–1.856** | YES |
| Involucre length / corolla length | **2.786** | **2.222–3.257** | YES |
| Leaf length / corolla length | **6.582** | **5.360–7.986** | YES |

Critically, **each of the six sex×species strata showed a strict decrease in length CV from leaf to involucre to floret**. Conditional bootstrap draws preserved this complete ordering in all six strata in **91%** of 3,000 draws. This is stronger than a test that averaged only three taxa without accounting for sex. The log-SD sensitivity comparisons also had ratios >1 for leaf length, width and aspect versus corresponding involucre features (2.424, 1.677, 1.461).

**What is proven:** in these plants, the distribution of relative *length* variation has an **outer→inner organ gradient**, and leaf size/shape is relatively more variable than corresponding involucral measures. **Not proven:** the rate of evolutionary change, genetic accessibility, fitness optimum, or whether pollen-/herbivore-driven stabilizing selection generated the pattern.

## Source 2: independent 2026 Ryukyu/Taiwan comparative morphometric table

Chang et al. (2026), DOI [10.1186/s12870-026-08097-6](https://doi.org/10.1186/s12870-026-08097-6), Tables 1 and 2 (https://pmc.ncbi.nlm.nih.gov/articles/PMC13020037/). Recompute `CV=SD/mean` from exact source summaries; do not treat six taxon/morphotype groups as independent evolutionary origins or assume individual-level paired sampling.

| Focal taxon/morph | Rosette-leaf-length CV ÷ head-length CV | Head-length CV ÷ floret-length CV | Leaf-length CV ÷ floret-length CV |
|---|---:|---:|---:|
| *C. japonicum* var. *albescens* | 1.894 | 2.802 | 5.306 |
| *C. brevicaule* | 6.594 | 1.798 | 11.853 |
| *C. irumtiense* | 2.969 | 9.903 | 29.405 |
| *C. japonicum* var. *takaoense*, white morph | 2.741 | 1.303 | 3.571 |
| *C. japonicum* var. *takaoense*, purple morph | **0.735** | 8.777 | 6.452 |
| *C. japonicum* var. *australe* | 1.399 | 1.794 | 2.510 |

Summary:
- Leaf-length CV > whole-capitulum-length CV in **5/6** groups.
- Leaf-width CV > whole-capitulum-width CV in **4/6** groups.
- Whole-capitulum-length CV > floret-length CV in **6/6** groups.
- Leaf-length CV > floret-length CV in **6/6** groups.

The purple *takaoense* morph is a **counterexample to universal leaf > whole-capitulum length variability**; head CV actually exceeds leaf CV. Yet floret CV remains low. In addition, the same paper documents visible between-lineage differences in floret length and colour: **low within-group variability is not a statement that florets never evolve.**

## Actual Azami reuse and its boundary

`zuizui0223/azami` already measures 22 visual endpoints and nine biological constructs. On the frozen common cohort (1,734 photos/heads, 42 taxa), an image-based **angle×head elongation** among-taxon RV is **0.339617**, but **angle×bract projection prominence** among-taxon RV is **0.003174**. Some capitulum trait combinations are integrated, others nearly independent on this image scale. Neither RV is a signed evolutionary transition, direct spine measurement, or selection coefficient.

The Azami Chapter 1 numerical archive intentionally excludes *leaf-only photographs* from head-measurement aggregates; it cannot directly yield a leaf↔head anatomical CV ratio or measure botanical corolla length. The two independent literature sources do that for limited subsets. Keep the original GEB paper frozen; this test is an **aza3 successor**.

## What could *cause* this gradient? Explicit competitors

1. **Organ-specific developmental canalization.** Repeated floret primordia share positional identity and coordinated pollen-presentation architecture, reducing within-plant variation in floret dimensions despite broader variation in leaf and outer protective bract development. Prediction: the core floret dimensions stay relatively stable under common-garden changes in light, water and herbivory; genetic variation and development need to be directly measured.
2. **Stabilizing reproductive selection.** Extreme floral dimensions disrupt anther/stigma contact and/or pollen/achene development more than comparable vegetative shape changes harm growth. Prediction: controlled floret-size extremes have lower viable seed production and/or pollen transfer, and natural pollinator regimes shift curvature or optimum. Our CV ratios do **not** already measure these fitness curves.
3. **An evolvable interface shielding a conserved core.** Involucral phyllaries, spines, orientation and stickiness change arthropod access/wetting while reproductive florets remain functionally stable. Prediction: safe manipulations of outer geometry show strong **guild-specific access and seed-fitness effects** without equivalent changes in pollen presentation; net fitness can reverse when pre-dispersal seed predators vs parasitoids change.
4. **Allometry / tissue stage / observation design.** Vegetative leaves may respond more strongly to growth opportunity, leaf position and maturity than adult heads or florets. Prediction: site- and plant-size-matched cohorts, phenophase alignment, repeatability controls and calibrated anatomy shrink or remove the leaf/head difference. The 2026 purple morph demonstrates why a universal rule is unsafe.
5. **Historical constraints versus present maintenance.** Contemporary stabilizing selection does not prove the original appearance of the organ pattern. If a floret program is developmentally conserved but ecological regimes favour alternative outer modules, the same head architecture may persist while the external modules evolve mosaic combinations. Requires ancestry, genomic segregation and timed fitness tests before coevolution/adaptation claims.

### Most discriminating hypothesis for aza3

> *The relative evolutionary freedom of a floral organ decreases toward the reproductive core, but functionally exposed outer modules can be reassembled in ways that depend on the temporal order and 3D access routes of pollinators, seed predators and their parasitoids.*

This is **only a prospective hypothesis** built on a newly established **phenotypic** outer→inner variation gradient. A Nature-caliber adaptive claim needs an independently predicted and tested fitness consequence of **the actual coupled organ architecture**, not another variance-ratio correlation.

## Reproducibility and next validation

- `analysis/retrieve_leaf_capitulum_source_2023.py`: original supplement fetch + SHA identity.
- `analysis/estimate_cirsium_cross_organ_variability_v1.py`: `artifact_tool` XLSX import + individual-level CV, log-SD sensitivity, sex × species bootstrap and leaf–involucre–corolla 3-tier gradient.
- `analysis/check_cirsium_ryukyu_cross_organ_cv_2026_v1.py`: independent external table summary CV and 5/6 vs 6/6 directional stress tests, runnable with standard Python.
- `data/evidence/cirsium_leaf_involucre_floret_variability_result_v1.json`: frozen quantitative headline, source identities and claim ceilings.

No original workbook bytes or copyrighted raw public images are silently committed; the publisher's openly licensed supplemental data remain source-linked.
