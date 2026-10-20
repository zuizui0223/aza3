# GloBI Cirsium temporal follow-up: design fixed before date recovery

Date: 2026-10-20. Follow-up to merged PR #54 / aza3 main ebf4b2e6c2d43f3fe7b16e0a75540ba472723a20.

## Question and controlling cohort
Does the descriptive host-versus-region visitor-composition contrast persist when comparing the SAME calendar month, rather than pooling every flowering month? Keep all 3,243 admitted observation IDs from the immutable 2026-10-08 GloBI archive. No new observations, host selection or post-result threshold changes. Use iNaturalist observation IDs only to recover dates and minimal observer/taxon metadata. The original species labels and archived coarse spatial mapping remain the primary source; present taxonomy agreement is a separate sensitivity.

## Fixed analyses
- Baseline: exact previous 2-degree grids, 10 records per host-region, at least three shared regions per host pair; equal-size rarefaction and 2x2 matching unchanged.
- Primary follow-up: same calendar month across both hosts and both regions; years pooled, so NOT contemporaneous communities. Average within host-pair/month, then equally across admitted months within pair, then equally across pairs.
- Strict sensitivity: same year AND month, same 10-record and three-region rules. Insufficient support must return null, not zero and not a relaxed primary analysis.
- Additional declared sensitivities: same quarter; same month after one record per observer/host/visitor/region/day; same month and same year-month at five records; same month retaining current visitor name/ancestor concordance only. Report every condition.
- Retain the ORIGINAL three-host/five-region core separately, whether or not temporal matching leaves enough comparisons. Changes in support are not causal adjustment effects.

## Recovery and privacy
Batch at most 200 explicit IDs, at least 1.1 seconds between API requests. Cache minimal source fields and per-response SHA256. No photos, exact coordinates, observer login or real name are retained in recovered metadata; observer numeric IDs become SHA256 pseudonyms. Do not use creation time as observation time. Admit only full valid ISO observation dates <=2026-10-08. Missing, withdrawn, future-edited, or malformed dates are reported and not imputed. Changed current taxa do not silently rewrite archived host/visitor labels.

## Inference limits
No p-values from repeated rectangles. Matching calendar month cannot remove interannual conditions, sampling effort, observer preferences or detectability; even year-month is not a matched sampling visit. A shift after matching may also reflect losing host-region units. No morphology, pollination efficacy, adaptive fitness, or re-expression is inferred. No field authorization is changed.

The dates and biological results have not been read when this design is committed. This is a declared exploratory follow-up, not prospective confirmatory preregistration of the original hypothesis.
