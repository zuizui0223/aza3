# Within-season Cirsium visitor association: declared v3 design

2026-10-20. Follow-up to merged PR #55 (cd94a991dcc86b5e109f4def22efc7755eaa22fb). Exploratory design, not independent preregistration: original visits and prior pooled results were already inspected. Before computing the new statistic, only host/time coverage has been checked.

## Question and fixed source

Within one original 2-degree region and calendar month, are source-recorded visitor species more similar within the same Cirsium host than across two hosts? Does any excess survive holding each observer's host allocation fixed?

Use the SAME 3,243 admitted adult visitsFlowersOf observation IDs; no new records and no imputed dates. Preserve the original GloBI ZIP SHA256 ce0513fde35d9c1be42cef87da7a6f4b262e9b4fe237ba0b78e9bf31771b58a4 and frozen minimal source-date SHA256 313142b1e4fcb1a179533c9b985ea8364b221122ed158b9e2d312ee213a731c9. Source names are not newly reidentified accepted taxa.

## Primary estimand, support and weights

For visitor counts a_s and b_s on two hosts within the same region-month, n_A,n_B >=10, compute

U = 0.5 * [sum_s a_s(a_s-1)/(n_A(n_A-1)) + sum_s b_s(b_s-1)/(n_B(n_B-1))] - sum_s a_s b_s/(n_A n_B).

This is within-host minus cross-host probability of matching recorded visitor identity, excluding self-comparisons. Under independent observations its population target is half the squared L2 distance between visitor probability vectors (half unbiased MMD-squared with the equality kernel). The estimator can be negative in finite samples; do not truncate it or interpret negative as a biological negative distance. Its exact expectation under unrestricted host-label reassignment within a region-month is zero. This fixes the positive small-sample baseline of raw dissimilarities without changing observation effort.

All admitted host pairs are reported; a pair must have at least THREE distinct geographic cells for the primary replicated summary. Average strata within each cell, then cells equally. Months may differ among cells because this is within-stratum host association, NOT the earlier two-region/two-host geographic-dominance estimand. Do not claim that changing estimands removes the previous empty year-month test.

## Competing observational explanation

Preserve each observer's recorded host totals inside every region-month and reassign the visitor identities only within that observer block. Compute the reference expectation analytically from finite-population hypergeometric means/variances. Report the observed U, this observer-conditioned expectation, and their difference. Single-host observers cannot exchange host labels: retain them and report the number of records/observers permitting exchange. A large conditional expectation means the observed host/visitor association is compatible with observer composition under this reference, not that observers causally explain nature.

Use 9,999 fixed-seed label permutations only as conditional-reference diagnostics for primary supported pairs; intervals are null reference ranges, not biological confidence intervals. The labels were not randomized in the field. Observer/day, fine-scale place, calendar year, flowering, visibility, abundance, and survey effort may still confound results.

## Fixed sensitivities

Retain monthly 10-record primary; daily observer deduplication; current visitor-name/ancestor concordance; monthly 5-record exploratory extension; strict year-month 10 and 5; and an insect-order version on the exact primary admitted records. Report all support losses and do not replace null/missing outcomes by zero. Include leave-one-region-out and leave-one-observer-out diagnostics on frozen primary units, stating any deletion that leaves fewer than two records for U.

## Verification and limits

Test the count formula against explicit ordered-pair enumeration, exact null expectations against all small allocations, observer-conditioned expectation against brute-force block allocations, label-swap symmetry, empty support, duplicate-ID rejection, unequal sample sizes and no-information observer blocks. Recompute actual outputs independently before claiming completion.

This analysis does not measure head angle, phyllary geometry, spines, pollen deposition, fitness, preference conditional on floral availability, or heritable re-expression. It does not authorize collection or alter the 298 planned Phase-A individuals. Keep original Q1/Q2/Q3 conclusions intact. Outputs must omit exact coordinates and observer identifiers.

Methodological context: Gretton et al. (2012), A Kernel Two-Sample Test, Journal of Machine Learning Research 13:723-773, https://jmlr.org/papers/v13/gretton12a.html. The equality-kernel specialization and finite-block expectations above are derived and tested here, not a new statistical theorem.
