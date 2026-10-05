# R1B own-data WGS transferability pilot v1

**Status:** DESIGN FROZEN — NOT COLLECTION AUTHORIZATION  
**Date:** 2026-10-05

## Purpose

This pilot is not a GWAS and not the Nature Figure-3 experiment.

Its sole purpose is to determine whether dense genomic inference in focal *Cirsium sieboldii* individuals is stable enough across existing public references to justify proceeding without first building a focal chromosome-scale reference.

The pilot opens only after the bounded public R1B-0 preflight has been classified, or earlier only if own samples already exist under valid authorization.

## Target sample size

Primary target: **8 individuals**.

Allowed range: 6–10 only when availability makes the exact eight-individual design impossible.

Eight is chosen for reference-bias discrimination, not statistical power for genotype–phenotype association.

## Biological sampling logic

Prefer two polymorphic or closely paired populations, four individuals per population.

Within each population, maximize directly observed contrast across capitulum modules, prioritizing:

1. nodding versus upright presentation;
2. colour state/continuous colour contrast where present;
3. phyllary state;
4. stickiness.

If one population does not contain both presentation states, use matched nearby populations and retain population as an explicit factor.

No authority-level species phenotype may substitute for the phenotype of the sequenced individual.

Each pilot individual must link:
- immutable individual ID;
- population ID;
- voucher or diagnostic image;
- direct colour measurement;
- direct orientation measurement;
- direct phyllary score;
- direct stickiness score;
- 2C genome-size measurement or explicit failure code;
- cytotype or explicit failure code;
- DNA extraction ID.

## Sequencing design

Target **~8× mean nuclear depth per individual** on a common short-read WGS platform.

Acceptable pilot range: 6–10× after QC.

Reason:
- sufficient aggregate information for genotype-likelihood PCA/relatedness;
- enough heterozygous sites to compare allele-balance distributions;
- enough window-level depth/callability to expose reference-specific dropout;
- deliberately below the depth required for the final fine-scale haplotype/SV claims.

Use the same library protocol across all individuals. PCR-free is preferred when DNA quantity permits.

Retain high-molecular-weight tissue/DNA from at least one fully phenotyped focal individual so that the same biological system can be promoted to HiFi + Hi-C if the pilot returns AMBER or RED.

## Reference panel

Map the identical reads independently to:

1. *C. heterophyllum* GCA_965225835.1;
2. *C. dissectum* GCA_965276805.1;
3. public *C. nipponicum* assembly.

Do not choose the “best” reference before the cross-reference comparison is complete.

## Primary metrics

Genome-wide:
- mapped-pair fraction;
- unique/high-confidence mapped fraction;
- depth distribution;
- callable genome fraction;
- heterozygous allele-balance distribution;
- genotype-likelihood PCA;
- kinship/relatedness matrix.

Window-level:
- callable fraction;
- missingness;
- normalized depth;
- heterozygosity;
- local PCA or local relatedness;
- reference-specific signal gain/loss.

Bias tests:
- sample × reference interaction;
- population × reference interaction;
- phenotype-state × reference interaction opened only after reference-blind QC is frozen;
- 2C/cytotype × reference interaction.

## Primary reference-invariance principle

The pilot passes only if the **biological relationship structure and the window-level regions that would enter later module-history analyses are materially stable across references**.

A high global mapping rate alone is not a pass.

## Decision classes

### OWN_WGS_GREEN

Genome-wide ancestry/relatedness and candidate local-history structure are reference-invariant enough that existing references can support Phase-A/Tier-2 discovery, with reference sensitivity retained throughout.

This still does not authorize absolute haplotype age or difficult structural-variant claims.

### OWN_WGS_AMBER

Genome-wide ancestry is stable, but local window inference or genotype callability is materially reference-dependent.

Action: build a chromosome-scale focal *C. sieboldii* reference before the Nature Figure-3 local-history test.

### OWN_WGS_RED

Reference effects correlate with major phenotype/ancestry classes or remain severe across all three public references.

Action: build at least one focal HiFi + Hi-C reference and anticipate a two-haplotype/mini-pangenome route.

### NOT_IDENTIFIABLE

Depth, biological sampling, cytotype uncertainty or DNA quality prevents reference bias from being separated from biological signal.

Action: redesign the pilot; do not open module-specific genomic history.

## Hard boundary

This pilot does **not** test whether colour, orientation, phyllary and stickiness map to different loci.

It only determines whether the genomic coordinate system is safe enough to ask that question later.
