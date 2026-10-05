# aza3 R1B-0 — public-data reference-transfer preflight v1

**Status:** EXECUTABLE PUBLIC-DATA PREFLIGHT — DOES NOT AUTHORIZE LOCAL-ANCESTRY CLAIMS  
**Date:** 2026-10-05

## Why this preflight exists

Gate R1 established that usable public Cirsium resources exist, but the central aza3 genomic claim requires local-history inference in focal Japanese populations. Before new C. sieboldii WGS is collected, the public target-capture data can answer a narrower and cheaper question:

> do the candidate reference backbones behave consistently for homologous captured nuclear loci in C. sieboldii and nearby Japanese Cirsium taxa?

This is a **reference-transferability preflight**, not a substitute for the planned 6–10-individual focal WGS pilot.

## New public-data finding

NCBI SRA query of Moreyra et al. BioProject `PRJNA957074` identifies one public `Cirsium sieboldii` target-capture run:

- run: `SRR30887308`
- experiment: `SRX26291290`
- BioSample: `SAMN44017917`
- library: `Cirsium_sieboldii_366`
- strategy: Targeted-Capture / Hybrid Selection / GENOMIC
- layout: paired
- platform: Illumina HiSeq X Ten
- reported bases: 656,395,250
- average read length: 250 bp

The same BioProject also contains public target-capture runs for several aza3-relevant or nearby Japanese taxa. The exact frozen panel is:

`data/planning/reference_public_preflight_runs_v1.csv`

This matters because reference behaviour can be compared across a small phylogenetic gradient instead of judging C. sieboldii in isolation.

## Candidate references

Use the same three resources frozen in Gate R1:

1. `C. heterophyllum` — `GCA_965225835.1`
2. `C. dissectum` — `GCA_965276805.1`
3. `C. nipponicum` — public assembly from PRJNA1127082 / Figshare DOI `10.6084/m9.figshare.26927092`

## Exact R1B-0 estimands

Because PRJNA957074 is target capture, R1B-0 is restricted to the captured-locus compartment.

For each run × reference:

1. fraction of read pairs mapped;
2. fraction mapped uniquely under the same mapping policy;
3. target-locus breadth covered;
4. target-locus depth distribution;
5. mismatch/edit-distance distribution where the mapper exposes it;
6. allele-balance distribution at callable heterozygous sites;
7. number/fraction of target loci passing frozen callability criteria.

Across references:

8. rank correlation of per-locus coverage;
9. concordance of genotype-likelihood PCA or distance matrix;
10. concordance of per-locus divergence ordering among Japanese taxa;
11. reference-specific locus dropout.

## Focal decision for C. sieboldii

R1B-0 asks only:

> Is C. sieboldii an outlier in reference-specific mapping or locus dropout relative to other Japanese Cirsium samples?

A useful result would be:
- one or more references retain the same broad locus set;
- pairwise genetic relationships among the Japanese comparison panel are stable;
- C. sieboldii does not show exceptional reference-specific dropout.

An adverse result would be:
- C. sieboldii loses a large or nonrandom subset of loci under one or more congener references;
- inferred pairwise relationships shift materially with the reference;
- one reference produces systematic allele-balance or callability distortion.

## What R1B-0 can authorize

If stable:
- proceed confidently to the 6–10-individual C. sieboldii low-cost WGS transferability pilot;
- prioritize the best one or two chromosome backbones plus C. nipponicum as sensitivity reference;
- use public target-capture loci to aid orthology and annotation.

If unstable:
- do not spend on a broad 6–10-individual cross-reference WGS pilot first;
- move directly toward a focal C. sieboldii long-read reference design, while retaining a small WGS check for phenotype-dependent bias.

## What R1B-0 cannot authorize

Even a perfect R1B-0 result does **not** authorize:
- local ancestry inference;
- recombination-block inference;
- haplotype-age inference;
- structural-variant reuse;
- module-specific local genealogy;
- genotype–phenotype association.

Target-capture data interrogate a sparse, bait-defined subset of the genome. The Nature Fig.3 claim still requires dense focal population data.

## Predeclared interpretation classes

### PUBLIC_TRANSFER_GREEN

Captured-locus mapping/callability is broadly stable across references and the Japanese comparison-panel relationship matrix is reference-invariant.

Next:
run the 6–10-individual focal WGS R1B pilot.

### PUBLIC_TRANSFER_AMBER

Global captured-locus signal is usable, but substantial reference-specific locus dropout or distance distortion remains.

Next:
retain public resources for orthology/topology; design focal C. sieboldii reference in parallel with a small WGS bias audit.

### PUBLIC_TRANSFER_RED

C. sieboldii or the Japanese radiation is strongly and inconsistently represented across congener references.

Next:
focal C. sieboldii HiFi + chromosome scaffolding becomes a prerequisite for dense local-history work.

### NOT_IDENTIFIABLE

Public target-capture structure, target definitions or recoverable metadata are insufficient for a comparable three-reference test.

Next:
skip further public-data fishing and proceed to the own-data R1B design.

## Stop rule

R1B-0 is one bounded use of the public target-capture panel.

Do not expand it into a new phylogenomics project. It exists only to lower uncertainty about the reference strategy before own-data sequencing.
