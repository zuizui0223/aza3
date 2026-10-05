# aza3 R1B-0 author-matrix recovery result v1

**Status:** TARGET RECOVERY GREEN — GENOME-REFERENCE LOCALIZATION STILL PENDING  
**Date:** 2026-10-05  
**Scope:** public-data preflight only

## Question

Before spending on a focal *Cirsium sieboldii* genome or new WGS, does the existing Moreyra et al. target-capture dataset show that the focal species is poorly recoverable at homologous nuclear loci?

## Direct public result

The Moreyra author repository contains the actual HybPiper summary table and per-locus sequence-length matrix used for the published analysis.

For the focal sample `Cirsium-sieboldii_366` (public run `SRR30887308`), the author HybPiper output reports:

| Metric | C. sieboldii |
|---|---:|
| Reads after the author workflow | 4,109,610 |
| Reads mapped | 1,317,250 |
| On-target | 32.1% |
| Genes mapped | 1,043 |
| Genes with contigs | 1,021 |
| Genes with sequences | **1,017** |
| Genes >=25% target length | 1,012 |
| Genes >=50% target length | **955** |
| Genes >=75% target length | **777** |
| Long paralog warnings | 237 |
| Depth paralog warnings | 275 |

Using the reproducible public 1,061 named-locus universe as a compatibility denominator, `GenesWithSeqs=1,017` corresponds to **95.9% broad recovery**.

This denominator is explicitly a compatibility layer: the paper reports 1,064 mapped target loci, and the exact historical 1,064-target intermediate file has not been recovered.

## High-stringency result

EAzami had already audited the Moreyra author repository and froze a conservative locus layer:

- no paralog warning;
- raw sequence occupancy >=0.80;
- **241 loci**.

The per-locus author sequence-length matrix shows:

> **C. sieboldii recovered 241 / 241 of these loci.**

Recovered sequence:
- total = 99,720 bp;
- median locus = 402 bp;
- minimum = 81 bp;
- maximum = 1,059 bp.

Thus the strongest currently reproducible public orthology layer has **zero focal-locus dropout** in the author output.

## Frozen ten-sample Japanese comparison panel

High-stringency 241-locus recovery is:

| Sample | recovered / 241 |
|---|---:|
| **C. sieboldii_366** | **241 / 241** |
| C. japonicum_horridumFJ315 | 241 / 241 |
| C. lineare_PE229 | 241 / 241 |
| C. nipponicum_276 | 241 / 241 |
| C. yoshinoi_FJ405 | 237 / 241 |
| C. tonense_FJ403 | 240 / 241 |
| C. pendulum_MW157 | 238 / 241 |
| C. dipsacolepis_FJ363 | 241 / 241 |
| C. tanakae_271 | 241 / 241 |
| C. maackii_MW159 | 241 / 241 |

Seven of ten samples recover all 241 loci. *C. sieboldii* belongs to the complete-recovery group.

For the broad HybPiper metrics, the focal sample is also ordinary rather than poor:

- on-target 32.1%; panel median 36.0%;
- genes with sequences 1,017; panel median 1,016.5;
- genes >=50% length 955; panel median 946.5;
- genes >=75% length 777; panel median 765.5.

## Paralog-warning check

Within the frozen ten-sample panel, *C. sieboldii* has the largest long/depth warning counts.

That initially looks concerning, but the full author matrix contains 265 *Cirsium* samples. Relative to those 265:

- long-warning count 237 is around the **63rd percentile**;
- depth-warning count 275 is around the **46th percentile**.

Therefore the focal warning burden is **not exceptional at genus scale**. It is not currently evidence for a *C. sieboldii*-specific duplication/homeology problem.

## Decision

### Target-recovery stage: GREEN

The current public data do **not** support a focal target-capture transfer failure.

Specifically:
- broad locus recovery is high;
- the conservative 241-locus layer is complete;
- focal sequence recovery is near or above the frozen Japanese-panel median;
- paralog warnings are not exceptional across the genus.

### Whole-genome reference-transfer stage: PENDING

This result does **not** answer whether the recovered *C. sieboldii* loci have stable placement, divergence ranking or local genomic context when mapped to:

1. *C. heterophyllum*;
2. *C. dissectum*;
3. *C. nipponicum*.

Therefore the overall `PUBLIC_TRANSFER_GREEN / AMBER / RED` decision remains closed.

## Consequence for the reference strategy

The evidence now argues **against immediately building a focal genome merely because C. sieboldii sequence is too divergent for Cirsium target capture**.

The next discriminating test is reference localization, not additional proof of target recovery.

A focal HiFi + chromosome assembly becomes justified if:
- homologous loci localize inconsistently across existing references;
- dense WGS produces reference-specific callability/local-history distortions;
- or structural variation becomes central to the Nature claim.

## Claim ceiling

This result is read from the authors' public output matrices, not yet from an independent aza3 rerun.

It does not establish:
- local ancestry;
- haplotype age;
- recombination blocks;
- structural-variant reuse;
- genomic modularity;
- genotype–phenotype association.

The independent focal HybPiper workflow remains useful as a reproducibility check, but it is no longer needed to answer the narrow question of whether the authors themselves recovered the focal loci.
