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

The NCBI Run Browser experiment design adds an important technical boundary: **unenriched libraries were mixed with target-enriched libraries at a 3:2 ratio before sequencing**. Therefore raw read-to-target mapping percentage contains a designed off-target component and is **technical QC only**, not a reference-transferability estimand. The primary Stage-1 quantities are HybPiper locus recovery, recovered length, paralogy and locus dropout.

The same BioProject also contains public target-capture runs for several aza3-relevant or nearby Japanese taxa. The exact frozen panel is:

`data/planning/reference_public_preflight_runs_v1.csv`

This matters because reference behaviour can be compared across a small phylogenetic gradient instead of judging C. sieboldii in isolation.

## Moreyra target-file reconstruction is possible from public assets

Moreyra et al. did not use only the generic Compositae1061 reference. Their workflow used:

1. the original Compositae1061 target reference; plus
2. the corresponding exons recovered from the highest-coverage `Cirsium tioganum` sample.

The original HybPiper reference is publicly available as:

`github:carol-siniscalchi/Comp1061-Angio353/comp1061_hybpiper_reference.fasta`

The public SRA contains the reported high-coverage augmentation source under the current NCBI taxon name `Cirsium scariosum var. americanum`:

- historical/sample name: `Cirsium tioganum`
- run: `SRR25265669`
- experiment: `SRX21011548`
- BioSample: `SAMN34240347`
- library: `Cirsium-tioganum_21`
- reported bases: 36,939,020,978
- average run read length: 302 bp

Therefore a Cirsium-adapted target can be reconstructed from public materials if the exact intermediate target file cannot be recovered.

The reconstructed file must be labelled `RECONSTRUCTED_CIRSIUM_TARGET`, not treated as byte-identical to the authors' original intermediate file.

## Candidate whole-genome references

Use the same three resources frozen in Gate R1:

1. `C. heterophyllum` — `GCA_965225835.1`
2. `C. dissectum` — `GCA_965276805.1`
3. `C. nipponicum` — public assembly from PRJNA1127082 / Figshare DOI `10.6084/m9.figshare.26927092`

## Two-stage R1B-0 pipeline

### Stage 1 — target-file sensitivity

Run the frozen public comparison panel twice:

A. `ORIGINAL_COMP1061`  
B. `RECONSTRUCTED_CIRSIUM_TARGET`

Use the same HybPiper version and trimming policy across all samples.

For each sample × target file record:
- number of loci recovered;
- exon length recovered per locus;
- fraction of expected target length;
- depth or read support where recoverable;
- paralog warnings;
- loci failing recovery.

Primary question:

> Does adding Cirsium-specific target sequences materially change locus recovery for C. sieboldii relative to the other Japanese samples?

If yes, target-file choice itself is a measurable transfer-bias source and must be propagated into later reference comparisons.

### Stage 2 — genome-reference localization

Use only a frozen clean-locus subset after Stage 1.

For each recovered orthologous locus, localize the sequence independently against:
- `GCA_965225835.1`;
- `GCA_965276805.1`;
- the public `C. nipponicum` assembly.

Record:
- unique vs multiple placement;
- alignment span;
- sequence identity/edit distance;
- chromosome/contig coordinate;
- syntenic consistency between the two chromosome-scale references;
- locus dropout by reference.

This stage asks whether homologous Cirsium loci have stable genomic placements and divergence rankings across available references.

It does **not** infer focal-population local ancestry.

## Clean-locus rule

Paralogy is a major confound in Asteraceae target capture and in the Moreyra workflow.

The published Moreyra et al. (2025) orthology procedure is now recoverable from the article and is the primary historical replication layer:

1. discard genes with **more than 10 HybPiper paralog warnings**;
2. for genes with **1–10 warnings**, recover alternative copies with HybPiper paralog retriever;
3. align and inspect gene trees to retain loci for which ortholog/paralog designation is biologically coherent;
4. discard genes with more than one paralog or obvious ortholog/paralog misassignment;
5. retain loci with **<50% missing data** and **>=80% species presence**.

That process yielded the published 350-locus phylogenomic dataset.

R1B-0 therefore reports three layers:

1. **all 1,061 target loci** for transparent target-recovery diagnostics;
2. **PUBLISHED_RULE_COMPATIBLE** loci, using the published warning/missingness/presence logic as far as the frozen public comparison panel permits;
3. **AUTO_STRICT_CLEAN** loci, defined by a reproducible automated rule for the reference-transfer estimand, excluding unresolved multi-copy placement and problematic paralogy without using C. sieboldii-specific performance.

The automatic clean-locus rule must be frozen before looking at C. sieboldii-specific reference performance and applied identically to all samples and both target-reference versions.

The 350-locus number is **not** hard-coded as the expected R1B-0 output: it came from the full Moreyra taxon set and manual gene-tree curation. A smaller comparison panel may admit a different number while following the same logic.

## Exact R1B-0 estimands

### Target-file estimands
1. loci recovered per sample;
2. recovered target length per locus;
3. Cirsium-target gain over original-target recovery;
4. paralog-warning burden;
5. sample × target-file locus dropout.

### Whole-genome reference estimands
6. fraction of clean loci uniquely localized;
7. alignment-span and identity distributions;
8. chromosome/contig placement concordance;
9. cross-reference rank correlation of per-locus divergence;
10. reference-specific locus dropout;
11. pairwise sample-distance concordance across reference choices.

## Focal decision for C. sieboldii

R1B-0 asks only:

> Is C. sieboldii an outlier in target-file or genome-reference transferability relative to other Japanese Cirsium samples?

A useful result would be:
- the same broad clean-locus set is recovered under both target files;
- one or more genome references uniquely localize most clean loci;
- pairwise relationships among the Japanese comparison panel are stable;
- C. sieboldii does not show exceptional target- or reference-specific dropout.

An adverse result would be:
- C. sieboldii depends strongly on the Cirsium-adapted target for locus recovery;
- C. sieboldii loses a large or nonrandom subset of loci under one or more congener genome references;
- inferred sample-distance structure changes materially with reference;
- many loci have reference-specific multiple placements or strong sequence-identity distortion.

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

Target-file sensitivity is small or well-bounded, captured-locus localization is broadly stable across references, and the Japanese comparison-panel relationship matrix is reference-invariant.

Next:
run the 6–10-individual focal WGS R1B pilot.

### PUBLIC_TRANSFER_AMBER

Public loci are usable for topology/orthology, but target-file dependence, reference-specific dropout or distance distortion remains.

Next:
retain public resources for orthology/topology; design focal C. sieboldii reference in parallel with a small WGS bias audit.

### PUBLIC_TRANSFER_RED

C. sieboldii or the Japanese radiation is strongly and inconsistently represented across target/reference choices.

Next:
focal C. sieboldii HiFi + chromosome scaffolding becomes a prerequisite for dense local-history work.

### NOT_IDENTIFIABLE

Public target-capture structure, target reconstruction, or recoverable metadata are insufficient for a comparable test.

Next:
skip further public-data fishing and proceed to the own-data R1B design.

## Stop rule

R1B-0 is one bounded use of the public target-capture panel.

Do not expand it into a new phylogenomics project. It exists only to lower uncertainty about the reference strategy before own-data sequencing.
