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

## Public compatibility target versus the unrecovered author target

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

The public materials support a **compatibility reconstruction**: Cirsium tioganum homologs can be appended to the public 1,061-locus Compositae1061 reference. This is useful for testing whether a closer Cirsium source sequence changes locus recovery.

However, the author repository and both the 2023 and 2025 Moreyra papers consistently report **1,064 mapped target loci**, whereas the public Compositae1061 FASTA and the named-locus columns recover **1,061 loci**. The identities of the extra three historical target loci are not publicly resolved.

Therefore:

- `PUBLIC_COMP1061_1061` = the reproducible public 1,061-locus compatibility reference;
- `PUBLIC_1061_PLUS_TIOGANUM_COMPATIBILITY` = the same 1,061 locus IDs with added Cirsium-source sequences where recovered;
- neither is the exact historical author 1,064-target file;
- the 1,064-versus-1,061 difference remains an explicit unresolved provenance item, not three inferred loci.

The EAzami public-repository audit is handed off in `data/contracts/aza3_eazami_moreyra_locus_handoff_v1.json`.

## Candidate whole-genome references

Use the same three resources frozen in Gate R1:

1. `C. heterophyllum` — `GCA_965225835.1`
2. `C. dissectum` — `GCA_965276805.1`
3. `C. nipponicum` — public assembly from PRJNA1127082 / Figshare DOI `10.6084/m9.figshare.26927092`

## Two-stage R1B-0 pipeline

### Stage 1 — target-file sensitivity

Run the frozen public comparison panel twice:

A. `PUBLIC_COMP1061_1061`  
B. `PUBLIC_1061_PLUS_TIOGANUM_COMPATIBILITY`

This is a **compatibility sensitivity**, not an attempt to manufacture the unrecovered historical 1,064-target file.

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

EAzami already reconstructed reproducible public locus sets from the authors' released HybPiper stats, sequence-length matrix and paralog report. aza3 reuses those frozen sets rather than rebuilding them after seeing the focal result:

1. **PUBLIC_1061_BROAD_RECOVERY** — all 1,061 reproducible named loci; recovery/dropout diagnostics only.
2. **MOREYRA_COMPATIBILITY_241** — 241 loci with zero public paralog warnings and raw sequence occupancy >=0.80; this is the **primary high-stringency reference-transfer layer**. Frozen file: `data/evidence/moreyra_conservative_241_no_warning_loci_v1.txt`; SHA256 `d561c6e393b1964fdd4b3acf14fda8b10f2f43923b1074cd35f86bfed07ebf73`.
3. **AUTO_STRICT_CLEAN** — optional additional aza3-specific automatic layer, frozen before opening focal transferability results.

For sensitivity only, the reproducible 531-locus warning<=10/high-occupancy set can also be used.

The published final 350-locus matrix is **not** any of these sets. Its identities depend on manual gene-tree decisions and remain unrecovered. Likewise, the paper-reported 1,064-target universe must not be silently collapsed to 1,061.

## Exact R1B-0 estimands

### Target-file estimands
1. loci recovered per sample;
2. recovered target length per locus;
3. Cirsium-source compatibility gain over public-1061 recovery;
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
- the broad 1,061 recovery pattern and the frozen 241-locus high-stringency layer are stable across compatibility target versions;
- one or more genome references uniquely localize most clean loci;
- pairwise relationships among the Japanese comparison panel are stable;
- C. sieboldii does not show exceptional target- or reference-specific dropout.

An adverse result would be:
- C. sieboldii depends strongly on Cirsium-source augmentation for recovery within the public 1,061 named-locus universe;
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
