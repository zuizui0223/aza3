# aza3 R1B-0 execution runbook v1

**Status:** EXECUTION CONTRACT — PUBLIC DATA ONLY  
**Date:** 2026-10-05  
**Parent:** `REFERENCE_TRANSFER_PUBLIC_PREFLIGHT_V1.md`

## Purpose

Execute the bounded public-data preflight without turning it into a new phylogenomics project.

The workflow has three compute phases:

1. reconstruct a Cirsium-adapted Compositae1061 target from public data;
2. run the frozen ten-sample Japanese comparison panel against the original and reconstructed target files;
3. localize a pre-frozen clean-locus subset against the three candidate genome references.

No phenotype or module-specific genomic claim is opened by this workflow.

## Software environment

Use a dedicated environment and record exact resolved versions.

Required core tools:
- HybPiper **2.3.4**;
- SRA Toolkit;
- BWA and samtools;
- MAFFT;
- TrimAl;
- a sequence aligner for genome-reference localization (minimap2 is preferred for assembled locus sequences);
- Python 3.12 for aza3 helper/receipt scripts.

Example environment creation:

```bash
mamba create -n aza3-r1b0 -c conda-forge -c bioconda \
  python=3.12 hybpiper=2.3.4 sra-tools bwa samtools minimap2 mafft trimal fastp
mamba activate aza3-r1b0
hybpiper --version
```

Do not silently change HybPiper major/minor version after the first run. If a version change is required, rerun all compared samples/targets or classify the comparison as non-identical.

## Inputs frozen in the repository

Run panel:
`data/planning/reference_public_preflight_runs_v1.csv`

Original target:
`carol-siniscalchi/Comp1061-Angio353/comp1061_hybpiper_reference.fasta`

Public target audit:
`data/evidence/comp1061_public_target_audit_v1.json`

Expected original-target structure:
- 2,597 FASTA records;
- 1,061 distinct loci;
- source prefixes `lett`, `saff`, `sunf`.

Cirsium augmentation source:
- run `SRR25265669`;
- library `Cirsium-tioganum_21`;
- current NCBI taxon `Cirsium scariosum var. americanum`.

Focal public sample:
- `C. sieboldii` run `SRR30887308`.

## Phase 0 — integrity check

Fetch the public original target and confirm its SHA/content structure before running HybPiper.

The GitHub blob SHA frozen by aza3 is:

`4f89e234007f367ffa8aa5e2be536bc44f31f445`

Do not substitute another Compositae1061 target without recording it as a separate target version.

## Phase 1 — reconstruct the Cirsium-adapted target

### 1A. Download tioganum reads

```bash
prefetch SRR25265669
fasterq-dump --split-files --threads 8 SRR25265669
gzip SRR25265669_1.fastq SRR25265669_2.fastq
```

Record checksums of the downloaded FASTQ files.

### 1B. Trim under one frozen policy

Use the same trimming policy for every run in R1B-0.

The historical Moreyra workflow used Trimmomatic. R1B-0 may use a modern tool such as fastp because the estimand is target/reference sensitivity, not byte-identical historical reconstruction; however **one policy must be used for every sample and both target versions**.

If historical-replication fidelity is prioritized, use a Trimmomatic-compatible policy and record all parameters.

### 1C. Run HybPiper with the original target

For paired reads:

```bash
hybpiper assemble \
  -t_dna comp1061_hybpiper_reference.fasta \
  -r SRR25265669_R1.trim.fastq.gz SRR25265669_R2.trim.fastq.gz \
  --bwa \
  --prefix tiog_original \
  --cpu 16
```

Do not use the future reconstructed target to build the augmentation source; that would be circular.

### 1D. Build the reconstructed target

```bash
python analysis/build_reconstructed_cirsium_target_v1.py \
  --original-target comp1061_hybpiper_reference.fasta \
  --tioganum-root tiog_original \
  --output reconstructed_cirsium_target_v1.fasta \
  --report reconstructed_cirsium_target_v1_receipt.json
```

Required receipt fields include:
- original sequence records;
- original distinct loci;
- number of tioganum genes found;
- number matching target loci;
- number appended;
- target loci without a recovered tioganum sequence.

Classification boundary:
`RECONSTRUCTED_CIRSIUM_TARGET` is a public-data reconstruction of the published logic, **not** a byte-identical claim about the authors' intermediate target file.

## Phase 2 — target-file sensitivity on the frozen public panel

For each run in `reference_public_preflight_runs_v1.csv`:

1. download paired reads;
2. apply the frozen trimming policy;
3. run HybPiper once with `ORIGINAL_COMP1061`;
4. run HybPiper once with `RECONSTRUCTED_CIRSIUM_TARGET`.

Example:

```bash
hybpiper assemble \
  -t_dna TARGET.fasta \
  -r SAMPLE_R1.trim.fastq.gz SAMPLE_R2.trim.fastq.gz \
  --bwa \
  --prefix SAMPLE__TARGET_ID \
  --cpu 16
```

Retain the standard HybPiper directory hierarchy. Current HybPiper post-processing assumes it.

For each sample × target version, record:
- `genes_with_seqs.txt`;
- sequence length/recovery statistics;
- BAM flagstat / mapping summary where generated;
- paralog warnings;
- extracted FNA sequences;
- HybPiper log;
- software/version receipt.

Primary target-sensitivity output:
`sample × 1061-locus` recovery matrix.

Do **not** begin from the published 350 loci. The 350 are a downstream curated subset of the 1,061-locus target universe.

## Phase 3 — orthology layers

### Layer A — ALL_1061_TARGET_LOCI

Used only for recovery/dropout diagnostics.

### Layer B — PUBLISHED_RULE_COMPATIBLE

Reproduce the published Moreyra logic as far as the frozen ten-sample panel permits:

- >10 paralog warnings: discard;
- 1–10 warnings: recover candidate paralogs and inspect gene-tree placement;
- discard loci with multiple unresolved paralogs or ortholog/paralog misassignment;
- require <50% missing data;
- require >=80% sample presence.

The published full dataset yielded 350 loci. The ten-sample preflight is **not required to yield 350**.

### Layer C — AUTO_STRICT_CLEAN

Before viewing focal `C. sieboldii` transferability:
- freeze an automatic rule excluding multi-copy / ambiguous loci;
- apply it identically to every sample and both target versions;
- use this layer for the primary reference-concordance metric.

## Phase 4 — three-reference locus localization

Use assembled clean-locus sequences, not raw target-capture reads, for the primary cross-reference localization test.

Candidate genome references:
1. `C. heterophyllum GCA_965225835.1`;
2. `C. dissectum GCA_965276805.1`;
3. public `C. nipponicum` assembly.

For each locus sequence, align independently to each reference. Example:

```bash
minimap2 -x asm20 -a reference.fa clean_loci.fasta > loci.sam
samtools view -bS loci.sam | samtools sort -o loci.sorted.bam
samtools index loci.sorted.bam
```

The final preset may be revised after a small blind benchmark, but it must then be frozen for all samples/references.

Record:
- unique/multiple placement;
- aligned fraction;
- sequence identity/edit distance;
- chromosome/contig coordinate;
- relative ordering/synteny between chromosome-scale references;
- locus dropout.

## Phase 5 — decision

Do not inspect capitulum phenotypes in R1B-0.

Classify only reference transferability:

- `PUBLIC_TRANSFER_GREEN`
- `PUBLIC_TRANSFER_AMBER`
- `PUBLIC_TRANSFER_RED`
- `NOT_IDENTIFIABLE`

If GREEN:
proceed to the 6–10-individual own-data `C. sieboldii` WGS transferability pilot.

If AMBER:
retain public references for topology/orthology and design the focal long-read reference in parallel with a small WGS bias audit.

If RED:
make focal `C. sieboldii` HiFi + chromosome scaffolding a prerequisite for dense local-history inference.

## Hard stop

R1B-0 does not test:
- local ancestry;
- haplotype age;
- recombination blocks;
- structural-variant reuse;
- genomic modularity;
- genotype–phenotype association.

Those require dense focal population genomic data.
