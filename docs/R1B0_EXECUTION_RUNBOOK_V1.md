# aza3 R1B-0 execution runbook v1

**Status:** EXECUTION CONTRACT — PUBLIC DATA ONLY  
**Date:** 2026-10-05  
**Parent:** `REFERENCE_TRANSFER_PUBLIC_PREFLIGHT_V1.md`

## Purpose

Execute the bounded public-data preflight without turning it into a new phylogenomics project.

The workflow has three compute phases:

1. establish focal recovery against the reproducible public 1,061-locus Compositae1061 compatibility reference;
2. test sensitivity to adding public Cirsium tioganum source sequences **without claiming reconstruction of the unrecovered historical 1,064-target file**;
3. localize the frozen 241-locus high-stringency compatibility layer (plus prespecified sensitivities) against the three candidate genome references.

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

Expected public-target structure:
- 2,597 FASTA records;
- 1,061 reproducible named loci;
- source prefixes `lett`, `saff`, `sunf`.

Historical provenance boundary:
- Moreyra 2023/2025 report 1,064 mapped target loci;
- the public author matrices and public FASTA expose 1,061 named loci;
- the three-locus difference is unresolved;
- do not call the public 1,061 FASTA the exact author target.

Frozen high-stringency compatibility layer:
- `data/evidence/moreyra_conservative_241_no_warning_loci_v1.txt`
- 241 loci
- SHA256 `d561c6e393b1964fdd4b3acf14fda8b10f2f43923b1074cd35f86bfed07ebf73`.

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

## Phase 1 — optional public Cirsium-source compatibility target

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
the output is `PUBLIC_1061_PLUS_TIOGANUM_COMPATIBILITY`. It adds source sequences only to the 1,061 public locus IDs. It **does not reconstruct or explain the unresolved paper-reported 1,064-target universe** and is not a byte-identical author-file claim.

## Phase 2 — target-file sensitivity on the frozen public panel

For each run in `reference_public_preflight_runs_v1.csv`:

1. download paired reads;
2. apply the frozen trimming policy;
3. run HybPiper once with `PUBLIC_COMP1061_1061`;
4. run HybPiper once with `PUBLIC_1061_PLUS_TIOGANUM_COMPATIBILITY`.

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

Primary target-sensitivity outputs:
- `sample × 1,061 named-locus` broad recovery matrix;
- recovery of the frozen `MOREYRA_COMPATIBILITY_241` layer.

Do **not** equate either with the published final 350, and do not treat 1,061 as the exact historical target count.

## Phase 3 — orthology layers

### Layer A — PUBLIC_1061_BROAD_RECOVERY

Used only for recovery/dropout diagnostics.

### Layer B — MOREYRA_COMPATIBILITY_241

Primary high-stringency transfer layer imported from the EAzami public-author-repository audit. It contains 241 loci with zero public paralog warnings and raw sequence occupancy >=0.80.

This is **not** the published 350-locus matrix.

### Layer C — AUTO_STRICT_CLEAN

Before viewing focal `C. sieboldii` transferability:
- freeze an automatic rule excluding multi-copy / ambiguous loci;
- apply it identically to every sample and both target versions;
- use this layer for the primary reference-concordance metric.

## Phase 4 — three-reference locus localization

Use assembled clean-locus sequences, not raw target-capture reads, for the primary cross-reference localization test.

Primary candidate genome references:
1. `C. heterophyllum GCA_965225835.1` (hap1; 17 chromosome-scale scaffolds);
2. `C. dissectum GCA_965276805.1` (hap1; 17 chromosome-scale scaffolds);
3. public `C. nipponicum` assembly — Figshare file `48979489`, `C.nipponicum_softmasked_genome.fa`, MD5 `e9390e23ffd0dc5e3da8271db4d1d3ca`.

Reference-internal controls:
4. `C. heterophyllum GCA_965225975.1` (hap2);
5. `C. dissectum GCA_965276745.1` (hap2).

The two hap2 controls answer a different question from the three primary references: how much can locus placement change merely by switching haplotypes from the **same diploid individual**? This is a negative-control scale for interpreting among-species reference differences. NCBI reports explicit haplotypic inversions in the `C. dissectum` pair on chromosomes 1, 6 and 7.

Do not infer that a cross-reference placement difference is biologically meaningful unless it exceeds or is qualitatively distinct from this within-individual haplotype sensitivity.

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
- locus dropout;
- hap1-vs-hap2 placement stability within `C. heterophyllum`;
- hap1-vs-hap2 placement stability within `C. dissectum`.

The primary three-reference comparison is interpreted only after the within-individual haplotype baseline is known.

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
