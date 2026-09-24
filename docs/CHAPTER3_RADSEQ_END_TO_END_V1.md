# Chapter 3 RAD-seq end-to-end design v1

Status date: 2026-09-02

## Why this document exists

The sampling ledgers define **who** should be collected. This document defines **how RAD-seq is allowed to turn those individuals into ancestry evidence**. The central rule is that RAD-seq is not one universal answer for every Chapter 3 question.

RAD-seq has four bounded roles:

1. **focal population ancestry** for P01/P02 and M01;
2. a possible **all-Japan topology/network sensitivity** for P03/P04;
3. **within-cytotype or ploidy-aware ancestry** in mixed-ploidy systems;
4. the **population-structure background** for M01 selection work.

It is not automatically the definitive Japanese species tree, and it is not allowed to turn a restriction-site outlier into evidence of adaptation.

The existing Moreyra/Comp1061 target-capture framework remains the primary cross-species scaffold. RAD asks the shallower questions that the species-tip tree cannot answer well: individual/population ancestry, reticulation, geographically structured lineages and whether an apparent trait transition survives when species-tip compression is reduced.

## R0 — DNA and identity before library construction

Every production individual must already satisfy the Chapter 3 individual identity contract. RAD material is linked to the same immutable `individual_id` used for phenotype, voucher/diagnostic images, population, cytotype status and authorization records.

Preferred production material is fresh or rapidly silica-dried leaf. Old herbarium DNA is not mixed into the primary production batch unless a separate degraded-DNA protocol is validated first.

Before digestion record:

- fluorometric DNA concentration;
- purity metrics;
- DNA fragment integrity;
- input mass used for the library;
- extraction batch;
- tissue and extraction identifiers.

Whole-genome amplification is not part of the primary protocol. DNA input is normalized before digestion.

## R1 — an empirical RAD assay pilot is mandatory

The production enzyme pair and size window are **not chosen from habit**. They are selected from an empirical pilot because the Chapter 3 panel spans phylogenetic distance, different cytotypes and a large reported genome-size contrast in the M01 case.

The pilot anchor ledger contains eight systems, two DNA templates each:

- `JPN_03 C. alpicola` — high-ploidy anchor;
- `JPN_05 C. aomorense` — diploid benchmark;
- `JPN_06 C. dipsacolepis` — P02 focal;
- `JPN_15 C. lineare` — P02 focal;
- `JPN_35 C. nipponicum` — tetraploid anchor;
- `JPN_36 C. sieboldii` — P01 focal;
- `C. brevicaule` — M01 larger-genome anchor;
- `C. irumtiense` — M01 smaller-genome anchor.

Confirmed cultivated or other assay-only material may be used for this **technical pilot**, but such material does not count toward biological Chapter 3 n.

### R1A — enzyme/fragment screen

At least three candidate rare/common-cutter combinations are screened on one qualified template from each anchor. Candidate pairs may be suggested by in-silico digestion when a suitable nuclear assembly exists, but the final pair must survive the empirical panel.

The pilot asks:

- how many usable fragments/loci are recovered at equal read depth;
- whether one genome-size or ploidy class dominates or collapses;
- whether short-fragment/adaptor contamination is excessive;
- whether the size-selection window is reproducible;
- whether locus recovery is stable across focal systems.

The enzyme pair and size window are frozen **before production and before trait-history outcomes are inspected**.

### R1B — reproducibility screen

The best two candidate protocols are then evaluated across all 16 anchor templates, with at least eight preregistered cross-batch technical repeats.

Metrics include retained reads, depth, usable loci, shared-locus recovery, replicate concordance and signatures of paralog/homeolog collapse.

### Single-protocol viability gate

A single production protocol is used across Japan38 and M01 only if it works across every mandatory anchor and does not create a shared-locus matrix dominated by one genome-size/ploidy class.

If this fails, the correct response is **protocol stratification**, not forcing one bad protocol. For example, the Japan38 core and M01 Ryukyu panel may be processed as separate RAD strata. Different RAD protocols are never concatenated into one SNP matrix. Cross-stratum species history remains anchored by target capture.

## R2 — production library design

The production design is randomized against the biological questions.

Taxa, populations, phenotype states, cytotypes and M01 lineages are distributed across:

- extraction batches;
- digestion/ligation batches;
- PCR/indexing batches;
- sequencing lanes or equivalent run partitions.

No focal population can occur exclusively in one library batch or lane.

Use unique sample indexing and discard ambiguous index assignment. Paired-end sequencing is required within a production stratum because it improves RAD-locus reconstruction and gives additional QC information. The exact read length is frozen after the pilot and kept constant within the stratum.

PCR cycle count is minimized and frozen after the pilot. If a UMI or degenerate-base duplicate-identification scheme is adopted, it must be applied to every library in the stratum; protocols with and without that feature are not mixed.

Technical library duplicates comprise 5–10% of the production libraries. They are selected before biological outcomes are known and deliberately placed in different processing batches. They never increase biological sample size.

## R3 — sequencing depth is pilot-derived, not guessed

No fixed reads-per-sample target is declared before the pilot.

Pilot reads are downsampled to generate saturation curves for:

- number of usable loci;
- median locus depth;
- genotype/allele-depth concordance between technical repeats;
- missingness and locus sharing.

The production target is frozen from these curves before P01–P05 or M01 outcomes are inspected. Libraries are pooled equimolarly after QC, and read imbalance is audited against taxon, population, genome size and sequencing batch.

A low-read sample may be resequenced only when its library complexity and identity show a technical yield failure. Weak or inconvenient biological results are not grounds for recollection or selective resequencing.

## R4 — primary bioinformatics

### Raw processing

Record and report, by sample and batch:

- raw reads;
- reads retained after demultiplexing/trimming;
- adapter contamination;
- base-quality summaries;
- fragment-size/insert diagnostics where available;
- depth and gross contamination outliers.

### Primary assembly: Stacks 2 de novo

The primary assembly is Stacks 2 in de novo mode because no single lineage-matched nuclear reference is assumed across Japan38.

Do not accept default `m/M/n` values without evaluation. Optimize them on outcome-blind technical/pilot metrics. The initial grid is:

- `m = 2–6`;
- `M = 1–6`;
- test catalog `n` around the chosen `M` (`M-1`, `M`, `M+1` where valid).

The optimization target is a stable set of broadly replicated loci (r80-style logic), technical-replicate concordance and low apparent paralog burden. This follows the rationale of Paris et al. 2017 rather than choosing the parameter set that happens to give the preferred biological topology.

Freeze a primary parameter set before trait-history inspection. Keep adjacent parameter sets as assembly sensitivities.

### Reference-aligned sensitivity

Reference-aligned `gstacks` is allowed only when a sufficiently close, lineage-appropriate nuclear reference exists. Mapping all focal taxa to one distant reference and interpreting unequal mapping success as evolutionary divergence is prohibited.

## R5 — missing data are an analysis dimension, not a cleanup nuisance

Restriction-site polymorphism means RAD missingness is partly biological and partly technical. Therefore one universal missing-data threshold is not tuned after seeing results.

Produce predeclared, purpose-specific matrices.

### Focal population matrix

For P01/P02 and within-lineage M01 ancestry, the primary matrix emphasizes loci broadly present within populations (r80-style occupancy), with at least one less-stringent prespecified sensitivity.

Maintain:

- a one-SNP-per-locus matrix for methods requiring approximate marker independence;
- a haplotype/locus-sequence matrix for distance/network and sequence-based sensitivities.

### All-Japan strict matrix

The Japan38 RAD topology/network is built only if the pilot confirms a viable common protocol and shared-locus recovery is not collapsing with phylogenetic distance, genome size or ploidy.

This is a **secondary sensitivity product**, never the definitive Japanese species tree.

If shared-locus information collapses, cancel this product. Do not loosen filters until a desired tree appears. Use the target-capture scaffold and overlay focal RAD population clusters instead.

## R6 — mixed ploidy is handled before SNP pooling

Taxon names are not accepted as perfect proxies for cytotype.

- Diploid panels may use conventional Stacks genotype calls after QC.
- Polyploid panels should retain allele read depths and use a ploidy-aware probabilistic approach such as polyRAD, or another frozen genotype-likelihood/dosage method when its inheritance assumptions are supportable.
- Individuals with unresolved cytotype do not enter dosage-dependent pooled analyses.
- Mixed cytotypes are never forced through one diploid SNP caller.
- Cross-ploidy organismal history remains anchored by target capture; RAD provides within-cytotype or otherwise qualified population information.

Loci are screened for extreme depth, excess heterozygosity relative to the relevant cytotype/population, implausible allele counts and technical-replicate discordance as indicators of paralog/homeolog collapse.

## R7 — ancestry and history outputs

At minimum, focal population RAD analyses include:

- genotype-space visualization such as PCA;
- ancestry/admixture sensitivity;
- pairwise differentiation and within-population diversity with missingness sensitivity;
- a distance/network representation where reticulation is plausible;
- an explicit audit of batch, read depth and missingness against inferred clusters.

The crucial ordering is:

**RAD QC/assembly → ancestry/network frozen → trait-history test.**

Orientation, phyllary, stickiness or flower colour cannot be used to choose `m/M/n`, missingness filters, MAF filters or the preferred ancestry model.

## R8 — what RAD can and cannot do for M01 selection

This is the largest correction to the earlier design.

RAD is useful for M01 E3 to establish:

- within-lineage population structure;
- ancestry and differentiation;
- the neutral demographic/background stratification against which a candidate mechanism must be interpreted.

RAD alone is **not sufficient** for:

- a two-species FST outlier claim;
- a one-statistic outlier scan;
- a historical selective sweep claim at an incompletely qualified RAD locus;
- identifying butterflies, flies, climate or any ecological agent.

A recent RAD-versus-WGS comparison showed that population-specific restriction-site dropout can create false selection signals and that reduced-representation data can miss true selected regions. Therefore, if E2 identifies a replicated pigment/expression mechanism, final E3 candidate-region evidence should preferentially move to a lineage-aware non-RAD confirmation route: targeted capture, amplicon resequencing or qualified low-pass/whole-genome resequencing. RAD remains the population-background layer.

Before any RAD outlier is treated as biologically interesting, test whether locus presence/missingness itself differs among populations or phenotype states. A candidate signal coupled to cut-site dropout is rejected as selection evidence.

## R9 — explicit failure routes

- **No pilot:** no production RAD.
- **No universal protocol across anchors:** stratify protocols or drop the all-Japan RAD matrix.
- **Population confounded with batch:** no population inference until balance is restored.
- **Technical repeats fail:** rebuild/resequence the affected batch before biological inference.
- **Shared loci collapse with divergence:** no all-Japan RAD topology claim.
- **Cytotype unresolved:** no pooled dosage-dependent analysis.
- **RAD conflicts with target capture:** investigate dropout, paralogy, introgression and sampling; do not silently replace the target-capture scaffold.
- **Cut-site dropout correlates with a candidate phenotype/population:** no RAD-only selection claim.
- **Trait outcome suggests a nicer filter:** filter retuning is prohibited.

## Method anchors

This design uses the current Stacks 2 workflow and filtering model; the r80-style parameter-optimization logic of Paris et al. (2017, DOI `10.1111/2041-210X.12775`); polyRAD for ploidy-aware probabilistic genotyping where appropriate (Clark et al. 2019); and the RAD missing-data / restriction-site dropout literature. The 2026 RAD-versus-WGS selection comparison (PMID `41562436`) is treated as a specific warning against using reduced-representation outliers as stand-alone evidence of selection.
