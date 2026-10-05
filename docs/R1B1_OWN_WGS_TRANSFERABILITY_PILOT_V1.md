# aza3 R1B-1 — own-data WGS reference-transferability pilot v1

**Status:** DESIGN FROZEN — NOT COLLECTION AUTHORIZATION
**Date:** 2026-10-05
**Parent gate:** R1B-0 public-data reference preflight

## Purpose

R1B-0 supports two points: public C. sieboldii target-capture recovery is not a failure, and recovered coding loci can be placed reproducibly across available Cirsium references.

R1B-1 asks:

> can dense focal C. sieboldii population data preserve the same ancestry and window-level signal when the reference genome is changed?

This is a reference-transferability pilot, not a genotype–phenotype study and not a final local-ancestry analysis.

## Sample target

Core design: **8 individuals**.
Preferred design: **10 individuals**.

Selection hierarchy:

1. mixed populations containing directly observed variation in one or more capitulum modules;
2. within-population phenotype contrasts;
3. geographically close paired populations with contrasting states;
4. broader geographic endpoints only after the first three priorities.

The core eight should, where biologically available, include at least two populations, at least two individuals per population, at least two within-population phenotype contrasts, and more than one directly measured head-module combination.

If mixed populations cannot be found, record that limitation explicitly rather than treating geography as phenotype replication.

## Same-individual metadata

Every WGS individual must link to immutable individual ID, population ID, voucher or diagnostic image, standardized colour, gravity-referenced presentation angle, direct phyllary state, direct stickiness state, genome size / 2C or explicit missingness, cytotype or explicit missingness, leaf-DNA sample and authorization record.

## Sequencing target

Use paired-end short-read WGS.

Preferred:
- PCR-free where DNA quality permits;
- PE150;
- nominal **8–10×** per individual against an approximately 0.9–1.0 Gb Cirsium genome.

This depth is for genotype-likelihood, callability and reference-bias diagnostics, not final causal-haplotype or structural-variant work.

## Reference panel

Primary coordinate backbone:
1. C. heterophyllum hap1 — GCA_965225835.1

Within-individual control:
2. C. heterophyllum hap2 — GCA_965225975.1

Independent coordinate sensitivity:
3. C. dissectum hap1 — GCA_965276805.1

Within-individual control:
4. C. dissectum hap2 — GCA_965276745.1

East-Asian sequence/annotation anchor:
5. C. nipponicum public assembly — Figshare DOI 10.6084/m9.figshare.26927092

Use identical preprocessing and mapping policy across all references.

## Primary logic

The same WGS reads are mapped to hap1 and hap2 of the same reference individual. Those differences define an empirical reference-noise floor.

Then ask whether switching among heterophyllum, dissectum and nipponicum changes the biological conclusion more than this within-individual haplotype-reference noise.

## Stage A — phenotype-blind

Before phenotype labels are opened, freeze:
- mapped-pair fraction;
- uniquely mapped fraction;
- depth and callable fraction;
- heterozygous allele balance;
- genotype-likelihood missingness;
- genome-wide PCA;
- pairwise genetic-distance and kinship matrices;
- window-level callable fraction;
- local PCA/local-distance summaries where SNP density permits.

## Stage B — population labels

Open population labels only after Stage A is frozen.

Test reference × population effects on callability, allele balance and population-structure summaries. A reference that creates or erases population structure is unacceptable for Nature Fig.3.

## Stage C — phenotype labels

Only after Stage A/B are frozen, open colour, orientation, phyllary and stickiness.

Test whether any state predicts reference-specific missingness, callability or allele-balance artefacts. This is a bias test, not an association scan.

## Decision classes

### WGS_TRANSFER_GREEN

Genome-wide ancestry/relatedness is stable, window-level diagnostics are stable relative to same-individual haplotype noise, and no phenotype class shows reference-linked missingness.

Action: keep heterophyllum hap1 as primary coordinate backbone, dissectum as independent sensitivity reference, nipponicum as East-Asian sequence/annotation anchor, and proceed without first building a focal C. sieboldii reference.

### WGS_TRANSFER_AMBER

Genome-wide structure is stable but a meaningful subset of windows is more reference-sensitive than the haplotype-noise baseline.

Action: public references remain valid for ancestry/triage, but build a focal C. sieboldii long-read reference before Nature Fig.3 local-history claims.

### WGS_TRANSFER_RED

Reference substitution materially changes population structure, relatedness, or creates phenotype-linked callability/missingness.

Action: focal C. sieboldii HiFi + chromosome scaffolding becomes mandatory; if bias persists against one focal linear reference, move to a mini-pangenome/graph route.

### NOT_IDENTIFIABLE

Coverage or sample design cannot distinguish reference artefact from biological structure.

## Claim ceiling

Even GREEN authorizes reference-framework adequacy only. It does not establish adaptive loci, causal variants, haplotype age, trait introgression, recombination-block reuse, SV reuse, genomic modularity or genotype–phenotype association.

## Stop rule

Do not scale this pilot into the full Phase-A panel before the reference-transferability decision is made.
