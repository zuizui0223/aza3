# R1B-0 focal public-data smoke result v1

**Status:** FOCAL_ORIGINAL_TARGET_TECHNICAL_SMOKE_PASS  
**Date:** 2026-10-05

The first empirical R1B-0 execution used the public *Cirsium sieboldii* target-capture run `SRR30887308` against the untouched public Compositae1061 HybPiper nucleotide target.

## Result

- 2,625,581 paired spots = 5,251,162 primary reads (2 × 125 bp).
- Public target: 2,597 reference sequences representing **1,061 loci**.
- **1,033 / 1,061 loci = 97.36%** received at least one mapped read.
- **1,006 / 1,061 = 94.82%** had >=10 mapped reads.
- **977 / 1,061 = 92.08%** had >=50 mapped reads.
- **936 / 1,061 = 88.22%** had >=100 mapped reads.
- Median mapped reads per locus = **683** (Q25 273; Q75 1,421).
- Only **28 / 1,061 = 2.64%** were zero-hit loci.
- 1,280,356 reads mapped to the target (**24.37%**).
- 1,052,392 reads were properly paired (**20.04%**).

The overall mapping fraction must be interpreted in light of the NCBI experimental-design metadata: **unenriched libraries were mixed with target-enriched libraries at a 3:2 ratio** before sequencing. Therefore total read mapping percentage is not the primary transferability metric.

## Immediate consequence

The original generic Compositae1061 target is already broadly reachable by this *C. sieboldii* library.

At an any-hit definition, a Cirsium-specific target augmentation can rescue at most the 28 currently zero-hit loci, so its maximum possible improvement in **gross locus presence** is only **2.64 percentage points**.

That does **not** mean the augmentation is useless. It may still improve:
- assembled exon length;
- sequence identity;
- allele balance;
- paralog discrimination;
- intron-flank recovery;
- per-locus depth.

Those remain the reason for the frozen ORIGINAL_COMP1061 vs RECONSTRUCTED_CIRSIUM_TARGET comparison.

## Reference-record detail

The public target contains two or three source reference sequences per locus.

- 586 loci have two source records: 432 hit both records, 126 hit one, 28 hit none.
- 475 loci have three source records: 328 hit all three, 126 hit two, 21 hit one, **0 hit none**.

Thus the zero-hit problem is confined to a small subset of the two-reference loci.

## NCBI STAT context

Independent NCBI STAT analysis places this sample at:
- identified spots: 23.36%;
- Viridiplantae: 16.50%;
- Asteraceae: 13.97%;
- Bacteria: 3.60%.

Across the frozen ten-run Japanese comparison panel, the *C. sieboldii* Asteraceae value is low but not extreme (third-lowest of ten), while its bacterial fraction is the second-lowest. There is no simple indication that gross contamination explains its lower NCBI identification fraction.

## Claim boundary

This result is **not** HybPiper gene recovery and is **not** the R1B-0 reference decision.

It does not authorize:
- `PUBLIC_TRANSFER_GREEN`;
- local ancestry;
- haplotype age;
- recombination blocks;
- structural-variant reuse;
- module-specific local genealogy;
- genotype-phenotype association.

The next inferential step remains the frozen target-file sensitivity and cross-reference test.
