# aza3 v7 — Does the whole 104,890-transcript Cirsium background invalidate selected DFR/ANS-like read matches?

**Status:** Source-backed hypothesis/analysis code. Biological conclusions are admitted only after an actual six-source CI job completes. All active Cirsium work is in aza3; anto is the independent general mathematical accessibility programme.

## Motivation: our strongest current positive signal has an essential denominator problem

The validated 2026 PRJNA1311153 six-source leaf-RNA pilot used 16 MiB of each original ENA R1 compressed FASTQ. Its 1,499,192 complete filtered reads produced 636 strict mappings to a preselected panel of **26** 2020 Cirsium japonicum DFR/ANS-like BLAST-labelled transcripts and zero mappings to 26 composition-shuffled negative controls. Source candidate TR39766 (DFR-like) had 161 strict fragments across all six source plants, TR41106 (ANS-related publisher annotation) had 270 across all six.

These are **real public young-leaf sequence observations**, NOT biochemical pathway function, orthology, absence/retention of complete genomic alleles or ancestral colour re-expression.

**Crucial bias:** MAPQ measures uniqueness *within the reference panel used for alignment*. The previous panel contained only 26 selected transcripts. If other transcripts in the 2020 source transcriptome share their sequences, a MAPQ 20+ alignment to the small panel might not be uniquely attributable to an annotated DFR or ANS family. A composition-shuffled 26-transcript decoy is not a substitute for real homologous paralogs.

## Harder discrimination: full-source competition

1. Obtain Park et al. (2020) publicly deposited author supplementary FASTA (104,890 unigene entries, uncompressed 90,263,821 bytes). Require exact member SHA256 **0c7b7d61a8e2ad5e89b9199347a50545be707329715cdbe0a8a351b0eef2815f** and **104,890 independent FASTA headers**; reject changed data rather than silently substituting another reference.
2. Independently validate that every previously selected 2020 candidate's nucleotide sequence appears **exactly** under the same transcript ID in this original full source. Add 26 composition-shuffled controls. Reference now has 104,890 real sequences + 26 decoys rather than only 26+26.
3. Index the entire reference once in minimap2; align the **first 8 MiB** of registered young-leaf R1 compressed FASTQ for each of the six public Ryukyu samples, with HTTP 206 and exact declared original source size, memory and per-run read caps. No R2, complete FASTQ or gene expression estimate is created.
4. Allow the rest of the 104,864 unselected original source transcripts to **compete** for each read. A read supports a prior DFR/ANS-like source ID only if it still maps with MAPQ>=20, at least 90 aligned bases and NM/alignment-length <=0.20. Record strict other-transcript and decoy counts, and keep low/ambiguous alignment counts.
5. Output the complete per-source aggregate and source provenance, **not** raw FASTQ or 104,890 source nucleotide sequences to Git. The source FASTA is kept in the temporary CI directory.

This experiment is more conservative than the previous 26-transcript mapping because the candidate reference now competes with many more real proteins' transcripts. However it is not a whole genome, contains many possible splice-isoform/assembly duplicates and cannot resolve function even if a candidate remains unique.

## What outcomes would mean

| Result under full-source competition | What can be said | What cannot be said |
|---|---|---|
| TR39766 or TR41106 still receives strict mapped fragments | Young-leaf RNA carries sequence whose best mappable full-2020-source transcript is a previously BLAST-labelled candidate | Functional DFR/ANS, intact gene, genetic ancestry, flower re-expression |
| Candidate mappings become MAPQ-low or change to other 2020 transcripts | The earlier candidate-specific signal is **not robust to real source-reference competition** | The focal plants have lost pigment genes |
| Composition-shuffled decoy attracts strict reads | The mapping scheme has unspecific sequence matches; controls/inference must be revisited | Decoy RNA exists in plant |
| Any source run cannot be fetched through strict Range or source FASTA differs | DATA_ACCESS_OR_SOURCE_PROVENANCE_BLOCKED, not a molecular phenotype conclusion | Any genetic absence or ecological interpretation |

**Absolutely do not compare the 8 MiB whole-reference count directly to 16 MiB selected-panel count as a biological expression fold change.** The source prefix length and reference competition both changed. A qualitative question of whether previously detected candidate IDs retain *any* uniquely mapping supporting reads is the permissible finding; a quantitative abundance/selection comparison is not.

## Why this is particularly important for the aza3 × anto bridge

anto predicts that plant lineages can carry hidden functional resource inventories. aza3 must determine if the inventory exists **functionally** in Cirsium. A selected source transcript's R1 fragments are *homology-compatible evidence*, not the retained genome state C or flower-specific regulatory state F. The current 2020 DFR vs ANR and ANS vs FLS gene-family trees are unresolved under NJ and FastTree. Only independent full-coverage focal alleles, functionally characterized enzyme family paralogs, multi-organ pigments and floral expression can classify those latent states.

The longer-term head LEGO test must separately resolve capitulum orientation/protective phyllary modules and viable-achene fitness; head phenotypic correlation is not a molecular mechanism.

## Reproduction

Offline biological-claim negative controls:

    python analysis/test_ryukyu_full_unigene_competition_v7.py -v

Live capped public-source competitive mapping (needs minimap2, samtools and Biopython):

    python analysis/anto_empirical/ryukyu_full_unigene_competition_v7.py --limit-bytes 8388608 --output /tmp/aza3_all_unigene_Ryukyu_v7.json

The GitHub Actions job uses a bounded source-archive download and six 8 MiB R1 ENA prefixes. It stores only JSON summaries and logs. Unvalidated 2020 source transcripts, gene function and floral reactivation remain separately unknown until independent experimental assays.
