# aza3 v7 — Does the whole 104,890-transcript Cirsium background invalidate selected DFR/ANS-like read matches?


## Completed source result — reference competition does NOT erase all candidate mappings

**Result completed 19 October 2026:** [GitHub Actions run 39231345705](https://github.com/zuizui0223/aza3/actions/runs/39231345705), real public 2020 author FASTA source verified against original **104,890 transcript sequences**, with six source-verified 2026 Ryukyu leaf R1 prefixes (8 MiB each). Previously selected DFR/ANS-like candidates remain supported by **215 strict unique selected-reference read fragments across the six samples** despite all other 2020 *Cirsium japonicum* source unigenes competing. The control panel of 26 nucleotide-composition-shuffled sequences received **zero strict mappings**. Full source data remain temporary on the runner, not committed.

| 2026 focal young-leaf run | Taxon in original study | Complete R1 prefix reads | Strict mapped to original selected 26 after 104,890-reference competition | TR39766 | TR41106 |
|---|---|---:|---:|---:|---:|
| SRR35152724 YGN-4296 | *C. irumtiense* | 130,790 | 37 | 3 | 24 |
| SRR35152725 OKN-4295 | *C. brevicaule* | 130,243 | 31 | 17 | 11 |
| SRR35152729 YRN-4166 | *C. brevicaule* | 117,106 | 26 | 5 | 14 |
| SRR35152730 AMM-4163 | *C. brevicaule* | 116,686 | 29 | 11 | 13 |
| SRR35152731 ISG-4111 | *C. irumtiense* | 129,883 | 62 | 17 | 37 |
| SRR35152732 MMY-4078 | *C. irumtiense* | 122,526 | 30 | 13 | 7 |
| **Total** | **6 public leaf libraries** | **747,234** | **215** | **66** | **106** |

Thus **both** TR39766 (DFR-family homology candidate) and TR41106 (author best-BLAST-hit ANS-family homology candidate) retain strict mappings in **all six original young-leaf samples**, not only when competing with a prescreened panel of 26. This is a qualitatively stronger source-specific read-level result than the original limited-panel observation. However it remains **2020 source-reference nucleotide similarity**, not 2026 individual gene-family orthology or functional anthocyanin competence.

**Do not calculate fold decreases from the original 636/378 counts:** those prior runs had different first-R1 compressed prefix sizes and different reference/identity filters; the samples are the SAME original six biological libraries and cannot be pooled into independent replication. The original paper's young-leaf tissue does not show floral pigment synthesis, cis/trans reactivation, regulatory ancestry or trait reversal. Protein-family function remains unresolved by the 2020 NJ/FastTree gene-tree tests. Full 2020 source transcriptome is also not a complete genome or pan-transcriptome of the Ryukyu lineages.

**Frozen source evidence:** [full-reference competitive mapping across six Ryukyu young-leaf samples](../data/anto_source_handoff/public_sra/ryukyu_six_full_104890_unigene_competitive_mapping_v7.json). [Independent source-count, BioSample and biological claim audit](../analysis/validate_ryukyu_full_104890_unigene_mapping_v7.py) checks the exact report, its 104,890 source count/hash and that *all six* instances remain detectable without any white-to-coloured re-evolution claim.

**Decision:** retain TR39766 and TR41106 as priority **gene-family homology** candidates for focused full-length focal-genomic/transcript and catalytic testing in aza3. Do not assign them a "functional genetic capability present" binary state yet, and do not invoke them as evidence that the southern coloured lineage re-evolved pigmentation. Avoid repeated large public full-source downloads for identical inputs: CI now runs the verification and negative-control tests offline; fresh full RNA/reference acquisition requires explicit manual dispatch.



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
