# Real Ryukyu Cirsium leaf-RNA evidence: candidate fragments persist under 104,890-transcript competition

**Actual publicly sourced RNA sequence analysis, 19 October 2026.** Workflow [39231345705](https://github.com/zuizui0223/aza3/actions/runs/39231345705) completed successfully, including full-source checksum control. Numerical results are in [the 6-specimen data freeze](../data/anto_source_handoff/public_sra/ryukyu_full_104890_reference_competition_8MiB_v7.json) and checked by [an independent, standard-library-only evidence validator](../analysis/validate_ryukyu_104890_full_reference_v7.py).

## Main source-backed outcome

The original 2026 Ryukyu comparison includes **three white-flowered taxon C. brevicaule and three coloured-flowered taxon C. irumtiense public YOUNG-LEAF RNA libraries** (PRJNA1311153). Neither taxon's published corolla colours are individual measured pigment concentrations in these six RNA donors.

All six R1 libraries were examined using only the first **8 MiB compressed bytes per sample**, as explicitly registered and verified by 206 Content-Range. The read data are an ORDERED gzip file prefix, not an independent random transcriptomic subsample or a complete RNA library. Across the six prefixes, **747,234 valid R1 reads** passed source processing.

The biological advance was not another selected-family BLAST: all previous 26 Park (2020) *Cirsium japonicum* colour-pathway-labelled transcript candidates were forced to compete with the **entire published 104,890 author unigene FASTA**, validated by original source-member SHA256 **0c7b7d61a8e2ad5e89b9199347a50545be707329715cdbe0a8a351b0eef2815f**. Twenty-six composition-shuffled control sequences were also admitted. Reproducible criteria required MAPQ >= 20, >=90 aligned bases and NM/aligned length <=0.20.

| Submitted 2026 taxon | Source specimen | ENA R1 | Complete R1 reads | Strict fragments to preselected 26 source genes | DFR-like TR39766 | ANS-annotated TR41106 |
|---|---|---|---:|---:|---:|---:|
| *C. irumtiense* | YGN-4296 | SRR35152724 | 130,790 | 37 | 3 | 24 |
| *C. brevicaule* | OKN-4295 | SRR35152725 | 130,243 | 31 | 17 | 11 |
| *C. brevicaule* | YRN-4166 | SRR35152729 | 117,106 | 26 | 5 | 14 |
| *C. brevicaule* | AMM-4163 | SRR35152730 | 116,686 | 29 | 11 | 13 |
| *C. irumtiense* | ISG-4111 | SRR35152731 | 129,883 | 62 | 17 | 37 |
| *C. irumtiense* | MMY-4078 | SRR35152732 | 122,526 | 30 | 13 | 7 |
| **All six** | | | **747,234** | **215** | **66** | **106** |

**Composition-shuffled references: zero strict supporting fragments**. This is one limited synthetic-negative panel, NOT a calibrated genomewide false discovery rate or proof of molecular orthology.

**Main result:** the TR39766 DFR-like and TR41106 ANS-annotated source-*C. japonicum* fragments remain a reproducible leaf-RNA homology-compatible signal across ALL SIX public Ryukyu samples even after the 104,890-transcript source background competes for read alignments. Their presence is not an artefact of mapping solely to a reference panel of 26 selected genes.

This is a real within-genus comparative SOURCE-LEVEL positive pattern: public leaf R1 fragments homologous to previously colour-associated *Cirsium* source transcript candidates can be found in both taxon-level white and coloured lineages.

## Two independent technical analysis histories are NOT biological replicates

Previous published-source read screens in this project found:
- First 8 MiB, smaller 26-transcript source reference and partly different alignment cutoffs: **378** strict candidate fragments in **747,995** accepted reads.
- First 16 MiB, smaller 26-transcript source reference: **636** strict candidate fragments in **1,499,192** accepted reads.

Under both limited-source references, TR39766 and TR41106 were independently visible in 6/6 run identities. These 8 and 16 MiB datasets are **nested early bytes of the SAME six libraries**, NOT distinct plants or independent evolutionary colour changes. Read filtering and reference denominator differ, so do NOT compute ratios (215/378 or 215/636) as biological expression, retained fraction, gene loss or selective reduction.

The full-source v7 assay is an additional specificity control. Its 215 strict fragments from the 26 labelled source genes include 66 TR39766-like and 106 TR41106-like reads; the remaining 43 belong to other source candidates. This is NOT a direct expression-count comparison between coloured and white plants.

## The most defensible biological interpretation

The white-flowered taxon *C. brevicaule* does not have an observed complete *absence of all homologous young-leaf transcripts* in this tiny capped source prefix. Instead, both Ryukyu taxon-title groups share fragments homologous to 2020 *Cirsium* colour-associated transcript candidates in young leaves.

This makes a total absence of related expressed sequences a **poorly supported default assumption**. It DOES NOT prove that the original coloured ancestor carried a functionally intact floral DFR/ANS pathway, nor that white florets retain the same active biochemical enzymes. Alternate paralogs, alternative flavonoid pathway reactions, partial source transcripts, regulation outside the flower, assembly ambiguity, residual history, introgression or entirely different pigments remain open.

In particular the prior [DFR/ANR and ANS/FLS gene-family study](CIRSIUM_FLAVONOID_GENE_FAMILY_UNCERTAINTY_V5.md) found ambiguous functional family placement under NJ and FastTree for selected 2020 protein fragments. TR39766's alignment with DFR catalytic-reference positions and TR41106's published ANS similarity are NOT enzymology or gene-function evidence. Leaf transcript fragments are not developmental capacity C/F measurements.

## Next experiment in aza3: from shared leaf fragments to true function

1. Reconstruct and sequence/validate genome-contextual DFR/ANR and ANS/FLS orthogroups from the ACTUAL six Ryukyu source samples using independent paired-source/genomic information and ploidy-aware paralog controls; do not transfer the 2020 *C. japonicum* transcript IDs as though they are known Ryukyu loci.
2. Establish tissue- and stage-specific floral versus involucral and nonfloral flavonoid/anthocyanin chemistry and regulatory expression **in the same new, individually identified and appropriately permitted plants**. A preserved coding module and an off floral switch are separable testable states.
3. Test catalysis with appropriate DFR/ANR/ANS/FLS substrate/product and controls, then rooted allelic history and standing-variation/introgression/de novo alternatives before claiming re-expression.
4. Integrate validated head modules with orientation/protective phyllary/stickiness and direct viable-achene fitness in the existing [genetic LEGO experiment](ANTO_HIGHER_ORDER_LEGO_MEASUREMENT_PLAN_V2.md). Leaf morphology remains peripheral to the head-centred question.
5. Export only validated non-sensitive uncertainty-bearing source-gene module summaries to anto for broad evolutionary history accessibility tests; all Ryukyu molecular data live in aza3.

**No new field sampling was performed or authorized by this analysis. No ancestral white-to-coloured reacquisition, enzyme activity or ecological adaptation was experimentally established.**

## Execution sources

- [Six R1 official source/MD5 manifest](../data/anto_source_handoff/public_sra/ryukyu_six_R1_prefix_source_v6.json).
- [Full reference competitive mapping implementation](../analysis/anto_empirical/ryukyu_full_unigene_competition_v7.py).
- [Frozen six-source result](../data/anto_source_handoff/public_sra/ryukyu_full_104890_reference_competition_8MiB_v7.json).
- [Validation and adversarial tests](../analysis/validate_ryukyu_104890_full_reference_v7.py).
- Chang et al. 2026 DOI 10.1186/s12870-026-08097-6; Park et al. 2020 DOI 10.1021/acsomega.0c04001.
