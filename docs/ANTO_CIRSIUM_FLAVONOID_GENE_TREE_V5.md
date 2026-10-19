# Cirsium v5 — DFR/ANR and ANS/FLS source protein family trees

**Status:** exploratory source-backed molecular family analysis. It does NOT establish Ryukyu flower-colour reacquisition or functional orthology.

## Why a protein gene-family tree is necessary

Published 2020 Cirsium japonicum transcript data yielded 26 author BLAST-annotated DFR/ANS-related sequences. Five passed the original local protein-alignment thresholds, across four source assembly groups. Only two DFR-like transcripts retained stable nearest-family similarity under reference-taxon omissions (TR35526|c0_g1_i1 and TR39766|c0_g1_i2).

Related flavonoid enzyme families overlap in sequence, structure and sometimes substrate/product selectivity. DFR/ANR/FNR family relationships: https://pmc.ncbi.nlm.nih.gov/articles/PMC10531346/; ANS/FLS related Fe/2OG oxygenases: https://pmc.ncbi.nlm.nih.gov/articles/PMC3907857/. A BLAST protein family name does not prove pigment function.

## Biological inputs (all verifiably published)

- Park 2020 Cirsium japonicum transcript sequences, NOT the 2026 C. brevicaule or C. irumtiense young-leaf experiments.
- Candidates: TR35526|c0_g1_i1, TR38875|c0_g4_i5, TR38875|c0_g4_i7, TR39766|c0_g1_i2, TR58302|c0_g1_i1. The two TR38875 entries share an assembly group and are treated as alternative transcript choices, NOT independent genetic origins.
- Nine independent NCBI-accessioned reference proteins: three DFR (Saussurea, Gerbera, Arabidopsis), two ANR (Arabidopsis, Vitis), two ANS (Saussurea, Arabidopsis), two FLS (Arabidopsis, Vitis). See data/anto_source_handoff/molecular/DFR_ANR_ANS_FLS_reference_provenance_v5.json.
- Two distinct enzyme superfamilies: SDR reductase DFR+ANR and Fe/2OG oxygenase ANS+FLS. **Never infer a single mixed-enzyme evolutionary tree.**

## Three preregistered exploratory gene tree tests

1. DFR_ANR_primary: three distinct Cirsium assembly groups (including TR38875 isoform i5), plus 3 DFR and 2 ANR source-verified reference proteins.
2. DFR_ANR_alt_isoform: same samples except substitute TR38875 isoform i7 for i5; this is a sensitivity analysis, not a biological replicate.
3. ANS_FLS: source TR58302 and 2 ANS + 2 FLS reference proteins.

Original source coding orientation is unknown. Stop-delimited protein fragments from all six reading frames are compared to **both related enzyme families** before selecting a segment. This has family-reference bias and does not reconstruct a complete genomic coding locus.

MAFFT --auto aligns each family separately. Unrooted neighbor joining uses shared-nongap amino acid p-distances, requiring at least 60 shared positions per pair. Site bootstrap resamples alignment columns 200 times with a fixed seed. Results report the fraction of candidate placements on the positive, related-enzyme, intermediate, or interleaved reference sides, plus alignment overlap and unrooted Newick.

This is a preliminary **gene-family placement screen**. Neighbor joining, sparse taxon sampling, fragment alignment and site-bootstrap support do not establish single-copy orthology, biochemical substrate preference, enzyme activity, paralog dosage, allelic age, epigenetic floral suppression, colour re-expression or adaptive function.

## Interpretation gates

- Reference families separate and a candidate is consistently attached to the DFR or ANS side: advance for deeper orthology/gene-tree and biochemical analysis, **not** a conclusion of functional floral competence.
- Positive and related-enzyme references interleave, candidate is on an intermediate stem, or alignments lack shared residues: family classification NOT_IDENTIFIABLE. Do not force a gain/loss label or infer absence of the corresponding gene.
- Same-id fragments / isoforms from a single assembly group: never count as two independent inherited resource variants.

For the biological next step, compare genuine Ryukyu genomes with floret and phyllary expression and pigment chemistry in the same individual, rooted allele genealogies, cytotype and introgression controls, then controlled gain/loss/rescue experiments where authorized. The six currently registered Ryukyu libraries were made from young leaves and do not measure floral re-expression.

## Reproduce and source ownership

Run with MAFFT and Biopython:

    sudo apt-get install mafft
    python -m pip install "biopython>=1.80,<2"
    python analysis/test_anto_cirsium_gene_tree_v5.py -v
    python analysis/anto_empirical/cirsium_flavonoid_gene_tree_v5.py --bootstraps 200 --output results/cirsium_gene_tree_v5.json --alignment-directory results/cirsium_gene_tree_alignments_v5

Source/analysis owner: aza3. Cross-angiosperm mathematical accessibility and plant dispersal theory: anto.
