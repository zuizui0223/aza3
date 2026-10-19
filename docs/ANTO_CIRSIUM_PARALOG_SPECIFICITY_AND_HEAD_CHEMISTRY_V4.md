# Cirsium empirical v4 — Are DFR/ANS-like BLAST hits specific to pigment synthesis?

**Source and ownership:** aza3 is the canonical Cirsium empirical home. The generic history-dependent accessibility theory stays in anto. This is a real-source *sequence-family specificity screen*, not Ryukyu-specific genotype analysis, functional orthology, anthocyanin reacquisition, head LEGO genetics or adaptation.

## The question raised by v3
From 26 author-annotated *Cirsium japonicum* 2020 transcript sequences, only 5 passed the original two-protein similarity screen (4 DFR and 1 ANS, across 4 source assembly groups). However, DFR belongs to a broader short-chain reductase family, while ANS belongs to a broader 2-oxoglutarate/Fe(II)-dependent oxygenase family. The closest protein by sequence may catalyze a different step. Therefore an unqualified DFR/ANS BLAST annotation is not a pigment-retention assay.

## Seven accession-verified protein references (NCBI EFetch)

| Class in this experiment | Verified reference | Organism | Role |
|---|---|---|---|
| DFR | ABQ97018.1 | *Saussurea medusa* | Closely related Cardueae pigment-pathway enzyme reference |
| DFR | CAA78930.1 | *Gerbera hybrida* | Independent Asteraceae DFR reference |
| DFR | BAD95233.1 | *Arabidopsis thaliana* | DFR across order boundary |
| ANS/LDOX | AAS48200.1 | *Saussurea medusa* | Closely related Cardueae ANS |
| ANS/LDOX | NP_194019.1 | *Arabidopsis thaliana* | Independent ANS and experimentally characterized LDOX |
| ANR | NP_176365.1 | *Arabidopsis thaliana* | DFR-like **related-enzyme negative control**, BANYULS |
| FLS | NP_196481.1 | *Arabidopsis thaliana* | ANS-like **related-enzyme negative control**, flavonol synthase |

The names "negative control" here mean *a different annotated biochemical role*, not proof that the homolog cannot participate in floral pigmentation regulation or another metabolic pathway. The native Arabidopsis BAN/ANR protein functions as an anthocyanidin reductase (see [Arabidopsis gene BAN, NCBI](https://www.ncbi.nlm.nih.gov/gene/842469), and [Xie et al. 2003, Science](https://pubmed.ncbi.nlm.nih.gov/12532018/)). ANS and FLS have related oxygenase folds and different canonical biochemical products; [Arabidopsis LDOX](https://www.ncbi.nlm.nih.gov/gene/828387) and [FLS1](https://www.ncbi.nlm.nih.gov/gene/830765) are annotated examples. Also see the empirical overlap/structural similarities of related flavonoid reductases ([2023 analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC10531346/)).

[Download-free curated FASTA in aza3](../data/anto_source_handoff/molecular/extended_DFR_ANS_ANR_FLS_reference_proteins_v4.fasta), with [source accessions and protein lengths](../data/anto_source_handoff/molecular/extended_reference_protein_provenance_v4.json).

## Experimental screening, not a certified gene tree

For each of the 26 existing *C. japonicum* transcript sequences:
1. Translate all six reading frames and retain stop-free peptides ≥35 amino acids, without claiming full CDS integrity.
2. Perform local BLOSUM62 protein alignment to 7 reference proteins, keeping per-family strongest alignment, identity, length and reference coverage. Alignment normalized by reference protein length is **only a similarity metric**.
3. Require at least 100 aligned amino acids, 30% exact identity among aligned positions, and 50% reference coverage for a **provisional strong similarity screen**. These thresholds are NOT functional or phylogenetic calibration.
4. Compare each best DFR to ANR and each best ANS to FLS, retaining ambiguous/negative-family results. Do not force all 26 proteins into a colour-making category.
5. Freeze individual sequence IDs, checksums and contradictory annotations; next test is homologous gene-tree/coalescent paralogy analysis, biochemical activity and same-individual flower-stage expression rather than another opportunistic BLAST.

The prior five-reference-positive transcripts are an *initial set*, not a guarantee that all five pass a stricter related-family control.

## Reproducible result: reference sensitivity reduces the robust shortlist

The genuine *C. japonicum* 2020 public sequences were analysed with all seven positive/related-family negative reference proteins. **The previous five strong candidate transcript records remained positive-family nearest matches under the full panel**, but they did not all survive reference-taxon omission.

| Source transcript | Full seven-reference best family | Omit Saussurea positive refs | Omit Arabidopsis positive refs | Across-scenario stability |
|---|---|---|---|---|
| TR39766|c0_g1_i2 | DFR | DFR | DFR | **Stable** |
| TR35526|c0_g1_i1 | DFR | DFR | DFR | **Stable, but borderline ~30.6% amino-acid identity** |
| TR38875|c0_g4_i5 | DFR | **ANR** | DFR | Reference-sensitive |
| TR38875|c0_g4_i7 | DFR | **ANR** | DFR | Reference-sensitive |
| TR58302|c0_g1_i1 | ANS | ANS | **FLS** | Reference-sensitive |

Precisely **2/5** retained the original positive-family best match and the minimum similarity screen across all tested reference removal scenarios. **3/5** remained positive when the Saussurea references alone were excluded. Both TR38875 entries belong to one source assembly group; they cannot be counted as independent genes or independent evolution events.

This is a genuine **reference robustness result**, not evidence that the others actually encode ANR or FLS enzymes. Their classification is *reference-dependent*. It demonstrates why generic DFR/ANS protein-family names should not be used to claim that the full floral anthocyanin toolkit remains functional or was re-evolved.

For targeted gene-family functional follow-up, prioritize TR39766 for the stronger, longer DFR-like alignment; inspect TR35526 as a weaker but robust DFR-like candidate, and explicitly **retain rather than discard** TR38875 and TR58302 for gene-tree and catalytic-family ambiguity resolution. The first two are not functionally proven either. An existing 2020 root:leaf CjDFR expression contrast is **not linked to these exact source transcript IDs** and should not be silently attributed to them.

Machine-readable results with all per-scenario reference scores and the independent source provenance are in [frozen seven-reference family report](../data/anto_source_handoff/molecular/cirsium_2020_seven_reference_family_control_v4.json) and [reference omission sensitivity](../data/anto_source_handoff/molecular/cirsium_2020_reference_omission_sensitivity_v4.json). The [jackknife runner](../analysis/run_anto_cirsium_reference_jackknife_v4.py) rechecks the results against its own saved evidence fingerprint, and fails if source-supported conclusions silently change.

## New Cirsium and genetic LEGO interpretation

This screen could support competing mechanisms, but cannot decide among them alone:
- **Retained DFR/ANS route:** expression capacity in a white floret is maintained but floral regulation switches off. Requires functional structural gene and floral regulation assays.
- **Alternative metabolic partition:** ANR-related proanthocyanidin synthesis or FLS-related flavonols might use overlapping precursors without accumulating visible anthocyanin. This is an **unverified testable candidate explanation**, not measured in white Ryukyu specimens.
- **Degraded pathway:** genuine coding/regulatory disruption in the white floral lineage, to be validated with genomic DNA, ploidy and developmental tissue controls.
- **Imported/introgressed allele:** local ancestry and pigment locus genealogies disagree; a coloured plant may appear with no intrinsic resident white-lineage reacquisition.

Within aza3's capitulum LEGO design, the functional units should distinguish **floret pigment amount, phyllary chemistry/structure, organ-specific flux allocation and regulatory activation**. A colour/defense trade-off requires a *directly measured allocation and ecological consequence*, not merely different enzyme names. Do not let leaf morphology replace the head-centred experiment. Same-individual florets and phyllaries, VIS/UV and LC-MS/HPLC plus direct relevant products and viable achenes are the correct ecological endpoint.

## What would qualify as a positive biological result?

A strong sequence hit is only a research input. A mechanistic claim about *Ryukyu Cirsium* requires an independently validated candidate coding/functional module actually present in those plants, organ/stage-specific RNA and pigment chemistry, ancestry/genetic origin, and differentiation between ancestral retention, gene flow, immigration and novel regulatory reconstruction. Current 2026 Ryukyu transcriptomes are young **leaf** material, not flower-specific expression.

## Reproduce

    python -m pip install "biopython>=1.80,<2"
    python analysis/test_anto_cirsium_paralog_panel_v4.py -v
    python analysis/run_anto_cirsium_paralog_panel_v4.py --out results/cirsium_paralog_controls_v4.json

No new field material is collected and no regulatory or evolutionary re-expression is inferred by this sequence comparison.
