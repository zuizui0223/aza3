# aza3 v6 — source-capped direct Ryukyu young-leaf RNA pilot

**As of 19 October 2026:** This is a first test designed to examine actual public 2026 *Cirsium brevicaule / irumtiense* reads rather than only the 2020 *C. japonicum* candidate-reference transcriptome.

## Why change the input data now?
The Park 2020 *C. japonicum* 26 DFR/ANS-like transcript candidates have ambiguous family attribution even after independent MAFFT+NJ and WAG FastTree analyses. Those are reference candidates, **not** actual alleles from the Ryukyu 2026 samples. Repeated re-ranking does not identify whether either Ryukyu lineage carries a matched transcript fragment.

The original Chang et al. (2026) **PRJNA1311153** project offers real young-leaf RNA reads from three *C. brevicaule* and three *C. irumtiense* specimens. Full paired FASTQs for the six are >10 GB for R1 alone and ~20 GB for both mates. We cap this pilot to **8 MiB of the first mate (R1) compressed gzip prefix per sample**, beginning with a single registered sample. No full reads, R2 or new field material are acquired.

Source manifest (exact original SRR/SAMN and specimen code, from official ENA read-run source):
[data/anto_source_handoff/public_sra/ryukyu_six_R1_prefix_source_v6.json](../data/anto_source_handoff/public_sra/ryukyu_six_R1_prefix_source_v6.json)

| Registered species | Source accession | Source specimen code | Source tissue |
|---|---|---|---|
| *C. brevicaule* | SRR35152725 | OKN-4295 | young leaf |
| *C. brevicaule* | SRR35152729 | YRN-4166 | young leaf |
| *C. brevicaule* | SRR35152730 | AMM-4163 | young leaf |
| *C. irumtiense* | SRR35152724 | YGN-4296 | young leaf |
| *C. irumtiense* | SRR35152731 | ISG-4111 | young leaf |
| *C. irumtiense* | SRR35152732 | MMY-4078 | young leaf |

**Actual tissue disclaimer:** These are registered young leaves, not floral transcriptomes. We do not infer individual flower colour from the source taxon name, and do not conclude flower-stage pigment expression from leaf mapped fragments.

## First actually executed source result — *C. brevicaule* AMM-4163

The first real public young-leaf-RNA test **succeeded** in [GitHub Actions run 39219960099](https://github.com/zuizui0223/aza3/actions/runs/39219960099), processing exactly 8 MiB of compressed **R1** from SAMN50798032 / SRR35152730. It recovered **116,686 complete high-length reads** and **53 high-confidence matched R1 fragments** distributed among ten source-2020 candidate transcripts. **Zero** reads passed the corresponding source-shuffled-negative reference mapping filter.

| Published Park 2020 candidate transcript ID | Source-AMM-4163 matching R1 fragments |
|---|---:|
| TR39766\|c0_g1_i2 (previous DFR-like source alignment) | **14** |
| TR41106\|c0_g1_i5 (source ANS-related annotation, family unresolved) | 19 |
| TR52731\|c6_g7_i2 | 9 |
| All other seven individual matched source transcript IDs | 11 |
| **Total filtered source-transcript R1 fragments** | **53** |

This is the first **direct 2026 white-taxon C. brevicaule young-leaf read-to-source-Cirsium homolog candidate link** in the project. It neither proves a whole gene is intact nor that the homolog is actually catalytic DFR or ANS, nor that the source specimen has functional anthocyanin capability in the flower. Furthermore, we aligned against a **selected 26-transcript panel**, not all 104,890 original Park 2020 transcripts; MAPQ uniqueness is only within that limited comparison database. Later full-transcriptome reference competition is required before treating any count as transcript-specific identification.

Frozen, sourced experimental report: [AMM-4163 8 MiB positive mapping record](../data/anto_source_handoff/public_sra/ryukyu_AMM4163_first_R1_8MiB_fragment_screen_v6.json). No raw read IDs, DNA strings or sensitive collection localities were committed.

### Input-quality fault found when scaling to all six specimens

An early six-sample expansion **failed in the FASTQ record parser**, not because any gene was absent: real *C. irumtiense* ISG-4111 SRR35152731 contains reads of **54–73 nucleotides** whose quality strings have the correct same lengths. Initial strict parser improperly considered any record shorter than the analytic 75-nt alignment threshold to be structurally malformed. The corrected method counts these as **valid-but-filtered-short**, while retaining hard failure for invalid headers, absent  markers and mismatched sequence/quality lengths. A synthetic regression reproduces short 54–74nt reads. No claim about RNA library structural corruption is supported by the observed source.

Six-source quantitative source mapping is **PENDING SUCCESSFUL REPLAY** until an independent CI run completes under the revised parser. Even after successful replay, prefix-only read counts are not expression estimates and **no non-detection can be interpreted as absence or loss**.

## Pilot source and quality safeguards

The [Python-only bounded prefix downloader and FASTQ parser](../analysis/anto_empirical/ryukyu_ena_prefix_homology_v6.py) uses strict HTTPS and an exact official ENA source allowlist; requires HTTP 206 Content-Range before reading any compressed data; refuses unbounded 200 responses; caps compressed input per R1 (8 MiB initially; maximum 32 MiB); verifies specimen/sample/taxon accession against the frozen 25-run manifest; and recovers only complete FASTQ four-line records from an inevitably truncated gzip stream.

The original provider MD5 is recorded for provenance, **not** misleadingly recomputed over the truncated prefix. A gzip stream which lacks an intact footer is expected; only complete read records are considered. The first reads in a file can be sorted by acquisition/technical factors, so the prefix is **not a random, representative sample of biological RNA**.

Mapped evidence is further constrained:

- Use a published, checksum-matched *C. japonicum* 2020 nucleotide FASTA with 26 candidate fragments; they are **all still functional gene-family candidates** and not full-locus orthology.
- Add 26 independently shuffled same-composition negative-decoy nucleotide sequences as a computational nonspecific-alignment control.
- Short-read minimap2 and SAM-level filtering (MAPQ ≥20; ≥90 aligned nucleotides; mismatch/indel edit rate ≤20%) under a fixed reference/read parser. Report actual unique alignment read names and biological transcript IDs.
- Retain all ambiguous/low-quality alignments in separate QC counts rather than force a molecular gain/loss call.

### What a positive or negative result means

**A read mapping from a 2026 Ryukyu young-leaf prefix** supports the limited observation that some reads are similar to a published 2020 *C. japonicum* transcript candidate within the admitted alignment filters. It does not establish that the 2026 individual has a complete coding locus, retains a catalytic anthocyanin enzyme, expresses it in its florets, contains a pigment, or has re-evolved it.

**No mapped read** is even more limited: NOT_DETECTED_IN_BIASED_PREFIX. It is not a gene deletion, absent expressed transcript in the entire RNA library, developmental incapacity, white-flower pathway loss or a comparison of expression rates between species.

If ENA does not support correctly bounded HTTP 206 retrieval, the pilot stops before download and records a source-access blocker; it must not fetch 20 GB as a fallback.

## Criteria for next research step

1. Confirm source retrieval of actual R1 prefix for at least one sample; preserve read count/verified R1 source and negative-decoy background.
2. Assess whether any *source-frozen* 2020 Cirsium homologous transcript reference has ≥2 uniquely mapping R1 fragments; require full case-by-case reference-traceability before calling a candidate.
3. Only then, potentially run all six at equal bounded prefix sizes and compare **detection**, with no inference of differential expression, gene absence or any fixed loss/regain route.
4. For a biological pathway conclusion, independently recover taxon/population-specific genomic loci, functional family identity and matched flower-stage RNA/chemistry; this is the downstream experimental aza3 plan.
5. Keep anto's general dispersal/retained-capacity theory separate and no Cirsium sample data in anto.

## Reproduce

Requires Python 3 standard library plus external minimap2 and samtools. Source read retrieval uses the official public ENA FASTQ endpoint and a prespecified bandwidth cap.

    python analysis/test_ryukyu_ena_prefix_homology_v6.py -v
    python analysis/anto_empirical/ryukyu_ena_prefix_homology_v6.py \
        --runs SRR35152730 --bytes-per-run 8388608 --output results/ryukyu_one_run_prefix_v6.json

After the first run is validated, run all six by omitting --runs and specifying a new fixed data budget. Neither download is a probabilistic sample of whole RNA expression.

**This is a real sequence-homology screening workflow, not a published ecological result or an established ancestral re-expression.**
