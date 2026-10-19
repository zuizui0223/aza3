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

### Completed six-library public source result — 16 MiB first-R1 prefixes

The **corrected six-sample run succeeded** in [GitHub Actions run 39228212947](https://github.com/zuizui0223/aza3/actions/runs/39228212947). Each library contributed only the first **16 MiB of compressed R1**; no R2 or full source file was retrieved. Among **1,499,192 complete R1 reads**, **636 primary read fragments** satisfied the stated source-transcript mapping filters (MAPQ >=20, >=90 nt aligned, edit rate <=20%) against 26 *C. japonicum* 2020 author-annotated flavonoid-related transcript references. **0** strict matches were assigned to the shuffled reference controls.

| Published source taxon | Original specimen code | Complete R1 prefix reads | Strict 2020 source-transcript fragment matches |
|---|---|---:|---:|
| C. brevicaule | OKN-4295 | 261,137 | 85 |
| C. brevicaule | YRN-4166 | 235,060 | 92 |
| C. brevicaule | AMM-4163 | 234,011 | 103 |
| C. irumtiense | YGN-4296 | 262,388 | 119 |
| C. irumtiense | ISG-4111 | 260,579 | 149 |
| C. irumtiense | MMY-4078 | 246,017 | 88 |
| **Total** | **Six source libraries** | **1,499,192** | **636** |

**18 different 2020 source transcript IDs** attracted one or more strict matches across the six original young-leaf libraries. Two are notable because they had matching fragments in **all six**:

- TR39766|c0_g1_i2: **161** matched read fragments total. Its 2020 published protein resembles DFR, but DFR vs ANR catalytic family assignment remains unresolved under NJ/FastTree and the associated site-bootstrap CI.
- TR41106|c0_g1_i5: **270** matched fragments total. The 2020 author BLAST annotation matches *Saussurea medusa* ANS, but catalytic ANS specificity is **not** established by a BLAST label or these short read mappings.

**Secondary bounded method:** [GitHub Actions run 39228426774](https://github.com/zuizui0223/aza3/actions/runs/39228426774) ran the first **8 MiB** of R1 per identical source library with a separately implemented (MAPQ >=20, >=85 aligned read bases, identity >=85%) mapping filter. It recovered **747,995 complete reads, 378 strict source-transcript fragments and zero strict shuffled decoy fragments**. TR39766 and TR41106 also had strict matches in all six libraries. This is a nested read-prefix / methodological sensitivity analysis, **not another six independent biological replicate experiment**; the different counts and thresholds must not be pooled.

**Scientifically permitted positive result:** in two lineages with contrasting *published taxon-level floral colour*, real 2026 young-leaf RNA contains sequences similar to a shared set of source-annotated flavonoid-related transcripts. This offers an empirical starting point for anto's possibility that molecular resources persist outside the visibly pigmented floral state. However **even a white-flowered taxon with shared leaf transcripts has not been shown to retain a catalytically functional flower-pigment pathway, express one in its florets, or re-evolve colour**. At the level of the six source individuals, floral colour was not individually measured. The full original transcriptome was NOT used as competing alignment reference: MAPQ uniqueness within the selected 26 source sequences and 26 shuffled nulls does not resolve all paralogy or contamination.

Two immutable, source-audited result records are now stored in aza3:
- [Six R1 libraries, 16 MiB primary source run](../data/anto_source_handoff/public_sra/ryukyu_six_source_leaf_R1_mapping_16MiB_v6.json)
- [The same six libraries, 8 MiB alternative mapping sensitivity](../data/anto_source_handoff/public_sra/ryukyu_six_source_leaf_R1_mapping_8MiB_alternate_v6.json)

The [cross-run biological-evidence validator](../analysis/validate_ryukyu_real_leaf_prefix_evidence_v6.py) requires exact sample-to-BioSample mapping, original read totals and read-to-transcript count conservation, disallows any reinterpreted floral pigment event, and refuses to count the nested prefix samples as independent lineages.

### Third, stricter within-read validation — two exact 31-mers

The independent code route [GitHub Actions run 39229918762](https://github.com/zuizui0223/aza3/actions/runs/39229918762) analysed the same six R1 files using 8 MiB from the beginning of each file. Its different screen required **two nonoverlapping exact 31-base matches within the same quality-filtered R1 read**, unique to the same original 2020 *C. japonicum* source assembly group within the preselected 26-transcript panel. It also used **dinucleotide-composition-preserving shuffled decoys**, a stronger but still artificial negative control.

- **747,995** quality-validated source R1 prefix reads from the same six samples;
- **268** fragments pass the double-exact-31-mer source-group screen, while **zero** pass the dinucleotide-shuffled decoy screen;
- the DFR-like TR39766 source assembly group has **70** matching reads, and the 2020 BLAST ANS-annotated TR41106 group has **90**, with **at least one match in each of the same six biological libraries**;
- prefix byte ranges, BioSample IDs and SHA256 fingerprints are preserved separately for each public R1 source.

[Complete source-anchored two-31-mer result records](../data/anto_source_handoff/public_sra/ryukyu_six_R1_two31mer_source_group_results_v6.json) are retained in aza3. These matches are more restrictive **within this candidate reference panel**. They do NOT establish unique origin against all ~104,890 transcripts in the original *C. japonicum* assembly, do not prove paralogy/orthology, and do not prove catalytic activity or flower-stage anthocyanin expression. The third pipeline is **not independent biological replication**: all analyses repeatedly use the same six first-in-file R1 libraries.

The three read-level pipelines jointly justify the *limited positive observation* that flavonoid-related **source-assembly-like RNA sequence fragments** are found in young leaves of both taxa. The current highest-level biological answer remains **LEAF_HOMOLOGY_DETECTED / FLORAL_REEXPRESSION_UNTESTED**.

### Source parser issue found and resolved

An earlier six-sample execution failed **in the original FASTQ parsing rule**, not because any gene was absent: real *C. irumtiense* ISG-4111 SRR35152731 contains valid short read/quality records of **54–73 nt**. The earlier parser incorrectly treated read lengths below the analytic 75-nt mapping threshold as invalid FASTQ. The corrected method retains short valid FASTQ records for QC and excludes them from the predeclared alignment filter. The fix was followed by the **successful six-sample run above**, independently checked by CI.

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

1. **Completed:** source retrieval of bounded R1 prefixes for all six original libraries, with complete read counts and zero strict shuffled-decoy matches. No full-library read coverage established.
2. **Completed as a selected-reference screen only:** TR39766 and TR41106 had candidate read fragment matches in all six libraries; do not label an enzyme orthologue or functional pathway from selected-reference mapping.
3. **Next discriminator:** compare the same bounded public R1 reads against the larger, independently curated Cirsium transcriptome and/or closely related ANR/FLS gene family references; selected-panel MAPQ does not establish unique gene identity. A broader reference analysis may change the current per-source allocations.
4. For a biological pathway conclusion, independently recover taxon/population-specific genomic loci, functional family identity and matched flower-stage RNA/chemistry; this is the downstream experimental aza3 plan.
5. Keep anto's general dispersal/retained-capacity theory separate and no Cirsium sample data in anto.

## Reproduce

Requires Python 3 standard library plus external minimap2 and samtools. Source read retrieval uses the official public ENA FASTQ endpoint and a prespecified bandwidth cap.

    python analysis/test_ryukyu_ena_prefix_homology_v6.py -v
    python analysis/anto_empirical/ryukyu_ena_prefix_homology_v6.py \
        --runs SRR35152730 --bytes-per-run 8388608 --output results/ryukyu_one_run_prefix_v6.json

After the first run is validated, run all six by omitting --runs and specifying a new fixed data budget. Neither download is a probabilistic sample of whole RNA expression.

**This is a real sequence-homology screening workflow, not a published ecological result or an established ancestral re-expression.**
