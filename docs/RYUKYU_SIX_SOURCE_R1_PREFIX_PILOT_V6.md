# aza3 v6 — Bounded public RNA-seq prefix pilot for six Ryukyu thistles

**Protocol date:** 2026-10-19. **Evidence status:** A live result is admitted only after the Github Actions job downloads the bounded source bytes. A planned scan is not an observed gene detection.

## Objective and genuine source

The 2026 Chang et al. public project PRJNA1311153 holds three young-leaf RNA libraries attributed to Cirsium brevicaule and three attributed to Cirsium irumtiense. The accessions and experiment-title sample names have been checked against the frozen ENA 25-run manifest. All active Cirsium biological analyses belong to aza3; anto retains genus-independent dispersal/history-dependent accessibility theory.

| Submitted taxon | R1 run | Original source specimen code |
|---|---|---|
| C. irumtiense | SRR35152724 | YGN-4296 |
| C. brevicaule | SRR35152725 | OKN-4295 |
| C. brevicaule | SRR35152729 | YRN-4166 |
| C. brevicaule | SRR35152730 | AMM-4163 |
| C. irumtiense | SRR35152731 | ISG-4111 |
| C. irumtiense | SRR35152732 | MMY-4078 |

The six entire compressed R1 files together are ~9.66 GB. The pilot requests **the initial 8 MiB of each file, at most 48 MiB of original compressed data**. This is NOT a random or representative downsample of each sequencing library.

The full [six-source ENA R1/MD5/size manifest](../data/anto_source_handoff/public_sra/ryukyu_six_R1_prefix_source_v6.json) is the source of immutable run IDs. The original study's NovaSeq Methods wording and ENA's MiniSeq metadata remain unresolved, not proof of invalid data.

## Mandatory source and read limits

1. Only HTTPS from the exact ftp.sra.ebi.ac.uk public R1 path is allowed. Request byte positions 0 through 8,388,607; require HTTP status 206 and Content-Range matching the declared full R1 byte size. A server returning HTTP 200 is **rejected before data transfer** to prevent accidentally downloading 1.5+ GB.
2. Decompress only the initial gzip segment (up to 96 MiB decoded bytes), parse complete valid FASTQ records, discard incomplete tail, verify base/quality consistency. An incomplete gzip cannot be compared with a provider MD5 of the entire file.
3. Use all 26 Park (2020) Cirsium japonicum BLAST-annotated transcript sequences (23 source assembly groups), with original publication sequence checksums. They are from **another Cirsium species/study** and have NOT been certified as functional DFR or ANS enzyme orthologues.
4. Count candidate evidence only if the **same R1 sequence** contains at least TWO quality-filtered, nonoverlapping, exact canonical 31-base matches UNIQUE to the same 2020 assembled transcript group. A single seed or conflicting groups is nonconfirming. Low-quality positions are excluded.
5. Scan a matched **dinucleotide-preserving shuffled reference source control**, excluding any k-mers shared with the original 26 transcripts. This is only a limited null/background, not a genomewide specificity estimate or certified false-discovery rate.
6. Upload per-run statistics and SHA256 of the obtained bytes only; never store the raw FASTQ in Git, log sequences, or report personal/sensitive locality coordinates.

## Allowed interpretations

A genuine read hit means only: a PARTIAL registered **young-leaf** R1 read contains nucleotide sequence homologous/compatible with a previously assembled C. japonicum source transcript. This does not establish orthology, catalytic DFR/ANS functionality, complete coding regions, an allele in the focal individual, ancestral pathway age, pigment in its florets, or flower reactivation.

If no hit meets the strict two-seed requirement, the biological result is **NOT_DETECTED_IN_ORDERED_R1_PREFIX**. It is NEVER evidence of gene deletion, permanent white-flower incapacity, or a closed re-expression window. A library inaccessible through Range is **NOT_ASSESSED**.

Do not compare read hit counts among taxa as expression differences: original read order, R1-only observations, differing library input and reference divergence are uncontrolled.

## What becomes decisive after this pilot

When source access works, progress from gene-family k-mers to complete/source-controlled paired reads and broader DFR/ANR/ANS/FLS orthology/paralogy analysis. Then assay the SAME Cirsium plant for mature and pigmentation-onset floret RNA and pigment chemistry, phyllary and nonfloral tissues, genomic ancestry/cytotype and heritable function/segregation. Only a rooted, validated historical loss/retention genealogy can discriminate southern pigment regain from ancestral coloured retention or immigration/introgression.

The capitulum LEGO hypothesis requires independent genetically backed head orientation/protective envelope modules and viable-achene fitness outcomes. Observing a pigment-related leaf RNA read cannot substitute for these measurements.

## Commands

Run offline source tests:

    python analysis/test_ryukyu_R1_prefix_kmer_pilot_v6.py -v

Perform explicit capped network pilot:

    python analysis/anto_empirical/ryukyu_R1_prefix_kmer_pilot_v6.py --limit-bytes 8388608 --output /tmp/aza3_ryukyu_R1_prefix_v6.json

The GitHub Actions workflow separates offline tests and network access, preserving a source-status record even if the ENA server does not support bounded Range requests. Only six per-run records with status REAL_R1_PREFIX_SCANNED count as complete pilot execution.

Primary source: Chang et al. 2026 DOI 10.1186/s12870-026-08097-6, PRJNA1311153. Independent reference transcript source: Park et al. 2020 DOI 10.1021/acsomega.0c04001.
