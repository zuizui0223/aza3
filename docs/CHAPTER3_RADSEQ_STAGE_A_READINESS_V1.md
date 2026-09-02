# Chapter 3 RAD-seq Stage A readiness v1

Status date: 2026-09-02

This checklist opens no sequencing or collection authorization. It defines what must exist before the 24-library Stage-A enzyme screen can start.

## A. Anchor material

For each `RAD_A01`–`RAD_A08`:

- [ ] one Stage-A DNA template identified and linked to an immutable sample ID;
- [ ] one independent reserve/Stage-B template identified where possible;
- [ ] voucher or diagnostic image linked;
- [ ] taxon-concept identity reviewed;
- [ ] cytotype/genome-size role recorded;
- [ ] tissue-use authorization recorded where required;
- [ ] DNA source is fresh/rapid silica-dried material or explicitly qualified equivalent.

Stage-A assay templates do **not** count toward the 167/193 Chapter 3 biological sample size unless a later contract separately admits them as biological samples with all same-individual phenotype fields.

## B. DNA QC before digestion

For every Stage-A template:

- [ ] fluorometric concentration recorded;
- [ ] purity metrics recorded;
- [ ] fragment integrity checked;
- [ ] normalized input mass prepared;
- [ ] no whole-genome amplification;
- [ ] remaining DNA banked for Stage-B reproducibility or troubleshooting.

## C. Candidate protocol preparation

The candidate registry is frozen as:

- `RAD_C01_AVAII_MSPI` — AvaII + MspI;
- `RAD_C02_ECORI_MSPI` — EcoRI + MspI;
- `RAD_C03_PSTI_MSPI` — PstI + MspI.

Before wet-lab execution:

- [ ] enzyme compatibility and supplier reaction conditions checked;
- [ ] adapter overhang compatibility documented;
- [ ] indexing scheme checked for uniqueness and balanced base composition;
- [ ] provisional 400–700 bp size-selection workflow tested on noncritical DNA if needed;
- [ ] PCR cycle policy frozen for Stage A;
- [ ] paired-end sequencing configuration selected for the pilot run;
- [ ] negative/library blank control positions defined separately from biological library slots.

The 400–700 bp interval remains a **screening window**, not a production commitment.

## D. Stage-A library allocation

Use `data/planning/chapter3_radseq_pilot_library_allocation_v1.csv` exactly for the 24 Stage-A biological assay slots:

- 8 anchor systems;
- one Stage-A template per anchor;
- all three candidate protocols per template;
- 24 shallow libraries total.

Protocol identity must not be confounded with sequencing lane/run. If more than one lane is used, distribute all three candidates across lanes rather than assigning one protocol to one lane.

## E. Stage-A technical outputs required

Before candidate ranking, produce one frozen table containing for every library:

- raw read pairs;
- demultiplexed retained read pairs;
- adapter/off-window fraction;
- fragment-size summary;
- equal-depth usable locus count;
- median locus depth;
- high-depth tail metrics;
- 6/8, 7/8 and 8/8 anchor shared-locus summaries by protocol.

All protocol comparisons use equal-read downsampling.

## F. Advancement decision

A protocol cannot advance if it has an unexplained mandatory-anchor library failure or a protocol-specific loss of either M01 genome-size anchor.

Advance the best two technically viable candidates using the frozen technical metrics. Do **not** inspect P01/P02 trait-history placement, M01 colour or any topology that would make one protocol scientifically more convenient.

If fewer than two candidates are viable, stop and redesign enzyme/size-window candidates. Do not lower the biological QC gates.

## G. Stage B remains closed until Stage A decision is frozen

Only after Stage A is complete may `TOP1` and `TOP2` be substituted into the 40-slot Stage-B allocation.

Stage B then contains:

- 32 primary reproducibility libraries;
- 8 cross-batch technical repeats;
- genotype concordance gate >=0.95;
- core-locus recovery gate >=0.90.

Production RAD remains closed after Stage B unless the full end-to-end production-opening conditions pass.
