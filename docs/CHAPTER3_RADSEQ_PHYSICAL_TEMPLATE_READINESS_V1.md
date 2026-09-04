# Chapter 3 RAD-seq physical-template readiness v1

Status date: 2026-09-02

This layer separates **a public sequence/reference exists** from **a physical tissue or DNA template is actually available for the pilot**. NCBI BioSamples, published target-capture accessions, herbarium metadata and literature reports are not treated as DNA in hand.

## Current result

Repository search found no admitted inventory demonstrating physical possession of RAD_A01–RAD_A08 material. Therefore all 16 planned biological templates start as `UNKNOWN_NOT_VERIFIED` in `data/planning/chapter3_radseq_physical_template_inventory_v1.csv`.

This is not a statement that the material does not exist. It means the repository cannot presently verify a tube, tissue packet, DNA extract or qualified external material transfer that can be assigned to a pilot slot.

## Inventory structure

Each of the eight anchors has two physical slots:

- `TEMPLATE_1`: required before Stage A can open;
- `TEMPLATE_2`: required before Stage B can open and must be an independent biological individual, not an aliquot or library duplicate of TEMPLATE_1.

Every row must eventually carry:

- immutable physical sample ID;
- material type and current holder/source;
- linked voucher or diagnostic image;
- identity determination;
- sample-level cytotype or the genome-size role needed by the anchor;
- tissue-use/transfer/collection authorization status;
- DNA extraction status;
- fluorometric concentration and total DNA amount;
- purity and fragment-integrity state;
- explicit Stage-A/Stage-B eligibility.

Sensitive exact localities and permit documents do not belong in this public ledger.

## Candidate acquisition routes

The routes are candidates, not evidence of possession.

- `RAD_A01 C. alpicola` and `RAD_A02 C. aomorense`: distributed collaborator acquisition is acceptable for assay-only pilot material after identity and cytotype confirmation.
- `RAD_A03 C. dipsacolepis` and `RAD_A04 C. lineare`: investigator/close-collaborator targeted acquisition from currently verified populations because these are focal P02 systems and conservation/identity gates matter.
- `RAD_A05 C. nipponicum` and `RAD_A06 C. sieboldii`: new wild-Japan material is preferred. Existing cultivated or public sequence references are not silently promoted into pilot DNA.
- `RAD_A07 C. brevicaule` and `RAD_A08 C. irumtiense`: assay-only qualified fresh/silica material may be used for the enzyme feasibility test, or material can be nested inside the preregistered M01 discovery acquisition. Published genome-size/sequence results do not establish physical possession.

## Stage-A physical gate

Stage A remains unauthorized until all eight `TEMPLATE_1` rows satisfy all of the following:

1. `physical_material_status = VERIFIED_IN_HAND_OR_TRANSFER_CONFIRMED`;
2. a nonempty immutable `physical_sample_id`;
3. voucher/diagnostic identity is verified;
4. the required cytotype/genome-size anchor role is qualified at the sample level or explicitly accepted for assay-only use;
5. tissue-use authorization/transfer status permits the assay;
6. DNA has been extracted or a qualified DNA extract is received;
7. fluorometric concentration and total DNA are recorded;
8. DNA purity and integrity pass the lab-frozen pilot QC;
9. `stage_a_eligible = true`.

Having seven of eight anchors is not sufficient for the preregistered eight-anchor comparison. A missing anchor triggers acquisition or a formally versioned redesign; it is not replaced because another taxon is easier to obtain.

## Stage-B physical gate

Stage B additionally requires the eight `TEMPLATE_2` rows to pass identity and DNA QC. TEMPLATE_2 must be biologically independent of TEMPLATE_1 for the same anchor. Technical repeats in the Stage-B library allocation are generated from TEMPLATE_1 DNA and are a different concept from TEMPLATE_2.

## Immediate acquisition order

The practical order is:

1. inventory any already-existing lab/freezer material without assuming it exists;
2. verify provenance and voucher links before extracting DNA;
3. prioritize RAD_A03/RAD_A04/RAD_A06 because they overlap focal Chapter 3 systems;
4. acquire RAD_A07/RAD_A08 together with M01 discovery material where logistics allow;
5. use collaborator assay material for RAD_A01/RAD_A02 if identity/cytotype documentation is adequate;
6. secure wild-Japan replacement material for RAD_A05 where possible;
7. run DNA QC only after a row has a real physical ID.

## Claim boundary

This inventory does not authorize collection, material transfer, DNA extraction, Stage A, Stage B or production RAD. It only turns an abstract 16-template pilot into a falsifiable physical-readiness checklist.
