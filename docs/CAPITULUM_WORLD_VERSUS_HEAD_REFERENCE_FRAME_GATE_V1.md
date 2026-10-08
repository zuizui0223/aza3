# アザミ頭花の向きは「重力座標の昆虫」を選別するのか：頭花座標と世界座標の識別可能性

2026-10-08. Supplemental PR #44. Scope = **Cirsium capitulum only** (floral face, true involucre/phyllary armature and spines, colour/UV, stigma/pollen, phenological stage, antagonist guild). No leaves, no approved new intervention, no experimental result, and no change to frozen Azami/EAzami lineages.

## New null competition, beyond an ordinary attraction–defence story

A hanging head rotates the *same head-relative structures* against a gravity-fixed world. Its role cannot be identified from one image angle and a visitor category.

- **W: world-/gravity-anchored approach**. Given insect taxon and arrival opportunity, arthropods tend to approach from similar real-world bearings (e.g., above or from a stem) despite different head orientation. World-bearing distribution is more stable across orientation than head-relative bearing; the set of phyllaries actually encountered changes. Prediction: angle × real access gap × independently identified insect approach guild changes first-contact failure and later oviposition (not necessarily total visit number).
- **H: head-anchored approach**. Arthropods adjust to the disc normal or head wall and use similar head-relative bearings. Real-world bearing rotates with the head; angle alone need not change which anatomical portals are touched. Under H, the observed image-angle × head-elongation correlation need not imply a selective spine barrier.
- **A: abiotic/visibility alternative**. Rain, pollen wetting, light/UV visibility and flower display change arrivals, pollen function or head movement directly, without changing route-specific mechanical entry. Weather can be a common cause of measured head angle and insect activity; no behavioural frame is resolved from raw correlation.
- **P: phenology/portal-switch alternative**. A bud's floret portal is closed or inaccessible, while anthesis may expose another route. A seed feeder blocked at the wall may switch; only continuous within-episode time-stamped video discriminates first failure from final success. Taxon-specific host/age effects are essential, because insect species do not share one immutable capitulum entry route.
- **C: community context alternative**. Co-flowering *Cirsium* hosts, plant display, head rank and parasitoid availability change arrival/oviposition pressure independently of focal spine geometry.

These are *competing hypotheses*, not ranked causal discoveries.

## Mathematically distinct reference frames

At a standardized **pre-contact** sphere (radius = 1.5 measured head diameters from head centre), measure world-frame unit vectors:
- `n`: head floral-disc normal directed outward from the receptacle, not vertical in a cropped image;
- `r`: head-centre-to-insect approach ray at first crossing of the sphere;
- `z=(0,0,+1)`: calibrated world vertical, opposite gravity.

Compute `world_up = r·z` and `head_front = r·n`. This pair measures **relative approach direction only**, not contact location, insect intent, whether a spine was touched or pollen/seed fitness. A calibrated synchronized multiview/stereo view plus a gravity reference and a botanical head axis are required. A single uncalibrated Raspberry Pi side image cannot recover both 3D vectors. If both n and r are copied or derived from the same segmented insect silhouette, the experiment is circular.

**Synthetic counterexamples** with two orientation bins (upright `n=z`; lateral `n=(1,0,0)`):
- World-fixed insect ray `r=z`: mean `world_up=1` in both bins, `head_front` changes 1 to 0.
- Head-fixed ray `r=n`: mean `head_front=1` in both bins, `world_up` changes 1 to 0.
- If only upright heads are sampled, `r·n = r·z` exactly: both explanations predict the same observed directional score, so **non-identifiable even with infinitely many contacts**. This is a rank/positivity boundary, not a hypothesis test.

The above are **geometry examples only**. No Cirsium data in PR #44 currently have 3D pre-entry ray measurements.

## Real analysis gate

Additional **empty** CSV: `data/intake/capitulum_world_head_3d_approach_v1.csv`. Reproducible strict audit: `analysis/check_capitulum_reference_frame_gate_v1.py`. A measured row has source evidence for gravity calibration, independent anatomical axis and insect ray, verified before portal access and without using oviposition/pollination outcome as its entry label.

A *descriptive, preregistration-candidate* frame comparison is allowed only within an exact stratum of **population × taxon × head stage × head rank × independently identified pre-entry guild**, across at least two natural orientation bins with **three different plants in each**. Multiple contacts per plant are averaged, not treated as independent plants. Report mean score differences in both world and head frames, without selecting whichever wins after looking at fitness. If there is no frame variation or no cross-frame geometry, print `NO_REAL_STEREO_GEOMETRY` or `GEOMETRY_PRESENT_BUT_FRAME_NOT_COMPARABLE`, **never zero effect**.

These minimums are *logical overlap checks*, **not** sufficient power or causal identification. Plant/head co-variation, weather, first-arrival selection, sexual function, and camera occlusion remain. Future causal inference requires exogenous head reorientation with botanical controls/permissions, or tightly matched natural variation and independent replication, as well as both early contact and full head fitness. Raw `azami` photo `presentation_angle` must never be silently relabelled calibrated head elevation.

Keep two different observation tiers:
1. One continuous, reviewed Pi video can describe **coarse visible above/side/below and contact-zone sequence** with the existing two ledgers, provided false-negative and viewing coverage are audited. It **cannot** prove world/head-frame preference or depth-specific physical barrier effects.
2. A source-calibrated **two-view/stereo pilot** can test the 3D reference-frame prediction if real natural orientation contrast and the same guild/stage are available. This script requires a fixed approach shell to prevent cherry-picking the frame classification point.

Read resulting access failures with actual phyllary/spine anatomy, and only then test confirmed pollen deposition, seed feeder oviposition, early effective parasitoid host killing and filled/viable achenes. Effective pollination must not be inferred from visit or head position; parasitoid appearance after seed damage cannot be called seed rescue. A phylogenetic adaptive/non-synchrony interpretation additionally requires replicated independently derived Cirsium lineages and ancestry/traits, not 42 taxon image RVs or source-based minimum transition counts.

## Mandatory novelty boundary from primary literature

- **Shibata et al. (2025), Nature Communications 16:4132**, DOI [10.1038/s41467-025-59337-6](https://doi.org/10.1038/s41467-025-59337-6): in *Arabidopsis halleri*, weather-driven orientation switches between sunny pollinator attraction and rainy floral protection, with functional experiments. Therefore **flower orientation balances attraction and protection** is already demonstrated elsewhere; it is not a new theorem, mechanism or Nature-level novelty by itself.
- **Fenster et al. (2009), New Phytologist**, DOI [10.1111/j.1469-8137.2009.02852.x](https://doi.org/10.1111/j.1469-8137.2009.02852.x): floral orientation and visitor movement were an established question.
- **Yu et al. (2021), Plant Biology**, DOI [10.1111/plb.13197](https://doi.org/10.1111/plb.13197): in *Abelia*, orientation modifies rain/wetting, visitation and pollen-load outcomes.
- **Cirsium nutans/acanthodes seed-release orientation study**, [PMC2274966](https://pmc.ncbi.nlm.nih.gov/articles/PMC2274966/): a thistle orientation difference did not explain seed-release differences in its specific aerodynamic experiment. This result concerns **dispersal**, not pollen, egg-entry or mechanical anti-herbivore selection.
- Earlier *Cirsium discolor* sticky-surface null/avoidance observations and *C. palustre* multi-antagonist timing mean a visible barrier need not create a pollination–defence trade-off; consult the original PR #44 source audit.

**Claim ceiling:** A new, potentially biologically informative *three-way conditional interface test* (head orientation × authentic armature × insect guild/entry route, stratified by phenophase) is now mathematically and observationally falsifiable. It has **not** demonstrated selective access, orientation-dependent fitness gain, adaptive module reassembly, or repeated evolution. Scientific value will depend on whether direct differential routes and mature seeds actually disagree with attraction-only / abiotic alternatives.
