# Cirsium sieboldii CS0 same-individual image capture v1

Status: prospective acquisition contract for Gate CS observations that may later enter the Fig.3 F-M-R genomic test.

This contract does not authorize collection, manipulation, off-path access or tissue sampling.

## Purpose

Gate CS itself only needs directly observed colour and anthesis orientation.

The Nature Fig.3 test is stronger: it asks whether the **frozen Azami phenotypic-module partition** is also the unit of genomic inheritance.

The frozen Azami constructs are image-derived. A small addition to the CS0 photo protocol can therefore preserve the same individual for later measurement of:

- presentation_angle;
- floral_lightness, floral_chroma, floral_hue;
- head_elongation, head_compactness;
- involucre_form, projection_prominence, projection_pattern where image quality permits.

This prevents a successful Gate CS census from becoming unusable for the later F-M-R test.

## Source measurement boundaries

The acquisition rules below are derived from the frozen Azami contracts:

- presentation_angle: head and floral end visible; source minimum head dimension 96 px;
- colour endpoints: visible corolla with usable exposure; source minimum 96 px;
- head outline endpoints: usable head or side/oblique outline; source minimum 96 px;
- involucre architecture and projection endpoints: high-resolution side/oblique head or involucre contour; source minimum 150 px;
- projection asymmetry is not admissible without pose control.

These are **image-derived constructs**, not direct 3D morphology, reflectance, pigment concentration, spine length or developmental asymmetry.

## Capture A — whole-head side profile

Required for every CS0 individual if visible from permitted access.

Frame:
- same flowering head used for the A1/A2 orientation score;
- whole capitulum outline visible;
- floral end and head base visible;
- a peduncle segment visible;
- gravity reference or standardized vertical frame visible;
- natural posture only; do not reposition the plant.

Later QC:
- usable head outline;
- source-compatible minimum head dimension >=96 px;
- gravity reference visible;
- EXIF/orientation metadata retained.

Potential constructs:
- presentation_angle;
- head_elongation;
- head_compactness.

## Capture B — corolla colour

Required for every CS0 individual if corolla is visible.

Frame:
- fresh visible florets from the same scored head;
- avoid clipping/highlight saturation and deep cast shadow where possible;
- use the same device/camera mode within a census;
- include a colour reference in-frame where this can be done without touching the plant or violating access rules;
- if a colour reference cannot be placed, record that explicitly rather than fabricating calibration.

Later QC:
- visible corolla with usable exposure;
- chromatic corolla pixels available;
- source-compatible minimum floral/head target dimension >=96 px.

Potential constructs:
- floral_lightness;
- floral_chroma;
- floral_hue.

Boundary:
these remain image-derived colour values unless a separate spectral/pigment protocol is opened.

## Capture C — high-resolution involucre contour

Required prospectively for individuals that may enter the later multi-module genomic test, where visibility permits.

Frame:
- side or controlled oblique view of the same head;
- involucre outline and outer phyllary projections as unobscured as possible;
- no physical rotation or repositioning of the head;
- high enough resolution that the involucre target can pass the frozen 150 px minimum;
- retain natural background if required by access restrictions, but maximize contour contrast by changing observer position only.

Optional:
- a second approximately orthogonal observer view, obtained by moving the observer rather than the plant, when public access permits.

Later QC:
- involucre target dimension >=150 px;
- outline usable;
- side/oblique pose class recorded;
- projection-pattern/asymmetry constructs require a frozen pose-control pass.

Potential constructs:
- involucre_form;
- projection_prominence;
- projection_pattern.

## Strong Fig.3 admission route

The full R-versus-M Nature test requires at least:
- two predeclared multi-trait modules;
- at least two constructs within each admitted multi-trait module;
- at least five admitted constructs total.

A minimal feasible image route is therefore:

**colour module**
- floral_lightness;
- floral_chroma;
- floral_hue;

plus

**head-form module**
- head_elongation;
- head_compactness.

That is five constructs across two multi-trait modules.

The involucre-armature module is additional leverage, not a post-hoc rescue. It is admitted only when the prospective high-resolution/pose QC passes.

## Individual linkage

All three capture roles must use the same:
- observation_id;
- genet_or_individual_id;
- head_id where the same head is visible;
- date;
- population_id;
- flower_stage.

Do not substitute a better photograph from another individual.

## Required image-QC fields

For each individual record:

- whole_head_side_image_id;
- corolla_colour_image_id;
- involucre_highres_image_id;
- involucre_second_view_image_id if available;
- gravity_reference_visible;
- colour_reference_visible;
- whole_head_min_dimension_px after QC;
- colour_target_min_dimension_px after QC;
- involucre_target_min_dimension_px after QC;
- whole_head_outline_usable;
- colour_exposure_usable;
- involucre_outline_usable;
- involucre_pose_control_pass;
- capture_device_id or device class;
- image_qc_notes.

Pixel dimensions may be filled during image QC rather than in the field.

## Fail-closed construct admission

A missing or poor image does not exclude an individual from Gate CS colour × orientation classification if the primary Gate measurements are otherwise valid.

However, that construct is unavailable for the later F-M-R test.

Do not impute missing constructs from:
- another head;
- another individual;
- species-level authority descriptions;
- a different developmental stage.

The construct-admission set is frozen before genotype-to-phenotype genomic discovery.

## Field priority

Do not spend so long acquiring perfect images of early individuals that the census becomes morph-biased.

Priority order:
1. unbiased population census;
2. primary colour + A1/A2 orientation;
3. Capture A and B;
4. Capture C where feasible.

Nature-scale imaging depth must never distort the Gate CS population sample.
