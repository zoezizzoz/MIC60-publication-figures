# Common-width figures — September 23, 2026

The current editable masters are in `/Users/picklejuice/Desktop/MIC60 Final Figs/`. `After_AI/` contains byte-identical saved copies of all eight final masters. `FINAL_MANIFEST.json` records their SHA-256 hashes. All were saved before delivery.

## Layout

- Every figure has a 612 pt (8.5 in / 215.9 mm) artboard width. The common artwork span is nominally 576 pt (8 in), with 18 pt side margins. Glyph outlines and strokes can extend by about 1 pt. Heights fit each figure's composition, with approximately 18 pt top/bottom margins.
- Original physical font sizes and font styles are retained. Artwork was first fitted proportionally; text frames were inverse-scaled to preserve their physical sizes.
- All 36 box bodies are exactly 24 pt wide. Compact Figure 3 B/E/H and C/F subgroup pairs have 37 pt center spacing (13 pt between bodies). Wider panels use larger separation to fill their available plot widths.
- Figure 2D uses three 154 pt plotting areas. Figure 3 J–L uses three 155 pt areas, and N uses 204 pt. Figure 4 TIMELESS uses 242 pt and MTT 199 pt. S4C uses 184 pt.
- Figure 3's TEM and fluorescence rows align to the common horizontal envelope. Microscopy panels and scale bars moved together; pixel content and image aspect ratios are preserved.
- Figure 4's quantitative row was widened, and its model enlarged proportionally. S4's blot was enlarged proportionally and balanced against the widened C plot. S1's GO charts extend to the common right edge.
- An old white backing rectangle that clipped Figure 3K after expansion was repositioned. The missing numeric label 2 was restored at the existing area-axis tick; no data or tick location was altered.
- Panel letters at the start of rows share the left margin. Existing text styling, statistical labels, and scientific content are retained.

## Files and provenance

- `Before_Row_Width_Changes/`: all eight original saved AI files, backed up and hash-verified before width edits, including the latest live Figure 3 overlap correction. This is the requested pre-change copy.
- `Before_AI/`: same starting point used by the automated verification.
- `Normalized_AI/`: intermediate snapshot after the common width and physical-font preservation stage.
- `After_AI/`: final editable Illustrator masters.
- `Previews/*_final.png`: final rendered figure previews.
- `Code/`: the native Illustrator layout scripts, plan preparation, inspection, and verification code used for this formatting pass.
- `Layout_Data/`: saved geometry, per-character typography, explicit per-object changes and path roles. These are the exact layout inputs; scientific analysis data and statistical source scripts remain in the existing code/data folders and were not recalculated or changed by this formatting pass.
- `Verification/`: recorded native geometry, box-width, PDF typography, image-pixel and image-proportion checks.

## Native layout sequence

These are state-specific Illustrator edits with UUID/geometry preconditions, not a plotting-data reanalysis. They must be applied to the matching baseline in an already-open document. They should not be reapplied to the final masters.

Set `Layout_Data/REQUEST.json` to `{ "name": "Fig3.ai" }` (or the relevant figure) before a native script call. The scripts target the current master path and guard saved state and original geometry where applicable.

1. `normalize_doc.jsx`: all eight baselines; 612 pt canvas / 576 pt artwork span, physical text sizes retained.
2. `apply_rows.jsx`: Fig3, Fig2, Fig4, FigS4 using the stored `*_ROW_PLAN.json` files.
3. `final_tidy.jsx`: Fig3, Fig1, Fig2, Fig4, FigS1, FigS4; backing shapes, panel-letter alignment, S1 box widths, final Figure 4 height.
4. `finish_details.jsx`: Fig3 (missing tick label) and FigS1 (GO row extension).
5. `verify_boxes.jsx` and `verify_global.py` record the final checks.

The original baselines, intermediate inventories, plans, and final masters are retained together so that the edits are reviewable and recoverable. Illustration and microscopy proportions are maintained; plots are widened through their axes and category spacing, rather than through typography.

## Verification results

- Eight figures saved at exactly 612 pt canvas width.
- All original non-whitespace text characters and physical font sizes preserved. One intentional addition: Figure 3K's missing 7 pt tick label 2.
- Embedded image pixel hashes and aspect ratios identical before/after.
- All 36 box bodies verified at 24 pt in the live Illustrator documents; compact Figure 3 pair spacing verified at 37 pt.
- Every affected quantitative path retained its y coordinate or, for resized point glyphs, its exact data center; 1,203 individual point centers explicitly checked in the row changes.
- Final renders reviewed for clipped plots, annotations, and labels. No statistical calculations or inferred results changed.
