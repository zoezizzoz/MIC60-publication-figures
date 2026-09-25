# Final Illustrator layout records — 23 September 2026

This directory records the last common-width formatting pass for all eight figures. It supplements the graph generators and frozen plotting tables; it does not recalculate measurements or statistics. The copied code, layout inputs, and verification records are byte-identical to the local completed layout package. SOURCE_FILES.json records each source and checksum.

## What is included

- Code/: nine original Illustrator inspection/editing and Python preparation/verification scripts.
- Layout_Data/: starting geometry, text and physical font sizes, normalized inventories, per-object row plans/results, and path-role mappings.
- Verification/: recorded canvas, typography, embedded-image, image-proportion, box-width and plotted-coordinate checks.
- BASELINE_MANIFEST.json and FINAL_MANIFEST.json: hashes identifying the exact native figure states before and after editing.
- SOURCE_README.md: the original local layout description and edit sequence. Its directory references describe the author's local archive, not files supplied in this repository.

The recorded final canvas width is 612 pt for all eight figures. Box bodies are 24 pt; compact Fig3 category pairs use 37 pt centers. Existing physical text sizes and embedded image pixels were checked before/after. A missing Fig3K area-axis tick label 2 was restored. Verification records report 1,203 individual data-point centers checked for preserved y positions during row changes. The local source package retains backups from before any common-width edits.

## Execution boundary

These scripts are an archival record of native edits, not portable graph entrypoints. They retain the original absolute macOS paths, Illustrator object UUIDs, and saved-state/geometry preconditions. Several scripts call save() on open documents. Do not run them on the final figures or an unrelated Illustrator session. Replaying them requires the exact baseline assemblies, the matching Illustrator environment, path configuration, and a separate backup. Those native AI files, embedded acquisition images, and render previews are outside the user's code/data deposit scope; their identifying hashes are supplied here.

The graph data and statistical generators remain in plotted_data/, reviewed_analysis/, and panels/. reproduce.py does not run these native editing scripts. This supplement improves the record of the final layout; it does not claim that a repository-only rebuild reproduces the complete final Illustrator composition. Scientific limitations in ../../LIMITATIONS.md still apply.
