# Proposed update to MIC60-publication-figures

The September 18 repository retained reviewed code but omitted the graph datasets needed to run most panels. This update restores the recoverable analysis inputs, adds a current-panel index of frozen plotted values, and includes the September 22–23 Fig3N, Fig4A/B, S4A and layout generators.

Changes preserve the existing `reviewed_analysis/` and `tables/` structure and Table S1 workflow. Portability edits replace machine-specific recovery paths, restore missing relative inputs and remove unrelated screenshot cleanup from the S4A SVG generator. No scientific measurements were altered. Sequencing reads, acquisition images, native figure assemblies and unrelated recovery alternatives are outside the update.

Validation includes successful reviewed analysis stages, panel generators, graph-input integrity checks and a byte-level manifest. Remaining source/censoring uncertainty for Fig3G and replicate limitations are explicit in LIMITATIONS.md.

Prepared locally on `codex/plotted-data-deposit-2026-09-23`, based on `0924a4e5256396266baa405587a617cfd20d0eed`. This preparation does not publish or push changes.
