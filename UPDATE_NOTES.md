# Proposed update to MIC60-publication-figures

The September 18 repository retained reviewed code but omitted the graph datasets needed to run most panels. This update restores the recoverable analysis inputs, adds a current-panel index of frozen plotted values, and includes the September 22–23 Fig3N, Fig4A/B, S4A and layout generators.

Changes preserve the existing `reviewed_analysis/` and `tables/` structure and Table S1 workflow. Portability edits replace machine-specific recovery paths, restore missing relative inputs and remove unrelated screenshot cleanup from the S4A SVG generator. No scientific measurements were altered. Sequencing reads, acquisition images, native figure assemblies and unrelated recovery alternatives are outside the update.

Validation includes successful reviewed analysis stages, panel generators, graph-input integrity checks and a byte-level manifest. Remaining source/censoring uncertainty for Fig3G and replicate limitations are explicit in LIMITATIONS.md.

Prepared locally on `codex/plotted-data-deposit-2026-09-23`, based on `0924a4e5256396266baa405587a617cfd20d0eed`. The September 23 deposit was published as draft pull request #1; the September 24 Figure S3 changes below extend that same draft branch.

The final common-width figure pass is recorded in provenance/final_layout_2026-09-23/, including the native edit scripts, per-object geometry plans, recorded preservation checks, and hashes for all eight baseline/final assemblies. No plotted table or scientific analysis result changed in this supplement. Native AI files and the manuscript remain outside the deposit.

## Figure S3, 24 September 2026

Replaces the former 12-gene AMPK list with 131 unique fly homologs of the complete human KEGG AMPK pathway (122 human genes, DIOPT 9.1 high/moderate-confidence mapping), arranged in 13 explicitly defined display groups. Includes frozen source records, all retained/excluded mappings, exact archived expression values, matching caption/Methods/Results wording, source references and full/standalone vector-PDF and 600-dpi PNG exports. The five other selected gene lists and archived DESeq2 estimates are unchanged. The new offline `figs3` rebuild and deposit verifier check membership, grouping, missing estimates and numeric agreement. Four S3 tables are indexed; the full deposit now has 43 tables covering 29 graph panels. The Methods wording and CSV availability statement match the updated gene-set workflow.

Figure S3 display clarification: the 18 rows without fold-change estimates now use an em dash, explained once in the figure key as “— Fold change unavailable.” All gene memberships, numeric values and point positions are unchanged. Matching caption/Methods wording and 600-dpi exports are updated.

Figure S3 now follows the surrounding figure style: no overall title or shared-context/footer lines, plain module/component headings, separate upper-left panel letters and open left/bottom axes. Axis/tick/outline width matches the measured S4 axes. The canvas is cropped without scaling; every remaining native character retains its font, size and scale, and all plotted data and symbol geometry are preserved. The caption retains sample size, contrast, grouping and cutoff information. The generator and full/standalone exports are synchronized.
