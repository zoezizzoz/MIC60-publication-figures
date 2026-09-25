# Proposed update to MIC60-publication-figures

The September 18 repository retained reviewed code but omitted the graph datasets needed to run most panels. This update restores the recoverable analysis inputs, adds a current-panel index of frozen plotted values, and includes the September 22–23 Fig3N, Fig4A/B, S4A and layout generators.

Changes preserve the existing `reviewed_analysis/` and `tables/` structure and Table S1 workflow. Portability edits replace machine-specific recovery paths, restore missing relative inputs and remove unrelated screenshot cleanup from the S4A SVG generator. No scientific measurements were altered. Sequencing reads, acquisition images, native figure assemblies and unrelated recovery alternatives are outside the update.

Validation includes successful reviewed analysis stages, panel generators, graph-input integrity checks and a byte-level manifest. Remaining source/censoring uncertainty for Fig3G and replicate limitations are explicit in LIMITATIONS.md.

Prepared locally on `codex/plotted-data-deposit-2026-09-23`, based on `0924a4e5256396266baa405587a617cfd20d0eed`. The September 23 deposit and September 24 Figure S3 updates were prepared together in pull request #1.

The final common-width figure pass is recorded in provenance/final_layout_2026-09-23/, including the native edit scripts, per-object geometry plans, recorded preservation checks, and hashes for all eight baseline/final assemblies. No plotted table or scientific analysis result changed in this supplement. Native AI files and the manuscript remain outside the deposit.

## Figure S3, 24 September 2026

Replaces the former 12-gene AMPK list with 131 unique fly homologs of the complete human KEGG AMPK pathway (122 human genes, DIOPT 9.1 high/moderate-confidence mapping), arranged in 13 explicitly defined display groups. Includes frozen source records, all retained/excluded mappings, exact archived expression values, matching caption/Methods/Results wording, source references and full/standalone vector-PDF and 600-dpi PNG exports. The five other selected gene lists and archived DESeq2 estimates are unchanged. The new offline `figs3` rebuild and deposit verifier check membership, grouping, missing estimates and numeric agreement. Four S3 tables are indexed; the full deposit now has 43 tables covering 29 graph panels. The Methods wording and CSV availability statement match the updated gene-set workflow.

Figure S3 display clarification: the 18 rows without fold-change estimates now use an em dash, explained once in the figure key as “— Fold change unavailable.” All gene memberships, numeric values and point positions are unchanged. Matching caption/Methods wording and 600-dpi exports are updated.

The initial Figure S3 formatting pass used: no overall title or shared-context/footer lines, plain module/component headings, separate upper-left panel letters and open left/bottom axes. Axis/tick/outline width matches the measured S4 axes. The canvas is cropped without scaling; all plotted data and symbol geometry are preserved. Typography is standardized as described below. The caption retains sample size, contrast, grouping and cutoff information. The generator and full/standalone exports are synchronized.

Figure S3 legend labels now read “Higher in CS” and “Higher in WT,” matching the surrounding figures. Legend colors and internal expression-status codes are unchanged.

Figure S3 typography now matches the current Fig4, S4 and S2 artwork: Arial at 7 pt for gene labels, ticks and legends; 8 pt for the common-axis title; 9 pt for module/component headings; and 13 pt for panel letters. Both native Illustrator copies use the same profile. All 302 text frames per copy were checked; all text content, plotted values and path geometry are preserved. Source PDFs and 600-dpi PNGs use the same sizes.

Figure S3 now uses light-gray header strips and black outlines following the supplied visual reference. Six adjusted-P labels are drawn for every estimable gene with adjusted P < 0.05, including gray points that fall below the fold-change threshold. Values below 0.001 are displayed as “adj. p < 0.001”; other labels use up to five decimal places. Data, point positions and the 7/8/9/13 pt typography are unchanged. Both native copies, full/standalone exports and caption wording are synchronized.

## Focused publication-supported AMPK panel

Figure S3D now displays nine genes supported by primary Drosophila experiments: AMPKalpha, alc, SNF4Agamma, Lkb1, ACC, gig, S6k, Thor and Atg1. Gene-specific evidence and limitations are recorded for Pan and Hardie (2002), Castanieto et al. (2014), Kim and Lee (2015), and Ulgherait et al. (2014). This non-exhaustive selection was assembled after RNA-seq analysis using biological evidence rather than significance. The complete 131-gene survey, including Takl1 and ninaD, remains unchanged in supporting data; an explicit flag identifies displayed genes. No DESeq2 analysis or P adjustment was rerun.

The current six-panel figure contains 86 entries and has a 612 × 669 pt canvas. Gray headers, black outlines, upper-left letters, the 7/8/9/13 pt typography and adjusted-P labeling convention are preserved from the side chat. Four labels remain in the compact figure; the full-survey supporting plot retains its two labels. Corresponding Results, Methods and caption text and four references are supplied. Five displayed AMPK genes have adjusted P values and none meets both thresholds; four have unavailable estimates, so no conclusion about AMPK activity is drawn.

## Consolidation of current S3 code

Removed the unused S3 style helper, four superseded layout reports and a duplicate report. The full-survey renderer now renders only the supporting plot; its unused former six-panel layout branch was removed. Both the standalone S3 command and the reviewed full-panel workflow use the complete current rebuild. Current code and QA indexes distinguish the nine-gene figure from the intentionally retained 131-gene supporting survey. No numerical input or plotted value was changed.
