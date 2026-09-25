# Figure S3: literature-reviewed selected genes

The current figure contains 81 group–gene entries (79 unique genes): FOXO-associated (10), DNA replication/genome maintenance (28), mitochondrial/metabolic (7), AMPK-associated (11), chromatin regulation (15), and DNA repair/checkpoints (10). These are non-exhaustive descriptive groups assembled and revised after the RNA-seq analysis. They are not direct-target signatures, enrichment tests or measurements of pathway activity.

## Revision following individual review

All 89 original entries were checked before this revision. Four unverified FOXO assignments (Lsp1alpha, Lsp1gamma, Gadd45 and dap) and three replication entries supported only by annotation in the checked evidence (DNAlig1, DNAlig3 and RfC3) were omitted. This does not establish absence of biological relevance. Panel C retains its seven genes under a mitochondrial/metabolic heading, without implying that they are validated Spargel targets. Panels B and F now encompass their published repair/genome-maintenance functions. Chromatin membership is unchanged.

Sesn and Atg8a were restored to AMPK after verifying fly experiments, bringing that panel from nine to 11 genes. The focused selection retains AMPKalpha, alc, SNF4Agamma, Lkb1, ACC, gig, S6k, Thor and Atg1. Seven displayed AMPK genes have adjusted P values; none meets both differential-expression thresholds. Three have no fold-change estimate and Thor has no adjusted P value. Missing estimates were not a selection criterion.

## Evidence and original data

- Supporting_Data/S3_original_gene_display_decisions.csv: all 89 original entries, inclusion decisions, reasons, publication evidence and interpretation limits.
- Supporting_Data/S3_original_gene_audit.csv and Documentation/S3_original_gene_audit.md: the dated audit of the original lists, including uncertainties and full references.
- Supporting_Data/selected_gene_modules_previous.csv: unchanged original memberships and expression values.
- Supporting_Data/AMPK_literature_selection.csv and AMPK_literature_references.csv: evidence for the 11 displayed AMPK genes from five primary publications.
- Rebuilt_Output/FigS3_supporting_data.csv: the exact 81 displayed entries, with unchanged archived expression values and evidence.
- Documentation/FigS3_figure_legend.txt and Manuscript_wording.md: matching caption and suggested Results/Methods.

Biological functions and named-regulator relationships are distinguished. Evidence can be indirect or tissue-specific; those limits remain explicit. Present-day publication support does not recover the original selection rationale or prove that the lists were prespecified. Selection decisions did not use dMIC60 effect sizes or P values.

## Rebuild

Run `python3 Code/rebuild_FigS3.py --rscript /path/to/Rscript`. Dependencies are R with grid and Quartz/Arial, Python with pypdfium2 and Pillow. One offline entrypoint checks 127 frozen inputs, reconstructs the full survey, applies the explicit reviewed selections and renders the current six-panel figure plus supporting plots. DESeq2 is not rerun. See Code/README.md.

## Preserved full survey

The complete 131-gene survey is unchanged in AMPK_gene_set_and_expression.csv and AMPK_grouped_gene_set_and_expression.csv. It maps 122 human KEGG hsa04152 genes using DIOPT 9.1 high/moderate-confidence matches, retaining 200 edges across 115 mapped human genes and deduplicating by FlyBase identifier. Full mapping records remain frozen under Supporting_Data/AMPK_sources. The 13 display groups describe human-source components, not demonstrated fly functions.

AMPK_full_survey_with_display_selection.csv flags nine members displayed in the focused panel. Sesn and Atg8a are independently supported fly genes outside that 131-gene orthology survey; they are retained in the focused-panel table without altering the survey. FigS3_full_orthology_survey.csv retains the prior 208-row full-figure table. The complete standalone survey plot remains available.

## Plot and artwork

Values use the same female dMIC60-CS versus dMIC60-WT contrast, with three biological libraries per genotype and both transgenes in a dMIC60-null background. Filled red/blue points meet BH-adjusted P <0.05 and |log2 fold change| >=0.58; gray points with adjusted P values fail one or both criteria. Labels indicate adjusted P <0.05 regardless of the fold-change cutoff. Size is -log10(adjusted P), capped at 50. Open circles indicate unavailable adjusted P; em dashes indicate unavailable fold change.

The figure is 612 × 669 pt and preserves the saved embedded artwork's staggered columns. Fonts remain Arial: 7 pt labels/legends/ticks, 8 pt common axis title, 9 pt panel titles and 13 pt upper-left panel letters. Header fill is #EFEFEF. Axis, tick and open-point outline widths remain 0.426791 pt. Gene rows are redistributed within each column (14.53125 pt left, 9.4898 pt right), with equal 28 pt panel gaps and common column bounds from 47 to 628 pt. Physical text and point sizes are unchanged. Current native AI artwork is editable vector content without external links or raster plots.

QA/rebuild_validation.json checks membership, frozen inputs and unchanged values. QA/compact_publication_validation.json records the final visual and native-artwork inspection. Earlier versions remain in the local backup and Git history.
