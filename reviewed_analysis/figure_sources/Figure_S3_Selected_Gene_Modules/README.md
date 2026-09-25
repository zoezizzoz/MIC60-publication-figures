# Figure S3 selected genes with a publication-supported AMPK panel

Figure S3D now displays nine AMPK-associated fly genes: AMPKalpha, alc, SNF4Agamma, Lkb1, ACC, gig, S6k, Thor and Atg1. Each has specific experimental support from a primary publication in Drosophila. The panel is a focused, non-exhaustive selection assembled after the RNA-seq analysis, without selecting for differential-expression significance. It is not a published target-gene signature, a systematic literature review, an enrichment test, or an AMPK activity assay.

The five other selected lists retain their memberships and values: FOXO (14), DNA replication (31), Spargel (7), chromatin (15), and checkpoint (10). The current figure has 86 gene-by-panel entries. This AMPK revision does not establish the original published provenance of those five lists.

## Rebuild

Run `python3 Code/rebuild_FigS3.py --rscript /path/to/Rscript`. Dependencies are R with grid and Quartz/Arial, Python with pypdfium2 and Pillow. The rebuild checks 127 frozen KEGG/DIOPT/DESeq2 inputs, reconstructs the complete homolog survey, selects the nine documented genes, and renders the compact figure and standalone supporting plots. It runs offline without fitting DESeq2 again.

## Publication support

Supporting_Data/AMPK_literature_selection.csv records each selected gene, its role, the exact supporting publication, experiment or figure, and interpretation limits. Supporting_Data/AMPK_literature_references.csv contains the four source references. Documentation/AMPK_published_selection.md explains the scope and selection process.

The AMPK subunits and ACC are supported by Pan and Hardie (2002), Lkb1 by Castanieto et al. (2014), gig and the S6k/Thor translation outputs by Kim and Lee (2015), and Atg1 by Ulgherait et al. (2014). Downstream connections are not presented as proof of direct AMPK phosphorylation or transcriptional regulation.

## Current plotted data

- Rebuilt_Output/FigS3_supporting_data.csv: 86 gene-by-panel entries, with unchanged DESeq2 values and publication evidence for the nine AMPK entries.
- Rebuilt_Output/AMPK_literature_selected_gene_set_and_expression.csv: nine selected genes with expression values and evidence.
- Rebuilt_Output/FigS3_selected_gene_modules.pdf and .png: compact six-panel figure.
- Rebuilt_Output/AMPK_literature_selected_panel.pdf and .png: standalone focused panel.
- Documentation/FigS3_figure_legend.txt and Documentation/Manuscript_wording.md: matching caption, Results and Methods.

## Preserved comprehensive supporting survey

The complete 131-gene survey is retained unchanged in AMPK_gene_set_and_expression.csv and AMPK_grouped_gene_set_and_expression.csv. These are fly homologs of all 122 human KEGG hsa04152 genes, mapped with DIOPT 9.1, retaining high- and moderate-confidence matches and deduplicating by FlyBase identifier. There are 200 mapping edges covering 115 human genes. Full mappings, confidence definitions and frozen inputs remain in Supporting_Data/AMPK_sources and AMPK_component_group_mapping.csv.

AMPK_full_survey_with_display_selection.csv contains all 131 genes plus a flag identifying the nine in the current panel. FigS3_full_orthology_survey.csv preserves the former 208-row figure data. AMPK_signaling_components_grouped_panel.pdf and .png retain the full 13-group survey as a supporting plot. Those categories describe the mapped human components; they do not establish the fly genes’ functions. Takl1 and ninaD remain in these files with their original values.

## Plot interpretation

All values use the existing genome-wide female dMIC60-CS versus dMIC60-WT DESeq2 contrast (three biological libraries per genotype; both transgenes in a dMIC60-null background). Filled red/blue points meet BH-adjusted P <0.05 and absolute log2 fold change >=0.58; gray points have adjusted P values but fail one or both criteria. Labels show adjusted P values below 0.05, including genes below the fold-change threshold. Point size represents -log10(adjusted P), capped at 50. Open circles indicate missing adjusted P values; em dashes indicate missing fold-change estimates.

For the nine selected genes, five have adjusted P values and none meets both thresholds. Thor has a fold-change estimate without an adjusted P value; AMPKalpha, Lkb1 and gig lack fold-change estimates. Missing estimates do not establish absence of expression or lack of pathway activity. In the complete survey, 113 genes have fold changes and 82 have adjusted P values; Takl1 and ninaD meet both thresholds for decreased expression. The original missing-estimate rows have zero counts in the supplied female count matrix; upstream counting/identifier handling remains unresolved.

## Artwork

The compact figure is 612 x 669 pt. Gene labels and legend/tick text remain 7 pt, panel titles 9 pt, panel letters 13 pt, and the common axis title 8 pt, matching the current saved native figure at the start of this revision. Axis, tick and open-point outline widths remain 0.426791 pt. Panel letters are upper left. Native AI artwork is editable vector content with no external links or rasterized plots.

The gray headers (#EFEFEF), black plot outlines and adjusted-P annotations incorporate the side chat’s saved formatting (repository commit 9b4be42). The compact figure has four annotations; the full survey retains two. Labels are shown for adjusted P <0.05 irrespective of the fold-change cutoff, preserving the distinction between annotated gray points and points meeting both thresholds.

QA/compact_publication_validation.json and QA/rebuild_validation.json describe the current compact figure. Other older QA records describe earlier survey/layout stages and are retained as history.
