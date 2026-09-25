# Plotted graph data

Data are organized using the current figure panel letters. Old source directories retain historical numbering so the existing code still resolves its input paths. Sequencing reads and acquisition images are excluded at the author’s request.

| Panel | Plotted table | Rows | Observation unit |
|---|---|---:|---|
| Fig1B | [PCA_coordinates.csv](Fig1B/PCA_coordinates.csv) | 12 | RNA-seq library |
| Fig1C | [volcano_plotted_genes.csv](Fig1C/volcano_plotted_genes.csv) | 7456 | gene with non-missing adjusted P |
| FigS1C | [volcano_plotted_genes.csv](FigS1C/volcano_plotted_genes.csv) | 9695 | gene with non-missing adjusted P |
| Fig1D | [heatmap_row_z_scores.csv](Fig1D/heatmap_row_z_scores.csv) | 60 | gene (columns are RNA-seq libraries for the matrix) |
| Fig1D | [heatmap_selected_genes.csv](Fig1D/heatmap_selected_genes.csv) | 60 | gene (columns are RNA-seq libraries for the matrix) |
| FigS1D | [heatmap_row_z_scores.csv](FigS1D/heatmap_row_z_scores.csv) | 60 | gene (columns are RNA-seq libraries for the matrix) |
| FigS1D | [heatmap_selected_genes.csv](FigS1D/heatmap_selected_genes.csv) | 60 | gene (columns are RNA-seq libraries for the matrix) |
| Fig2A | [up_displayed_GO_terms.csv](Fig2A/up_displayed_GO_terms.csv) | 16 | GO term |
| Fig2A | [down_displayed_GO_terms.csv](Fig2A/down_displayed_GO_terms.csv) | 21 | GO term |
| FigS1E | [up_displayed_GO_terms.csv](FigS1E/up_displayed_GO_terms.csv) | 3 | GO term |
| FigS1E | [down_displayed_GO_terms.csv](FigS1E/down_displayed_GO_terms.csv) | 17 | GO term |
| Fig2B | [V1_up_28_GO_gene_memberships.csv](Fig2B/V1_up_28_GO_gene_memberships.csv) | 116 | gene membership or STRING edge |
| Fig2B | [V1_up_STRING_edges.csv](Fig2B/V1_up_STRING_edges.csv) | 13 | gene membership or STRING edge |
| Fig2C | [V1_up_28_GO_gene_memberships.csv](Fig2C/V1_up_28_GO_gene_memberships.csv) | 116 | gene membership or STRING edge |
| Fig2C | [V1_up_STRING_edges.csv](Fig2C/V1_up_STRING_edges.csv) | 13 | gene membership or STRING edge |
| Fig2D | [female_normalized_counts.csv](Fig2D/female_normalized_counts.csv) | 18 | gene × female RNA-seq library |
| Fig2D | [adjusted_p_annotations.csv](Fig2D/adjusted_p_annotations.csv) | 6 | gene × sex |
| Fig3A | [sleep_plotted_values.csv](Fig3A/sleep_plotted_values.csv) | 192 | 30-minute bin × genotype summary |
| Fig3B | [sleep_plotted_values.csv](Fig3B/sleep_plotted_values.csv) | 177 | fly |
| Fig3C | [sleep_plotted_values.csv](Fig3C/sleep_plotted_values.csv) | 322 | fly × light phase |
| Fig3D | [activity_plotted_values.csv](Fig3D/activity_plotted_values.csv) | 192 | 30-minute bin × genotype summary |
| Fig3E | [activity_plotted_values.csv](Fig3E/activity_plotted_values.csv) | 170 | fly |
| Fig3F | [activity_plotted_values.csv](Fig3F/activity_plotted_values.csv) | 304 | fly × light phase |
| Fig3G | [recovered_survival_observations_UNVERIFIED.csv](Fig3G/recovered_survival_observations_UNVERIFIED.csv) | 91 | fly, including missing survival times |
| Fig3H | [locomotion_plotted_values.csv](Fig3H/locomotion_plotted_values.csv) | 8 | experimental population |
| Fig3J | [TEM_field_means.csv](Fig3J/TEM_field_means.csv) | 32 | microscopy field |
| Fig3K | [TEM_field_means.csv](Fig3K/TEM_field_means.csv) | 32 | microscopy field |
| Fig3L | [TEM_field_means.csv](Fig3L/TEM_field_means.csv) | 24 | microscopy field |
| Fig3N | [TMRM_MTG_field_means.csv](Fig3N/TMRM_MTG_field_means.csv) | 32 | microscopy field |
| Fig4A | [TIMELESS_DAPI_values.csv](Fig4A/TIMELESS_DAPI_values.csv) | 64 | provided cell/image measurement without biological replicate IDs |
| Fig4B | [mtt_plot_values.csv](Fig4B/mtt_plot_values.csv) | 78 | technical well |
| Fig4B | [mtt_experiment_means.csv](Fig4B/mtt_experiment_means.csv) | 16 | experiment × genotype × dose |
| Fig4B | [mtt_biological_summary.csv](Fig4B/mtt_biological_summary.csv) | 10 | genotype × dose summary |
| FigS1B | [western_blot_quantification.csv](FigS1B/western_blot_quantification.csv) | 6 | paired blot measurement |
| FigS2 | [KEGG_GSEA_plotted_pathways.csv](FigS2/KEGG_GSEA_plotted_pathways.csv) | 20 | pathway × sex |
| FigS3 | [selected_gene_modules.csv](FigS3/selected_gene_modules.csv) | 208 | gene × figure panel |
| FigS3 | [AMPK_grouped_gene_set_and_expression.csv](FigS3/AMPK_grouped_gene_set_and_expression.csv) | 131 | unique fly gene in the AMPK panel |
| FigS3 | [AMPK_component_group_mapping.csv](FigS3/AMPK_component_group_mapping.csv) | 200 | retained human-to-fly mapping edge |
| FigS3 | [AMPK_component_group_definitions.csv](FigS3/AMPK_component_group_definitions.csv) | 13 | AMPK display group |
| FigS4A | [figS4_plotted_values.csv](FigS4A/figS4_plotted_values.csv) | 64 | gene × curated pathway |
| FigS4A | [figS4_enrichment.csv](FigS4A/figS4_enrichment.csv) | 6 | curated pathway |
| FigS4A | [figS4_stress_pathways.csv](FigS4A/figS4_stress_pathways.csv) | 64 | gene × curated pathway |
| FigS4C | [transfection_efficiency.csv](FigS4C/transfection_efficiency.csv) | 20 | image × experiment |

Each row of [INDEX.csv](INDEX.csv) records its source, generator, selection, caveat and SHA-256. The tables retain native column names; gene identifiers are Drosophila symbols, `padj`/`padj_F` are gene-level DESeq2 adjusted P, `log2FoldChange`/`lfc_F` are CS-versus-WT log2 fold changes, and missing adjusted P is not zero.

**Fig3G is not certified as the current plotted dataset.** Historical survival records/code are included for reconciliation; blank times are not assigned a censoring status.

Schematic, workflow, model, blot-photo and micrograph panels do not have numeric plotting tables. Fig1A, Fig3I/M, Fig4C, FigS1A and FigS4B are therefore not treated as missing graph tables.

Figure S3 was updated on 24 September 2026. Only the AMPK membership was redefined, using high/moderate-confidence fly homologs of the complete human KEGG AMPK pathway. The full table retains the five other selected lists unchanged. The 13 AMPK categories describe human-source roles/protein families; they are display groups, not enrichment tests or evidence of AMPK activation. Missing fold changes and adjusted P values are retained explicitly.
