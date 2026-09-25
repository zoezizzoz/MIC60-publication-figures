# Figure S3: selected genes and fly homologs of human AMPK-pathway genes

This is the current grouped Figure S3, updated 24 September 2026. It replaces the former 12-gene AMPK panel with 131 unique fly genes mapped from the complete human KEGG AMPK signaling pathway, hsa04152. The other five selected lists are unchanged: FOXO (14), DNA replication (31), spargel-associated (7), chromatin regulation (15) and checkpoint (10). The complete figure contains 208 gene-by-panel entries, including genes shared between panels.

## Rebuild

Requirements: Python 3 with `pypdfium2` and Pillow; R with the standard `grid` package. The PDF renderer uses macOS Quartz and Arial to reproduce the reviewed typography. The archived `org.Dm.eg.db` 3.19.1 identifier check is supplied; no annotation database download is needed for the frozen rebuild.

```sh
python3 Code/rebuild_FigS3.py --rscript /path/to/Rscript
```

All 127 frozen inputs are hash-checked first. The command rebuilds membership from the archived DIOPT responses, joins the archived female DESeq2 contrast, assigns display groups, checks the data, and creates editable vector PDFs plus 600-dpi PNGs in `Rebuilt_Output/`. It runs offline, does not fit DESeq2 again, and does not change any Illustrator file. `Code/generate_FigS3_selected_modules.R` can also redraw the PDFs directly from the included plotting tables.

## Membership and display groups

The complete KEGG hsa04152 entry was retrieved on 24 September 2026 (122 human genes). All high- and moderate-confidence DIOPT 9.1 matches were retained and deduplicated by FlyBase ID: 115 human genes have retained matches, giving 200 mapping edges and 131 unique fly genes. High confidence requires best matches in both directions and score >=2; moderate confidence requires score >=4 or a best match in either direction and score >=2. Gene selection did not depend on RNA-seq results.

Panel D groups these genes into 13 human-source component/function categories. AMPK subunits are ordered alpha, beta, gamma; remaining genes are alphabetical within groups. These are manually specified display categories, not official KEGG subpathway sets. Broad groups such as PP2A-, CFTR-, CD36- and CIDEA-related proteins denote mapping relationships and do not establish the corresponding function in flies. Only AMPK membership was redefined; the other five lists remain selected exploratory lists.

## Files

- `Rebuilt_Output/FigS3_supporting_data.csv`: all 208 gene-by-panel entries, six unchanged DESeq2 numeric columns, gene-set provenance and AMPK display groups.
- `Rebuilt_Output/AMPK_grouped_gene_set_and_expression.csv`: all 131 AMPK-panel rows in display order.
- `Rebuilt_Output/AMPK_component_group_mapping.csv`: all 200 retained human-to-fly edges, confidence scores, KEGG nodes, group assignments and rationale.
- `Rebuilt_Output/AMPK_component_group_definitions.csv`: 13 category definitions and counts.
- `Supporting_Data/AMPK_sources/`: frozen KEGG entry/KGML, all 122 DIOPT responses, retained and excluded mappings, mapping audit, identifier check and source documentation.
- `Supporting_Data/female_DE_source.csv`: unchanged archived female DESeq2 results, with three biological libraries per genotype.
- `Supporting_Data/selected_gene_modules_previous.csv`: previous membership preserved solely to verify the unchanged five lists and document the AMPK membership replacement.
- `Documentation/Manuscript_wording.md`: matching caption, Methods, Results and source references.
- `QA/`: membership, grouping and reconstruction checks.

## Plot interpretation

Red/blue points meet both BH-adjusted P <0.05 and absolute log2 fold change >=0.58. Gray points have adjusted P values but fail one or both thresholds. Open circles have fold-change estimates without adjusted P values. Rows with no fold-change estimate are marked with an em dash (—). Point size represents -log10(adjusted P), capped at 50. All values use the existing genome-wide female dMIC60-CS-versus-dMIC60-WT contrast; both transgenes are in a dMIC60-null background.

Of 131 AMPK-panel genes, 113 have fold-change estimates and 82 have adjusted P values. Takl1 and ninaD meet both thresholds for decreased expression; none meets both for increased expression. The 18 missing fold-change rows have zero counts in the supplied female count matrix; the upstream counting/identifier handling has not been independently resolved. Missing estimates are not evidence of absent gene expression or an unchanged pathway. This is a descriptive gene-level display, not enrichment or a measurement of AMPK activity.

## Artwork

Full figure: 612 x 981 pt. Standalone AMPK panel: 612 x 500 pt. Typography matches the current Fig4, S4 and S2 reference artwork: Arial, with 7 pt gene labels, ticks and legends; 8 pt common-axis title; 9 pt module and component headings; and 13 pt panel letters. Gene labels remain italic, and headings and panel letters remain bold. Both native Illustrator copies now use the same typography. The overall figure title and contextual notes are in the caption; plain headings, standalone upper-left panel letters and open left/bottom axes match the other figures. Axes, ticks and point outlines are 0.426791 pt, measured from the current S4 axes. All text content, plotted values, data-symbol geometry and canvas dimensions are preserved in this typography update. High-resolution PNGs are intended for viewing; vector PDFs and native AI retain editable artwork.

## References

- KEGG hsa04152, https://www.kegg.jp/entry/hsa04152 (accessed 24 September 2026). The exact membership is the archived entry, not a gene list taken from the general database paper.
- Kanehisa et al. (2025), *Nucleic Acids Research* 53, D672–D677. https://doi.org/10.1093/nar/gkae909
- Hu et al. (2011), *BMC Bioinformatics* 12, 357. https://doi.org/10.1186/1471-2105-12-357
- DIOPT confidence definitions: https://www.flyrnai.org/DIOPT_help.html
