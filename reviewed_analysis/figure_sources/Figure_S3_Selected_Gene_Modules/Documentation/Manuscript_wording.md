# Figure S3 manuscript wording

## Caption

Figure S3 Expression changes in selected genes and fly homologs of human AMPK-pathway genes.

Panels show selected FOXO-associated genes (A), DNA-replication genes (B), spargel-associated genes (C), fly homologs of human AMPK-pathway genes (D), chromatin-regulation genes (E), and DNA-damage-checkpoint genes (F). Panel D contains 131 fly genes mapped from the human KEGG AMPK signaling pathway (hsa04152) using DIOPT 9.1, retaining high- and moderate-confidence matches (see Methods). Genes are grouped by the pathway roles or protein families of their mapped human components; these display categories do not establish each fly gene’s AMPK-related function. Each point shows the DESeq2 log₂ fold change for female dMIC60-CS relative to female dMIC60-WT; both genotypes express the indicated transgene in a dMIC60-null background. Red and blue points indicate increased and decreased expression, respectively, meeting both thresholds: Benjamini-Hochberg-adjusted P < 0.05 and |log₂ fold change| ≥ 0.58. Gray points have adjusted P values but do not meet both thresholds. Labels report adjusted P values below 0.05, including genes that do not meet the fold-change threshold. Point size represents −log₁₀(adjusted P), capped at 50. Open circles indicate fold-change estimates without adjusted P values; rows without fold-change estimates are marked with an em dash (—). Values come from the same female DESeq2 contrast used in the main analyses (n = 3 independent biological libraries per genotype). These panels are descriptive gene-level summaries, not pathway-enrichment tests or measurements of AMPK activity.

## Methods

**AMPK-pathway gene selection and ortholog mapping**

The 122 human genes in the KEGG AMPK signaling pathway (hsa04152; accessed 24 September 2026; Kanehisa et al., 2025) were mapped to *Drosophila melanogaster* using DIOPT 9.1 (Hu et al., 2011). High- and moderate-confidence matches were retained, including multiple fly matches per human gene, and deduplicated by FlyBase identifier. This yielded 131 fly genes corresponding to 115 human pathway genes. High-confidence matches required best matches in both directions and a DIOPT score ≥2; moderate-confidence matches required a score ≥4 or a best match in either direction with a score ≥2. Gene selection was independent of differential-expression results.

Fly gene symbols were checked against org.Dm.eg.db 3.19.1 and matched to the existing female dMIC60-CS-versus-dMIC60-WT DESeq2 results. Fold changes and genome-wide Benjamini–Hochberg-adjusted P values were retained without rerunning DESeq2 or adjusting P values within the subset. Genes were arranged into 13 display groups according to the pathway roles or protein families of their mapped human counterparts; these categories do not establish the corresponding functions in flies. Of the 131 genes, 113 had fold-change estimates, including 82 with adjusted P values. An em dash (—) indicates an unavailable fold-change estimate; open circles indicate estimates without adjusted P values. Gene membership, expression values, mapping scores and display-group assignments are provided in AMPK_component_group_mapping.csv and AMPK_grouped_gene_set_and_expression.csv in the accompanying GitHub repository. This analysis summarizes gene expression; it does not test pathway enrichment or measure AMPK activity.

## Results

We also examined 131 fly homologs mapped from the human KEGG AMPK signaling pathway (Fig. S3D). Takl1 and ninaD showed decreased expression in dMIC60-CS females (log₂ fold changes −0.88 and −0.84; adjusted P = 0.0023 and 0.0128, respectively); no gene in this set met both thresholds for increased expression. These gene-level differences do not establish a change in AMPK activity.

## References

Hu, Y., Flockhart, I., Vinayagam, A., Bergwitz, C., Berger, B., Perrimon, N. and Mohr, S. E. (2011). An integrative approach to ortholog prediction for disease-focused and other functional studies. BMC Bioinformatics 12, 357. https://doi.org/10.1186/1471-2105-12-357

Kanehisa, M., Furumichi, M., Sato, Y., Matsuura, Y. and Ishiguro-Watanabe, M. (2025). KEGG: biological systems database as a model of the real world. Nucleic Acids Res 53, D672-D677. https://doi.org/10.1093/nar/gkae909
