# Published basis for the focused AMPK panel

The revised Figure S3D is a non-exhaustive, publication-supported selection of nine fly genes. The scope is the AMPK heterotrimer, the Lkb1 upstream connection, ACC metabolism, the TSC2–TOR translation branch, and Atg1-mediated autophagy. Genes were chosen for experimentally documented connections in Drosophila in the four primary studies below. The panel does not claim that all nine genes are direct AMPK targets or transcriptionally regulated by AMPK.

This curation was performed after inspection of the RNA-seq data and is not a prospectively specified gene set. Neither effect size nor statistical significance was an inclusion criterion. It is not a systematic literature review or an exhaustive list of published AMPK-associated genes. The nine genes are a subset of the archived 131-gene homolog survey; family similarity alone was insufficient for the focused display. Other members remain in the full survey, so exclusion from this display is not evidence that a gene has no AMPK function.

The table AMPK_literature_selection.csv records each gene, its experimental evidence, the source location, and the limits of interpretation. Core subunits and direct substrates are distinguished from downstream pathway outputs. The selection is encoded explicitly and joined to the unchanged female DESeq2 results. Missing estimates are displayed rather than used to exclude genes.

## Primary publications

- Pan, D. A. and Hardie, D. G. (2002). A homologue of AMP-activated protein kinase in Drosophila melanogaster is sensitive to AMP and is activated by ATP depletion. Biochem. J. 367, 179–186. https://doi.org/10.1042/BJ20020703

- Castanieto, A., Johnston, M. J. and Nystul, T. G. (2014). EGFR signaling promotes self-renewal through the establishment of cell polarity in Drosophila follicle stem cells. eLife 3, e04437. https://doi.org/10.7554/eLife.04437

- Kim, M. and Lee, J. H. (2015). Identification of an AMPK phosphorylation site in Drosophila TSC2 (gigas) that regulate cell growth. Int. J. Mol. Sci. 16, 7015–7026. https://doi.org/10.3390/ijms16047015

- Ulgherait, M., Rana, A., Rera, M., Graniel, J. and Walker, D. W. (2014). AMPK modulates tissue and organismal aging in a non-cell-autonomous manner. Cell Rep. 8, 1767–1780. https://doi.org/10.1016/j.celrep.2014.08.006

## Preserved broad survey

AMPK_gene_set_and_expression.csv and AMPK_grouped_gene_set_and_expression.csv retain all 131 genes and their original values. AMPK_component_group_mapping.csv retains every human-to-fly mapping. The 13 human-source display categories describe this broad survey only; they are not the basis for the focused panel. AMPK_full_survey_with_display_selection.csv adds an explicit flag for membership in the focused panel.

The prior 208-row full-figure table is retained as FigS3_full_orthology_survey.csv. The current FigS3_supporting_data.csv contains the five unchanged lists (77 entries) and nine AMPK entries, totaling 86 gene-by-panel entries.

## Evidence checks

The retracted Mirouse et al. (2007) polarity paper was not used. Lkb1 inclusion is supported by the independent Castanieto et al. (2014) experiment. Primary sources support the specific roles recorded here; they do not establish AMPK activation in the present RNA-seq experiment or validate the memberships of the other five selected panels.
