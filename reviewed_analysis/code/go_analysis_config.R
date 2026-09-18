# One analysis specification for every sex-by-direction GO comparison.
# These are the previously reviewed settings; centralization does not change them.
# Preserve selection during presentation edits; do not change it to fit the layout.
# The top-eight display rule was introduced in the earlier review rerun; the
# original artwork selection code was unavailable, so this is not a verified
# reconstruction of the paper's original term-selection rule.
GO_CONFIG <- list(
  contrast = "dMIC60-CS versus dMIC60-WT",
  sexes = c("female", "male"),
  directions = c("up", "down"),
  ontologies = c("BP", "CC", "MF"),
  de_adjusted_p_cutoff = 0.05,
  absolute_log2_fold_change_cutoff = 0.58,
  background_rule = "Non-missing DESeq2 adjusted P, separately within each sex",
  annotation_package = "org.Dm.eg.db",
  key_type = "SYMBOL",
  minimum_gene_set_size = 10L,
  maximum_gene_set_size = 500L,
  p_adjust_method = "BH",
  pool_ontologies = FALSE,
  go_adjusted_p_cutoff = 0.05,
  top_terms_per_ontology = 8L,
  selection_order = c("p.adjust", "pvalue", "ID"),
  adjustment_family = "Separate ontology, sex, and direction",
  redundancy_filter = "None"
)
