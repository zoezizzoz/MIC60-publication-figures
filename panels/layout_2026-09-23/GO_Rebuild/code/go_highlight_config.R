# Display-only examples tied to the manuscript Results paragraph beginning
# "Gene Ontology (GO) analysis of dMIC60-CS females identified enrichment...".
# Include all three members of the specific ATP-synthesis and NADH-dehydrogenase
# terms. Repeat the same NADH genes in the displayed male oxidoreductase term
# for a matched sex/direction comparison. Do not rank examples by P or fold change.
# No analysis, term-selection, count, or enrichment parameter is changed.
GO_HIGHLIGHTS <- data.frame(
  sex=c("female","female","male"),
  direction=c("up","up","down"),
  ID=c("GO:0005753","GO:0003954","GO:0016491"),
  genes=c("knon/ATPsynG/ATP8","ND2/ND4L/ND6","ND2/ND4L/ND6"),
  reason=c(
    "All three contributing genes in the specific mitochondrial ATP-synthase term; ATP synthesis is part of the manuscript oxidative-metabolism interpretation.",
    "All three contributing genes in the specific NADH-dehydrogenase term; mitochondrial electron transport and NADH metabolism are discussed in the manuscript.",
    "The same three NADH-dehydrogenase genes highlighted in females, all confirmed members of this male downregulated term; matched comparison of the manuscript oxidative-metabolism theme."),
  stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(paste(GO_HIGHLIGHTS$sex,GO_HIGHLIGHTS$direction,GO_HIGHLIGHTS$ID)))
