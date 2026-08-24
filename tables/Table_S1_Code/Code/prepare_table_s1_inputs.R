#!/usr/bin/env Rscript

# Reproduce the three analysis products assembled into Supplementary Table S1.
# Female and male libraries are analyzed separately with WT as the reference.

suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  normalizePath("Code")
}
table_dir <- normalizePath(file.path(script_dir, ".."))
input_path <- file.path(table_dir, "Original_Data", "count_matrix_symbol.csv")
output_dir <- file.path(table_dir, "Rebuilt_Output", "generated_inputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

count_mat <- read.csv(input_path, row.names = 1, check.names = FALSE)
count_mat <- as.matrix(count_mat)
storage.mode(count_mat) <- "integer"

stopifnot(
  nrow(count_mat) == 17173L,
  ncol(count_mat) == 12L,
  !anyDuplicated(rownames(count_mat)),
  all(count_mat >= 0L)
)

fit_sex <- function(sample_names) {
  sex_counts <- count_mat[, sample_names, drop = FALSE]
  genotype <- factor(
    ifelse(grepl("^CS", sample_names), "CS", "WT"),
    levels = c("WT", "CS")
  )
  col_data <- data.frame(row.names = sample_names, genotype = genotype)
  dds <- DESeqDataSetFromMatrix(
    countData = sex_counts,
    colData = col_data,
    design = ~ genotype
  )
  dds <- DESeq(dds, quiet = TRUE)
  result <- as.data.frame(results(dds, contrast = c("genotype", "CS", "WT")))
  result$gene <- rownames(result)
  # The deposited publication inputs are ordered from smallest raw P value,
  # with untestable rows retained at the end.
  result <- result[order(result$pvalue, na.last = TRUE), ]
  list(result = result, normalized = counts(dds, normalized = TRUE))
}

female_samples <- colnames(count_mat)[grepl("F[123]$", colnames(count_mat))]
male_samples <- colnames(count_mat)[grepl("M[123]$", colnames(count_mat))]
stopifnot(length(female_samples) == 6L, length(male_samples) == 6L)

female <- fit_sex(female_samples)
male <- fit_sex(male_samples)

write.csv(
  female$result,
  file.path(output_dir, "DEG_CSF_vs_WRF.csv"),
  row.names = TRUE
)
write.csv(
  male$result,
  file.path(output_dir, "DEG_CSM_vs_WRM.csv"),
  row.names = TRUE
)

normalized <- cbind(
  female$normalized[, female_samples, drop = FALSE],
  male$normalized[, male_samples, drop = FALSE]
)
colnames(normalized) <- sub("^WRF", "WTF", colnames(normalized))
colnames(normalized) <- sub("^WRM", "WTM", colnames(normalized))
write.csv(normalized, file.path(output_dir, "normalized_counts.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

cat("Generated Table S1 inputs in:", output_dir, "\n")
cat("Female rows:", nrow(female$result), "\n")
cat("Male rows:", nrow(male$result), "\n")
cat("Normalized-count dimensions:", nrow(normalized), "x", ncol(normalized), "\n")
cat("R session details:", file.path(output_dir, "R_sessionInfo.txt"), "\n")
