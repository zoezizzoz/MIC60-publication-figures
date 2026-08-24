# Table S1 code and source data

This folder contains the code and data used to produce `../Table_S1.xlsx`.
Table S1 is an output table rather than a plotted figure, but its analysis and
workbook-assembly code are included for reproducibility.

## Contents

- `Code/prepare_table_s1_inputs.R`: reruns the female and male DESeq2 analyses
  from the 12-library raw count matrix and exports normalized counts.
- `Code/build_table_s1.mjs`: assembles the four-sheet Excel workbook and applies
  the publication formatting, formulas, tables, and classification thresholds.
- `Original_Data/count_matrix_symbol.csv`: raw gene-by-library count matrix.
- `Original_Data/DEG_CSF_vs_WRF.csv`: female DESeq2 results used in the
  publication workbook.
- `Original_Data/DEG_CSM_vs_WRM.csv`: male DESeq2 results used in the
  publication workbook.
- `Original_Data/normalized_counts.csv`: normalized counts used in the
  publication workbook.
- `Rebuilt_Output/`: regenerated files and verification previews.

## Rebuild from the archived analysis products

From this `Table_S1_Code` directory, run:

```bash
node Code/build_table_s1.mjs
```

This writes `Rebuilt_Output/Table_S1.xlsx`.

## Rerun DESeq2 from the raw count matrix, then rebuild

The R analysis requires `DESeq2`. The workbook builder requires Node.js and
`@oai/artifact-tool`.

```bash
Rscript Code/prepare_table_s1_inputs.R
TABLE_S1_INPUT_DIR=Rebuilt_Output/generated_inputs node Code/build_table_s1.mjs
```

The analysis is sex-specific (`design = ~ genotype`), uses WT as the reference,
retains all 17,173 genes without a manual count prefilter, and includes all 12
libraries. Workbook direction labels use adjusted P < 0.05 and absolute log2
fold change >= 0.58. The script also records `R_sessionInfo.txt`; retain that
file when rerunning because DESeq2 and R version differences can produce small
numerical differences. The archived CSV inputs are the exact values assembled
into the publication workbook.
