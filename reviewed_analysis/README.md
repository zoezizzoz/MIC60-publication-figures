# Reviewed analysis code

This directory contains the code-only snapshot used for the final MIC60 figure review and regeneration workflow. It intentionally excludes manuscript files, generated figures, review reports, and duplicated source datasets.

## Contents

- `code/`: shared analysis, audit, plotting, styling, and Illustrator-notation scripts.
- `figure_sources/`: the edited panel generators called by `code/08_full_figure_panels.py`.
- `run_all.py`: workflow runner.
- `publish_figures.py`: controlled publication-output copier.
- `validate_outputs.py`: output and numerical consistency checks.
- `WORKFLOW.md`: detailed methods, requirements, limitations, and execution order.

The scripts expect the corresponding data files and output folders described in `WORKFLOW.md`. Figure-specific source code and inputs are also retained under `../figures/`.

## Table S1

Table S1 code is maintained separately at `../tables/Table_S1_Code/`:

- `Code/prepare_table_s1_inputs.R` reruns the sex-specific DESeq2 analyses and normalized-count export.
- `Code/build_table_s1.mjs` builds and formats the four-sheet Excel workbook.

See `../tables/Table_S1_Code/README.md` for the rebuild commands and thresholds.
