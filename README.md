# MIC60 publication analysis code

This repository contains the current code used for the final MIC60 publication analysis and the current Table S1 workflow. Superseded figure folders, old generated graphs, duplicated input files, and obsolete panel-number versions were removed from the active branch on 18 September 2026. They remain recoverable from Git history.

## Current contents

- `reviewed_analysis/`: final reviewed R, Python, and Illustrator scripts; shared figure styling; edited panel generators; workflow runner; publishing helper; and validation code.
- `tables/Table_S1_Code/`: the current Table S1 DESeq2 preparation and Excel-workbook generation scripts, with the exact archived inputs used for the workbook.
- `tables/Table_S1.xlsx`: current supplementary Table S1 output.
- `R_PACKAGES.md`: R dependencies used by the analysis and figure code.

## Reviewed figure-analysis code

Start with `reviewed_analysis/README.md`, then consult `reviewed_analysis/WORKFLOW.md` for the analysis sequence, methods, requirements, and limitations.

The repository intentionally keeps this as a code-only reviewed snapshot rather than duplicating large experimental datasets and generated figure exports. The scripts preserve their expected relative input/output paths so they can be reunited with the separately archived final source package.

## Table S1

Table S1 can be rebuilt by following `tables/Table_S1_Code/README.md`.

- `prepare_table_s1_inputs.R` reruns the female and male DESeq2 analyses and exports normalized counts.
- `build_table_s1.mjs` builds and formats the four-sheet Excel workbook.

The Table S1 scripts and archived CSV inputs were verified against the active local project before this cleanup.

## Version policy

The `main` branch contains only the current code snapshot and Table S1 workflow. Do not add alternate, previous, or manually renamed figure versions to `main`; use Git history or a clearly named archival branch when an older state must be retained.

## License

No reuse license was supplied with the project. Reuse permission should be clarified by the authors before code or data are reused.
