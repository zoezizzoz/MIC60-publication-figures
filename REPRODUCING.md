# Reproducing the graphs

## Setup

Use Python 3 and R with the dependencies in `reviewed_analysis/requirements.txt` and `R_PACKAGES.md`. `environment/` records the preparation runtime and available package versions. Bioconductor packages require a compatible R/Bioconductor installation. The original Table S1 formatter additionally requires Node.js and `@oai/artifact-tool`; no replacement formatter or new workbook layout was introduced in this update.

The core table/checksum verifier uses only Python's standard library. The Figure 4 data extractor uses `openpyxl`. Optional PDF verification needs `pdfplumber` and Poppler (`pdftoppm` on PATH). Current Figure 4 PDF output uses the macOS Quartz device and Arial. Shared plot helpers choose the available PDF device; graphical appearance can vary across systems or font installations.

## Exact deposited plotting values

`plotted_data/` is the frozen table set for review. Every table has a source and selection in `INDEX.csv`. A copied table is byte-identical to its recorded source; filtered tables preserve values as text and only select rows (for example non-missing volcano adjusted P or included sleep observations). Missing P values remain missing.

Do not replace these tables silently with a new statistical fit. The archived DESeq2 results supplied the plotted fold changes/P values. A fresh fit during preparation preserved all DEG classifications and non-missing-P counts but showed numerical differences in low-count/untested genes; `reviewed_analysis/results/rnaseq/archive_comparison.csv` records those differences. `plotted_data/Fig2D/` retains the recovered normalized plotting values independently of fresh fits.

## Commands

From the repository root:

```sh
python3 tools/verify_deposit.py --checksums  # initial snapshot integrity
python3 reproduce.py --list
python3 reproduce.py --step fig4b          # pooled normalized MTT wells, mean +/- SEM
python3 reproduce.py --step recent-panels  # current 3N, 4A/B, S4A
python3 reproduce.py --step current-layout # September 23 GO and sleep/activity
python3 reproduce.py --step figs3          # current 11-gene AMPK panel and supporting survey; offline
python3 reproduce.py --step networks       # recovered ring and STRING inputs
python3 reproduce.py --step reviewed       # full reviewed analysis, including DESeq2 refit
```

Choose `--rscript /path/to/Rscript` if Rscript is not on PATH. `reproduce.py` records exit codes and output in ignored `verification/local_runs/`. It stops at the first failure. `reviewed` uses the existing nine-stage `run_all.py`; the `--publish` flag is deliberately not passed, so Illustrator publication files are never overwritten.

Individual generators are listed in `plotted_data/INDEX.csv`. Script directory names preserve old panel numbering: `Figure_2E_Individual_Genes` corresponds to current Fig2D; old `Figure_3O_TMRM` corresponds to current Fig3N; old `Figure_4C_TIMELESS` corresponds to current Fig4A; `Figure_S4B_Transfection_Efficiency` corresponds to current FigS4C. Use the explicit index rather than inferring panel letters from historical folders.

Rerunning changes generated outputs and session records; full-snapshot checksum verification should be performed before rerunning. The default verifier continues to check the frozen plotted tables and key numerical source relationships after a rerun.

## Validation performed during preparation

The reviewed RNA-seq fit, sleep/activity processing, assay audit, review panels, GO, RNA-seq displays and eight full-panel generators were run successfully from the deposit paths. Ring and STRING-network scripts, the Fig3N generator, Figure 4 extraction/vector generators, the S4A renderer, TEM field displays, additional assays, compact PCA and female heatmap were also exercised. See `verification/PREPARATION_REPORT.json` for actual recorded outcomes and any limits. This validates executable paths and recoverable numerical inputs, not pixel identity to every manually edited Illustrator panel.

## Final native layout provenance

[provenance/final_layout_2026-09-23/README.md](provenance/final_layout_2026-09-23/README.md) describes the archived code and geometry used for the final eight-figure common-width pass. These state-specific Illustrator scripts are excluded from reproduce.py and retain the original absolute paths and object UUIDs. Replaying the native edits requires the separately retained baseline AI assemblies and backups; the public code/data deposit does not contain those assemblies. Recorded verification results and exact before/after hashes are supplied for audit.

## Figure S3 update, 25 September 2026

The `figs3` step rebuilds from 127 hash-checked frozen inputs: the complete KEGG hsa04152 entry/KGML, DIOPT 9.1 responses, identifier mapping and archived female DESeq2 table. It requires Python with `pypdfium2` and Pillow, plus R with `grid`; PDF output uses macOS Quartz and Arial. It performs no network requests or DESeq2 refit. The current figure has 81 gene-by-panel entries, including 11 AMPK genes selected from the explicit publication-evidence table. The full 131-gene survey remains in supporting data with its 13 display groups and a standalone plot. An individual audit of all 89 original entries documents seven display omissions and broader functional headings. Sesn and Atg8a are supported independently of the full orthology survey. The original lists and every numeric RNA-seq value are preserved. See the Figure S3 package README for selection criteria, missing-value interpretation, source references and display limitations.

Both the `figs3` command and the reviewed full-panel workflow invoke the same current S3 rebuild, including the publication-based gene selection and supporting survey. The S3 step was rebuilt successfully from its installed deposit path. Checks verify all six numeric DESeq2 columns for the 81 displayed entries, the 89 original audit decisions, the 70 retained entries in the other five panels, and the complete survey with 200 retained mapping edges and 13 groups. A fresh temporary copy also reproduced all nine output CSVs byte for byte and all three PNGs pixel for pixel; see `verification/S3_CURRENT_REBUILD.json`. Included PDF/PNG outputs are deliberate exceptions to the general generated-artifact exclusion.

## Figure 4B update, 25 September 2026

The `fig4b` step rebuilds the pooled-well MTT panel. All 78 source-cell absorbances and matched-control normalizations are checked against the raw workbook before plotting. The 10 pooled means and s.e.m. values are independently checked by `tools/verify_deposit.py`. Both prior MTT generator entrypoints now delegate to this same implementation; the full assay audit also produces pooled-well summaries. No inferential tests are introduced. Historical culture-preparation summaries remain labelled as provenance, and the superseded active biological-summary files were removed.
