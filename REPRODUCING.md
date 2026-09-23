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
python3 reproduce.py --step recent-panels  # current 3N, 4A/B, S4A
python3 reproduce.py --step current-layout # September 23 GO and sleep/activity
python3 reproduce.py --step networks       # recovered ring and STRING inputs
python3 reproduce.py --step reviewed       # full reviewed analysis, including DESeq2 refit
```

Choose `--rscript /path/to/Rscript` if Rscript is not on PATH. `reproduce.py` records exit codes and output in ignored `verification/local_runs/`. It stops at the first failure. `reviewed` uses the existing nine-stage `run_all.py`; the `--publish` flag is deliberately not passed, so Illustrator publication files are never overwritten.

Individual generators are listed in `plotted_data/INDEX.csv`. Script directory names preserve old panel numbering: `Figure_2E_Individual_Genes` corresponds to current Fig2D; old `Figure_3O_TMRM` corresponds to current Fig3N; old `Figure_4C_TIMELESS` corresponds to current Fig4A; `Figure_S4B_Transfection_Efficiency` corresponds to current FigS4C. Use the explicit index rather than inferring panel letters from historical folders.

Rerunning changes generated outputs and session records; full-snapshot checksum verification should be performed before rerunning. The default verifier continues to check the frozen plotted tables and key numerical source relationships after a rerun.

## Validation performed during preparation

The reviewed RNA-seq fit, sleep/activity processing, assay audit, review panels, GO, RNA-seq displays and eight full-panel generators were run successfully from the deposit paths. Ring and STRING-network scripts, the Fig3N generator, Figure 4 extraction/vector generators, the S4A renderer, TEM field displays, additional assays, compact PCA and female heatmap were also exercised. See `verification/PREPARATION_REPORT.json` for actual recorded outcomes and any limits. This validates executable paths and recoverable numerical inputs, not pixel identity to every manually edited Illustrator panel.
