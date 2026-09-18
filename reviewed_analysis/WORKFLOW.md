# Reviewed analysis workflow

Run from any working directory using Python with NumPy, pandas, and openpyxl, and R with DESeq2, damr, sleepr, data.table, ggplot2, jsonlite, clusterProfiler, org.Dm.eg.db, patchwork, and systemfonts:

```sh
python3 run_all.py
```

The validated command pattern is:

```sh
python3 run_all.py --rscript "$(command -v Rscript)"
```

The runner executes nine steps in dependency order, stops on errors, records logs, and checks that all inputs remain unchanged. Use `--step rnaseq`, `--step sleep`, `--step assays`, `--step panels`, or `--step go` to rerun one step. The panel step requires outputs from the assay step. Final Illustrator-refined exports are preserved unless `--publish` is supplied explicitly.

1. `01_rnaseq.R` fits each sex once, exports normalized counts and DE results, and compares the rerun with archived tables. All 12 libraries are retained. DESeq2's default independent-filtering alpha is 0.1 for reproduction; DEG classification still uses adjusted P < 0.05 and |log2FC| ≥ 0.58.
2. `02_sleep_activity.R` loads the nine DAM files once. It compares 0–60 h and 12–60 h windows, no IQR filtering, independent phase filtering, and paired phase filtering. Corrected panels use 12–60 h and paired phase filtering, matching the edited methods. It preserves every exclusion flag, documents five-minute acquisition intervals, and corrects the activity scale. It exports cohort summaries from recording start dates; the original `replicate` column is always 3 and does not identify independent batches. Pooled rank tests do not adjust for cohort effects.
3. `03_assay_audit.py` uses Python for extraction and R for numerical tests. It reconstructs all 78 MTT values from raw absorbance, checks them against the saved plotted workbook, reproduces TIMELESS ANOVA/Tukey values conditionally on independence, reproduces exact TEM field-level tests, enumerates the 70 locomotion rank permutations, and recalculates stress-set Fisher/BH results. It does not infer biological replication from the number of cells, fields, or wells.
4. `04_review_panels.R` creates a descriptive MTT replacement with experiment means and between-experiment SD and verifies the saved western-blot values. Transfection efficiency is summarized per experiment.

5. `05_go_pathways.R` reruns GO over-representation tests on the archived female and male DESeq2 results. It uses non-missing adjusted-P genes as the sex-specific background, DEG thresholds adjusted P < 0.05 and |log2FC| >= 0.58, gene-set size 10–500, and separate BH correction within each ontology and direction. It validates all 960 female-upregulated terms against the saved table, and plots up to eight significant terms per ontology/direction with gene counts. Full results and membership lists are retained. This independently specified rerun differs from the old assembled artwork; see `../Figures/GO_README.md`.

`data/PROVENANCE.json` records input provenance and SHA-256 hashes. `results/` holds generated numerical outputs, plots, and execution logs. Current editable Illustrator assemblies are deposited separately in `../Figures/Illustrator/`.

## Scope and limitations

This is an executable review workflow for the accessible data, not a replacement for the missing full publication pipeline. It does not reconstruct FASTQ processing, the original GO artwork selection, full GSEA, or STRING retrieval, survival analysis, or TMRM segmentation. Those steps need the original generating code and, in several cases, unavailable source files or specimen identities. See `../Review/statistical_review.md` and `../Review/figure_corrections.md`.

The original readable code is archived once per unique file hash in `../Archive/Original_Code/`, with every source path mapped in `INDEX.json`. Different style-file versions are retained separately. These archived scripts retain their original relative paths and are reference material; the supported execution entry point in this package is `run_all.py`.

New replacement panels are candidates for review and assembly. They do not overwrite the Illustrator figures. Sleep/activity P values include a six-comparison Holm sensitivity correction chosen during this review, not a retrospectively claimed prespecified analysis plan.

## Shared figure style and edited original generators

All 28 replacement panels use `code/figure_style.R`, copied and edited from the original sleep/TEM/TMRM style. Edit this one working file to keep typography, genotype colors, points, and axes consistent. `code/generate_sleep_activity_figures.R` and `code/generate_mtt_figure.R` are edited copies of the original generators. The numbered audit scripts call these generators; their earlier duplicate plotting blocks have been removed. GO had no supplied original barplot generator, so `05_go_pathways.R` uses the same shared style.

To redraw sleep/activity or MTT from existing reviewed tables, run the corresponding `generate_*.R` script with Rscript. To rerun analyses, use `run_all.py`. Publishing fresh code output over `../Figures/PDF`, `../Figures/PNG`, and the figure manifest is an explicit operation: use `run_all.py --publish` or run `python3 publish_figures.py`. See `../Review/code_and_typography_changes.md` and `../Review/code_diffs/` for source mapping and exact changes.

`code/go_analysis_config.R` centralizes all GO comparison settings; `results/go/analysis_configuration.json` records the parameters used. `../Review/go_consistency_validation.json` verifies that all four comparisons follow those rules and all 5,540 gene counts agree with their membership lists.

## Full-figure extension

6. `06_additional_assays.R`: locomotion exact rank permutation and descriptive TIMELESS ratios; original generators were unreadable.
7. `07_rnaseq_displays.R`: explicitly reconstructed design-blind rlog PCA/MDS and sex-specific top-30 heatmaps. Archived contrast statistics determine selection; one all-library rlog matrix supplies visualization only.
8. `08_full_figure_panels.py`: runs eight edited original generators under `figure_sources/*/Code` for volcanoes, individual genes, TEM, fly western-blot quantification, KEGG, S3 modules, stress sets and transfection efficiency. Each uses the shared style. Individual-gene plots reuse normalized counts instead of refitting DESeq2. The S3 checkpoint is retained as reference; it is no longer an active source of P values.
9. `09_schematics.R`: native vector replacements for Fig1A, FigS1A and Fig4C; original editable generators unavailable.

The descriptive Fig3N TMRM panel is republished from its deposited reviewed source. Its object-level segmentation is not rerun because the original raw-image segmentation workflow and field-to-fly assignments were unavailable.

Additional `--step` choices are `additional_assays`, `rnaseq_displays`, `full_panels`, and `schematics`. Run the full workflow for dependency order; `publish_figures.py` requires every configured output and fails if one is missing. It copies 29 PDF/PNG pairs and records their hashes.

Additional R packages: readxl, pheatmap, ggrepel, dplyr, tidyr, stringr, scales, grid and gridExtra (plus the package versions recorded with individual outputs). `figure_sources` retains readable source snapshots and reference scripts, including some inactive scripts whose inputs remain unavailable; only the eight generators called by step 08 are supported there. Original file names may contain old panel numbers; the publishing manifest maps them to current letters. Do not run every historical script recursively.

`Review/full_figure_source_inventory.json` records pre-edit source hashes and recovery locations. `Review/code_diffs` makes changes to original generators reviewable. Standardizing style does not imply that independent experimental units or all raw-image provenance have been recovered.

After regenerating the package, run `python3 validate_outputs.py` for manifest, figure text, heatmap selection/z-score, S3 contrast, assay-count and prior-table checks (requires pdfplumber, pandas and NumPy).
