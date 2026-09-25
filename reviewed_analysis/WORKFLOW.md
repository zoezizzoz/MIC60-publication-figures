# Reviewed analysis workflow

This is the reviewed analysis workflow retained from the September 18 baseline, with recoverable graph inputs restored by the September 23 update. Current plotted tables and newer panel/layout generators are indexed in ../plotted_data/INDEX.csv. Illustrator assemblies, acquisition images, sequencing reads and superseded recovery alternatives are maintained separately.

## Workflow entry point

The recoverable inputs are supplied under `data/` and `figure_sources/`. To refit and regenerate the reviewed analyses, run:

```sh
python3 run_all.py --rscript "$(command -v Rscript)"
```

The runner executes the stages in dependency order, stops on errors, and records logs. Use `--step` with `rnaseq`, `sleep`, `assays`, `panels`, `go`, `additional_assays`, `rnaseq_displays`, `full_panels`, or `schematics` to run one stage.

Fresh generated panels are not copied over publication exports unless `--publish` is supplied. This protects Illustrator-refined typography and layout.

## Analysis stages

1. `code/01_rnaseq.R`: sex-specific DESeq2 analysis, normalized counts, and differential-expression tables using all 12 libraries.
2. `code/02_sleep_activity.R`: sleep/activity processing, exclusion audit, phase summaries, and sensitivity comparisons.
3. `code/03_assay_audit.py`: numerical audits for MTT, TIMELESS, TEM, locomotion, and stress-set results.
4. `code/04_review_panels.R`: descriptive assay panels and western-blot verification.
5. `code/05_go_pathways.R`: female and male GO over-representation analysis and lollipop plots using the shared style and centralized GO settings.
6. `code/06_additional_assays.R`: additional locomotion and TIMELESS summaries.
7. `code/07_rnaseq_displays.R`: PCA/MDS, heatmaps, and RNA-seq display matrices.
8. `code/08_full_figure_panels.py`: invokes the current edited generators retained under `figure_sources/*/Code/`.
9. `code/09_schematics.R`: vector schematic replacements.

`publish_figures.py` performs the controlled publication-output copy. `validate_outputs.py` checks manifests, selected genes, numerical summaries, labels, and expected outputs after a complete run with the archived inputs restored.

## Current figure generators

The `figure_sources/` tree contains the reviewed generators and restored supporting inputs. It includes the female/male volcano plots, program ring, STRING network, individual-gene panels, locomotion, TEM, TMRM, MTT, TIMELESS, western-blot quantification, KEGG/GSEA, selected-gene modules, stress response, and transfection efficiency.

## Analysis conventions

- Differential-expression classification: adjusted P < 0.05 and absolute log2 fold change >= 0.58.
- GO analysis: non-missing adjusted-P genes as the sex-specific universe, gene-set size 10-500, and separate BH correction within ontology and direction.
- Genotype labels use the hyphenated forms `dMIC60-WT`, `dMIC60-CS`, `UAS-dMIC60-WT`, and `UAS-dMIC60-CS`.
- Final publication exports may contain manual Illustrator typography/layout refinements not recreated by the analysis scripts alone.

## Dependencies

Python requirements are listed in `requirements.txt`. R dependencies are summarized in the repository-level `R_PACKAGES.md`; individual regenerated outputs should retain their R `sessionInfo()` records alongside the generated outputs.

## Scope

This code supports the downstream analyses for which readable inputs were available. It does not reconstruct upstream FASTQ processing, unavailable raw-image segmentation, or the current survival curve from historical records with unresolved missing-time/censoring semantics. No missing source data or specimen metadata were inferred.
