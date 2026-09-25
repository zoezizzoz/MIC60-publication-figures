# Current Figure S3 code

Run `python3 Code/rebuild_FigS3.py` from this package, or `python3 reproduce.py --step figs3` from the repository root. Both build the current six-panel figure with 11 publication-supported AMPK genes, plus its supporting survey. The reviewed full-panel workflow uses the same rebuild.

| File | Purpose |
| --- | --- |
| `rebuild_FigS3.py` | Complete offline rebuild and validation; the main entrypoint. |
| `prepare_diopt.py` | Parse the frozen orthology records. |
| `assemble_ampk_data.py` | Reconstruct the complete supporting survey from archived expression values. |
| `group_ampk_components.py` | Assign the supporting survey's documented human-source display groups. |
| `generate_FigS3_selected_modules.R` | Render the current six-panel figure and standalone 11-gene AMPK panel. |
| `generate_AMPK_full_survey.R` | Render only the complete 131-gene supporting survey. |
| `export_png.py` | Export the vector PDFs as 600-dpi PNGs. |

The preparation scripts are rebuild dependencies, not alternative versions of the current figure. Historical module membership is a frozen input. S3_original_gene_display_decisions.csv records all decisions, and AMPK_literature_selection.csv records the focused AMPK selection. The rebuild verifies unchanged archived numeric values for every retained or restored gene. The complete 131-gene survey is retained intentionally as supporting data. Superseded figure-generation branches and the unused style helper are excluded from this folder.
