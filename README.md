# MIC60 publication figures: code and plotted data

This update adds the recoverable graph inputs and September 22–23 figure rebuilds to the existing [MIC60-publication-figures repository](https://github.com/zoezizzoz/MIC60-publication-figures). It is based on `main` commit `0924a4e5256396266baa405587a617cfd20d0eed` (18 September 2026).

## Find a graph's data

Start with **[plotted_data/README.md](plotted_data/README.md)**. Its 40 tables cover 29 graph panels using the current panel letters. [INDEX.csv](plotted_data/INDEX.csv) maps each table to its generator, original source, observation unit, row selection and checksum.

- `reviewed_analysis/`: existing reviewed pipeline, restored inputs, shared styles, generators and numerical outputs.
- `panels/Fig3N/`: current connected-object TMRM/MTG field-mean graph; 32 fields and 22,531 contributing objects.
- `panels/Fig4AB/`: current TIMELESS/DAPI and MTT graphs, with source-cell references for the 64 measurements and 78 wells.
- `panels/FigS4A/`: 64 stress-response points, enrichment results and the current SVG renderer with color-matched labels.
- `panels/layout_2026-09-23/`: latest recovered GO and sleep/activity formatting generators, with their calculation inputs.
- `data/original_inputs/`: supplied measurement workbooks, count matrices and analysis tables. The historical name `Original_Data` does **not** imply every file is raw instrument data.
- `tables/`: existing Table S1 workflow, archived inputs and workbook.
- `provenance/`, `verification/`: recovery records, checks, dependencies and limitations.

Sequencing reads are available separately under [GEO GSE346353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE346353), confirmed public on 23 September 2026. They are excluded from this graph-data update. Microscopy/blot acquisition images and Illustrator assemblies are also excluded from the graph-data repository. No data values were invented to fill a recovery gap.

## Check and rebuild

```sh
python3 tools/verify_deposit.py
python3 tools/verify_deposit.py --checksums
python3 reproduce.py --list
python3 reproduce.py --step recent-panels
```

See [REPRODUCING.md](REPRODUCING.md) for dependencies and the reviewed workflow. Builds write into the existing `results/`, `Rebuilt_Output/`, or `Final_Graphs/` directories; they do not modify `plotted_data/`. The package was exercised with R 4.4.1 and Python 3.12.14 on macOS; exact versions are recorded in `environment/`. Some PDF generators use macOS Quartz/Arial, and Table S1's original workbook formatter uses `@oai/artifact-tool`.

## Scope and remaining gaps

**Fig3G survival is not certified as the data behind the current curve.** Historical observations and code are provided in `provenance/survival_legacy/`; missing times/censoring, historical filtering and the displayed P value still need reconciliation. [LIMITATIONS.md](LIMITATIONS.md) also records missing biological replicate mappings, fresh DESeq2 numerical differences and Illustrator-only layout refinements.

Current generators and data belong on this branch. Older alternative code and failed recovery helpers are preserved in the local companion recovery folder, rather than mixed into this update. Earlier published versions remain available in Git history.

## License

The existing repository has no reuse license. This preparation does not assign a new license or change ownership.
