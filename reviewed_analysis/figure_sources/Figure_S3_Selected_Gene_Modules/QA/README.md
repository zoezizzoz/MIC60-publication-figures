# Current Figure S3 checks

- `rebuild_validation.json`: complete current rebuild, including 11 displayed AMPK genes, 81 gene-by-panel entries, 127 frozen inputs and unchanged archived expression values.
- `compact_publication_validation.json`: current figure typography, annotations and saved native-artwork inspection.
- `membership_validation.json`: construction of the complete 131-gene supporting survey; comparisons to the original AMPK list refer to that survey, not the current 11-gene display.
- `grouping_validation.json`: the supporting survey's 13 human-source display groups.

The repository-level `verification/S3_CURRENT_REBUILD.json` records a fresh rebuild of the current code: all nine output CSVs are byte-identical and all three plot PNGs are pixel-identical.

Superseded layout reports and duplicate reports have been removed. They remain recoverable from Git history.
