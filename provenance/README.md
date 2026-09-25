# Source provenance

The update starts from the existing repository main commit `0924a4e5256396266baa405587a617cfd20d0eed`. `source_files.json` lists direct copies still present in the update, their collection-relative origin, source hash and current deposit hash. Later verification can regenerate output tables/session records; a changed hash there is not represented as an unchanged original input. `FILE_MANIFEST.json` covers the final Git payload, including newly indexed/filtered plotting tables.

`cloud_graph_input_recovery.json` records the successful recovery of 103 graph-support/code files from the older iCloud project. Its `LOCAL_RECOVERY_REFERENCE/` paths identify copies retained in the local companion recovery folder; only the selected inputs needed by active generators were added to this Git update. `recovered_code_index.csv` catalogs older/recovery code retained locally and is not a list of current entrypoints. `portability_changes.json` records initial mechanical path/preview adjustments; the Git diff is the authoritative change record against the published baseline.

`survival_legacy/` contains original historical observations and a script whose correspondence to the current curve is unverified. No alternate recovery script is presented as an active replacement for this panel.

The author narrowed this deposit to the data used in plotted graphs. Sequencing reads, acquisition image files and the broad raw-data archive are not included. Source paths in historical metadata are provenance, not required runtime paths.

`final_layout_2026-09-23/` adds exact copies of the final Illustrator formatting code, layout plans, and recorded checks after the plotting-data snapshot. Its source index verifies the copied files, and baseline/final native hashes identify the separate local assemblies. Absolute source paths are preserved as historical provenance; these scripts are not portable entrypoints and are not executed by reproduce.py.
