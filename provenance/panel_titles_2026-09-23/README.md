# Panel title consistency

Shared sex and null-background context belongs in figure legends. Preserve headings that distinguish corresponding female/male analyses, up/down gene sets, gene modules, treatments, and sample identities. Preserve null-background labels within the experimental schematic and blot sample labels.

Removed the redundant dMIC60-null female heading from FigS4A and MIC60-null HeLa Cells heading from Fig4. Current Figure 1C, Figure 2A, Figure S1C/E, and Figure S2 Female/Male labels remain. The Figure 4 and S4 manuscript legends already identify the relevant cell system or sex. No manuscript edit was required.

All eight live figures were inspected. Saved PDF raster comparisons confirm that only the two removed title areas changed; embedded image bytes and all other text/positions/font attributes and artboards were preserved. Backups of the changed masters are in Before_AI locally. Scientific data and statistics were not modified.

The baseline for this title pass follows Axis_Label_Spacing_2026-09-23 and the user's live Figure 3 M/N edits. The existing repository provenance/final_layout_2026-09-23 records an earlier Common_Width stage; its hashes must not be interpreted as current master hashes. The intervening native axis-spacing package and native Illustrator baselines remain local. apply.jsx is a state-specific archival script, requiring its exact native UUIDs/content, and is not run by reproduce.py.

The active S4A renderer now omits its redundant heading; generated SVG comparison shows only that title element removed. Both deposited GO renderer copies now use Female/Male headings to match the existing native panels.
