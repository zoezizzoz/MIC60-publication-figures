# Figure 4A/B vector rebuild

These plots are regenerated from the latest recoverable analytical code, original measurement workbooks, and later September editing records. They are not image traces. Older code snapshots are retained in the local recovery companion, outside the current Git update.

## Inputs and calculations

- **4A TIMELESS/DAPI:** all 64 values (16 per condition) from `Fig4A.xlsx`, exactly matching the recovered plotting CSV. A continuous baseline, grouped treatment pairs, boxplot medians/IQR with 1.5-IQR whiskers, caps, and every observation are retained. The later September audit states that legacy P=.0477/.0218 could not be reproduced and were removed. Those unsupported annotations are not restored. Biological replicate identifiers are unavailable, so this panel remains descriptive.
- **4B MTT (25 September update):** all 78 wells from the July 15 and July 28 raw sheets. Blank-corrected absorbance is normalized to the corresponding genotype's 0 mM control mean within each culture preparation. Normalized wells are then pooled within genotype and dose. Lines show their arithmetic mean; capped error bars show sample SD / sqrt(number of wells), i.e. s.e.m. Small gray points are individual wells. There are 9 wells/genotype at 0, 20 and 40 mM (3 + 6 across two preparations), and 6 wells/genotype at 5 and 10 mM (one preparation). Each well is treated as a biological replicate according to the authors' requested presentation; preparation IDs remain in the source table. No inferential tests or significance annotations are added. Both raw workbooks contain identical cell contents. July 15's control is labelled DMSO; vehicle equivalence across dates is not established by these source files.

## Layout and typography

The MTT PDF uses the current Illustrator data-region dimensions (199 x 132.9134 pt), with Arial 8 pt axis titles, 7 pt tick/legend text and 5.6 pt chemical subscripts. It is imported at 100% without font or stroke scaling. Figure 4A has four equal 28 pt visible-width boxes and one continuous baseline. Figure 4B has a single genotype legend inside the upper center; culture-mean symbols and their legend were removed. Its long y-axis label wraps onto two lines to remain readable at the placed size.

Colors: WT #5AB4E5, CS #CB78A8, WT +H2O2 #47599E, CS +H2O2 #891740. No data values are changed to alter appearance.

## Reproduction

1. Run `Code/prepare_plot_data.py` with Python and openpyxl to read the immutable source copies and write traceable CSVs.
2. Run each panel's `Code/generate_*_vector.R` with R. Dependencies: ggplot2, jsonlite, grid, gtable; MTT also uses readxl. On this Mac, the Quartz PDF device embeds Arial without the missing X11/Cairo dependency.
3. Run `Code/verify_and_render_pdfs.py` with Python, pdfplumber and pdftoppm to verify zero raster image objects, retained text, all content inside page bounds, and 28 pt TIMELESS boxes; render the PDF previews.
4. Native Illustrator assembly files are retained locally; historical insertion helpers are archived under `provenance/fig4ab_vector_rebuild_2026-09-23/` and are not used for the current rebuild.

`Documentation/DATA_VERIFICATION.json` records raw-data provenance and checks. `Documentation/Relevant_Historical_Records.json` preserves the later P-value, spacing, and legend decisions. `Documentation/PDF_VECTOR_QA.json` records the vector export checks. Original source workbooks and earlier figures are preserved.

The plotted statistics are in `Supporting_Data/mtt_pooled_well_summary.csv`. `mtt_experiment_means.csv` is retained only for source-preparation provenance and is not the plotted mean or error bar. Both older MTT entrypoints delegate to the current vector generator.
