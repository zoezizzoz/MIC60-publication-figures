# Figure 4A/B vector rebuild

These plots are regenerated from the latest recoverable analytical code, original measurement workbooks, and later September editing records. They are not image traces. Older code snapshots are retained in the local recovery companion, outside the current Git update.

## Inputs and calculations

- **4A TIMELESS/DAPI:** all 64 values (16 per condition) from `Fig4A.xlsx`, exactly matching the recovered plotting CSV. A continuous baseline, grouped treatment pairs, boxplot medians/IQR with 1.5-IQR whiskers, caps, and every observation are retained. The later September audit states that legacy P=.0477/.0218 could not be reproduced and were removed. Those unsupported annotations are not restored. Biological replicate identifiers are unavailable, so this panel remains descriptive.
- **4B MTT:** all 78 technical wells from the July 15 and July 28 raw sheets. Both recovered workbooks contain identical cell contents. Blank-corrected absorbance is normalized to the same experiment/genotype 0 mM control mean. Lines show the equally weighted mean of experiment means. Capped error bars show between-experiment SD where there are two experiments; 5 and 10 mM have one experiment and no SD bar. Small gray/black points show wells. Colored circles are experiment 1 (July 15); triangles are experiment 2 (July 28). The archived descriptive policy is retained. July 15's control is labelled DMSO; the source does not establish vehicle equivalence across dates.

## Layout and typography

Both PDFs are exported at the dimensions of the current placed raster graphs, for 100% import without font scaling. Main type is Arial 8 pt; chemical subscripts use proportional 5.6 pt text. Figure 4A has four equal 28 pt visible-width boxes and one continuous baseline. Figure 4B's two-row legend is inside the upper center, following the recorded September placement. Its long y-axis label wraps onto two lines to remain readable at the placed size.

Colors: WT #5AB4E5, CS #CB78A8, WT +H2O2 #47599E, CS +H2O2 #891740. No data values are changed to alter appearance.

## Reproduction

1. Run `Code/prepare_plot_data.py` with Python and openpyxl to read the immutable source copies and write traceable CSVs.
2. Run each panel's `Code/generate_*_vector.R` with R. Dependencies: ggplot2, jsonlite, grid, gtable; MTT also uses readxl. On this Mac, the Quartz PDF device embeds Arial without the missing X11/Cairo dependency.
3. Run `Code/verify_and_render_pdfs.py` with Python, pdfplumber and pdftoppm to verify zero raster image objects, retained text, all content inside page bounds, and 28 pt TIMELESS boxes; render the PDF previews.
4. Native Illustrator assembly files and recovery insertion helpers are kept outside this code/data update.

`Documentation/DATA_VERIFICATION.json` records raw-data provenance and checks. `Documentation/Relevant_Historical_Records.json` preserves the later P-value, spacing, and legend decisions. `Documentation/PDF_VECTOR_QA.json` records the vector export checks. Original source workbooks and earlier figures are preserved.
