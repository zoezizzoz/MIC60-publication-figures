# Code for MIC60 publication figures

This repository contains figure-generation code and the associated inputs for
the manuscript **“Oxidation-Resistant dMIC60 Drives Mitochondrial–Nuclear
Crosstalk and Promotes Resilience.”** Folder names match the final figure and
panel labels in the assembled publication artwork.

The final RNA-seq supplementary table is included as `tables/Table_S1.xlsx`.
Its DESeq2 analysis and workbook-generation code are in
`tables/Table_S1_Code/`, together with the raw count matrix and archived inputs.

## Repository structure

Each panel folder is under `figures/` and contains the applicable subset of:

- `Code/`: R scripts and shared plotting style.
- `Original_Data/`: input data read by the scripts.
- `Supporting_Data/`: processed values or statistics used for verification.
- `Final_Graphs/`: the graph exported for assembly into the publication figure.
- `Source_Project/`: a non-script source project when the panel was made in
  GraphPad Prism.
- `tables/`: publication supplementary tables.

Table S1 can be rebuilt by following `tables/Table_S1_Code/README.md`.

The final code-only review workflow is available in `reviewed_analysis/`. It
contains the shared analysis scripts, edited panel generators, workflow runner,
publishing helper, and validation code without duplicating the underlying data
or generated outputs.

## Publication panel map

| Folder | Publication panel |
|---|---|
| `Figure_1C_Female_Volcano` | Female RNA-seq volcano plot |
| `Figure_2B_Program_Ring` | Mitochondrial-program ring |
| `Figure_2C_STRING_Network` | STRING interaction network |
| `Figure_2D_Selected_Genes` | Normalized counts for selected genes |
| `Figure_3A_Sleep_Profile` | Sleep profile |
| `Figure_3B_Total_Sleep` | Total sleep |
| `Figure_3C_Day_Night_Sleep` | Day/night sleep |
| `Figure_3D_Activity_Profile` | Activity profile |
| `Figure_3E_Total_Activity` | Total activity |
| `Figure_3F_Day_Night_Activity` | Day/night activity |
| `Figure_3H_Performance_Index` | Negative-geotaxis performance index |
| `Figure_3J_Mitochondrial_Perimeter` | Mitochondrial perimeter |
| `Figure_3K_Mitochondrial_Area` | Mitochondrial area |
| `Figure_3L_Mitochondrial_Aspect_Ratio` | Mitochondrial aspect ratio |
| `Figure_3N_TMRM_MTG_Ratio` | TMRM/MitoTracker Green ratio |
| `Figure_4A_TIMELESS` | TIMELESS/DAPI quantification |
| `Figure_4B_MTT_Viability` | MTT viability dose response |
| `Figure_S1B_Western_Blot` | dMIC60-Myc immunoblot quantification |
| `Figure_S1C_Male_Volcano` | Male RNA-seq volcano plot |
| `Figure_S2_KEGG_GSEA` | KEGG gene-set enrichment analysis |
| `Figure_S3_Selected_Gene_Modules` | Selected gene-module expression |
| `Figure_S4A_Stress_Response` | Stress-response gene sets |
| `Figure_S4C_Transfection_Efficiency` | Transfection efficiency |

Figure 3H now includes an R rebuild and statistical-audit export; its earlier
`.pzfx` project remains included for provenance. Panels 1A, 3I, 3M, 4C, S1A,
and S3B are schematics, representative images, a workflow graphic, or an
immunoblot image and therefore have no graph-generation code in this repository.

## Running the scripts

Run each R script with `Rscript figures/<panel>/Code/<script>.R`. Scripts resolve
their inputs relative to their own panel folder and write regenerated files to
`Rebuilt_Output/`; scripts that maintain the publication export also refresh
`Final_Graphs/`. RNA-seq scripts require the packages in `R_PACKAGES.md`; the
STRING panel and a first-time KEGG analysis may require access to their public
annotation services. Checked-in annotation snapshots keep the final S2 plotting
step offline-reproducible.

The repository does not currently include a software license. Reuse permission
should therefore be clarified by the authors before code or data are reused.
