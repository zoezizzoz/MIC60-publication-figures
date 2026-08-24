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

## Publication panel map

| Folder | Publication panel |
|---|---|
| `Figure_1B_RNAseq_QC` | PCA and MDS |
| `Figure_1C_Female_Volcano` | Female RNA-seq volcano plot |
| `Figure_1D_Female_Heatmap` | Female RNA-seq heatmap |
| `Figure_2A_GO_Enrichment` | Female GO enrichment |
| `Figure_2B_Program_Ring` | Mitochondrial-program ring |
| `Figure_2C_STRING_Network` | STRING interaction network |
| `Figure_2D_Selected_Genes` | Normalized counts for selected genes |
| `Figure_3A_Sleep_Profile` | Sleep profile |
| `Figure_3B_Total_Sleep` | Total sleep |
| `Figure_3C_Day_Night_Sleep` | Day/night sleep |
| `Figure_3D_Activity_Profile` | Activity profile |
| `Figure_3E_Total_Activity` | Total activity |
| `Figure_3F_Day_Night_Activity` | Day/night activity |
| `Figure_3G_Survival` | Kaplan–Meier survival |
| `Figure_3H_Performance_Index` | Negative-geotaxis performance index |
| `Figure_3J_Mitochondrial_Perimeter` | Mitochondrial perimeter |
| `Figure_3K_Mitochondrial_Area` | Mitochondrial area |
| `Figure_3L_Mitochondrial_Aspect_Ratio` | Mitochondrial aspect ratio |
| `Figure_3N_TMRM_MTG_Ratio` | TMRM/MitoTracker Green ratio |
| `Figure_4A_TIMELESS` | TIMELESS/DAPI quantification |
| `Figure_4B_MTT_Viability` | MTT viability dose response |
| `Figure_S1B_Western_Blot` | dMIC60-Myc immunoblot quantification |
| `Figure_S1C_Male_Volcano` | Male RNA-seq volcano plot |
| `Figure_S1D_Male_Heatmap` | Male RNA-seq heatmap |
| `Figure_S1E_Male_GO_Enrichment` | Male GO enrichment |
| `Figure_S2_Targeted_Modules` | Targeted gene modules |
| `Figure_S3A_Stress_Response` | Stress-response gene sets |
| `Figure_S3C_Transfection_Efficiency` | Transfection efficiency |

Figure 3H was created in GraphPad Prism; its `.pzfx` source project is included
instead of R/Python code. Panels 1A, 3I, 3M, 4C, S1A, and S3B are schematics,
representative images, a workflow graphic, or an immunoblot image and therefore
have no graph-generation code in this repository.

## Running the scripts

Run each R script from its panel directory unless its header says otherwise.
The scripts use relative paths and write rebuilt files into their panel's output
directory. RNA-seq scripts require the packages listed in `R_PACKAGES.md`; the
STRING panel also requires access to the remote STRING service.

The repository does not currently include a software license. Reuse permission
should therefore be clarified by the authors before code or data are reused.
