# R package dependencies

The figure scripts use R 4.4.x and the packages below. Bioconductor packages
should be installed with `BiocManager`; the remaining packages are available
from CRAN or their project repositories.

## Bioconductor

- `AnnotationDbi`
- `DESeq2`
- `KEGGREST`
- `clusterProfiler`
- `enrichplot`
- `org.Dm.eg.db`

## CRAN and other R packages

- `cowplot`
- `damr`
- `data.table`
- `geomtextpath`
- `ggetho`
- `ggforce`
- `ggplot2`
- `ggraph`
- `ggrepel`
- `igraph`
- `patchwork`
- `pheatmap`
- `rbioapi`
- `readxl`
- `Rtsne`
- `scales`
- `sleepr`
- `tidygraph`
- `tidyverse`
- `uwot`

## Additional dependencies exercised by the restored/current generators

- `BiocManager` (installation helper)
- `ggnewscale`
- `ggtext`
- `gtable`
- `jsonlite`
- `systemfonts`
- `grid` and `stats` (bundled with R)

See `environment/R_installed_packages.csv` and R session records for actual versions. The historical survival reference additionally calls `survival` and `survminer`; it is not a certified current analysis.
