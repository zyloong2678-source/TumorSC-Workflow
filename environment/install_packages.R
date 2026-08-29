# Install dependencies. Pass --core to install only the synthetic core workflow.
args <- commandArgs(trailingOnly = TRUE)
core_only <- "--core" %in% args
cran <- c("Seurat", "yaml", "here", "ggplot2", "ggrepel", "patchwork", "svglite", "ragg", "pheatmap", "dplyr", "tidyr")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

if (!core_only) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  bioc <- c("fgsea", "msigdbr", "UCell", "SingleCellExperiment", "slingshot", "decoupleR", "dorothea")
  missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message(if (core_only) "Core installation complete." else "Installation complete. Harmony and CellChat remain optional; follow their current official installation instructions before enabling those modules.")
