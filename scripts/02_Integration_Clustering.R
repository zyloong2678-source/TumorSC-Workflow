# Purpose: Perform PCA, optional batch-aware integration, clustering, and UMAP.
# Input:   results/objects/01_preprocessed.rds
# Output:  results/objects/02_clustered.rds and UMAP figures

source("functions/workflow_utils.R")
source("functions/plotting_functions.R")
cfg <- read_config()
assert_packages(c("Seurat", "ggplot2"))
object <- load_seurat("01_preprocessed", cfg)
validate_metadata(object, cfg$integration$batch_variable)

assay <- if (cfg$preprocessing$method == "SCTransform") "SCT" else "RNA"
Seurat::DefaultAssay(object) <- assay
if (assay == "RNA") object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE)
object <- Seurat::RunPCA(object, npcs = cfg$integration$n_pcs, verbose = FALSE)
dims <- seq_len(min(cfg$integration$dimensions, ncol(Seurat::Embeddings(object, "pca"))))
reduction <- "pca"

if (cfg$integration$method == "harmony") {
  assert_packages("harmony")
  object <- harmony::RunHarmony(object, group.by.vars = cfg$integration$batch_variable, reduction = "pca")
  reduction <- "harmony"
} else if (cfg$integration$method == "seurat_cca") {
  stop("seurat_cca requires sample-wise anchor integration before merging. Use harmony, none, or adapt this module to the assay design.", call. = FALSE)
} else if (cfg$integration$method != "none") stop("Unknown integration method.", call. = FALSE)

object <- Seurat::FindNeighbors(object, reduction = reduction, dims = dims, verbose = FALSE)
object <- Seurat::FindClusters(object, resolution = cfg$integration$resolution, random.seed = cfg$project$seed, verbose = FALSE)
object <- Seurat::RunUMAP(object, reduction = reduction, dims = dims, min.dist = cfg$integration$umap_min_dist,
                          seed.use = cfg$project$seed, verbose = FALSE)
p1 <- Seurat::DimPlot(object, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggplot2::ggtitle("Clusters")
p2 <- Seurat::DimPlot(object, group.by = cfg$integration$batch_variable) + ggplot2::ggtitle("Samples")
p <- (p1 | p2) & theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
save_figure(p, "02_umap_clusters_samples", cfg, width_mm = 180, height_mm = 85)
save_seurat(object, "02_clustered", cfg)
save_session_info("02_Integration_Clustering", cfg)

