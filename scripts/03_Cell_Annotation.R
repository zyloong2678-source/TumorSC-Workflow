# Purpose: Visualize canonical markers and apply transparent two-level annotations.
# Input:   results/objects/02_clustered.rds; user-reviewed mappings below
# Output:  results/objects/03_annotated.rds and marker plots

source("functions/workflow_utils.R")
source("functions/plotting_functions.R")
cfg <- read_config()
assert_packages(c("Seurat", "ggplot2"))
object <- load_seurat("02_clustered", cfg)

# Fill these named vectors after reviewing marker expression. Example:
# level1_map <- c(`0` = "T_NK", `1` = "Myeloid")
level1_map <- c()
level2_map <- c()

markers1 <- available_features(object, cfg$annotation$markers)
if (!length(markers1)) stop("No configured Level 1 markers occur in the object.", call. = FALSE)
dot1 <- Seurat::DotPlot(object, features = cfg$annotation$markers, group.by = "seurat_clusters") +
  ggplot2::coord_flip() + theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
save_figure(dot1, "03_level1_marker_dotplot", cfg, width_mm = 150, height_mm = 120)

feature_subset <- head(markers1, 12)
fp <- Seurat::FeaturePlot(object, features = feature_subset, order = TRUE, ncol = 4,
                          pt.size = cfg$plotting$point_size) &
  theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
save_figure(fp, "03_canonical_marker_featureplots", cfg, width_mm = 180, height_mm = 160)

clusters <- as.character(object$seurat_clusters)
if (length(level1_map)) {
  unmapped <- setdiff(unique(clusters), names(level1_map))
  if (length(unmapped)) stop("Level 1 mapping missing clusters: ", paste(unmapped, collapse = ", "), call. = FALSE)
  object[[cfg$annotation$level1_column]] <- unname(level1_map[clusters])
} else {
  object[[cfg$annotation$level1_column]] <- paste0("Unreviewed_cluster_", clusters)
  warning("No Level 1 map supplied. Labels remain explicitly unreviewed.")
}
if (length(level2_map)) {
  object[[cfg$annotation$level2_column]] <- unname(level2_map[clusters])
} else object[[cfg$annotation$level2_column]] <- "Not_assigned"

p <- Seurat::DimPlot(object, group.by = cfg$annotation$level1_column, label = TRUE, repel = TRUE) +
  ggplot2::ggtitle("Marker-informed Level 1 annotation") +
  theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
save_figure(p, "03_umap_level1_annotation", cfg)
save_seurat(object, "03_annotated", cfg)
save_session_info("03_Cell_Annotation", cfg)

