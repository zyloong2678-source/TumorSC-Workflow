# Purpose: Run cluster/cell-state differential expression and export tables/plots.
# Input:   results/objects/04_malignancy_candidates.rds and config contrasts
# Output:  marker tables, volcano plot, and heatmap

source("functions/workflow_utils.R")
source("functions/plotting_functions.R")
cfg <- read_config()
assert_packages(c("Seurat", "ggplot2", "ggrepel"))
object <- load_seurat("04_malignancy_candidates", cfg)
de <- cfg$differential_expression
validate_metadata(object, de$group_by)
Seurat::Idents(object) <- object[[de$group_by, drop = TRUE]]

if (!is.null(de$ident_1)) {
  markers <- Seurat::FindMarkers(object, ident.1 = de$ident_1, ident.2 = de$ident_2,
    test.use = de$test_use, min.pct = de$min_pct, logfc.threshold = de$logfc_threshold)
  markers$gene <- rownames(markers)
  markers$contrast <- paste(de$ident_1, "vs", de$ident_2 %||% "all_other")
  write_table(markers, "05_pairwise_differential_expression.csv", cfg)
  save_figure(plot_volcano(markers, cfg), "05_volcano", cfg)
} else {
  markers <- Seurat::FindAllMarkers(object, only.pos = TRUE, test.use = de$test_use,
    min.pct = de$min_pct, logfc.threshold = de$logfc_threshold)
  write_table(markers, "05_all_group_markers.csv", cfg)
}

fc <- if ("avg_log2FC" %in% names(markers)) "avg_log2FC" else "avg_logFC"
group_col <- intersect(c("cluster", "ident"), names(markers))[1]
if (!is.na(group_col)) {
  top <- do.call(rbind, lapply(split(markers, markers[[group_col]]), function(x) head(x[order(x[[fc]], decreasing = TRUE), ], 5)))
  if (nrow(top)) {
    heat <- Seurat::DoHeatmap(object, features = unique(top$gene), group.by = de$group_by, raster = TRUE) +
      theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
    save_figure(heat, "05_top_marker_heatmap", cfg, width_mm = 150, height_mm = 120)
  }
}
save_session_info("05_Differential_Expression", cfg)

