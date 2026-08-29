# Purpose: Infer Slingshot lineages and pseudotime within a biologically chosen subset.
# Input:   results/objects/04_malignancy_candidates.rds
# Output:  trajectory object, lineage definitions, pseudotime/weight tables,
#          curve coordinates, trajectory plot, and descriptive gene dynamics
# Caution: Pseudotime ordering does not establish temporal causality.

source("functions/workflow_utils.R")
source("functions/plotting_functions.R")
cfg <- read_config()
if (!isTRUE(cfg$trajectory$enabled)) stop("Optional module disabled. Set trajectory.enabled: true after defining a defensible subset.", call. = FALSE)
assert_packages(c("slingshot", "SingleCellExperiment", "Seurat", "ggplot2"), optional = TRUE)
object <- load_seurat("04_malignancy_candidates", cfg)
validate_metadata(object, c(cfg$trajectory$subset_column, cfg$trajectory$cluster_column))
keep <- object[[cfg$trajectory$subset_column, drop = TRUE]] %in% unlist(cfg$trajectory$subset_values)
sub <- subset(object, cells = colnames(object)[keep])
if (ncol(sub) < 50) stop("Trajectory subset has fewer than 50 cells; review subset definition.", call. = FALSE)
sce <- Seurat::as.SingleCellExperiment(sub)
SingleCellExperiment::reducedDim(sce, "UMAP") <- Seurat::Embeddings(sub, "umap")
start <- cfg$trajectory$start_cluster
sce <- slingshot::slingshot(sce, clusterLabels = sub[[cfg$trajectory$cluster_column, drop = TRUE]],
                            reducedDim = "UMAP", start.clus = start)
pt <- as.data.frame(slingshot::slingPseudotime(sce))
pt$cell <- rownames(pt)
write_table(pt, "08_pseudotime.csv", cfg)
weights <- as.data.frame(slingshot::slingCurveWeights(sce, as.probs = TRUE))
weights$cell <- rownames(weights)
write_table(weights, "08_lineage_weights.csv", cfg)

lineages <- slingshot::slingLineages(sce)
lineage_table <- do.call(rbind, lapply(seq_along(lineages), function(i) {
  data.frame(lineage = paste0("Lineage", i), order = seq_along(lineages[[i]]),
             cluster = lineages[[i]], stringsAsFactors = FALSE)
}))
write_table(lineage_table, "08_lineage_cluster_order.csv", cfg)

curves <- slingshot::slingCurves(sce)
curve_table <- do.call(rbind, lapply(seq_along(curves), function(i) {
  ord <- curves[[i]]$ord
  data.frame(lineage = paste0("Lineage", i), point_order = seq_along(ord),
             dim_1 = curves[[i]]$s[ord, 1], dim_2 = curves[[i]]$s[ord, 2])
}))
write_table(curve_table, "08_curve_coordinates.csv", cfg)

embedding <- as.data.frame(Seurat::Embeddings(sub, "umap"))
names(embedding)[1:2] <- c("dim_1", "dim_2")
embedding$cell <- rownames(embedding)
embedding$cluster <- sub[[cfg$trajectory$cluster_column, drop = TRUE]]
embedding$pseudotime <- pt[match(embedding$cell, pt$cell), 1]
p <- ggplot2::ggplot(embedding, ggplot2::aes(dim_1, dim_2, colour = pseudotime)) +
  ggplot2::geom_point(size = 0.6, alpha = 0.8) +
  ggplot2::geom_path(data = curve_table, ggplot2::aes(dim_1, dim_2, group = lineage),
                     inherit.aes = FALSE, colour = "black", linewidth = 0.6) +
  ggplot2::scale_colour_viridis_c(option = "magma", na.value = "grey85") +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2", colour = "Pseudotime",
                title = "Slingshot lineage and inferred pseudotime") +
  theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
save_figure(p, "08_slingshot_pseudotime", cfg)

dynamic_genes <- intersect(unlist(cfg$trajectory$dynamic_genes), rownames(sub))
if (length(dynamic_genes)) {
  expression <- as.matrix(Seurat::GetAssayData(sub, assay = "RNA", layer = "data")[dynamic_genes, , drop = FALSE])
  dynamics <- data.frame(
    cell = rep(colnames(sub), each = length(dynamic_genes)),
    gene = rep(dynamic_genes, times = ncol(sub)),
    expression = as.vector(expression),
    pseudotime = rep(embedding$pseudotime[match(colnames(sub), embedding$cell)], each = length(dynamic_genes))
  )
  dynamics <- dynamics[is.finite(dynamics$pseudotime), ]
  write_table(dynamics, "08_selected_gene_dynamics.csv", cfg)
  pd <- ggplot2::ggplot(dynamics, ggplot2::aes(pseudotime, expression)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.2, colour = "#777777") +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                         linewidth = 0.6, colour = "#2166AC") +
    ggplot2::facet_wrap(~gene, scales = "free_y") +
    ggplot2::labs(x = "Inferred pseudotime", y = "Normalized expression",
                  title = "Descriptive expression trends") +
    theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
  save_figure(pd, "08_selected_gene_dynamics", cfg, width_mm = 160, height_mm = 95)
}
saveRDS(sce, project_path(cfg$project$output_dir, "objects", "08_slingshot_sce.rds"))
save_session_info("08_Trajectory_Analysis", cfg)
