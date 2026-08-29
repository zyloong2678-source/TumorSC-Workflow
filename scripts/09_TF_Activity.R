# Purpose: Estimate regulator activity using a curated DoRothEA network and decoupleR.
# Input:   results/objects/04_malignancy_candidates.rds
# Output:  averaged expression matrix, filtered network summary, TF activity
#          table/matrix, within-state rankings, heatmap, and ranking plot
# Caution: Inferred activity is not direct proof of biochemical TF activity.

source("functions/workflow_utils.R")
source("functions/plotting_functions.R")
cfg <- read_config()
if (!isTRUE(cfg$tf_activity$enabled)) stop("Optional module disabled. Set tf_activity.enabled: true after installing dependencies.", call. = FALSE)
assert_packages(c("decoupleR", "dorothea", "Seurat", "pheatmap", "dplyr", "tidyr", "ggplot2"), optional = TRUE)
object <- load_seurat("04_malignancy_candidates", cfg)
validate_metadata(object, cfg$tf_activity$group_by)
expr <- as.matrix(Seurat::AverageExpression(object, group.by = cfg$tf_activity$group_by,
                                             assays = "RNA", layer = "data")$RNA)
utils::write.csv(data.frame(gene = rownames(expr), expr, check.names = FALSE),
                 project_path(cfg$project$output_dir, "tables", "09_average_expression_matrix.csv"),
                 row.names = FALSE)
network <- if (tolower(cfg$tf_activity$organism) == "human") dorothea::dorothea_hs else dorothea::dorothea_mm
network <- dplyr::filter(network, confidence %in% unlist(cfg$tf_activity$confidence_levels))
network <- dplyr::distinct(network, tf, target, .keep_all = TRUE)
network_summary <- dplyr::summarise(dplyr::group_by(network, tf, confidence),
                                    targets = dplyr::n_distinct(target), .groups = "drop")
write_table(as.data.frame(network_summary), "09_regulon_network_summary.csv", cfg)
activities <- decoupleR::run_ulm(mat = expr, network = network, .source = "tf",
                                 .target = "target", .mor = "mor", minsize = cfg$tf_activity$min_n)
activities <- as.data.frame(activities)
score_column <- intersect(c("score", "estimate"), names(activities))[1]
if (is.na(score_column)) stop("decoupleR output has no score/estimate column; check installed API.", call. = FALSE)
names(activities)[names(activities) == score_column] <- "activity_score"
write_table(activities, "09_tf_activities.csv", cfg)
ranking <- dplyr::group_by(activities, condition)
ranking <- dplyr::slice_max(ranking, order_by = abs(activity_score), n = cfg$tf_activity$top_tfs_per_group,
                            with_ties = FALSE)
ranking <- dplyr::ungroup(ranking)
write_table(as.data.frame(ranking), "09_tf_activity_ranking.csv", cfg)
mat <- tidyr::pivot_wider(activities, names_from = condition, values_from = activity_score, values_fill = 0)
tf_names <- mat$source
mat <- as.matrix(mat[, setdiff(names(mat), "source"), drop = FALSE]); rownames(mat) <- tf_names
utils::write.csv(data.frame(tf = rownames(mat), mat, check.names = FALSE),
                 project_path(cfg$project$output_dir, "tables", "09_tf_activity_matrix.csv"), row.names = FALSE)
top_sources <- unique(ranking$source)
heat_mat <- mat[intersect(top_sources, rownames(mat)), , drop = FALSE]
grDevices::pdf(project_path(cfg$project$output_dir, "figures", "09_tf_activity_heatmap.pdf"), width = 6, height = 7)
pheatmap::pheatmap(heat_mat, scale = "row", border_color = NA, main = "Inferred TF activity")
grDevices::dev.off()
p <- ggplot2::ggplot(ranking, ggplot2::aes(x = stats::reorder(source, activity_score),
                                           y = activity_score, fill = activity_score > 0)) +
  ggplot2::geom_col(width = 0.7) + ggplot2::coord_flip() +
  ggplot2::facet_wrap(~condition, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"), guide = "none") +
  ggplot2::labs(x = "Transcription factor", y = "ULM activity score",
                title = "Top inferred TF activities by state") +
  theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
save_figure(p, "09_tf_activity_ranking", cfg, width_mm = 170, height_mm = 120)
save_session_info("09_TF_Activity", cfg)
