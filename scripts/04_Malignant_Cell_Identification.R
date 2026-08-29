# Purpose: Build an evidence table for candidate malignant epithelial cells.
# Input:   results/objects/03_annotated.rds; optional externally generated CNV calls
# Output:  results/objects/04_malignancy_candidates.rds and candidate evidence table
# Candidate status is supporting evidence, not malignant-cell ground truth.

source("functions/workflow_utils.R")
source("functions/plotting_functions.R")
cfg <- read_config()
assert_packages(c("Seurat", "ggplot2"))
object <- load_seurat("03_annotated", cfg)
level1 <- cfg$annotation$level1_column
validate_metadata(object, level1)

epi_cells <- colnames(object)[object[[level1, drop = TRUE]] == cfg$malignancy$epithelial_label]
object$malignancy_candidate <- "Not_evaluated"
if (!length(epi_cells)) {
  warning("No cells have the configured epithelial label; no candidate calls were made.")
} else {
  object$malignancy_candidate[epi_cells] <- "Epithelial_requires_CNV_context"
  marker_sets <- lapply(cfg$malignancy$candidate_marker_sets, function(x) intersect(x, rownames(object)))
  marker_sets <- marker_sets[lengths(marker_sets) > 0]
  if (length(marker_sets)) object <- Seurat::AddModuleScore(object, features = marker_sets, name = "candidate_marker_score")
}

if (cfg$malignancy$cnv_method != "none") {
  if (is.null(cfg$malignancy$cnv_result_path)) stop("Set malignancy.cnv_result_path for external CNV calls.", call. = FALSE)
  cnv_path <- project_path(cfg$malignancy$cnv_result_path)
  require_input(cnv_path, "External CNV result")
  cnv <- utils::read.csv(cnv_path, stringsAsFactors = FALSE)
  if (!all(c("cell", "cnv_call") %in% names(cnv))) stop("CNV result requires cell and cnv_call columns.")
  object$cnv_call <- cnv$cnv_call[match(colnames(object), cnv$cell)]
}

evidence_cols <- intersect(c(level1, "malignancy_candidate", "cnv_call",
                             grep("^candidate_marker_score", colnames(object[[]]), value = TRUE)), colnames(object[[]]))
evidence <- data.frame(cell = colnames(object), object[[]][, evidence_cols, drop = FALSE], check.names = FALSE)
write_table(evidence, "04_malignancy_candidate_evidence.csv", cfg)
p <- Seurat::DimPlot(object, group.by = "malignancy_candidate") +
  ggplot2::ggtitle("Candidate malignant-cell assessment status") +
  theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
save_figure(p, "04_malignancy_candidate_status", cfg)
save_seurat(object, "04_malignancy_candidates", cfg)
save_session_info("04_Malignant_Cell_Identification", cfg)
