# Purpose: Load counts, calculate QC metrics, filter cells, and normalize expression.
# Input:   Synthetic object, 10x directories, or Seurat RDS configured in config.yaml
# Output:  results/objects/01_preprocessed.rds and QC figures/tables

source("functions/workflow_utils.R")
source("functions/plotting_functions.R")
cfg <- read_config()
assert_packages(c("Seurat", "ggplot2", "patchwork"))
ensure_output_dirs(cfg)

load_input <- function(cfg) {
  mode <- cfg$project$input_mode
  if (mode == "synthetic") return(load_seurat("00_synthetic_raw", cfg))
  if (mode == "rds") {
    path <- project_path(cfg$project$input_path)
    require_input(path, "Configured RDS")
    return(readRDS(path))
  }
  if (mode == "10x_directory") {
    sheet_path <- project_path(cfg$project$sample_sheet)
    require_input(sheet_path, "Sample sheet")
    sheet <- utils::read.csv(sheet_path, stringsAsFactors = FALSE)
    if (!all(c("sample_id", "path") %in% names(sheet))) stop("Sample sheet requires sample_id and path.")
    objects <- lapply(seq_len(nrow(sheet)), function(i) {
      path <- project_path(sheet$path[i])
      require_input(path, paste("10x directory for", sheet$sample_id[i]))
      x <- Seurat::CreateSeuratObject(Seurat::Read10X(path), project = sheet$sample_id[i],
                                      min.cells = cfg$qc$min_cells_per_gene)
      x$sample_id <- sheet$sample_id[i]
      for (column in setdiff(names(sheet), "path")) x[[column]] <- sheet[[column]][i]
      x
    })
    return(if (length(objects) == 1L) objects[[1]] else merge(objects[[1]], y = objects[-1], add.cell.ids = sheet$sample_id))
  }
  stop("Unsupported project.input_mode: ", mode, call. = FALSE)
}

object <- load_input(cfg)
mt_pattern <- if (tolower(cfg$project$species) == "human") "^MT-" else "^mt-"
object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = mt_pattern)
qc_before <- plot_qc_violin(object, cfg)
save_figure(qc_before, "01_qc_before_filtering", cfg, width_mm = 180, height_mm = 70)

keep <- object$nFeature_RNA >= cfg$qc$min_features & object$nFeature_RNA <= cfg$qc$max_features &
  object$nCount_RNA >= cfg$qc$min_counts & object$nCount_RNA <= cfg$qc$max_counts &
  object$percent.mt <= cfg$qc$max_percent_mt
summary_table <- data.frame(stage = c("before", "after"), cells = c(ncol(object), sum(keep)))
object <- subset(object, cells = colnames(object)[keep])
if (ncol(object) < 20) stop("Fewer than 20 cells remain after QC; review config cutoffs.", call. = FALSE)

if (cfg$preprocessing$method == "SCTransform") {
  object <- Seurat::SCTransform(object, vars.to.regress = unlist(cfg$preprocessing$regress_variables),
                                variable.features.n = cfg$preprocessing$variable_features, verbose = FALSE)
} else if (cfg$preprocessing$method == "LogNormalize") {
  object <- Seurat::NormalizeData(object, scale.factor = cfg$preprocessing$scale_factor, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, nfeatures = cfg$preprocessing$variable_features, verbose = FALSE)
} else stop("preprocessing.method must be LogNormalize or SCTransform.", call. = FALSE)

write_table(summary_table, "01_qc_cell_counts.csv", cfg)
save_figure(plot_qc_violin(object, cfg), "01_qc_after_filtering", cfg, width_mm = 180, height_mm = 70)
save_seurat(object, "01_preprocessed", cfg)
save_session_info("01_QC_Preprocessing", cfg)

