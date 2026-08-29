# Shared validation, I/O, and reproducibility helpers.

assert_packages <- function(packages, optional = FALSE) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    message <- paste("Missing R package(s):", paste(missing, collapse = ", "))
    if (optional) stop(message, "\nInstall them before enabling this optional module.", call. = FALSE)
    stop(message, "\nRun: Rscript environment/install_packages.R", call. = FALSE)
  }
  invisible(TRUE)
}

read_config <- function(path = "config/config.yaml") {
  assert_packages("yaml")
  if (!file.exists(path)) stop("Config file not found: ", path, call. = FALSE)
  cfg <- yaml::read_yaml(path)
  required <- c("project", "qc", "preprocessing", "integration", "plotting")
  absent <- setdiff(required, names(cfg))
  if (length(absent)) stop("Missing config section(s): ", paste(absent, collapse = ", "), call. = FALSE)
  set.seed(cfg$project$seed %||% 2026)
  cfg
}

`%||%` <- function(x, y) if (is.null(x)) y else x

project_path <- function(...) {
  assert_packages("here")
  here::here(...)
}

ensure_output_dirs <- function(cfg) {
  root <- project_path(cfg$project$output_dir %||% "results")
  dirs <- file.path(root, c("objects", "tables", "figures", "logs", "example_figures"))
  vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

save_session_info <- function(script_name, cfg) {
  ensure_output_dirs(cfg)
  path <- project_path(cfg$project$output_dir, "logs", paste0(script_name, "_sessionInfo.txt"))
  capture.output(sessionInfo(), file = path)
  invisible(path)
}

require_input <- function(path, label = "Required input") {
  if (!file.exists(path) && !dir.exists(path)) {
    stop(label, " not found: ", path, "\nRun the preceding script or update config/config.yaml.", call. = FALSE)
  }
  invisible(path)
}

save_seurat <- function(object, name, cfg) {
  ensure_output_dirs(cfg)
  path <- project_path(cfg$project$output_dir, "objects", paste0(name, ".rds"))
  saveRDS(object, path)
  message("Saved: ", path)
  invisible(path)
}

load_seurat <- function(name, cfg) {
  path <- project_path(cfg$project$output_dir, "objects", paste0(name, ".rds"))
  require_input(path, "Seurat object")
  readRDS(path)
}

write_table <- function(x, filename, cfg) {
  ensure_output_dirs(cfg)
  path <- project_path(cfg$project$output_dir, "tables", filename)
  utils::write.csv(x, path, row.names = FALSE)
  message("Saved: ", path)
  invisible(path)
}

available_features <- function(object, requested) {
  intersect(unique(unlist(requested, use.names = FALSE)), rownames(object))
}

validate_metadata <- function(object, columns) {
  absent <- setdiff(columns, colnames(object[[]]))
  if (length(absent)) stop("Missing metadata column(s): ", paste(absent, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

