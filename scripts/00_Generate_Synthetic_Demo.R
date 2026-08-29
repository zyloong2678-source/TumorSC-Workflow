# Purpose: Create a small, explicitly synthetic Seurat object for smoke testing.
# Input:   config/config.yaml
# Output:  results/objects/00_synthetic_raw.rds
# The object tests software plumbing only and must not support biological claims.

source("functions/workflow_utils.R")
cfg <- read_config()
assert_packages(c("Seurat", "Matrix"))
ensure_output_dirs(cfg)

if (cfg$project$input_mode != "synthetic") {
  stop("Synthetic generation is allowed only when project.input_mode is 'synthetic'.", call. = FALSE)
}

n_genes <- cfg$synthetic$n_genes
cells_per_sample <- cfg$synthetic$cells_per_sample
sample_ids <- unlist(cfg$synthetic$sample_ids)
canonical <- unique(c(
  "EPCAM", "KRT8", "KRT18", "KRT19", "CD3D", "CD3E", "TRAC", "NKG7", "GNLY",
  "CD79A", "MS4A1", "CD74", "LYZ", "LST1", "TYROBP", "COL1A1", "COL1A2", "DCN",
  "PECAM1", "VWF", "KDR", "MKI67", "TOP2A", "C1QA", "C1QB", "FCER1A", "CD1C",
  paste0("MT-", c("CO1", "CO2", "CO3", "ND1", "ND2", "CYB"))
))
genes <- unique(c(canonical, paste0("GENE", seq_len(max(0, n_genes - length(canonical))))))
cell_names <- unlist(lapply(sample_ids, function(x) paste0(x, "_cell", seq_len(cells_per_sample))))
counts <- matrix(stats::rpois(length(genes) * length(cell_names), lambda = 1.2),
                 nrow = length(genes), dimnames = list(genes, cell_names))
counts[sample(seq_along(counts), floor(length(counts) * 0.65))] <- 0

# Add deterministic latent profiles so clustering and multi-group DE code paths
# are exercised. These labels are software fixtures, not inferred cell types.
profile_markers <- list(
  profile_1 = c("EPCAM", "KRT8", "KRT18", "KRT19", paste0("GENE", 1:30)),
  profile_2 = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", paste0("GENE", 31:60)),
  profile_3 = c("CD79A", "MS4A1", "CD74", paste0("GENE", 61:90)),
  profile_4 = c("LYZ", "LST1", "TYROBP", "C1QA", "C1QB", paste0("GENE", 91:120)),
  profile_5 = c("COL1A1", "COL1A2", "DCN", paste0("GENE", 121:150)),
  profile_6 = c("PECAM1", "VWF", "KDR", paste0("GENE", 151:180))
)
latent_profile <- rep(names(profile_markers), length.out = length(cell_names))
for (profile in names(profile_markers)) {
  rows <- intersect(profile_markers[[profile]], rownames(counts))
  cols <- which(latent_profile == profile)
  counts[rows, cols] <- counts[rows, cols] + matrix(
    stats::rpois(length(rows) * length(cols), lambda = 8), nrow = length(rows)
  )
}
counts <- Matrix::Matrix(counts, sparse = TRUE)
object <- Seurat::CreateSeuratObject(counts = counts, project = "SYNTHETIC_DEMO", min.cells = 0)
object$sample_id <- rep(sample_ids, each = cells_per_sample)
object$data_provenance <- "synthetic_smoke_test_only"
object$synthetic_latent_profile <- latent_profile
save_seurat(object, "00_synthetic_raw", cfg)
save_session_info("00_Generate_Synthetic_Demo", cfg)
