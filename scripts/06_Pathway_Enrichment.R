# Purpose: Run Hallmark GSEA and calculate configurable pathway/module scores.
# Input:   differential-expression table and results/objects/04_malignancy_candidates.rds
# Output:  GSEA table/plot and scored Seurat object

source("functions/workflow_utils.R")
source("functions/plotting_functions.R")
cfg <- read_config()
assert_packages(c("Seurat", "fgsea", "msigdbr", "ggplot2"), optional = TRUE)
object <- load_seurat("04_malignancy_candidates", cfg)
de_path <- project_path(cfg$project$output_dir, "tables", "05_pairwise_differential_expression.csv")
require_input(de_path, "Pairwise DE table (configure ident_1 and run script 05)")
de <- utils::read.csv(de_path, stringsAsFactors = FALSE)
fc <- if ("avg_log2FC" %in% names(de)) "avg_log2FC" else "avg_logFC"
ranks <- stats::setNames(de[[fc]], de$gene)
ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)

msig <- msigdbr::msigdbr(species = cfg$pathway$organism, category = cfg$pathway$msigdb_category)
pathways <- split(msig$gene_symbol, msig$gs_name)
gsea <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks,
  minSize = cfg$pathway$min_size, maxSize = cfg$pathway$max_size)
gsea <- as.data.frame(gsea)
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
write_table(gsea, "06_hallmark_gsea.csv", cfg)

top <- head(gsea[order(gsea$padj, -abs(gsea$NES)), ], 15)
p <- ggplot2::ggplot(top, ggplot2::aes(x = stats::reorder(pathway, NES), y = NES, fill = NES > 0)) +
  ggplot2::geom_col(width = 0.7) + ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"), guide = "none") +
  ggplot2::labs(x = NULL, y = "Normalized enrichment score") +
  theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
save_figure(p, "06_hallmark_gsea", cfg, width_mm = 135, height_mm = 100)

score_sets <- pathways[head(names(pathways), 5)]
if (cfg$pathway$score_method == "UCell") {
  assert_packages("UCell", optional = TRUE)
  object <- UCell::AddModuleScore_UCell(object, features = score_sets)
} else object <- Seurat::AddModuleScore(object, features = score_sets, name = "HallmarkScore")
save_seurat(object, "06_pathway_scored", cfg)
save_session_info("06_Pathway_Enrichment", cfg)

