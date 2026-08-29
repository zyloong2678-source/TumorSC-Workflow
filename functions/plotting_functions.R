# Unified, restrained plotting helpers for scientific output.

tumorsc_palette <- c(
  Epithelial = "#D55E00", T_NK = "#0072B2", B = "#56B4E9",
  Myeloid = "#009E73", Fibroblast = "#CC79A7", Endothelial = "#E69F00",
  Other = "#7A7A7A"
)

theme_tumorsc <- function(base_size = 8, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.title = ggplot2::element_text(face = "bold")
    )
}

save_figure <- function(plot, filename, cfg, width_mm = NULL, height_mm = NULL) {
  assert_packages(c("ggplot2", "svglite", "ragg"))
  ensure_output_dirs(cfg)
  width_mm <- width_mm %||% cfg$plotting$width_mm
  height_mm <- height_mm %||% cfg$plotting$height_mm
  if (identical(cfg$project$input_mode, "synthetic")) {
    disclaimer <- "Synthetic software test - no biological interpretation."
    if (inherits(plot, "patchwork")) {
      plot <- plot + patchwork::plot_annotation(caption = disclaimer)
    } else {
      plot <- plot + ggplot2::labs(caption = disclaimer)
    }
  }
  stem <- project_path(cfg$project$output_dir, "figures", filename)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width_mm, height = height_mm,
                  units = "mm", device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(stem, ".svg"), plot, width = width_mm, height = height_mm,
                  units = "mm", device = svglite::svglite)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width_mm, height = height_mm,
                  units = "mm", dpi = cfg$plotting$dpi, device = ragg::agg_png)
  invisible(stem)
}

plot_qc_violin <- function(object, cfg) {
  Seurat::VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3, pt.size = 0) &
    theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
}

plot_volcano <- function(markers, cfg) {
  assert_packages(c("ggplot2", "ggrepel"))
  fc <- if ("avg_log2FC" %in% names(markers)) "avg_log2FC" else "avg_logFC"
  required <- c(fc, "p_val_adj", "gene")
  if (!all(required %in% names(markers))) stop("Volcano input lacks required columns.", call. = FALSE)
  markers$status <- "Not significant"
  sig <- markers$p_val_adj < cfg$differential_expression$adjusted_p_cutoff &
    abs(markers[[fc]]) >= cfg$differential_expression$logfc_threshold
  markers$status[sig & markers[[fc]] > 0] <- "Higher"
  markers$status[sig & markers[[fc]] < 0] <- "Lower"
  markers$neg_log10_fdr <- -log10(pmax(markers$p_val_adj, .Machine$double.xmin))
  labels <- head(markers[order(markers$p_val_adj), ], 10)
  ggplot2::ggplot(markers, ggplot2::aes(x = .data[[fc]], y = neg_log10_fdr, colour = status)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.7) +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = gene), size = 2.2,
                             max.overlaps = Inf, min.segment.length = 0) +
    ggplot2::scale_colour_manual(values = c(Higher = "#B2182B", Lower = "#2166AC",
                                             `Not significant` = "#BDBDBD")) +
    ggplot2::labs(x = "Average log2 fold change", y = "-log10 adjusted P", colour = NULL) +
    theme_tumorsc(cfg$plotting$base_size, cfg$plotting$base_family)
}
