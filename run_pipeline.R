# Run the core workflow sequentially. Optional modules remain opt-in in config.
core <- sprintf("scripts/%02d_%s.R", 0:5, c(
  "Generate_Synthetic_Demo", "QC_Preprocessing", "Integration_Clustering",
  "Cell_Annotation", "Malignant_Cell_Identification", "Differential_Expression"
))

for (script in core) {
  message("\n>>> Running ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  if (!identical(status, 0L)) stop("Pipeline stopped at ", script, call. = FALSE)
}
