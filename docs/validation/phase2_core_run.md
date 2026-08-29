# Phase 2 local validation record

Date: 2026-08-29  
Platform: Windows 11 x64 (build 26200)  
R: 4.6.1 (2026-06-24 ucrt)  
Seurat: 5.5.1  
SeuratObject: 5.4.0

## Executed path

Scripts `00` through `05` were executed sequentially with an absolute
`Rscript.exe` path and a repository-local package library.

- `00`: generated 240 explicitly synthetic cells across two synthetic samples.
- `01`: retained 240/240 cells under the configured smoke-test QC thresholds.
- `02`: generated PCA, graph clusters, and UMAP; six clusters of 40 cells each
  matched the six deliberately injected latent software-test profiles.
- `03`: generated canonical-marker DotPlot/FeaturePlot output. Labels remained
  `Unreviewed_cluster_*` because no human-reviewed mapping was supplied.
- `04`: generated the candidate evidence table. No epithelial or malignant-cell
  candidate call was made because annotation remained unreviewed.
- `05`: tested six cluster contrasts and exported 238 marker rows plus a heatmap.

## Warnings retained

- The local R startup reported failures to set `C.UTF-8` locale categories and
  used locale `C`; plots were visually checked after switching the synthetic
  disclaimer to ASCII text.
- Seurat reported that four configured canonical markers were absent from the
  small synthetic fixture. Missing markers were reported, not imputed.
- Seurat suggested optional `presto` installation for faster Wilcoxon tests;
  the base implementation completed successfully.

## Not executed

Scripts `06`-`09`, Harmony integration, CNV result import, CellChat, Slingshot,
decoupleR/DoRothEA, inferCNV, and CopyKAT were not locally executed. Presence of
their code is not evidence of runtime validation.

