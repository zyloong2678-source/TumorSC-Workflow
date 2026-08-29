# Environment

Run `Rscript environment/install_packages.R --core` for the core synthetic
workflow, or omit `--core` to add the declared Bioconductor dependencies. Then
initialize a project-local lockfile:

On Windows, a repository-local library can be used without modifying the R
installation, for example by setting `R_LIBS_USER` to `environment/library`
before running the installer and workflow.

```r
install.packages("renv")
renv::init(bare = TRUE)
renv::snapshot()
```

Commit `renv.lock` after testing on your own machine. This repository does not
ship a fabricated lockfile because package versions were not resolved in the
creation environment. Optional CellChat, Harmony, inferCNV, and CopyKAT modules
should be installed from their current official sources and recorded by a new
`renv::snapshot()`.
