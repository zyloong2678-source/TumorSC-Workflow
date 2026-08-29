# Data directory

No private or unpublished data are distributed with this repository.

## Supported inputs

1. **Synthetic smoke test** (default): generated locally by
   `scripts/00_Generate_Synthetic_Demo.R`. It is only for checking software
   plumbing and must not be interpreted biologically.
2. **10x directories**: place one filtered feature-barcode matrix directory per
   sample under `data/raw/<sample_id>/` and populate `data/sample_metadata.csv`.
3. **Seurat RDS**: set `project.input_mode: rds` and `project.input_path` in
   `config/config.yaml`.

The sample sheet must contain `sample_id` and `path` columns. Optional columns
such as `patient_id`, `condition`, and `batch` are retained as metadata.
Do not commit controlled-access or identifiable data.

