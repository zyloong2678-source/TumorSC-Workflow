# Workflow diagram

The diagram is stored as Mermaid source so changes remain reviewable in Git.

```mermaid
flowchart TD
  A[Raw counts or synthetic smoke test] --> B[QC and preprocessing]
  B --> C[Batch-aware integration]
  C --> D[Clustering and UMAP]
  D --> E[Two-level, marker-informed annotation]
  E --> F[Epithelial-cell assessment]
  F --> G[Candidate malignant epithelial cells]
  G --> H{Downstream analyses}
  H --> I[Differential expression]
  H --> J[Pathway enrichment and scores]
  H --> K[Cell-cell communication]
  H --> L[Trajectory and pseudotime]
  H --> M[TF activity inference]
```

CNV inference is an optional evidence layer between epithelial-cell assessment
and candidate classification. It requires appropriate reference cells and
biological interpretation.

