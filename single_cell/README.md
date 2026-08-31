# Single-cell RNA-seq analysis

This module demonstrates an end-to-end single-cell RNA-seq workflow using
the public 10x Genomics PBMC 3k dataset.

## Workflow

The analysis was implemented in Python using Scanpy and AnnData.

Steps:

1. Load public 10x PBMC 3k count data
2. Calculate cell- and gene-level QC metrics
3. Assess mitochondrial transcript fraction
4. Filter low-quality cells and low-detection genes
5. Library-size normalisation
6. Log transformation
7. Highly variable gene selection
8. Regression and scaling
9. Principal component analysis
10. Nearest-neighbour graph construction
11. UMAP dimensionality reduction
12. Leiden community detection
13. Differential marker-gene analysis
14. Canonical-marker-based cell-type annotation
15. Export reproducible result tables and figures

## Dataset summary

| Metric | Result |
|---|---:|
| Initial cells | 2,700 |
| Initial genes | 32,738 |
| Cells retained after QC | 2,638 |
| Highly variable genes | 1,838 |
| Leiden clusters | 6 |
| Marker tests | 82,284 |

## Cell-type composition

| Cell type | Cells | Percentage |
|---|---:|---:|
| T cells | 1,185 | 44.92% |
| CD14 Monocytes | 636 | 24.11% |
| NK / Cytotoxic cells | 427 | 16.19% |
| B cells | 341 | 12.93% |
| Dendritic cells | 36 | 1.36% |
| Platelets | 13 | 0.49% |

Cell identities are provisional marker-based annotations based on canonical
PBMC marker programmes and differential marker-gene evidence.

Examples include:

- T cells: CD3D and related T-cell markers
- Myeloid/monocyte populations: LYZ, TYROBP and related markers
- NK/cytotoxic cells: NKG7/GNLY-associated programme
- B cells: MS4A1/CD79A-associated programme
- Platelets: PF4 and PPBP
- Dendritic cells: antigen-presentation/dendritic marker programme

## Main outputs

### Processed data

Processed AnnData objects are generated locally under:

single_cell/processed/

These .h5ad files are intentionally excluded from Git because they are
generated analysis artefacts.

### Result tables

single_cell/results/

Includes:

- QC summaries
- Leiden cluster counts
- differential marker genes
- top markers per cluster
- marker-panel scores
- cluster annotations
- cell-type composition

### Figures

single_cell/figures/

Includes:

- QC violin plots
- PCA variance plots
- Leiden UMAP
- canonical-marker UMAPs
- biologically annotated cell-type UMAP

## Reproducibility

Run:

```bash
python single_cell/01_scrna_pipeline.py
python single_cell/02_annotate_celltypes.py
pytest -q
