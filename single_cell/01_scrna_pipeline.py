"""
Public 10x PBMC 3k single-cell RNA-seq workflow.

QC -> normalisation -> HVGs -> PCA -> neighbours ->
UMAP -> Leiden clustering -> marker-gene analysis.
"""

from pathlib import Path

import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]

BASE = ROOT / "single_cell"
PROCESSED = BASE / "processed"
RESULTS = BASE / "results"
FIGURES = BASE / "figures"

for directory in (PROCESSED, RESULTS, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 2
sc.settings.figdir = FIGURES
sc.settings.set_figure_params(dpi=100, facecolor="white")


print("\n==========================================")
print("SINGLE-CELL RNA-SEQ ANALYSIS")
print("==========================================")

# ---------------------------------------------------------
# 1. Real public single-cell dataset
# ---------------------------------------------------------

print("\nLoading public 10x PBMC 3k dataset...")

adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()

initial_cells = adata.n_obs
initial_genes = adata.n_vars

print(f"Initial cells: {initial_cells:,}")
print(f"Initial genes: {initial_genes:,}")


# ---------------------------------------------------------
# 2. QC
# ---------------------------------------------------------

adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

qc_before = pd.DataFrame(
    {
        "metric": [
            "cells",
            "genes",
            "median_genes_per_cell",
            "median_counts_per_cell",
            "median_pct_mt",
        ],
        "value": [
            adata.n_obs,
            adata.n_vars,
            adata.obs["n_genes_by_counts"].median(),
            adata.obs["total_counts"].median(),
            adata.obs["pct_counts_mt"].median(),
        ],
    }
)

qc_before.to_csv(
    RESULTS / "qc_before_filtering.tsv",
    sep="\t",
    index=False,
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    show=False,
    save="_pbmc3k_qc_before.png",
)


# ---------------------------------------------------------
# 3. Filtering
# ---------------------------------------------------------

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)

adata = adata[
    (adata.obs["n_genes_by_counts"] < 2500)
    & (adata.obs["pct_counts_mt"] < 5)
].copy()

print("\nAfter QC:")
print(f"Cells retained: {adata.n_obs:,}")
print(f"Genes retained: {adata.n_vars:,}")

qc_after = pd.DataFrame(
    {
        "metric": [
            "cells_retained",
            "genes_retained",
            "median_genes_per_cell",
            "median_counts_per_cell",
            "median_pct_mt",
        ],
        "value": [
            adata.n_obs,
            adata.n_vars,
            adata.obs["n_genes_by_counts"].median(),
            adata.obs["total_counts"].median(),
            adata.obs["pct_counts_mt"].median(),
        ],
    }
)

qc_after.to_csv(
    RESULTS / "qc_after_filtering.tsv",
    sep="\t",
    index=False,
)


# Preserve counts
adata.layers["counts"] = adata.X.copy()


# ---------------------------------------------------------
# 4. Normalisation
# ---------------------------------------------------------

print("\nNormalising...")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Preserve full normalised expression matrix
adata.raw = adata


# ---------------------------------------------------------
# 5. Highly variable genes
# ---------------------------------------------------------

sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
)

n_hvg = int(adata.var["highly_variable"].sum())

print(f"Highly variable genes: {n_hvg:,}")

adata = adata[:, adata.var["highly_variable"]].copy()


# ---------------------------------------------------------
# 6. Scale + PCA
# ---------------------------------------------------------

sc.pp.regress_out(
    adata,
    ["total_counts", "pct_counts_mt"],
)

sc.pp.scale(
    adata,
    max_value=10,
)

print("Running PCA...")

sc.tl.pca(
    adata,
    svd_solver="arpack",
)

sc.pl.pca_variance_ratio(
    adata,
    log=True,
    show=False,
    save="_pbmc3k_pca_variance.png",
)


# ---------------------------------------------------------
# 7. Neighbours + UMAP
# ---------------------------------------------------------

print("Building neighbour graph...")

sc.pp.neighbors(
    adata,
    n_neighbors=10,
    n_pcs=40,
)

print("Running UMAP...")

sc.tl.umap(
    adata,
    random_state=42,
)


# ---------------------------------------------------------
# 8. Leiden clustering
# ---------------------------------------------------------

print("Running Leiden clustering...")

sc.tl.leiden(
    adata,
    resolution=0.5,
    key_added="leiden",
    random_state=42,
)

cluster_counts = (
    adata.obs["leiden"]
    .value_counts()
    .sort_index()
    .rename_axis("cluster")
    .reset_index(name="cells")
)

cluster_counts.to_csv(
    RESULTS / "leiden_cluster_counts.tsv",
    sep="\t",
    index=False,
)

print("\nClusters:")
print(cluster_counts.to_string(index=False))


# ---------------------------------------------------------
# 9. Marker-gene analysis
# ---------------------------------------------------------

print("\nCalculating marker genes...")

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
)

markers = sc.get.rank_genes_groups_df(
    adata,
    group=None,
)

markers.to_csv(
    RESULTS / "leiden_marker_genes.tsv",
    sep="\t",
    index=False,
)


# ---------------------------------------------------------
# 10. Figures
# ---------------------------------------------------------

sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    title="PBMC 3k - Leiden clusters",
    show=False,
    save="_pbmc3k_leiden.png",
)

canonical_markers = [
    "IL7R",
    "CCR7",
    "CD14",
    "LYZ",
    "MS4A1",
    "CD79A",
    "GNLY",
    "NKG7",
    "FCGR3A",
    "LST1",
    "PPBP",
]

available_markers = [
    gene
    for gene in canonical_markers
    if gene in adata.raw.var_names
]

if available_markers:
    sc.pl.umap(
        adata,
        color=available_markers,
        ncols=3,
        show=False,
        save="_pbmc3k_markers.png",
    )


# ---------------------------------------------------------
# 11. Save AnnData
# ---------------------------------------------------------

adata.write(
    PROCESSED / "pbmc3k_processed.h5ad"
)


# ---------------------------------------------------------
# 12. Summary
# ---------------------------------------------------------

summary = {
    "dataset": "10x PBMC 3k",
    "initial_cells": initial_cells,
    "initial_genes": initial_genes,
    "cells_after_qc": int(adata.n_obs),
    "highly_variable_genes": int(adata.n_vars),
    "leiden_clusters": int(adata.obs["leiden"].nunique()),
    "marker_tests": int(len(markers)),
}

pd.Series(summary).to_csv(
    RESULTS / "scrna_analysis_summary.tsv",
    sep="\t",
    header=False,
)

print("\n==========================================")
print("SINGLE-CELL RNA-SEQ SUMMARY")
print("==========================================")

for key, value in summary.items():
    print(f"{key}: {value}")

print("\nSaved:")
print(PROCESSED / "pbmc3k_processed.h5ad")
print(RESULTS / "leiden_marker_genes.tsv")
print(FIGURES)

print("\nSINGLE-CELL RNA-SEQ ANALYSIS PASSED")

