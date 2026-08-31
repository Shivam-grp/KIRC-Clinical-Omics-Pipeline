"""
Cell-type annotation of the PBMC 3k single-cell RNA-seq analysis.

Canonical immune-cell marker panels are scored in each cell.
Mean marker scores are calculated per Leiden cluster and used
to generate provisional biological cell-type annotations.

Annotations should be interpreted together with differential
marker-gene results.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]

INPUT = ROOT / "single_cell" / "processed" / "pbmc3k_processed.h5ad"
RESULTS = ROOT / "single_cell" / "results"
FIGURES = ROOT / "single_cell" / "figures"

RESULTS.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)


print("\n==========================================")
print("PBMC CELL-TYPE ANNOTATION")
print("==========================================")

adata = sc.read_h5ad(INPUT)

print(f"Cells: {adata.n_obs:,}")
print(f"Leiden clusters: {adata.obs['leiden'].nunique()}")


# ----------------------------------------------------------
# Canonical PBMC marker panels
# ----------------------------------------------------------

marker_panels = {

    "T cells": [
        "CD3D",
        "CD3E",
        "TRAC",
        "LCK",
        "IL7R",
        "CCR7",
    ],

    "NK / Cytotoxic cells": [
        "NKG7",
        "GNLY",
        "PRF1",
        "CTSW",
        "GZMB",
        "CCL5",
    ],

    "B cells": [
        "MS4A1",
        "CD79A",
        "CD79B",
        "CD37",
        "CD74",
        "HLA-DRA",
    ],

    "CD14 Monocytes": [
        "LST1",
        "LYZ",
        "TYROBP",
        "CTSS",
        "S100A8",
        "S100A9",
        "FCN1",
    ],

    "FCGR3A Monocytes": [
        "FCGR3A",
        "LST1",
        "TYROBP",
        "MS4A7",
        "IFITM3",
        "LGALS3",
        "LILRB1",
    ],

    "Dendritic cells": [
        "FCER1A",
        "CD1C",
        "CST3",
        "HLA-DPA1",
        "HLA-DRA",
    ],

    "Platelets": [
        "PPBP",
        "PF4",
        "GNG11",
        "RGS18",
        "NRGN",
        "GP9",
    ],
}


# ----------------------------------------------------------
# Check which markers exist
# ----------------------------------------------------------

gene_space = set(adata.raw.var_names)

print("\n===== MARKER PANELS =====")

available_panels = {}

for cell_type, genes in marker_panels.items():

    available = [
        gene
        for gene in genes
        if gene in gene_space
    ]

    available_panels[cell_type] = available

    print(
        f"{cell_type:22s}: "
        f"{len(available)}/{len(genes)} markers available"
    )


# ----------------------------------------------------------
# Score marker panels
# ----------------------------------------------------------

print("\nScoring cells...")

score_columns = []

for cell_type, genes in available_panels.items():

    if not genes:
        continue

    score_name = (
        "score_"
        + cell_type.lower()
        .replace(" ", "_")
        .replace("/", "")
    )

    sc.tl.score_genes(
        adata,
        gene_list=genes,
        score_name=score_name,
        use_raw=True,
        random_state=42,
    )

    score_columns.append(score_name)


# ----------------------------------------------------------
# Cluster-level marker scores
# ----------------------------------------------------------

cluster_scores = (
    adata.obs
    .groupby("leiden", observed=True)[score_columns]
    .mean()
)

cluster_scores.to_csv(
    RESULTS / "cluster_celltype_scores.tsv",
    sep="\t",
)


print("\n==========================================")
print("CELL-TYPE SCORES BY CLUSTER")
print("==========================================")

print(cluster_scores.round(3).to_string())


# ----------------------------------------------------------
# Select strongest marker programme
# ----------------------------------------------------------

score_to_label = {
    "score_" + name.lower().replace(" ", "_").replace("/", ""): name
    for name in available_panels
}

cluster_annotation = {}

for cluster in cluster_scores.index:

    winning_score = cluster_scores.loc[cluster].idxmax()

    cluster_annotation[str(cluster)] = score_to_label[winning_score]


annotation_table = pd.DataFrame(
    {
        "leiden_cluster": list(cluster_annotation.keys()),
        "cell_type": list(cluster_annotation.values()),
    }
)

annotation_table.to_csv(
    RESULTS / "cluster_celltype_annotations.tsv",
    sep="\t",
    index=False,
)


print("\n==========================================")
print("PROVISIONAL CLUSTER ANNOTATIONS")
print("==========================================")

print(annotation_table.to_string(index=False))


# ----------------------------------------------------------
# Add labels to AnnData
# ----------------------------------------------------------

adata.obs["cell_type"] = (
    adata.obs["leiden"]
    .astype(str)
    .map(cluster_annotation)
    .astype("category")
)


# Cell counts by annotation

cell_counts = (
    adata.obs["cell_type"]
    .value_counts()
    .rename_axis("cell_type")
    .reset_index(name="cells")
)

cell_counts["percent"] = (
    cell_counts["cells"]
    / adata.n_obs
    * 100
)

cell_counts.to_csv(
    RESULTS / "celltype_counts.tsv",
    sep="\t",
    index=False,
)


print("\n==========================================")
print("CELL-TYPE COMPOSITION")
print("==========================================")

print(cell_counts.to_string(index=False))


# ----------------------------------------------------------
# Labelled UMAP
# ----------------------------------------------------------

fig = sc.pl.umap(
    adata,
    color="cell_type",
    legend_loc="on data",
    legend_fontsize=8,
    title="PBMC 3k — canonical marker annotation",
    show=False,
    return_fig=True,
)

fig.savefig(
    FIGURES / "pbmc3k_celltype_umap.png",
    dpi=200,
    bbox_inches="tight",
)

plt.close(fig)


# ----------------------------------------------------------
# Marker expression figure
# ----------------------------------------------------------

selected_markers = [
    "CD3D",
    "IL7R",
    "NKG7",
    "GNLY",
    "MS4A1",
    "CD79A",
    "LYZ",
    "S100A8",
    "FCGR3A",
    "LST1",
    "FCER1A",
    "PPBP",
    "PF4",
]

selected_markers = [
    gene
    for gene in selected_markers
    if gene in gene_space
]

fig = sc.pl.umap(
    adata,
    color=selected_markers,
    ncols=4,
    show=False,
    return_fig=True,
)

fig.savefig(
    FIGURES / "pbmc3k_canonical_marker_umap.png",
    dpi=180,
    bbox_inches="tight",
)

plt.close(fig)


# ----------------------------------------------------------
# Save annotated object
# ----------------------------------------------------------

adata.write(
    ROOT
    / "single_cell"
    / "processed"
    / "pbmc3k_annotated.h5ad"
)


print("\nSaved:")

print(
    RESULTS
    / "cluster_celltype_annotations.tsv"
)

print(
    RESULTS
    / "cluster_celltype_scores.tsv"
)

print(
    FIGURES
    / "pbmc3k_celltype_umap.png"
)

print("\nCELL-TYPE ANNOTATION PASSED")
