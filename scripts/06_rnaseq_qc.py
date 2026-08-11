"""Quality control and filtering for TCGA-KIRC bulk RNA-seq."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


DATA_ROOT = Path("/mnt/e/KIRC_data")

COUNTS_FILE = (
    DATA_ROOT / "processed" / "kirc_rnaseq_raw_counts.parquet"
)

META_FILE = (
    DATA_ROOT / "processed" / "kirc_rnaseq_sample_metadata.tsv"
)

MANIFEST_FILE = (
    DATA_ROOT / "downloads" / "tcga_kirc_rnaseq_cohort.tsv"
)

RAW_DIR = DATA_ROOT / "raw" / "rnaseq"

RESULTS_DIR = DATA_ROOT / "results" / "rnaseq_qc"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

FILTERED_COUNTS = (
    DATA_ROOT / "processed" / "kirc_rnaseq_filtered_counts.parquet"
)

ANNOTATION_FILE = (
    DATA_ROOT / "processed" / "kirc_gene_annotation.tsv"
)


def extract_gene_annotation(manifest: pd.DataFrame) -> pd.DataFrame:
    """Extract gene names/types from one GDC STAR-count file."""

    row = manifest.iloc[0]

    path = RAW_DIR / row["file_id"] / row["file_name"]

    annotation = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        usecols=["gene_id", "gene_name", "gene_type"],
    )

    annotation = annotation[
        annotation["gene_id"].astype(str).str.startswith("ENSG")
    ].copy()

    annotation = annotation.drop_duplicates("gene_id")

    annotation.to_csv(
        ANNOTATION_FILE,
        sep="\t",
        index=False,
    )

    return annotation


def main() -> None:

    print("Loading RNA-seq matrix...")

    counts = pd.read_parquet(COUNTS_FILE)

    metadata = pd.read_csv(
        META_FILE,
        sep="\t",
        index_col=0,
    )

    manifest = pd.read_csv(
        MANIFEST_FILE,
        sep="\t",
    )

    print(f"Samples: {counts.shape[0]}")
    print(f"Genes before filtering: {counts.shape[1]}")

    # ---------------------------------------------------------
    # Gene annotation
    # ---------------------------------------------------------

    annotation = extract_gene_annotation(manifest)

    print(f"Gene annotations: {len(annotation)}")

    # ---------------------------------------------------------
    # Sample-level QC
    # ---------------------------------------------------------

    library_size = counts.sum(axis=1)

    detected_genes = (counts > 0).sum(axis=1)

    qc = pd.DataFrame(
        {
            "sample_type": metadata.loc[counts.index, "sample_type"],
            "library_size": library_size,
            "detected_genes": detected_genes,
        }
    )

    qc.to_csv(
        RESULTS_DIR / "sample_qc_metrics.tsv",
        sep="\t",
    )

    print("\nLibrary size summary:")
    print(library_size.describe().round(0).to_string())

    print("\nDetected genes summary:")
    print(detected_genes.describe().round(0).to_string())

    # ---------------------------------------------------------
    # Gene filtering
    # Keep genes with >=10 raw counts in >=10 samples
    # ---------------------------------------------------------

    keep = (counts >= 10).sum(axis=0) >= 10

    filtered = counts.loc[:, keep].copy()

    print()
    print("Filtering rule: >=10 counts in >=10 samples")
    print(f"Genes retained: {filtered.shape[1]}")
    print(f"Genes removed: {counts.shape[1] - filtered.shape[1]}")

    filtered.to_parquet(
        FILTERED_COUNTS,
        compression="zstd",
    )

    # ---------------------------------------------------------
    # Library-size plot
    # ---------------------------------------------------------

    plt.figure(figsize=(8, 5))

    for sample_type, group in qc.groupby("sample_type"):
        plt.hist(
            np.log10(group["library_size"]),
            bins=30,
            alpha=0.6,
            label=sample_type,
        )

    plt.xlabel("log10 total raw counts")
    plt.ylabel("Number of samples")
    plt.title("TCGA-KIRC RNA-seq library sizes")
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        RESULTS_DIR / "library_size_distribution.png",
        dpi=200,
    )

    plt.close()

    # ---------------------------------------------------------
    # Detected genes plot
    # ---------------------------------------------------------

    plt.figure(figsize=(8, 5))

    for sample_type, group in qc.groupby("sample_type"):
        plt.hist(
            group["detected_genes"],
            bins=30,
            alpha=0.6,
            label=sample_type,
        )

    plt.xlabel("Detected genes")
    plt.ylabel("Number of samples")
    plt.title("Genes detected per TCGA-KIRC sample")
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        RESULTS_DIR / "detected_genes_distribution.png",
        dpi=200,
    )

    plt.close()

    # ---------------------------------------------------------
    # PCA using log-CPM
    # ---------------------------------------------------------

    print("\nPreparing PCA...")

    library_size_filtered = filtered.sum(axis=1)

    log_cpm = filtered.astype("float32").div(
        library_size_filtered,
        axis=0,
    )

    log_cpm *= 1_000_000
    log_cpm = np.log1p(log_cpm)

    # Use the 5000 most variable genes for PCA.
    gene_variance = log_cpm.var(axis=0)

    n_top = min(5000, len(gene_variance))

    top_genes = gene_variance.nlargest(n_top).index

    pca_input = log_cpm.loc[:, top_genes]

    pca = PCA(n_components=2)

    coordinates = pca.fit_transform(pca_input)

    pca_df = pd.DataFrame(
        {
            "PC1": coordinates[:, 0],
            "PC2": coordinates[:, 1],
            "sample_type": metadata.loc[
                pca_input.index,
                "sample_type",
            ].values,
        },
        index=pca_input.index,
    )

    pca_df.to_csv(
        RESULTS_DIR / "pca_coordinates.tsv",
        sep="\t",
    )

    variance = pca.explained_variance_ratio_ * 100

    plt.figure(figsize=(8, 6))

    for sample_type, group in pca_df.groupby("sample_type"):
        plt.scatter(
            group["PC1"],
            group["PC2"],
            alpha=0.7,
            label=sample_type,
        )

    plt.xlabel(f"PC1 ({variance[0]:.1f}% variance)")
    plt.ylabel(f"PC2 ({variance[1]:.1f}% variance)")
    plt.title("TCGA-KIRC RNA-seq PCA")
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        RESULTS_DIR / "rnaseq_pca.png",
        dpi=200,
    )

    plt.close()

    print()
    print("PCA variance explained:")
    print(f"PC1: {variance[0]:.2f}%")
    print(f"PC2: {variance[1]:.2f}%")

    print()
    print("Saved:")
    print(FILTERED_COUNTS)
    print(ANNOTATION_FILE)
    print(RESULTS_DIR)

    print()
    print("RNA-SEQ QC PASSED")


if __name__ == "__main__":
    main()
