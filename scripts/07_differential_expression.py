"""Differential expression: TCGA-KIRC Primary Tumor vs Solid Tissue Normal."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


DATA_ROOT = Path("/mnt/e/KIRC_data")

COUNTS_FILE = (
    DATA_ROOT / "processed" / "kirc_rnaseq_filtered_counts.parquet"
)

META_FILE = (
    DATA_ROOT / "processed" / "kirc_rnaseq_sample_metadata.tsv"
)

ANNOTATION_FILE = (
    DATA_ROOT / "processed" / "kirc_gene_annotation.tsv"
)

RESULTS_DIR = DATA_ROOT / "results" / "differential_expression"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def main():

    print("Loading filtered counts...")

    counts = pd.read_parquet(COUNTS_FILE)

    metadata = pd.read_csv(
        META_FILE,
        sep="\t",
        index_col=0,
    )

    metadata = metadata.loc[counts.index].copy()

    # Rename for simpler modelling.
    metadata["condition"] = metadata["sample_type"].replace(
        {
            "Primary Tumor": "Tumor",
            "Solid Tissue Normal": "Normal",
        }
    )

    metadata["condition"] = pd.Categorical(
        metadata["condition"],
        categories=["Normal", "Tumor"],
    )

    print(f"Samples: {counts.shape[0]}")
    print(f"Genes: {counts.shape[1]}")

    print("\nConditions:")
    print(metadata["condition"].value_counts())

    # PyDESeq2 requires integer raw counts.
    counts = counts.astype("int64")

    print("\nInitialising PyDESeq2...")

    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design="~condition",
        refit_cooks=True,
        n_cpus=2,
        low_memory=True,
    )

    print("\nRunning DESeq2 model...")
    print("This may take some time on an 8 GB laptop.")

    dds.deseq2()

    print("\nRunning Tumor vs Normal statistical test...")

    stats = DeseqStats(
        dds,
        contrast=["condition", "Tumor", "Normal"],
        alpha=0.05,
        n_cpus=2,
    )

    stats.summary()

    results = stats.results_df.copy()

    results.index.name = "gene_id"
    results = results.reset_index()

    # Add gene annotation.
    annotation = pd.read_csv(
        ANNOTATION_FILE,
        sep="\t",
    )

    results = results.merge(
        annotation,
        on="gene_id",
        how="left",
    )

    # Define significance.
    results["significant"] = (
        (results["padj"] < 0.05)
        & (results["log2FoldChange"].abs() >= 1)
    )

    results["direction"] = "Not significant"

    results.loc[
        results["significant"]
        & (results["log2FoldChange"] >= 1),
        "direction",
    ] = "Upregulated"

    results.loc[
        results["significant"]
        & (results["log2FoldChange"] <= -1),
        "direction",
    ] = "Downregulated"

    results = results.sort_values(
        ["padj", "pvalue"],
        na_position="last",
    )

    results.to_csv(
        RESULTS_DIR / "tumor_vs_normal_deseq2.tsv",
        sep="\t",
        index=False,
    )

    # -----------------------------------------------------
    # Summary
    # -----------------------------------------------------

    up = (
        results["direction"] == "Upregulated"
    ).sum()

    down = (
        results["direction"] == "Downregulated"
    ).sum()

    print()
    print("Differential-expression summary")
    print("--------------------------------")
    print(f"Upregulated genes:   {up}")
    print(f"Downregulated genes: {down}")
    print(
        f"Total significant:   {up + down}"
    )

    # -----------------------------------------------------
    # Volcano plot
    # -----------------------------------------------------

    plot_df = results.dropna(
        subset=["padj", "log2FoldChange"]
    ).copy()

    plot_df["minus_log10_padj"] = (
        -np.log10(plot_df["padj"].clip(lower=1e-300))
    )

    plt.figure(figsize=(9, 6))

    groups = [
        ("Not significant", "grey"),
        ("Upregulated", "red"),
        ("Downregulated", "blue"),
    ]

    for label, color in groups:
        subset = plot_df[
            plot_df["direction"] == label
        ]

        plt.scatter(
            subset["log2FoldChange"],
            subset["minus_log10_padj"],
            s=8,
            alpha=0.5,
            label=label,
            color=color,
        )

    plt.axvline(1, linestyle="--", linewidth=1)
    plt.axvline(-1, linestyle="--", linewidth=1)

    plt.axhline(
        -np.log10(0.05),
        linestyle="--",
        linewidth=1,
    )

    plt.xlabel("log2 fold change (Tumor vs Normal)")
    plt.ylabel("-log10 adjusted p-value")
    plt.title("TCGA-KIRC Differential Expression")
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        RESULTS_DIR / "volcano_plot.png",
        dpi=200,
    )

    plt.close()

    # -----------------------------------------------------
    # MA plot
    # -----------------------------------------------------

    ma_df = results.dropna(
        subset=["baseMean", "log2FoldChange"]
    )

    plt.figure(figsize=(9, 6))

    plt.scatter(
        np.log10(ma_df["baseMean"] + 1),
        ma_df["log2FoldChange"],
        s=8,
        alpha=0.4,
    )

    plt.axhline(0, linestyle="--")
    plt.axhline(1, linestyle="--")
    plt.axhline(-1, linestyle="--")

    plt.xlabel("log10 mean normalized expression")
    plt.ylabel("log2 fold change")
    plt.title("TCGA-KIRC MA Plot")

    plt.tight_layout()

    plt.savefig(
        RESULTS_DIR / "ma_plot.png",
        dpi=200,
    )

    plt.close()

    # Save top DE genes.
    top = results[
        results["significant"]
    ].head(50)

    top.to_csv(
        RESULTS_DIR / "top_50_de_genes.tsv",
        sep="\t",
        index=False,
    )

    print()
    print("Saved results to:")
    print(RESULTS_DIR)

    print()
    print("DIFFERENTIAL EXPRESSION PASSED")


if __name__ == "__main__":
    main()
