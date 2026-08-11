"""Create publication-style differential-expression visualisations."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DATA_ROOT = Path("/mnt/e/KIRC_data")

DE_FILE = (
    DATA_ROOT
    / "results"
    / "differential_expression"
    / "tumor_vs_normal_deseq2.tsv"
)

COUNTS_FILE = (
    DATA_ROOT
    / "processed"
    / "kirc_rnaseq_filtered_counts.parquet"
)

META_FILE = (
    DATA_ROOT
    / "processed"
    / "kirc_rnaseq_sample_metadata.tsv"
)

OUT = DATA_ROOT / "results" / "differential_expression"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    de = pd.read_csv(DE_FILE, sep="\t")

    counts = pd.read_parquet(COUNTS_FILE)

    meta = pd.read_csv(
        META_FILE,
        sep="\t",
        index_col=0,
    )

    # --------------------------------------------
    # Annotated volcano plot
    # --------------------------------------------

    plot = de.dropna(
        subset=["padj", "log2FoldChange"]
    ).copy()

    plot["minus_log10_padj"] = -np.log10(
        plot["padj"].clip(lower=1e-300)
    )

    plot["direction"] = "Not significant"

    sig = (
        (plot["padj"] < 0.05)
        & (plot["log2FoldChange"].abs() >= 1)
    )

    plot.loc[
        sig & (plot["log2FoldChange"] > 0),
        "direction",
    ] = "Upregulated"

    plot.loc[
        sig & (plot["log2FoldChange"] < 0),
        "direction",
    ] = "Downregulated"

    plt.figure(figsize=(10, 7))

    for label in [
        "Not significant",
        "Upregulated",
        "Downregulated",
    ]:
        subset = plot[plot["direction"] == label]

        plt.scatter(
            subset["log2FoldChange"],
            subset["minus_log10_padj"],
            s=8,
            alpha=0.5,
            label=label,
        )

    # Label strongest protein-coding genes
    labels = (
        plot[
            sig
            & (plot["gene_type"] == "protein_coding")
        ]
        .assign(abs_fc=lambda x: x["log2FoldChange"].abs())
        .nlargest(12, "abs_fc")
    )

    for _, row in labels.iterrows():
        plt.text(
            row["log2FoldChange"],
            row["minus_log10_padj"],
            str(row["gene_name"]),
            fontsize=8,
        )

    plt.axvline(1, linestyle="--", linewidth=1)
    plt.axvline(-1, linestyle="--", linewidth=1)
    plt.axhline(
        -np.log10(0.05),
        linestyle="--",
        linewidth=1,
    )

    plt.xlabel("log2 fold change: Tumour vs Normal")
    plt.ylabel("-log10 adjusted p-value")
    plt.title("TCGA-KIRC Differential Expression")
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        OUT / "annotated_volcano_plot.png",
        dpi=250,
    )

    plt.close()

    # --------------------------------------------
    # Top protein-coding genes
    # --------------------------------------------

    protein = de[
        (de["padj"] < 0.05)
        & (de["gene_type"] == "protein_coding")
    ].copy()

    top_up = protein.nlargest(
        10,
        "log2FoldChange",
    )

    top_down = protein.nsmallest(
        10,
        "log2FoldChange",
    )

    top = pd.concat([top_up, top_down])

    gene_ids = top["gene_id"].tolist()

    available = [
        g for g in gene_ids
        if g in counts.columns
    ]

    expression = counts[available].copy()

    # log counts for visualisation only
    expression = np.log2(expression + 1)

    # Z-score by gene
    expression = (
        expression - expression.mean(axis=0)
    ) / expression.std(axis=0)

    gene_map = dict(
        zip(top["gene_id"], top["gene_name"])
    )

    expression = expression.rename(
        columns=gene_map
    )

    # order normal then tumour
    sample_order = meta.sort_values(
        "sample_type",
        ascending=False,
    ).index

    expression = expression.loc[sample_order]

    # transpose genes x samples
    heat = expression.T

    plt.figure(figsize=(14, 7))

    image = plt.imshow(
        heat,
        aspect="auto",
        interpolation="nearest",
    )

    plt.yticks(
        range(len(heat.index)),
        heat.index,
        fontsize=8,
    )

    plt.xlabel("Samples")
    plt.ylabel("Genes")
    plt.title(
        "Top Differentially Expressed Protein-Coding Genes"
    )

    plt.colorbar(
        image,
        label="Gene-wise Z-score",
    )

    plt.tight_layout()

    plt.savefig(
        OUT / "top_de_gene_heatmap.png",
        dpi=250,
    )

    plt.close()

    # --------------------------------------------
    # ccRCC biology check
    # --------------------------------------------

    genes_of_interest = [
        "VHL",
        "CA9",
        "VEGFA",
        "EGLN3",
        "SLC2A1",
        "FABP7",
        "CD70",
        "UMOD",
        "AQP2",
        "SLC12A1",
    ]

    check = de[
        de["gene_name"].isin(genes_of_interest)
    ][
        [
            "gene_name",
            "baseMean",
            "log2FoldChange",
            "padj",
        ]
    ].sort_values("log2FoldChange", ascending=False)

    print("\nSelected ccRCC / kidney biology genes")
    print("=" * 65)
    print(check.to_string(index=False))

    check.to_csv(
        OUT / "ccrcc_marker_check.tsv",
        sep="\t",
        index=False,
    )

    print()
    print("Saved:")
    print(OUT / "annotated_volcano_plot.png")
    print(OUT / "top_de_gene_heatmap.png")
    print(OUT / "ccrcc_marker_check.tsv")

    print()
    print("DE VISUALISATION PASSED")


if __name__ == "__main__":
    main()
