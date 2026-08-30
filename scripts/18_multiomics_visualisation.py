"""
Create publication/portfolio-quality figures for the
TCGA-KIRC RNA-seq × DNA methylation integration.

Outputs:
1. RNA vs promoter methylation integration scatter
2. Top concordant epigenetic candidates
3. Selected ccRCC/kidney marker figure
4. Integration category summary
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# PATHS
# ============================================================

BASE = Path(
    "/mnt/e/KIRC_data/results/multiomics_integration"
)

INTEGRATED = (
    BASE /
    "rna_methylation_integrated.tsv"
)

CANDIDATES = (
    BASE /
    "concordant_epigenetic_candidates.tsv"
)

STRICT = (
    BASE /
    "strict_epigenetic_candidates.tsv"
)

FIG_DIR = (
    BASE /
    "figures"
)

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# SETTINGS
# ============================================================

KNOWN_GENES = {
    "VHL",
    "CA9",
    "VEGFA",
    "SLC2A1",
    "AQP2",
    "UMOD",
    "EGLN3",
    "FABP7",
    "KNG1",
    "CD70",
}


# ============================================================
# LOAD
# ============================================================

def load():

    integrated = pd.read_csv(
        INTEGRATED,
        sep="\t"
    )

    candidates = pd.read_csv(
        CANDIDATES,
        sep="\t"
    )

    strict = pd.read_csv(
        STRICT,
        sep="\t"
    )

    return (
        integrated,
        candidates,
        strict
    )


# ============================================================
# FIGURE 1
# RNA EXPRESSION vs PROMOTER METHYLATION
# ============================================================

def integration_scatter(
    integrated,
    strict
):

    fig, ax = plt.subplots(
        figsize=(10, 8)
    )

    # Other genes
    other = integrated[
        integrated[
            "integration_class"
        ] == "Other"
    ]

    ax.scatter(
        other[
            "promoter_mean_delta_beta"
        ],
        other[
            "rna_log2fc"
        ],
        alpha=0.25,
        s=20,
        label="Other genes"
    )

    # Silencing candidates
    silencing = integrated[
        integrated[
            "integration_class"
        ]
        ==
        "Candidate epigenetic silencing"
    ]

    ax.scatter(
        silencing[
            "promoter_mean_delta_beta"
        ],
        silencing[
            "rna_log2fc"
        ],
        alpha=0.7,
        s=35,
        label="Candidate silencing"
    )

    # Activation candidates
    activation = integrated[
        integrated[
            "integration_class"
        ]
        ==
        "Candidate epigenetic activation"
    ]

    ax.scatter(
        activation[
            "promoter_mean_delta_beta"
        ],
        activation[
            "rna_log2fc"
        ],
        alpha=0.7,
        s=35,
        label="Candidate activation"
    )

    # Zero lines
    ax.axhline(
        0,
        linewidth=1
    )

    ax.axvline(
        0,
        linewidth=1
    )

    # Highlight top strict candidates
    top = (
        strict
        .sort_values(
            "ranking_score",
            ascending=False
        )
        .head(15)
    )

    for row in top.itertuples():

        ax.annotate(
            row.gene_name,
            (
                row.promoter_mean_delta_beta,
                row.rna_log2fc
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8
        )

    ax.set_xlabel(
        "Promoter methylation change (Δβ)\nTumour − Normal"
    )

    ax.set_ylabel(
        "RNA expression change (log2 fold change)\nTumour vs Normal"
    )

    ax.set_title(
        "TCGA-KIRC RNA Expression × Promoter DNA Methylation"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    out = (
        FIG_DIR /
        "rna_methylation_integration_scatter.png"
    )

    fig.savefig(
        out,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(out)


# ============================================================
# FIGURE 2
# TOP STRICT CANDIDATES
# ============================================================

def top_candidates_plot(
    strict
):

    if strict.empty:

        print(
            "No strict candidates available."
        )

        return

    top = (
        strict
        .sort_values(
            "ranking_score",
            ascending=False
        )
        .head(20)
        .copy()
    )

    top = top.sort_values(
        "ranking_score"
    )

    fig, ax = plt.subplots(
        figsize=(10, 9)
    )

    ax.barh(
        top["gene_name"],
        top["ranking_score"]
    )

    ax.set_xlabel(
        "Multi-omics candidate ranking score"
    )

    ax.set_ylabel(
        "Gene"
    )

    ax.set_title(
        "Top Concordant Epigenetic Candidates in TCGA-KIRC"
    )

    fig.tight_layout()

    out = (
        FIG_DIR /
        "top_epigenetic_candidates.png"
    )

    fig.savefig(
        out,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(out)


# ============================================================
# FIGURE 3
# KNOWN ccRCC / KIDNEY BIOLOGY
# ============================================================

def kidney_marker_plot(
    integrated
):

    marker = integrated[
        integrated[
            "gene_name"
        ].isin(
            KNOWN_GENES
        )
    ].copy()

    if marker.empty:

        print(
            "No selected kidney markers found."
        )

        return

    marker = marker.sort_values(
        "rna_log2fc"
    )

    genes = marker[
        "gene_name"
    ].tolist()

    y = np.arange(
        len(genes)
    )

    fig, ax = plt.subplots(
        figsize=(10, 7)
    )

    ax.scatter(
        marker[
            "promoter_mean_delta_beta"
        ],
        y,
        s=90,
        label="Promoter Δβ"
    )

    ax.scatter(
        marker[
            "rna_log2fc"
        ],
        y,
        s=90,
        marker="x",
        label="RNA log2FC"
    )

    ax.axvline(
        0,
        linewidth=1
    )

    ax.set_yticks(
        y
    )

    ax.set_yticklabels(
        genes
    )

    ax.set_xlabel(
        "Effect size"
    )

    ax.set_title(
        "Selected ccRCC / Kidney Biology Genes"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    out = (
        FIG_DIR /
        "ccrcc_marker_multiomics.png"
    )

    fig.savefig(
        out,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(out)


# ============================================================
# FIGURE 4
# INTEGRATION CATEGORY COUNTS
# ============================================================

def category_plot(
    integrated
):

    counts = (
        integrated[
            "integration_class"
        ]
        .value_counts()
    )

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    ax.bar(
        counts.index,
        counts.values
    )

    ax.set_ylabel(
        "Number of genes"
    )

    ax.set_title(
        "RNA–Methylation Integration Categories"
    )

    ax.tick_params(
        axis="x",
        rotation=20
    )

    for i, value in enumerate(
        counts.values
    ):

        ax.text(
            i,
            value,
            f"{value:,}",
            ha="center",
            va="bottom"
        )

    fig.tight_layout()

    out = (
        FIG_DIR /
        "integration_category_counts.png"
    )

    fig.savefig(
        out,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(out)


# ============================================================
# BIOLOGICAL SUMMARY TABLE
# ============================================================

def save_top_table(
    strict
):

    columns = [
        "gene_name",
        "integration_class",
        "promoter_direction",
        "promoter_dmps",
        "promoter_mean_delta_beta",
        "rna_direction",
        "rna_log2fc",
        "rna_fdr",
        "ranking_score"
    ]

    available = [
        c
        for c in columns
        if c in strict.columns
    ]

    top = (
        strict
        .sort_values(
            "ranking_score",
            ascending=False
        )
        .head(50)
    )

    out = (
        BASE /
        "top_50_multiomics_candidates.tsv"
    )

    top[
        available
    ].to_csv(
        out,
        sep="\t",
        index=False
    )

    print(out)


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print("=" * 60)

    print(
        "TCGA-KIRC MULTI-OMICS VISUALISATION"
    )

    print("=" * 60)

    (
        integrated,
        candidates,
        strict
    ) = load()

    print()
    print(
        f"Integrated genes: "
        f"{len(integrated):,}"
    )

    print(
        f"Concordant candidates: "
        f"{len(candidates):,}"
    )

    print(
        f"Strict candidates: "
        f"{len(strict):,}"
    )

    print()
    print("Creating figures...")

    integration_scatter(
        integrated,
        strict
    )

    top_candidates_plot(
        strict
    )

    kidney_marker_plot(
        integrated
    )

    category_plot(
        integrated
    )

    save_top_table(
        strict
    )

    print()
    print("=" * 60)

    print(
        "MULTI-OMICS VISUALISATION PASSED"
    )

    print("=" * 60)


if __name__ == "__main__":
    main()
