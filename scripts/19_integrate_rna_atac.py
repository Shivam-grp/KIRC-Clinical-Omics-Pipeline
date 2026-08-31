"""
TCGA-KIRC RNA-seq + ATAC-seq shared regulatory interpretation.

This analysis combines two independent evidence layers:

1. Tumour-versus-normal RNA differential expression
2. KIRC ATAC-seq regulatory accessibility

ATAC groups were discovered by unsupervised clustering and therefore
do NOT represent the same contrast as tumour-versus-normal RNA-seq.

Accordingly, ATAC and RNA effect directions are NOT interpreted as
directly concordant or discordant.

Instead, the workflow identifies genes with:
    - statistically significant RNA differential expression
    - significant ATAC accessibility variation
    - ATAC peaks near the gene transcription start site

This provides an exploratory cross-modal regulatory prioritisation.
"""

from pathlib import Path
import gzip
import re
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]

ATAC_FILE = (
    ROOT
    / "atac_seq"
    / "results"
    / "significant_cluster_associated_peaks.tsv"
)

GTF_FILE = (
    ROOT
    / "regulatory_integration"
    / "reference"
    / "gencode.v22.annotation.gtf.gz"
)

RESULTS = ROOT / "regulatory_integration" / "results"
FIGURES = ROOT / "regulatory_integration" / "figures"

RESULTS.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

PROMOTER_DISTANCE = 2000
PROXIMAL_DISTANCE = 50000


# ============================================================
# Utility functions
# ============================================================

def clean_name(value):
    return re.sub(
        r"[^a-z0-9]+",
        "",
        str(value).lower()
    )


def read_table(path, nrows=None):
    """Read TSV/CSV robustly."""

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            nrows=nrows,
            low_memory=False
        )

        if df.shape[1] > 1:
            return df

    except Exception:
        pass

    return pd.read_csv(
        path,
        sep=",",
        nrows=nrows,
        low_memory=False
    )


def find_column(columns, aliases):

    normalized = {
        clean_name(col): col
        for col in columns
    }

    for alias in aliases:

        key = clean_name(alias)

        if key in normalized:
            return normalized[key]

    return None


# ============================================================
# RNA differential-expression discovery
# ============================================================

def discover_rna_de_file():
    """
    Locate the existing KIRC RNA differential-expression table.

    Handles TSV/CSV/TXT and gzip-compressed tables, including DESeq2
    files where gene IDs are stored in an unnamed first column.
    """

    search_roots = [
        Path("/mnt/e/KIRC_data/results"),
        ROOT,
    ]

    aliases_gene = [
        "gene_name",
        "gene",
        "symbol",
        "gene_symbol",
        "gene_id",
        "ensembl_gene_id",
        "feature",
    ]

    aliases_fc = [
        "log2FoldChange",
        "log2fc",
        "log2_fold_change",
        "logFC",
        "log_fold_change",
    ]

    aliases_fdr = [
        "padj",
        "FDR",
        "qvalue",
        "q_value",
        "adj_pvalue",
        "adjusted_pvalue",
        "adj.P.Val",
        "p_adj",
    ]

    candidates = []

    suffixes = (
        ".tsv",
        ".csv",
        ".txt",
        ".tsv.gz",
        ".csv.gz",
        ".txt.gz",
    )

    print("\n===== SEARCHING FOR RNA DIFFERENTIAL EXPRESSION =====")

    for root in search_roots:

        if not root.exists():
            continue

        for path in root.rglob("*"):

            if not path.is_file():
                continue

            low = path.name.lower()

            if not low.endswith(suffixes):
                continue

            full_low = str(path).lower()

            # Avoid obviously unrelated result sets
            if any(
                x in full_low
                for x in [
                    "methyl",
                    "atac",
                    "single_cell",
                    "single-cell",
                    "regulatory_integration",
                    "multiomics_integration",
                ]
            ):
                continue

            try:
                preview = pd.read_csv(
                    path,
                    sep=None,
                    engine="python",
                    nrows=20,
                    compression="infer",
                )
            except Exception:
                continue

            if preview.shape[1] < 3:
                continue

            cols = list(preview.columns)

            gene_col = find_column(cols, aliases_gene)
            fc_col = find_column(cols, aliases_fc)
            fdr_col = find_column(cols, aliases_fdr)

            # Common DESeq2 case:
            # gene IDs are in an unnamed first column.
            if gene_col is None and len(cols) > 0:
                first = cols[0]

                values = (
                    preview[first]
                    .dropna()
                    .astype(str)
                    .head(10)
                )

                if (
                    str(first).lower().startswith("unnamed")
                    or values.str.startswith("ENSG").mean() >= 0.5
                ):
                    gene_col = first

            if fc_col and fdr_col and gene_col:

                score = 0

                filename = path.name.lower()
                pathname = str(path).lower()

                if "differential" in pathname:
                    score += 10

                if "expression" in pathname:
                    score += 10

                if "deseq" in pathname:
                    score += 10

                if "rna" in pathname:
                    score += 5

                if "result" in filename:
                    score += 3

                if "significant" in filename:
                    score += 1

                candidates.append(
                    (
                        score,
                        path,
                        gene_col,
                        fc_col,
                        fdr_col,
                    )
                )

                print("\nCandidate:")
                print(path)
                print(
                    f"  gene={gene_col} | "
                    f"log2FC={fc_col} | "
                    f"FDR={fdr_col}"
                )

    if not candidates:

        print("\nNo automatic candidate found.")
        print("\nPossible RNA-related files:")

        shown = 0

        for root in search_roots:
            if not root.exists():
                continue

            for path in root.rglob("*"):
                if not path.is_file():
                    continue

                low = str(path).lower()

                if any(
                    x in low
                    for x in [
                        "rna",
                        "deseq",
                        "expression",
                        "differential",
                    ]
                ):
                    print(path)
                    shown += 1

                    if shown >= 40:
                        break

        raise FileNotFoundError(
            "RNA differential-expression table still not identified."
        )

    candidates.sort(
        key=lambda x: x[0],
        reverse=True,
    )

    best = candidates[0]

    print("\n===== SELECTED RNA DE FILE =====")
    print(best[1])
    print(f"Gene column: {best[2]}")
    print(f"log2FC column: {best[3]}")
    print(f"FDR column: {best[4]}")

    return best[1], best[2], best[3], best[4]


# ============================================================
# GENCODE parsing
# ============================================================

def get_attribute(text, key):

    match = re.search(
        rf'{key} "([^"]+)"',
        text
    )

    return (
        match.group(1)
        if match
        else None
    )


def load_genes():

    print("\nReading GENCODE gene annotation...")

    records = []

    with gzip.open(
        GTF_FILE,
        "rt"
    ) as handle:

        for line in handle:

            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            if len(fields) != 9:
                continue

            chrom = fields[0]
            feature = fields[2]

            if feature != "gene":
                continue

            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = fields[8]

            gene_id = get_attribute(
                attrs,
                "gene_id"
            )

            gene_name = get_attribute(
                attrs,
                "gene_name"
            )

            gene_type = (
                get_attribute(
                    attrs,
                    "gene_type"
                )
                or
                get_attribute(
                    attrs,
                    "gene_biotype"
                )
            )

            if not gene_id:
                continue

            if not gene_name:
                gene_name = gene_id

            # TSS depends on strand
            if strand == "+":
                tss = start
            else:
                tss = end

            records.append(
                {
                    "chrom": chrom,
                    "gene_id": gene_id.split(".")[0],
                    "gene_name": gene_name,
                    "gene_type": gene_type,
                    "strand": strand,
                    "tss": tss,
                }
            )

    genes = pd.DataFrame(records)

    genes = genes.drop_duplicates(
        ["gene_id", "gene_name", "chrom", "tss"]
    )

    print(
        f"GENCODE genes loaded: "
        f"{len(genes):,}"
    )

    return genes


# ============================================================
# ATAC parsing
# ============================================================

def parse_peak_id(value):

    value = str(value)

    match = re.search(
        r"^(chr[^:]+):(\d+)-(\d+)",
        value
    )

    if not match:
        return None

    chrom = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))

    return (
        chrom,
        start,
        end,
        (start + end) // 2,
    )


def load_atac():

    if not ATAC_FILE.exists():
        raise FileNotFoundError(
            f"ATAC result missing: {ATAC_FILE}"
        )

    atac = pd.read_csv(
        ATAC_FILE,
        sep="\t",
        low_memory=False
    )

    required = [
        "peak_id",
        "delta_accessibility",
        "FDR",
    ]

    for col in required:

        if col not in atac.columns:
            raise ValueError(
                f"ATAC column missing: {col}"
            )

    atac["FDR"] = pd.to_numeric(
        atac["FDR"],
        errors="coerce"
    )

    atac["delta_accessibility"] = pd.to_numeric(
        atac["delta_accessibility"],
        errors="coerce"
    )

    atac = atac[
        (atac["FDR"] < 0.05)
        &
        (
            atac["delta_accessibility"]
            .abs()
            >= 1
        )
    ].copy()

    coords = atac[
        "peak_id"
    ].map(parse_peak_id)

    valid = coords.notna()

    atac = atac.loc[
        valid
    ].copy()

    coords = coords.loc[
        valid
    ]

    atac["chrom"] = [
        x[0]
        for x in coords
    ]

    atac["start"] = [
        x[1]
        for x in coords
    ]

    atac["end"] = [
        x[2]
        for x in coords
    ]

    atac["midpoint"] = [
        x[3]
        for x in coords
    ]

    print("\n===== ATAC INPUT =====")

    print(
        f"Significant ATAC peaks: "
        f"{len(atac):,}"
    )

    return atac


# ============================================================
# Link ATAC peaks to nearest TSS
# ============================================================

def link_peaks_to_genes(atac, genes):

    print(
        "\nLinking ATAC peaks to nearest "
        "GENCODE transcription start sites..."
    )

    chromosome_genes = {}

    for chrom, sub in genes.groupby(
        "chrom"
    ):

        sub = (
            sub.sort_values("tss")
            .reset_index(drop=True)
        )

        chromosome_genes[chrom] = sub

    records = []

    for row in atac.itertuples():

        chrom = row.chrom

        if chrom not in chromosome_genes:
            continue

        gene_table = chromosome_genes[
            chrom
        ]

        tss_values = gene_table[
            "tss"
        ].to_numpy()

        midpoint = int(
            row.midpoint
        )

        idx = np.searchsorted(
            tss_values,
            midpoint
        )

        candidates = []

        if idx < len(
            gene_table
        ):
            candidates.append(
                idx
            )

        if idx > 0:
            candidates.append(
                idx - 1
            )

        if not candidates:
            continue

        best_index = min(
            candidates,
            key=lambda i: abs(
                int(
                    gene_table.iloc[i]["tss"]
                )
                -
                midpoint
            )
        )

        gene = gene_table.iloc[
            best_index
        ]

        distance = (
            midpoint
            -
            int(gene["tss"])
        )

        abs_distance = abs(
            distance
        )

        if abs_distance <= PROMOTER_DISTANCE:

            link_class = "promoter"

        elif abs_distance <= PROXIMAL_DISTANCE:

            link_class = "proximal"

        else:

            link_class = "distal_unlinked"

        records.append(
            {
                "peak_id": row.peak_id,
                "chrom": chrom,
                "start": row.start,
                "end": row.end,
                "midpoint": midpoint,
                "gene_id": gene["gene_id"],
                "gene_name": gene["gene_name"],
                "gene_type": gene["gene_type"],
                "tss": gene["tss"],
                "distance_to_tss": distance,
                "abs_distance_to_tss": abs_distance,
                "link_class": link_class,
                "atac_delta": row.delta_accessibility,
                "atac_fdr": row.FDR,
            }
        )

    links = pd.DataFrame(
        records
    )

    links.to_csv(
        RESULTS
        / "atac_nearest_gene_links.tsv",
        sep="\t",
        index=False
    )

    print(
        f"Peak-gene links created: "
        f"{len(links):,}"
    )

    print("\nLink classes:")

    print(
        links["link_class"]
        .value_counts()
        .to_string()
    )

    return links


# ============================================================
# RNA processing
# ============================================================

def load_rna(genes):

    (
        path,
        gene_col,
        fc_col,
        fdr_col,
    ) = discover_rna_de_file()

    rna = read_table(
        path
    )

    rna = rna[
        [
            gene_col,
            fc_col,
            fdr_col,
        ]
    ].copy()

    rna.columns = [
        "gene_raw",
        "rna_log2fc",
        "rna_fdr",
    ]

    rna["gene_raw"] = (
        rna["gene_raw"]
        .astype(str)
        .str.strip()
    )

    rna["rna_log2fc"] = pd.to_numeric(
        rna["rna_log2fc"],
        errors="coerce"
    )

    rna["rna_fdr"] = pd.to_numeric(
        rna["rna_fdr"],
        errors="coerce"
    )

    # Map Ensembl IDs to gene symbols when necessary
    ensembl_fraction = (
        rna["gene_raw"]
        .str.startswith("ENSG")
        .mean()
    )

    if ensembl_fraction > 0.5:

        gene_map = (
            genes[
                [
                    "gene_id",
                    "gene_name",
                ]
            ]
            .drop_duplicates("gene_id")
        )

        rna["gene_id"] = (
            rna["gene_raw"]
            .str.split(".")
            .str[0]
        )

        rna = rna.merge(
            gene_map,
            on="gene_id",
            how="left"
        )

        rna["gene_name"] = (
            rna["gene_name"]
            .fillna(
                rna["gene_id"]
            )
        )

    else:

        rna["gene_name"] = (
            rna["gene_raw"]
        )

    rna = rna.dropna(
        subset=[
            "rna_log2fc",
            "rna_fdr",
        ]
    )

    rna = rna.sort_values(
        "rna_fdr"
    )

    rna = rna.drop_duplicates(
        "gene_name"
    )

    significant = rna[
        (rna["rna_fdr"] < 0.05)
        &
        (
            rna["rna_log2fc"]
            .abs()
            >= 1
        )
    ].copy()

    significant.to_csv(
        RESULTS
        / "significant_rna_genes.tsv",
        sep="\t",
        index=False
    )

    print("\n===== RNA INPUT =====")

    print(
        f"RNA genes tested: "
        f"{len(rna):,}"
    )

    print(
        f"Significant RNA genes: "
        f"{len(significant):,}"
    )

    return rna, significant, path


# ============================================================
# Cross-modal integration
# ============================================================

def integrate(rna_sig, links):

    regulatory = links[
        links["link_class"].isin(
            [
                "promoter",
                "proximal",
            ]
        )
    ].copy()

    print("\n===== REGULATORY ATAC LINKS =====")

    print(
        f"Promoter/proximal peaks: "
        f"{len(regulatory):,}"
    )

    # Pick strongest ATAC evidence per gene
    regulatory[
        "abs_atac_delta"
    ] = regulatory[
        "atac_delta"
    ].abs()

    regulatory[
        "atac_strength"
    ] = (
        -np.log10(
            regulatory[
                "atac_fdr"
            ].clip(lower=1e-300)
        )
        +
        regulatory[
            "abs_atac_delta"
        ]
    )

    regulatory = (
        regulatory.sort_values(
            [
                "gene_name",
                "atac_strength",
            ],
            ascending=[
                True,
                False,
            ]
        )
        .drop_duplicates(
            "gene_name"
        )
    )

    shared = rna_sig.merge(
        regulatory,
        on="gene_name",
        how="inner"
    )

    shared[
        "rna_strength"
    ] = (
        -np.log10(
            shared[
                "rna_fdr"
            ].clip(lower=1e-300)
        )
        +
        shared[
            "rna_log2fc"
        ].abs()
    )

    shared[
        "combined_evidence_score"
    ] = (
        shared[
            "rna_strength"
        ]
        +
        shared[
            "atac_strength"
        ]
    )

    shared = shared.sort_values(
        "combined_evidence_score",
        ascending=False
    )

    shared.to_csv(
        RESULTS
        / "shared_rna_atac_candidates.tsv",
        sep="\t",
        index=False
    )

    shared.head(
        50
    ).to_csv(
        RESULTS
        / "top_50_rna_atac_candidates.tsv",
        sep="\t",
        index=False
    )

    print("\n===== SHARED REGULATORY CANDIDATES =====")

    print(
        f"Genes with RNA + ATAC evidence: "
        f"{len(shared):,}"
    )

    if len(shared):

        display_columns = [
            "gene_name",
            "rna_log2fc",
            "rna_fdr",
            "link_class",
            "distance_to_tss",
            "atac_delta",
            "atac_fdr",
            "combined_evidence_score",
        ]

        print(
            "\nTop candidates:"
        )

        print(
            shared[
                display_columns
            ]
            .head(20)
            .to_string(
                index=False
            )
        )

    return shared, regulatory


# ============================================================
# Figures
# ============================================================

def make_figures(
    rna_sig,
    regulatory,
    shared,
):

    # Evidence counts
    counts = pd.Series(
        {
            "RNA DE genes":
                len(rna_sig),

            "ATAC-linked genes":
                regulatory[
                    "gene_name"
                ].nunique(),

            "Shared genes":
                shared[
                    "gene_name"
                ].nunique(),
        }
    )

    plt.figure(
        figsize=(7, 5)
    )

    counts.plot(
        kind="bar"
    )

    plt.ylabel(
        "Number of genes"
    )

    plt.title(
        "KIRC RNA–ATAC shared regulatory evidence"
    )

    plt.xticks(
        rotation=25,
        ha="right"
    )

    plt.tight_layout()

    plt.savefig(
        FIGURES
        / "rna_atac_evidence_counts.png",
        dpi=200
    )

    plt.close()

    # Top integrated genes
    if len(shared):

        top = (
            shared.head(20)
            .sort_values(
                "combined_evidence_score"
            )
        )

        plt.figure(
            figsize=(8, 7)
        )

        plt.barh(
            top["gene_name"],
            top[
                "combined_evidence_score"
            ]
        )

        plt.xlabel(
            "Combined evidence score"
        )

        plt.ylabel(
            "Gene"
        )

        plt.title(
            "Top KIRC genes supported by RNA and ATAC evidence"
        )

        plt.tight_layout()

        plt.savefig(
            FIGURES
            / "top_rna_atac_candidates.png",
            dpi=200
        )

        plt.close()

        # Distance to nearest TSS
        plt.figure(
            figsize=(8, 5)
        )

        plt.hist(
            shared[
                "distance_to_tss"
            ],
            bins=40
        )

        plt.xlabel(
            "ATAC peak distance to gene TSS (bp)"
        )

        plt.ylabel(
            "Shared candidate peaks"
        )

        plt.title(
            "Regulatory distance distribution"
        )

        plt.tight_layout()

        plt.savefig(
            FIGURES
            / "rna_atac_tss_distance.png",
            dpi=200
        )

        plt.close()


# ============================================================
# Main
# ============================================================

def main():

    print(
        "\n============================================================"
    )

    print(
        "TCGA-KIRC RNA + ATAC REGULATORY INTERPRETATION"
    )

    print(
        "============================================================"
    )

    genes = load_genes()

    atac = load_atac()

    links = link_peaks_to_genes(
        atac,
        genes
    )

    (
        rna,
        rna_sig,
        rna_path,
    ) = load_rna(
        genes
    )

    shared, regulatory = integrate(
        rna_sig,
        links
    )

    make_figures(
        rna_sig,
        regulatory,
        shared
    )

    summary = pd.DataFrame(
        [
            [
                "RNA_DE_file",
                str(rna_path),
            ],
            [
                "RNA_genes_tested",
                len(rna),
            ],
            [
                "significant_RNA_genes",
                len(rna_sig),
            ],
            [
                "significant_ATAC_peaks",
                len(atac),
            ],
            [
                "promoter_proximal_ATAC_genes",
                regulatory[
                    "gene_name"
                ].nunique(),
            ],
            [
                "shared_RNA_ATAC_genes",
                len(shared),
            ],
        ],
        columns=[
            "metric",
            "value",
        ]
    )

    summary.to_csv(
        RESULTS
        / "rna_atac_integration_summary.tsv",
        sep="\t",
        index=False
    )

    print(
        "\n===== SAVED RESULTS ====="
    )

    for path in sorted(
        RESULTS.glob("*")
    ):
        print(path)

    print(
        "\n===== SAVED FIGURES ====="
    )

    for path in sorted(
        FIGURES.glob("*")
    ):
        print(path)

    print(
        "\n============================================================"
    )

    print(
        "RNA + ATAC REGULATORY INTEGRATION PASSED"
    )

    print(
        "============================================================"
    )


if __name__ == "__main__":
    main()
