"""
Integrate TCGA-KIRC differential DNA methylation with RNA-seq
differential expression.

Main goals
----------
1. Extract promoter-associated significant methylation changes.
2. Map promoter CpGs to their corresponding genes.
3. Combine gene-level promoter methylation with RNA-seq DE.
4. Identify:
   - promoter hypermethylation + RNA downregulation
   - promoter hypomethylation + RNA upregulation
5. Save reproducible multi-omics candidate tables.
"""

from pathlib import Path
import math

import numpy as np
import pandas as pd


# ============================================================
# PATHS
# ============================================================

METH_FILE = Path(
    "/mnt/e/KIRC_data/results/"
    "differential_methylation/annotation/"
    "kirc_significant_dmps_annotated.tsv"
)

DE_DIR = Path(
    "/mnt/e/KIRC_data/results/differential_expression"
)

GENE_ANNOTATION = Path(
    "/mnt/e/KIRC_data/processed/"
    "kirc_gene_annotation.tsv"
)

OUT_DIR = Path(
    "/mnt/e/KIRC_data/results/"
    "multiomics_integration"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

PROMOTER_LINKS_OUT = (
    OUT_DIR /
    "promoter_dmp_gene_links.tsv"
)

METH_GENE_OUT = (
    OUT_DIR /
    "gene_promoter_methylation_summary.tsv"
)

INTEGRATED_OUT = (
    OUT_DIR /
    "rna_methylation_integrated.tsv"
)

CANDIDATE_OUT = (
    OUT_DIR /
    "concordant_epigenetic_candidates.tsv"
)

STRICT_OUT = (
    OUT_DIR /
    "strict_epigenetic_candidates.tsv"
)

SUMMARY_OUT = (
    OUT_DIR /
    "integration_summary.txt"
)


# ============================================================
# THRESHOLDS
# ============================================================

DE_FDR = 0.05
DE_LOG2FC = 1.0

PROMOTER_REGIONS = {
    "TSS1500",
    "TSS200",
    "5'UTR",
    "1stExon",
}


# ============================================================
# UTILITIES
# ============================================================

def delimiter(path):
    if path.suffix.lower() == ".csv":
        return ","

    return "\t"


def normalise_gene_id(series):
    """
    Remove Ensembl version suffix:
    ENSG000001234.15 -> ENSG000001234
    """

    return (
        series
        .astype(str)
        .str.strip()
        .str.replace(
            r"\.\d+$",
            "",
            regex=True
        )
    )


def find_column(columns, possibilities):

    lower_map = {
        str(c).lower():
        c
        for c in columns
    }

    for name in possibilities:

        if name.lower() in lower_map:
            return lower_map[
                name.lower()
            ]

    return None


# ============================================================
# FIND RNA-SEQ DE TABLE
# ============================================================

def find_de_table():

    print()
    print("Searching for RNA-seq differential-expression table...")

    if not DE_DIR.exists():
        raise RuntimeError(
            f"Directory not found:\n{DE_DIR}"
        )

    files = (
        list(DE_DIR.rglob("*.tsv"))
        +
        list(DE_DIR.rglob("*.csv"))
    )

    if not files:
        raise RuntimeError(
            "No TSV/CSV files found in "
            f"{DE_DIR}"
        )

    valid = []

    for path in files:

        try:

            header = pd.read_csv(
                path,
                sep=delimiter(path),
                nrows=3
            )

        except Exception:
            continue

        cols = header.columns

        lfc = find_column(
            cols,
            [
                "log2FoldChange",
                "log2_fold_change",
                "log2fc",
                "logFC",
            ]
        )

        padj = find_column(
            cols,
            [
                "padj",
                "fdr",
                "adj_pvalue",
                "adjusted_pvalue",
            ]
        )

        gene_name = find_column(
            cols,
            [
                "gene_name",
                "gene_symbol",
                "symbol",
            ]
        )

        gene_id = find_column(
            cols,
            [
                "gene_id",
                "ensembl_gene_id",
            ]
        )

        if (
            lfc is not None
            and padj is not None
            and (
                gene_name is not None
                or gene_id is not None
            )
        ):

            score = 0

            text = path.name.lower()

            for keyword in [
                "differential",
                "deseq",
                "result",
                "all",
            ]:
                if keyword in text:
                    score += 1

            valid.append(
                (
                    score,
                    path,
                    lfc,
                    padj,
                    gene_name,
                    gene_id,
                )
            )

    if not valid:

        print()
        print("Available files:")

        for path in files:
            print(path)

        raise RuntimeError(
            "Could not automatically identify "
            "the RNA-seq DE results table."
        )

    valid.sort(
        key=lambda x: x[0],
        reverse=True
    )

    selected = valid[0]

    (
        _,
        path,
        lfc_col,
        padj_col,
        gene_name_col,
        gene_id_col,
    ) = selected

    print(f"Selected RNA-seq file:")
    print(path)

    return (
        path,
        lfc_col,
        padj_col,
        gene_name_col,
        gene_id_col,
    )


# ============================================================
# LOAD RNA-SEQ RESULTS
# ============================================================

def load_rna():

    (
        path,
        lfc_col,
        padj_col,
        gene_name_col,
        gene_id_col,
    ) = find_de_table()

    rna = pd.read_csv(
        path,
        sep=delimiter(path)
    )

    rename = {
        lfc_col:
            "rna_log2fc",

        padj_col:
            "rna_fdr",
    }

    if gene_name_col is not None:
        rename[
            gene_name_col
        ] = "gene_name"

    if gene_id_col is not None:
        rename[
            gene_id_col
        ] = "gene_id"

    rna = rna.rename(
        columns=rename
    )

    # --------------------------------------------------------
    # Attach gene symbols if RNA table contains only Ensembl IDs
    # --------------------------------------------------------

    if (
        "gene_name" not in rna.columns
        and
        "gene_id" in rna.columns
    ):

        if not GENE_ANNOTATION.exists():

            raise RuntimeError(
                "RNA results contain no gene_name "
                "and gene annotation file was not found."
            )

        annotation = pd.read_csv(
            GENE_ANNOTATION,
            sep="\t"
        )

        ann_gene_id = find_column(
            annotation.columns,
            [
                "gene_id",
                "ensembl_gene_id",
            ]
        )

        ann_gene_name = find_column(
            annotation.columns,
            [
                "gene_name",
                "gene_symbol",
                "symbol",
            ]
        )

        if (
            ann_gene_id is None
            or ann_gene_name is None
        ):

            raise RuntimeError(
                "Could not identify gene ID/name "
                "columns in RNA gene annotation."
            )

        annotation = annotation.rename(
            columns={
                ann_gene_id:
                    "gene_id_annotation",

                ann_gene_name:
                    "gene_name",
            }
        )

        rna[
            "_gene_key"
        ] = normalise_gene_id(
            rna["gene_id"]
        )

        annotation[
            "_gene_key"
        ] = normalise_gene_id(
            annotation[
                "gene_id_annotation"
            ]
        )

        annotation = (
            annotation[
                [
                    "_gene_key",
                    "gene_name",
                ]
            ]
            .drop_duplicates(
                "_gene_key"
            )
        )

        rna = rna.merge(
            annotation,
            on="_gene_key",
            how="left"
        )

    if "gene_name" not in rna.columns:

        raise RuntimeError(
            "RNA-seq gene symbols could not be obtained."
        )

    rna["gene_name"] = (
        rna["gene_name"]
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

    rna = rna[
        rna["gene_name"].notna()
    ].copy()

    rna = rna[
        rna["gene_name"] != ""
    ].copy()

    rna = rna[
        rna["gene_name"].str.lower()
        != "nan"
    ].copy()

    # One gene symbol per row for integration
    rna = (
        rna
        .sort_values(
            "rna_fdr"
        )
        .drop_duplicates(
            "gene_name"
        )
    )

    rna[
        "rna_significant"
    ] = (
        (rna["rna_fdr"] < DE_FDR)
        &
        (
            rna["rna_log2fc"].abs()
            >= DE_LOG2FC
        )
    )

    rna["rna_direction"] = (
        "Not significant"
    )

    rna.loc[
        rna["rna_significant"]
        &
        (
            rna["rna_log2fc"] > 0
        ),
        "rna_direction"
    ] = "Upregulated"

    rna.loc[
        rna["rna_significant"]
        &
        (
            rna["rna_log2fc"] < 0
        ),
        "rna_direction"
    ] = "Downregulated"

    print()
    print("RNA-SEQ SUMMARY")
    print("=" * 35)

    print(
        f"Genes available: "
        f"{len(rna):,}"
    )

    print(
        f"Significant RNA genes "
        f"(FDR<{DE_FDR}, |log2FC|>={DE_LOG2FC}): "
        f"{rna['rna_significant'].sum():,}"
    )

    print()
    print(
        rna["rna_direction"]
        .value_counts()
        .to_string()
    )

    return rna


# ============================================================
# CORRECTLY PAIR GENE NAMES AND GENE REGIONS
# ============================================================

def build_promoter_links(meth):

    records = []

    for row in meth.itertuples(
        index=False
    ):

        gene_names = getattr(
            row,
            "gene_names",
            None
        )

        gene_regions = getattr(
            row,
            "gene_region",
            None
        )

        if pd.isna(gene_names):
            continue

        names = [
            x.strip()
            for x in str(
                gene_names
            ).split(";")
        ]

        if pd.isna(gene_regions):

            regions = [
                ""
            ] * len(names)

        else:

            regions = [
                x.strip()
                for x in str(
                    gene_regions
                ).split(";")
            ]

        # Manifest normally stores corresponding
        # gene/region entries positionally.
        # Pad if lengths differ.

        if len(regions) < len(names):

            regions += (
                [""]
                *
                (
                    len(names)
                    -
                    len(regions)
                )
            )

        if len(names) < len(regions):

            names += (
                [""]
                *
                (
                    len(regions)
                    -
                    len(names)
                )
            )

        for name, region in zip(
            names,
            regions
        ):

            if not name:
                continue

            region_parts = {
                x.strip()
                for x in str(
                    region
                ).split(";")
                if x.strip()
            }

            is_promoter = (
                bool(
                    region_parts
                    &
                    PROMOTER_REGIONS
                )
            )

            # Also handle joined annotations
            if not is_promoter:

                is_promoter = any(
                    promoter
                    in str(region)
                    for promoter
                    in PROMOTER_REGIONS
                )

            if not is_promoter:
                continue

            records.append(
                {
                    "cpg_id":
                        row.cpg_id,

                    "gene_name":
                        name,

                    "gene_region":
                        region,

                    "direction":
                        row.direction,

                    "delta_beta":
                        row.delta_beta,

                    "methylation_fdr":
                        row.fdr,

                    "mean_beta_tumour":
                        row.mean_beta_tumour,

                    "mean_beta_normal":
                        row.mean_beta_normal,
                }
            )

    links = pd.DataFrame(
        records
    )

    if links.empty:

        raise RuntimeError(
            "No promoter-associated "
            "CpG-gene links were identified."
        )

    links["gene_name"] = (
        links["gene_name"]
        .astype(str)
        .str.strip()
    )

    links = links.drop_duplicates(
        [
            "cpg_id",
            "gene_name",
            "gene_region",
        ]
    )

    links.to_csv(
        PROMOTER_LINKS_OUT,
        sep="\t",
        index=False
    )

    print()
    print("PROMOTER METHYLATION")
    print("=" * 35)

    print(
        f"Promoter CpG-gene links: "
        f"{len(links):,}"
    )

    print(
        f"Genes with promoter DMPs: "
        f"{links['gene_name'].nunique():,}"
    )

    return links


# ============================================================
# GENE-LEVEL METHYLATION
# ============================================================

def summarise_promoter_methylation(
    links
):

    def classify(group):

        hyper = int(
            (
                group["direction"]
                ==
                "Hypermethylated"
            ).sum()
        )

        hypo = int(
            (
                group["direction"]
                ==
                "Hypomethylated"
            ).sum()
        )

        mean_delta = (
            group["delta_beta"]
            .mean()
        )

        if hyper > 0 and hypo == 0:
            direction = "Hypermethylated"

        elif hypo > 0 and hyper == 0:
            direction = "Hypomethylated"

        elif mean_delta > 0:
            direction = "Mixed_hyper_dominant"

        elif mean_delta < 0:
            direction = "Mixed_hypo_dominant"

        else:
            direction = "Mixed"

        return pd.Series(
            {
                "promoter_dmps":
                    group[
                        "cpg_id"
                    ].nunique(),

                "promoter_hyper_dmps":
                    hyper,

                "promoter_hypo_dmps":
                    hypo,

                "promoter_mean_delta_beta":
                    mean_delta,

                "promoter_max_abs_delta_beta":
                    group[
                        "delta_beta"
                    ].abs().max(),

                "promoter_min_fdr":
                    group[
                        "methylation_fdr"
                    ].min(),

                "promoter_direction":
                    direction,
            }
        )

    summary = (
        links
        .groupby(
            "gene_name",
            group_keys=False
        )
        .apply(
            classify,
            include_groups=False
        )
        .reset_index()
    )

    summary.to_csv(
        METH_GENE_OUT,
        sep="\t",
        index=False
    )

    return summary


# ============================================================
# RNA + METHYLATION INTEGRATION
# ============================================================

def integrate(
    rna,
    methylation
):

    integrated = methylation.merge(
        rna[
            [
                "gene_name",
                "rna_log2fc",
                "rna_fdr",
                "rna_significant",
                "rna_direction",
            ]
        ],
        on="gene_name",
        how="inner"
    )

    print()
    print("MULTI-OMICS OVERLAP")
    print("=" * 35)

    print(
        f"Genes with both promoter methylation "
        f"and RNA data: {len(integrated):,}"
    )

    # --------------------------------------------------------
    # Broad biological classification
    # --------------------------------------------------------

    integrated[
        "integration_class"
    ] = "Other"

    # Promoter hypermethylation + RNA down
    silencing = (
        (
            integrated[
                "promoter_mean_delta_beta"
            ] > 0
        )
        &
        (
            integrated[
                "rna_direction"
            ] == "Downregulated"
        )
    )

    # Promoter hypomethylation + RNA up
    activation = (
        (
            integrated[
                "promoter_mean_delta_beta"
            ] < 0
        )
        &
        (
            integrated[
                "rna_direction"
            ] == "Upregulated"
        )
    )

    integrated.loc[
        silencing,
        "integration_class"
    ] = (
        "Candidate epigenetic silencing"
    )

    integrated.loc[
        activation,
        "integration_class"
    ] = (
        "Candidate epigenetic activation"
    )

    # --------------------------------------------------------
    # Strict candidates:
    # promoter CpGs all move in a single direction
    # --------------------------------------------------------

    strict_silencing = (
        (
            integrated[
                "promoter_direction"
            ] == "Hypermethylated"
        )
        &
        (
            integrated[
                "rna_direction"
            ] == "Downregulated"
        )
    )

    strict_activation = (
        (
            integrated[
                "promoter_direction"
            ] == "Hypomethylated"
        )
        &
        (
            integrated[
                "rna_direction"
            ] == "Upregulated"
        )
    )

    integrated[
        "strict_concordant"
    ] = (
        strict_silencing
        |
        strict_activation
    )

    # --------------------------------------------------------
    # Ranking score
    #
    # Used only to prioritise candidates,
    # not as a statistical test.
    # --------------------------------------------------------

    safe_fdr = (
        integrated["rna_fdr"]
        .clip(lower=1e-300)
    )

    integrated[
        "ranking_score"
    ] = (
        integrated[
            "rna_log2fc"
        ].abs()
        *
        integrated[
            "promoter_mean_delta_beta"
        ].abs()
        *
        (
            -np.log10(
                safe_fdr
            )
        ).clip(
            upper=50
        )
    )

    integrated = integrated.sort_values(
        [
            "strict_concordant",
            "ranking_score",
        ],
        ascending=[
            False,
            False,
        ]
    )

    integrated.to_csv(
        INTEGRATED_OUT,
        sep="\t",
        index=False
    )

    candidates = integrated[
        integrated[
            "integration_class"
        ].isin(
            [
                "Candidate epigenetic silencing",
                "Candidate epigenetic activation",
            ]
        )
    ].copy()

    candidates.to_csv(
        CANDIDATE_OUT,
        sep="\t",
        index=False
    )

    strict = integrated[
        integrated[
            "strict_concordant"
        ]
    ].copy()

    strict.to_csv(
        STRICT_OUT,
        sep="\t",
        index=False
    )

    return (
        integrated,
        candidates,
        strict
    )


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print("=" * 60)
    print(
        "TCGA-KIRC RNA-SEQ × DNA METHYLATION INTEGRATION"
    )
    print("=" * 60)

    if not METH_FILE.exists():

        raise RuntimeError(
            f"Methylation annotation not found:\n"
            f"{METH_FILE}"
        )

    # --------------------------------------------------------
    # RNA
    # --------------------------------------------------------

    rna = load_rna()

    # --------------------------------------------------------
    # Methylation
    # --------------------------------------------------------

    print()
    print("Loading significant annotated DMPs...")

    meth = pd.read_csv(
        METH_FILE,
        sep="\t"
    )

    print(
        f"Significant methylation CpGs: "
        f"{len(meth):,}"
    )

    required = {
        "cpg_id",
        "gene_names",
        "gene_region",
        "direction",
        "delta_beta",
        "fdr",
        "mean_beta_tumour",
        "mean_beta_normal",
    }

    missing = (
        required
        -
        set(
            meth.columns
        )
    )

    if missing:

        raise RuntimeError(
            "Missing methylation columns: "
            +
            ", ".join(
                sorted(
                    missing
                )
            )
        )

    # --------------------------------------------------------
    # Correct gene-region pairing
    # --------------------------------------------------------

    promoter_links = (
        build_promoter_links(
            meth
        )
    )

    methylation_gene = (
        summarise_promoter_methylation(
            promoter_links
        )
    )

    # --------------------------------------------------------
    # Integrate
    # --------------------------------------------------------

    (
        integrated,
        candidates,
        strict,
    ) = integrate(
        rna,
        methylation_gene
    )

    silencing = (
        candidates[
            "integration_class"
        ]
        ==
        "Candidate epigenetic silencing"
    ).sum()

    activation = (
        candidates[
            "integration_class"
        ]
        ==
        "Candidate epigenetic activation"
    ).sum()

    # --------------------------------------------------------
    # FINAL SUMMARY
    # --------------------------------------------------------

    print()
    print("=" * 60)
    print("MULTI-OMICS INTEGRATION SUMMARY")
    print("=" * 60)

    print(
        f"Genes with promoter methylation "
        f"and RNA-seq data: "
        f"{len(integrated):,}"
    )

    print(
        f"Concordant epigenetic candidates: "
        f"{len(candidates):,}"
    )

    print(
        f"  Silencing candidates: "
        f"{silencing:,}"
    )

    print(
        f"  Activation candidates: "
        f"{activation:,}"
    )

    print(
        f"Strict concordant candidates: "
        f"{len(strict):,}"
    )

    print()
    print(
        "TOP 25 STRICT CONCORDANT CANDIDATES"
    )
    print("=" * 60)

    columns = [
        "gene_name",
        "integration_class",
        "promoter_direction",
        "promoter_dmps",
        "promoter_mean_delta_beta",
        "rna_log2fc",
        "rna_fdr",
        "ranking_score",
    ]

    if len(strict):

        print(
            strict[
                columns
            ]
            .head(25)
            .to_string(
                index=False
            )
        )

    else:

        print(
            "No strict candidates found "
            "under current thresholds."
        )

    # --------------------------------------------------------
    # Known ccRCC biology check
    # --------------------------------------------------------

    kidney_genes = {
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

    biology = integrated[
        integrated[
            "gene_name"
        ].isin(
            kidney_genes
        )
    ]

    print()
    print("SELECTED ccRCC / KIDNEY BIOLOGY CHECK")
    print("=" * 60)

    if len(biology):

        print(
            biology[
                [
                    "gene_name",
                    "promoter_direction",
                    "promoter_mean_delta_beta",
                    "rna_direction",
                    "rna_log2fc",
                    "rna_fdr",
                    "integration_class",
                ]
            ]
            .sort_values(
                "gene_name"
            )
            .to_string(
                index=False
            )
        )

    else:
        print(
            "No selected marker genes had "
            "significant promoter DMPs."
        )

    summary = f"""
TCGA-KIRC RNA-seq × DNA methylation integration

Differential expression thresholds:
FDR < {DE_FDR}
absolute log2 fold change >= {DE_LOG2FC}

Methylation thresholds inherited from Step 15:
FDR < 0.05
absolute delta-beta >= 0.20

Promoter regions:
TSS1500
TSS200
5'UTR
1stExon

Genes with both promoter methylation and RNA data:
{len(integrated)}

Concordant epigenetic candidates:
{len(candidates)}

Candidate epigenetic silencing:
{silencing}

Candidate epigenetic activation:
{activation}

Strict concordant candidates:
{len(strict)}

Interpretation:
Promoter hypermethylation combined with RNA downregulation
is considered compatible with candidate epigenetic silencing.

Promoter hypomethylation combined with RNA upregulation
is considered compatible with candidate epigenetic activation.

These relationships are associations and do not by themselves
demonstrate causal regulation.
""".strip()

    SUMMARY_OUT.write_text(
        summary + "\n",
        encoding="utf-8"
    )

    print()
    print("Saved:")
    print(PROMOTER_LINKS_OUT)
    print(METH_GENE_OUT)
    print(INTEGRATED_OUT)
    print(CANDIDATE_OUT)
    print(STRICT_OUT)
    print(SUMMARY_OUT)

    print()
    print("=" * 60)
    print("RNA-METHYLATION INTEGRATION PASSED")
    print("=" * 60)


if __name__ == "__main__":
    main()
