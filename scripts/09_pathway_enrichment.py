"""Pathway enrichment of TCGA-KIRC differentially expressed genes."""

from pathlib import Path

import gseapy as gp
import pandas as pd


DATA_ROOT = Path("/mnt/e/KIRC_data")

DE_FILE = (
    DATA_ROOT
    / "results"
    / "differential_expression"
    / "tumor_vs_normal_deseq2.tsv"
)

OUT = DATA_ROOT / "results" / "pathway_enrichment"
OUT.mkdir(parents=True, exist_ok=True)


def clean_gene_list(df):
    genes = (
        df["gene_name"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    genes = genes[genes != ""]

    return sorted(set(genes))


def run_enrichment(genes, label):

    libraries = [
        "GO_Biological_Process_2023",
        "KEGG_2021_Human",
        "Reactome_2022",
    ]

    print(f"\nRunning enrichment for {label}")
    print(f"Genes: {len(genes)}")

    result = gp.enrichr(
        gene_list=genes,
        gene_sets=libraries,
        organism="human",
        outdir=None,
        cutoff=0.05,
    )

    table = result.results.copy()

    table = table.sort_values(
        "Adjusted P-value"
    )

    table.to_csv(
        OUT / f"{label}_enrichment.tsv",
        sep="\t",
        index=False,
    )

    print("\nTop enriched pathways:")
    print(
        table[
            [
                "Gene_set",
                "Term",
                "Adjusted P-value",
                "Combined Score",
            ]
        ]
        .head(15)
        .to_string(index=False)
    )

    return table


def main():

    de = pd.read_csv(DE_FILE, sep="\t")

    sig = de[
        (de["padj"] < 0.05)
        & (de["log2FoldChange"].abs() >= 1)
        & (de["gene_type"] == "protein_coding")
    ].copy()

    up = sig[
        sig["log2FoldChange"] >= 1
    ]

    down = sig[
        sig["log2FoldChange"] <= -1
    ]

    up_genes = clean_gene_list(up)
    down_genes = clean_gene_list(down)

    print("Protein-coding DE genes")
    print("----------------------")
    print(f"Upregulated:   {len(up_genes)}")
    print(f"Downregulated: {len(down_genes)}")

    run_enrichment(
        up_genes,
        "upregulated",
    )

    run_enrichment(
        down_genes,
        "downregulated",
    )

    print()
    print("PATHWAY ENRICHMENT PASSED")
    print(f"Results: {OUT}")


if __name__ == "__main__":
    main()
