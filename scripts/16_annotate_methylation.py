"""
Annotate TCGA-KIRC differential methylation CpGs using the
Illumina HumanMethylation450 v1.2 manifest.

Input:
    HumanMethylation450_15017482_v1-2.csv
    kirc_differential_methylation.tsv
    kirc_significant_dmps.tsv

Outputs:
    annotated all CpGs
    annotated significant DMPs
    CpG-to-gene links
    gene-level methylation summary
    annotation provenance
"""

from pathlib import Path
import pandas as pd
import numpy as np


# ============================================================
# PATHS
# ============================================================

MANIFEST = Path(
    "/mnt/e/KIRC_data/reference/"
    "HumanMethylation450_15017482_v1-2.csv"
)

DM_DIR = Path(
    "/mnt/e/KIRC_data/results/differential_methylation"
)

ALL_DMPS = (
    DM_DIR /
    "kirc_differential_methylation.tsv"
)

SIG_DMPS = (
    DM_DIR /
    "kirc_significant_dmps.tsv"
)

OUT_DIR = (
    DM_DIR /
    "annotation"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

ALL_ANNOTATED = (
    OUT_DIR /
    "kirc_differential_methylation_annotated.tsv"
)

SIG_ANNOTATED = (
    OUT_DIR /
    "kirc_significant_dmps_annotated.tsv"
)

GENE_LINKS_OUT = (
    OUT_DIR /
    "kirc_dmp_gene_links.tsv"
)

GENE_SUMMARY_OUT = (
    OUT_DIR /
    "kirc_gene_methylation_summary.tsv"
)

PROVENANCE_OUT = (
    OUT_DIR /
    "annotation_provenance.txt"
)


# ============================================================
# LOAD ILLUMINA ANNOTATION
# ============================================================

def load_annotation():

    print()
    print("Reading Illumina HumanMethylation450 v1.2 manifest...")
    print(f"Manifest: {MANIFEST}")

    if not MANIFEST.exists():
        raise RuntimeError(
            f"Manifest not found:\n{MANIFEST}"
        )

    print(
        f"Manifest size: "
        f"{MANIFEST.stat().st_size / (1024**2):.1f} MB"
    )

    # --------------------------------------------------------
    # Official 450K CSV has a 7-line metadata preamble.
    # The actual column header begins after those 7 lines.
    # --------------------------------------------------------

    header = pd.read_csv(
        MANIFEST,
        skiprows=7,
        nrows=0,
        encoding="utf-8-sig"
    )

    header.columns = [
        str(c).strip()
        for c in header.columns
    ]

    print()
    print(
        f"Columns detected: {len(header.columns)}"
    )

    # Prefer IlmnID, which is the standard CpG identifier.
    if "IlmnID" in header.columns:
        id_col = "IlmnID"

    elif "Name" in header.columns:
        id_col = "Name"

    else:

        print("\nFirst detected columns:")
        print(
            "\n".join(
                header.columns[:20]
            )
        )

        raise RuntimeError(
            "Could not identify IlmnID/Name column."
        )

    wanted = [
        id_col,
        "CHR",
        "MAPINFO",
        "Strand",
        "UCSC_RefGene_Name",
        "UCSC_RefGene_Accession",
        "UCSC_RefGene_Group",
        "UCSC_CpG_Islands_Name",
        "Relation_to_UCSC_CpG_Island",
        "Regulatory_Feature_Name",
        "Regulatory_Feature_Group",
        "DHS",
        "Enhancer",
    ]

    usecols = [
        c
        for c in wanted
        if c in header.columns
    ]

    print()
    print("Annotation columns selected:")

    for column in usecols:
        print(f"  {column}")

    annotation = pd.read_csv(
        MANIFEST,
        skiprows=7,
        usecols=usecols,
        dtype=str,
        encoding="utf-8-sig",
        low_memory=False,
    )

    annotation.columns = [
        str(c).strip()
        for c in annotation.columns
    ]

    # --------------------------------------------------------
    # Standardise column names
    # --------------------------------------------------------

    annotation = annotation.rename(
        columns={
            id_col:
                "cpg_id",

            "CHR":
                "chromosome",

            "MAPINFO":
                "position_hg19",

            "Strand":
                "strand",

            "UCSC_RefGene_Name":
                "gene_names",

            "UCSC_RefGene_Accession":
                "gene_accessions",

            "UCSC_RefGene_Group":
                "gene_region",

            "UCSC_CpG_Islands_Name":
                "cpg_island",

            "Relation_to_UCSC_CpG_Island":
                "cpg_island_relation",

            "Regulatory_Feature_Name":
                "regulatory_feature",

            "Regulatory_Feature_Group":
                "regulatory_feature_group",

            "DHS":
                "dhs",

            "Enhancer":
                "enhancer",
        }
    )

    annotation["cpg_id"] = (
        annotation["cpg_id"]
        .astype(str)
        .str.strip()
    )

    # Keep actual methylation CpG probes
    annotation = annotation[
        annotation["cpg_id"]
        .str.startswith("cg")
    ].copy()

    # Remove duplicated CpG identifiers
    annotation = (
        annotation
        .drop_duplicates(
            subset="cpg_id"
        )
        .reset_index(drop=True)
    )

    print()
    print(
        f"450K CpG annotation records: "
        f"{len(annotation):,}"
    )

    print()
    print("Example annotation:")
    print(
        annotation
        .head(5)
        .to_string(index=False)
    )

    return annotation


# ============================================================
# ANNOTATE RESULTS
# ============================================================

def annotate_results(results, annotation, label):

    merged = results.merge(
        annotation,
        on="cpg_id",
        how="left",
        validate="many_to_one"
    )

    if "chromosome" in merged.columns:

        matched = (
            merged["chromosome"]
            .notna()
            .sum()
        )

    else:
        matched = 0

    percentage = (
        matched /
        len(merged) *
        100
        if len(merged)
        else 0
    )

    print()
    print(
        f"{label}: "
        f"{matched:,}/{len(merged):,} "
        f"CpGs annotated "
        f"({percentage:.2f}%)"
    )

    return merged


# ============================================================
# CREATE CpG -> GENE LINKS
# ============================================================

def build_gene_links(sig):

    if "gene_names" not in sig.columns:

        raise RuntimeError(
            "gene_names annotation unavailable."
        )

    links = sig[
        sig["gene_names"].notna()
    ].copy()

    # --------------------------------------------------------
    # One row can refer to multiple genes separated by ;
    # --------------------------------------------------------

    links["gene_name"] = (
        links["gene_names"]
        .astype(str)
        .str.split(";")
    )

    links = links.explode(
        "gene_name"
    )

    links["gene_name"] = (
        links["gene_name"]
        .astype(str)
        .str.strip()
    )

    links = links[
        (links["gene_name"] != "")
        &
        (links["gene_name"].str.lower() != "nan")
    ].copy()

    links = links.drop_duplicates(
        subset=[
            "cpg_id",
            "gene_name",
        ]
    )

    # --------------------------------------------------------
    # Promoter classification
    #
    # TSS1500
    # TSS200
    # 5'UTR
    # 1stExon
    # --------------------------------------------------------

    if "gene_region" in links.columns:

        promoter_pattern = (
            r"TSS1500|TSS200|5'UTR|1stExon"
        )

        links[
            "promoter_associated"
        ] = (
            links["gene_region"]
            .fillna("")
            .str.contains(
                promoter_pattern,
                case=False,
                regex=True
            )
        )

    else:

        links[
            "promoter_associated"
        ] = False

    links[
        "abs_delta_beta"
    ] = (
        links["delta_beta"]
        .abs()
    )

    return links


# ============================================================
# GENE SUMMARY
# ============================================================

def build_gene_summary(links):

    summary = (
        links
        .groupby("gene_name")
        .agg(

            n_dmps=(
                "cpg_id",
                "nunique"
            ),

            hypermethylated_dmps=(
                "direction",
                lambda x:
                    int(
                        (
                            x ==
                            "Hypermethylated"
                        ).sum()
                    )
            ),

            hypomethylated_dmps=(
                "direction",
                lambda x:
                    int(
                        (
                            x ==
                            "Hypomethylated"
                        ).sum()
                    )
            ),

            promoter_dmps=(
                "promoter_associated",
                "sum"
            ),

            mean_delta_beta=(
                "delta_beta",
                "mean"
            ),

            strongest_delta_beta=(
                "abs_delta_beta",
                "max"
            ),

            minimum_fdr=(
                "fdr",
                "min"
            ),
        )
        .reset_index()
    )

    summary = summary.sort_values(
        [
            "n_dmps",
            "strongest_delta_beta"
        ],
        ascending=[
            False,
            False
        ]
    )

    return summary


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print("=" * 55)
    print("TCGA-KIRC CpG -> GENE / GENOMIC ANNOTATION")
    print("=" * 55)

    # --------------------------------------------------------
    # Check files
    # --------------------------------------------------------

    for path in [
        MANIFEST,
        ALL_DMPS,
        SIG_DMPS
    ]:

        if not path.exists():

            raise RuntimeError(
                f"Required file missing:\n{path}"
            )

    # --------------------------------------------------------
    # Load annotation
    # --------------------------------------------------------

    annotation = load_annotation()

    # --------------------------------------------------------
    # Load differential methylation
    # --------------------------------------------------------

    print()
    print("Loading differential methylation results...")

    all_dm = pd.read_csv(
        ALL_DMPS,
        sep="\t"
    )

    significant = pd.read_csv(
        SIG_DMPS,
        sep="\t"
    )

    print(
        f"All tested CpGs: "
        f"{len(all_dm):,}"
    )

    print(
        f"Significant DMPs: "
        f"{len(significant):,}"
    )

    # --------------------------------------------------------
    # Annotate ALL tested probes
    # --------------------------------------------------------

    all_annotated = annotate_results(
        all_dm,
        annotation,
        "All tested CpGs"
    )

    all_annotated.to_csv(
        ALL_ANNOTATED,
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Annotate SIGNIFICANT probes
    # --------------------------------------------------------

    sig_annotated = annotate_results(
        significant,
        annotation,
        "Significant DMPs"
    )

    sig_annotated.to_csv(
        SIG_ANNOTATED,
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # CpG -> Gene
    # --------------------------------------------------------

    print()
    print("Building CpG-to-gene links...")

    gene_links = build_gene_links(
        sig_annotated
    )

    gene_links.to_csv(
        GENE_LINKS_OUT,
        sep="\t",
        index=False
    )

    n_genes = (
        gene_links[
            "gene_name"
        ]
        .nunique()
    )

    n_promoter = int(
        gene_links[
            "promoter_associated"
        ].sum()
    )

    print(
        f"CpG-gene associations: "
        f"{len(gene_links):,}"
    )

    print(
        f"Unique genes: "
        f"{n_genes:,}"
    )

    print(
        f"Promoter-associated DMP links: "
        f"{n_promoter:,}"
    )

    # --------------------------------------------------------
    # Gene summary
    # --------------------------------------------------------

    gene_summary = build_gene_summary(
        gene_links
    )

    gene_summary.to_csv(
        GENE_SUMMARY_OUT,
        sep="\t",
        index=False
    )

    print()
    print("=" * 55)
    print(
        "TOP 20 GENES BY SIGNIFICANT "
        "METHYLATION SITES"
    )
    print("=" * 55)

    print(
        gene_summary
        .head(20)
        .to_string(index=False)
    )

    # --------------------------------------------------------
    # CpG island context
    # --------------------------------------------------------

    if (
        "cpg_island_relation"
        in sig_annotated.columns
    ):

        print()
        print("=" * 35)
        print("CpG ISLAND CONTEXT")
        print("=" * 35)

        print(
            sig_annotated[
                "cpg_island_relation"
            ]
            .fillna("Unannotated")
            .value_counts()
            .to_string()
        )

    # --------------------------------------------------------
    # Gene-region context
    # --------------------------------------------------------

    if "gene_region" in sig_annotated.columns:

        print()
        print("=" * 35)
        print("GENE REGION CONTEXT")
        print("=" * 35)

        print(
            sig_annotated[
                "gene_region"
            ]
            .fillna("Unannotated")
            .value_counts()
            .head(20)
            .to_string()
        )

    # --------------------------------------------------------
    # Chromosome distribution
    # --------------------------------------------------------

    if "chromosome" in sig_annotated.columns:

        print()
        print("=" * 35)
        print("CHROMOSOME DISTRIBUTION")
        print("=" * 35)

        print(
            sig_annotated[
                "chromosome"
            ]
            .fillna("Unannotated")
            .value_counts()
            .sort_index()
            .to_string()
        )

    # --------------------------------------------------------
    # Hyper / Hypo summary
    # --------------------------------------------------------

    print()
    print("=" * 35)
    print("DMP DIRECTION")
    print("=" * 35)

    print(
        significant[
            "direction"
        ]
        .value_counts()
        .to_string()
    )

    # --------------------------------------------------------
    # Provenance
    # --------------------------------------------------------

    provenance = f"""
TCGA-KIRC methylation CpG annotation

Array:
Illumina HumanMethylation450 BeadChip

Manifest:
HumanMethylation450 v1.2

Manifest file:
{MANIFEST}

Reference coordinate system:
Illumina HumanMethylation450 v1.2 / hg19 (GRCh37)

Differential methylation input:
{SIG_DMPS}

Significant DMPs:
{len(significant)}

CpG-to-gene associations:
{len(gene_links)}

Unique genes:
{n_genes}

Promoter-associated CpG-gene links:
{n_promoter}

Significance thresholds inherited from Step 15:
FDR < 0.05
absolute delta-beta >= 0.20
""".strip()

    PROVENANCE_OUT.write_text(
        provenance + "\n",
        encoding="utf-8"
    )

    # --------------------------------------------------------
    # FINAL
    # --------------------------------------------------------

    print()
    print("Saved:")
    print(ALL_ANNOTATED)
    print(SIG_ANNOTATED)
    print(GENE_LINKS_OUT)
    print(GENE_SUMMARY_OUT)
    print(PROVENANCE_OUT)

    print()
    print("=" * 55)
    print("CpG ANNOTATION PASSED")
    print("=" * 55)


if __name__ == "__main__":
    main()
