"""Prepare TCGA-KIRC 450K methylation cohort overlapping RNA-seq samples."""

from pathlib import Path

import pandas as pd


DATA_ROOT = Path("/mnt/e/KIRC_data")

METH_FILE = (
    DATA_ROOT
    / "downloads"
    / "tcga_kirc_methylation_manifest.tsv"
)

RNA_META = (
    DATA_ROOT
    / "processed"
    / "kirc_rnaseq_sample_metadata.tsv"
)

OUTPUT = (
    DATA_ROOT
    / "downloads"
    / "tcga_kirc_methylation_450_rna_overlap.tsv"
)


PLATFORM = "Illumina Human Methylation 450"


def main():

    meth = pd.read_csv(
        METH_FILE,
        sep="\t",
    )

    rna = pd.read_csv(
        RNA_META,
        sep="\t",
        index_col=0,
    )

    # RNA biological sample IDs are stored in the index.
    rna_samples = set(rna.index.astype(str))

    # Select 450K platform and tumour/normal only.
    cohort = meth[
        (meth["platform"] == PLATFORM)
        & (
            meth["sample_type"].isin(
                [
                    "Primary Tumor",
                    "Solid Tissue Normal",
                ]
            )
        )
    ].copy()

    # Retain only samples with exact RNA-seq overlap.
    cohort = cohort[
        cohort["sample_submitter_id"].isin(
            rna_samples
        )
    ].copy()

    print("450K RNA-overlap cohort")
    print("=======================")

    print(f"\nFiles: {len(cohort)}")

    print("\nSample types:")
    print(
        cohort["sample_type"]
        .value_counts()
        .to_string()
    )

    print(
        "\nUnique biological samples:",
        cohort["sample_submitter_id"].nunique(),
    )

    print(
        "Unique patients:",
        cohort["case_submitter_id"].nunique(),
    )

    duplicated = cohort[
        cohort.duplicated(
            "sample_submitter_id",
            keep=False,
        )
    ].sort_values("sample_submitter_id")

    print(
        "Duplicated sample IDs:",
        cohort["sample_submitter_id"]
        .duplicated()
        .sum(),
    )

    total_gib = (
        cohort["file_size"].sum()
        / (1024 ** 3)
    )

    print(
        f"\nDownload size for RNA-overlap cohort: "
        f"{total_gib:.2f} GiB"
    )

    if not duplicated.empty:

        print("\nDUPLICATED METHYLATION SAMPLES")
        print("==============================")

        cols = [
            "sample_submitter_id",
            "case_submitter_id",
            "sample_type",
            "file_id",
            "file_name",
            "file_size",
        ]

        print(
            duplicated[cols]
            .to_string(index=False)
        )

    cohort = cohort.sort_values(
        [
            "sample_type",
            "case_submitter_id",
            "sample_submitter_id",
        ]
    )

    cohort.to_csv(
        OUTPUT,
        sep="\t",
        index=False,
    )

    print("\nSaved:")
    print(OUTPUT)

    print("\nMETHYLATION COHORT PREPARATION PASSED")


if __name__ == "__main__":
    main()
