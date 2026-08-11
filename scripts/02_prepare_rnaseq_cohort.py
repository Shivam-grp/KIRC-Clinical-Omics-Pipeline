"""Prepare the TCGA-KIRC tumour/normal RNA-seq analysis cohort."""

from pathlib import Path

import pandas as pd


DATA_ROOT = Path("/mnt/e/KIRC_data")
INPUT = DATA_ROOT / "downloads" / "tcga_kirc_rnaseq_manifest.tsv"
OUTPUT = DATA_ROOT / "downloads" / "tcga_kirc_rnaseq_cohort.tsv"


def main() -> None:
    df = pd.read_csv(INPUT, sep="\t")

    print("Original files:", len(df))

    keep_types = [
        "Primary Tumor",
        "Solid Tissue Normal",
    ]

    cohort = df[df["sample_type"].isin(keep_types)].copy()

    print("\nSelected sample types:")
    print(cohort["sample_type"].value_counts())

    print("\nUnique patients:")
    print(cohort["case_submitter_id"].nunique())

    print("\nUnique samples:")
    print(cohort["sample_submitter_id"].nunique())

    duplicated_samples = cohort["sample_submitter_id"].duplicated().sum()
    print("\nDuplicated sample IDs:", duplicated_samples)

    total_bytes = cohort["file_size"].sum()
    total_gb = total_bytes / (1024 ** 3)

    print(f"\nEstimated download size: {total_gb:.2f} GiB")

    cohort = cohort.sort_values(
        ["sample_type", "case_submitter_id"]
    ).reset_index(drop=True)

    cohort.to_csv(OUTPUT, sep="\t", index=False)

    print(f"\nCohort manifest saved to:")
    print(OUTPUT)

    print("\nCOHORT PREPARATION PASSED")


if __name__ == "__main__":
    main()
