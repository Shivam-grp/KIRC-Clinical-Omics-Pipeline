"""Compare TCGA-KIRC methylation platforms and overlap with RNA-seq."""

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

OUT = (
    DATA_ROOT
    / "downloads"
    / "tcga_kirc_methylation_platform_comparison.tsv"
)


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

    # Keep only tumour and normal samples.
    meth = meth[
        meth["sample_type"].isin(
            [
                "Primary Tumor",
                "Solid Tissue Normal",
            ]
        )
    ].copy()

    print("Methylation files:", len(meth))

    print("\nRNA-seq biological samples:")
    print(len(rna))

    print("\nRNA-seq sample types:")
    print(
        rna["sample_type"]
        .value_counts()
        .to_string()
    )

    rna_samples = set(
        rna.index.dropna()
    )

    rna_patients = set(
        rna["patient_id"].dropna()
    )

    results = []

    print()
    print("METHYLATION PLATFORM COMPARISON")
    print("=" * 72)

    for platform, group in meth.groupby("platform"):

        unique_samples = (
            group["sample_submitter_id"]
            .dropna()
            .nunique()
        )

        unique_patients = (
            group["case_submitter_id"]
            .dropna()
            .nunique()
        )

        duplicated_samples = (
            group["sample_submitter_id"]
            .duplicated()
            .sum()
        )

        meth_samples = set(
            group["sample_submitter_id"]
            .dropna()
        )

        meth_patients = set(
            group["case_submitter_id"]
            .dropna()
        )

        sample_overlap = (
            meth_samples & rna_samples
        )

        patient_overlap = (
            meth_patients & rna_patients
        )

        total_gib = (
            group["file_size"].sum()
            / (1024 ** 3)
        )

        tumor_files = (
            group["sample_type"]
            == "Primary Tumor"
        ).sum()

        normal_files = (
            group["sample_type"]
            == "Solid Tissue Normal"
        ).sum()

        print()
        print(platform)
        print("-" * len(str(platform)))

        print(f"Files: {len(group)}")
        print(
            f"Primary Tumor files: "
            f"{tumor_files}"
        )
        print(
            f"Solid Tissue Normal files: "
            f"{normal_files}"
        )

        print(
            f"Unique methylation samples: "
            f"{unique_samples}"
        )

        print(
            f"Unique methylation patients: "
            f"{unique_patients}"
        )

        print(
            f"Duplicated sample IDs: "
            f"{duplicated_samples}"
        )

        print(
            f"Exact RNA-seq sample overlap: "
            f"{len(sample_overlap)}"
        )

        print(
            f"RNA-seq patient overlap: "
            f"{len(patient_overlap)}"
        )

        print(
            f"Estimated download: "
            f"{total_gib:.2f} GiB"
        )

        results.append(
            {
                "platform": platform,
                "files": len(group),
                "tumor_files": tumor_files,
                "normal_files": normal_files,
                "unique_samples": unique_samples,
                "unique_patients": unique_patients,
                "duplicated_samples": duplicated_samples,
                "rna_sample_overlap": len(sample_overlap),
                "rna_patient_overlap": len(patient_overlap),
                "download_gib": total_gib,
            }
        )

    result = pd.DataFrame(results)

    result.to_csv(
        OUT,
        sep="\t",
        index=False,
    )

    print()
    print("Saved:")
    print(OUT)

    print()
    print("PLATFORM COMPARISON PASSED")


if __name__ == "__main__":
    main()
