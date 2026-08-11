"""Build a TCGA-KIRC raw RNA-seq count matrix from GDC STAR-count files."""

from pathlib import Path

import pandas as pd


DATA_ROOT = Path("/mnt/e/KIRC_data")

MANIFEST = (
    DATA_ROOT
    / "downloads"
    / "tcga_kirc_rnaseq_cohort.tsv"
)

RAW_DIR = DATA_ROOT / "raw" / "rnaseq"

OUT_DIR = DATA_ROOT / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COUNTS_OUT = OUT_DIR / "kirc_rnaseq_raw_counts.parquet"
META_OUT = OUT_DIR / "kirc_rnaseq_sample_metadata.tsv"


def read_counts(file_path: Path, sample_id: str) -> pd.Series:
    df = pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        usecols=["gene_id", "unstranded"],
    )

    # Remove STAR summary rows such as N_unmapped.
    df = df[df["gene_id"].astype(str).str.startswith("ENSG")]

    counts = pd.to_numeric(
        df["unstranded"],
        errors="raise",
    ).astype("int64")

    counts.index = df["gene_id"].astype(str)

    counts.name = sample_id

    return counts


def main() -> None:
    manifest = pd.read_csv(MANIFEST, sep="\t")

    print(f"RNA-seq files in manifest: {len(manifest)}")

    sample_series = []

    for number, row in enumerate(
        manifest.itertuples(index=False),
        start=1,
    ):
        file_path = (
            RAW_DIR
            / row.file_id
            / row.file_name
        )

        if not file_path.exists():
            raise FileNotFoundError(file_path)

        counts = read_counts(
            file_path,
            row.sample_submitter_id,
        )

        sample_series.append(counts)

        if number % 50 == 0 or number == len(manifest):
            print(
                f"Processed {number}/{len(manifest)} files"
            )

    print("\nCombining count vectors...")

    # genes x sequencing files
    matrix = pd.concat(sample_series, axis=1)

    print(f"Before technical-replicate collapse: {matrix.shape}")

    # Duplicate column names correspond to technical aliquots
    # from the same biological TCGA sample.
    matrix = matrix.T.groupby(level=0).sum().T

    print(f"After technical-replicate collapse: {matrix.shape}")

    # PyDESeq2 expects samples x genes.
    matrix = matrix.T

    # Ensure integer raw counts.
    matrix = matrix.astype("int64")

    print(f"\nFinal matrix: {matrix.shape[0]} samples x {matrix.shape[1]} genes")

    matrix.to_parquet(
        COUNTS_OUT,
        compression="zstd",
    )

    # One metadata record per biological sample.
    metadata = (
        manifest[
            [
                "sample_submitter_id",
                "case_submitter_id",
                "sample_type",
            ]
        ]
        .drop_duplicates(subset="sample_submitter_id")
        .rename(
            columns={
                "sample_submitter_id": "sample_id",
                "case_submitter_id": "patient_id",
            }
        )
        .set_index("sample_id")
    )

    # Match metadata order to count matrix.
    metadata = metadata.loc[matrix.index]

    metadata.to_csv(
        META_OUT,
        sep="\t",
    )

    print("\nSample types:")
    print(metadata["sample_type"].value_counts().to_string())

    print("\nSaved:")
    print(COUNTS_OUT)
    print(META_OUT)

    print("\nCOUNT MATRIX BUILD PASSED")


if __name__ == "__main__":
    main()
