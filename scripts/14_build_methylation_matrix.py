from pathlib import Path

import numpy as np
import pandas as pd


MANIFEST = Path(
    "/mnt/e/KIRC_data/downloads/tcga_kirc_methylation_450_rna_overlap.tsv"
)

RAW_DIR = Path("/mnt/e/KIRC_data/raw/methylation_450")

OUT_DIR = Path("/mnt/e/KIRC_data/processed/methylation_450")
OUT_DIR.mkdir(parents=True, exist_ok=True)

MATRIX_OUT = OUT_DIR / "kirc_methylation_beta_matrix.npy"
PROBE_OUT = OUT_DIR / "kirc_methylation_probe_ids.tsv"
META_OUT = OUT_DIR / "kirc_methylation_sample_metadata.tsv"
QC_OUT = OUT_DIR / "kirc_methylation_qc.tsv"


def read_beta_file(path):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["cpg_id", "beta"],
        dtype={
            "cpg_id": "string",
            "beta": "float32",
        },
        na_values=["NA", "NaN", ""],
    )


def main():

    print("\nTCGA-KIRC 450K METHYLATION MATRIX BUILD")
    print("=" * 45)

    manifest = pd.read_csv(MANIFEST, sep="\t")

    required = {
        "file_id",
        "file_name",
        "sample_submitter_id",
        "sample_type",
    }

    missing = required - set(manifest.columns)

    if missing:
        raise SystemExit(
            "Missing manifest columns: "
            + ", ".join(sorted(missing))
        )

    print(f"Manifest files: {len(manifest)}")

    # Biological samples
    samples = (
        manifest[
            ["sample_submitter_id", "sample_type"]
        ]
        .drop_duplicates("sample_submitter_id")
        .reset_index(drop=True)
    )

    print(
        f"Unique biological samples: {len(samples)}"
    )

    duplicates = (
        manifest["sample_submitter_id"]
        .value_counts()
    )

    duplicates = duplicates[duplicates > 1]

    print(
        f"Samples with technical duplicates: {len(duplicates)}"
    )

    # --------------------------------------------------
    # Establish canonical probe order
    # --------------------------------------------------

    first = manifest.iloc[0]

    first_path = (
        RAW_DIR
        / str(first["file_id"])
        / str(first["file_name"])
    )

    if not first_path.exists():
        raise SystemExit(
            f"File not found:\n{first_path}"
        )

    first_df = read_beta_file(first_path)

    probes = first_df["cpg_id"].astype(str).to_numpy()
    n_probes = len(probes)
    n_samples = len(samples)

    print(f"CpG probes: {n_probes:,}")
    print(
        f"Matrix: {n_probes:,} probes x "
        f"{n_samples:,} samples"
    )

    if len(set(probes)) != n_probes:
        raise RuntimeError(
            "Duplicate CpG IDs detected."
        )

    pd.DataFrame(
        {"cpg_id": probes}
    ).to_csv(
        PROBE_OUT,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # Disk-backed matrix — suitable for 8 GB RAM
    # --------------------------------------------------

    matrix = np.lib.format.open_memmap(
        MATRIX_OUT,
        mode="w+",
        dtype=np.float32,
        shape=(n_probes, n_samples),
    )

    qc_records = []

    # --------------------------------------------------
    # Process one biological sample at a time
    # --------------------------------------------------

    for sample_index, sample_row in samples.iterrows():

        sample_id = sample_row["sample_submitter_id"]
        sample_type = sample_row["sample_type"]

        rows = manifest[
            manifest["sample_submitter_id"] == sample_id
        ]

        vectors = []

        for _, row in rows.iterrows():

            path = (
                RAW_DIR
                / str(row["file_id"])
                / str(row["file_name"])
            )

            if not path.exists():
                raise RuntimeError(
                    f"Missing file:\n{path}"
                )

            df = read_beta_file(path)

            if len(df) != n_probes:
                raise RuntimeError(
                    f"Probe count mismatch:\n{path}"
                )

            current_probes = (
                df["cpg_id"]
                .astype(str)
                .to_numpy()
            )

            if not np.array_equal(
                current_probes,
                probes,
            ):
                raise RuntimeError(
                    "CpG ordering mismatch in:\n"
                    f"{path}"
                )

            beta = df["beta"].to_numpy(
                dtype=np.float32
            )

            vectors.append(beta)

        # Collapse technical duplicates by mean
        if len(vectors) == 1:
            beta = vectors[0]
        else:
            stacked = np.vstack(vectors)

            with np.errstate(invalid="ignore"):
                beta = np.nanmean(
                    stacked,
                    axis=0,
                ).astype(np.float32)

        valid = beta[~np.isnan(beta)]

        outside_range = int(
            ((valid < 0) | (valid > 1)).sum()
        )

        if outside_range:
            raise RuntimeError(
                f"Invalid beta values in {sample_id}"
            )

        matrix[:, sample_index] = beta

        missing_count = int(
            np.isnan(beta).sum()
        )

        qc_records.append(
            {
                "sample_id": sample_id,
                "sample_type": sample_type,
                "files_collapsed": len(vectors),
                "n_probes": n_probes,
                "missing_beta": missing_count,
                "missing_percent":
                    100 * missing_count / n_probes,
                "mean_beta": float(np.nanmean(beta)),
                "median_beta": float(np.nanmedian(beta)),
                "std_beta": float(np.nanstd(beta)),
                "min_beta": float(np.nanmin(beta)),
                "max_beta": float(np.nanmax(beta)),
            }
        )

        if (
            (sample_index + 1) % 25 == 0
            or sample_index + 1 == n_samples
        ):
            print(
                f"Processed "
                f"{sample_index + 1}/{n_samples} samples"
            )

    matrix.flush()

    # --------------------------------------------------
    # Save metadata and QC
    # --------------------------------------------------

    samples.rename(
        columns={
            "sample_submitter_id": "sample_id"
        }
    ).to_csv(
        META_OUT,
        sep="\t",
        index=False,
    )

    qc = pd.DataFrame(qc_records)

    qc.to_csv(
        QC_OUT,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # Summary
    # --------------------------------------------------

    print("\nMETHYLATION QC SUMMARY")
    print("=" * 30)

    print(f"Samples: {n_samples}")
    print(f"CpG probes: {n_probes:,}")

    print("\nSample types:")
    print(
        samples["sample_type"]
        .value_counts()
        .to_string()
    )

    print(
        "\nTechnical duplicate samples collapsed:",
        int((qc["files_collapsed"] > 1).sum()),
    )

    print(
        "Median missingness:",
        f"{qc['missing_percent'].median():.4f}%",
    )

    print(
        "Maximum missingness:",
        f"{qc['missing_percent'].max():.4f}%",
    )

    print("\nSaved:")
    print(MATRIX_OUT)
    print(PROBE_OUT)
    print(META_OUT)
    print(QC_OUT)

    print(
        "\nMETHYLATION MATRIX BUILD PASSED"
    )


if __name__ == "__main__":
    main()
