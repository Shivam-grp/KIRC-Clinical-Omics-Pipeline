from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy import stats


BASE = Path("/mnt/e/KIRC_data/processed/methylation_450")
OUT = Path("/mnt/e/KIRC_data/results/differential_methylation")
OUT.mkdir(parents=True, exist_ok=True)

MATRIX_FILE = BASE / "kirc_methylation_beta_matrix.npy"
PROBE_FILE = BASE / "kirc_methylation_probe_ids.tsv"
META_FILE = BASE / "kirc_methylation_sample_metadata.tsv"
QC_FILE = BASE / "kirc_methylation_qc.tsv"

SAMPLE_FILTER_OUT = OUT / "methylation_sample_filter.tsv"
PROBE_FILTER_OUT = OUT / "methylation_probe_filter.tsv"
RESULTS_OUT = OUT / "kirc_differential_methylation.tsv"
SIGNIFICANT_OUT = OUT / "kirc_significant_dmps.tsv"


SAMPLE_MISSING_MAX = 20.0
PROBE_MISSING_MAX = 10.0

FDR_THRESHOLD = 0.05
DELTA_BETA_THRESHOLD = 0.20

BLOCK_SIZE = 20000
EPS = 1e-5


def bh_adjust(pvalues):
    """
    Benjamini-Hochberg FDR correction.
    """
    pvalues = np.asarray(pvalues, dtype=float)

    adjusted = np.full_like(pvalues, np.nan)

    valid = np.isfinite(pvalues)

    p = pvalues[valid]

    if len(p) == 0:
        return adjusted

    order = np.argsort(p)
    ranked = p[order]

    n = len(ranked)

    q = ranked * n / np.arange(1, n + 1)

    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    restored = np.empty_like(q)
    restored[order] = q

    adjusted[valid] = restored

    return adjusted


def beta_to_m(beta):
    """
    Convert beta values to M-values.

    Statistical testing is performed on M-values,
    while delta-beta is retained as biological
    effect size.
    """
    beta = np.clip(beta, EPS, 1 - EPS)

    return np.log2(
        beta / (1 - beta)
    )


def main():

    print()
    print("TCGA-KIRC DIFFERENTIAL METHYLATION")
    print("=" * 42)

    # --------------------------------------------------
    # Load data
    # --------------------------------------------------

    matrix = np.load(
        MATRIX_FILE,
        mmap_mode="r"
    )

    probes = pd.read_csv(
        PROBE_FILE,
        sep="\t"
    )

    metadata = pd.read_csv(
        META_FILE,
        sep="\t"
    )

    qc = pd.read_csv(
        QC_FILE,
        sep="\t"
    )

    print(f"Original matrix: {matrix.shape}")
    print(f"Metadata samples: {len(metadata)}")
    print(f"QC samples: {len(qc)}")

    if matrix.shape[1] != len(metadata):
        raise RuntimeError(
            "Matrix columns do not match metadata rows."
        )

    # --------------------------------------------------
    # SAMPLE QC
    # --------------------------------------------------

    sample_table = metadata.merge(
        qc[
            [
                "sample_id",
                "missing_percent",
            ]
        ],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )

    sample_table["qc_pass"] = (
        sample_table["missing_percent"]
        <= SAMPLE_MISSING_MAX
    )

    sample_table.to_csv(
        SAMPLE_FILTER_OUT,
        sep="\t",
        index=False,
    )

    keep_samples = (
        sample_table["qc_pass"]
        .to_numpy()
    )

    keep_indices = np.where(
        keep_samples
    )[0]

    kept_meta = (
        sample_table.loc[
            keep_samples
        ]
        .reset_index(drop=True)
    )

    print()
    print("SAMPLE QC")
    print("=" * 20)

    print(
        f"Threshold: <= "
        f"{SAMPLE_MISSING_MAX:.0f}% missing"
    )

    print(
        f"Samples retained: "
        f"{len(keep_indices)}/{len(sample_table)}"
    )

    print()
    print(
        kept_meta["sample_type"]
        .value_counts()
        .to_string()
    )

    tumour_mask = (
        kept_meta["sample_type"]
        == "Primary Tumor"
    ).to_numpy()

    normal_mask = (
        kept_meta["sample_type"]
        == "Solid Tissue Normal"
    ).to_numpy()

    n_tumour = int(
        tumour_mask.sum()
    )

    n_normal = int(
        normal_mask.sum()
    )

    print()
    print(f"Primary Tumor: {n_tumour}")
    print(f"Solid Tissue Normal: {n_normal}")

    if n_tumour == 0 or n_normal == 0:
        raise RuntimeError(
            "Tumour or normal group is empty."
        )

    # --------------------------------------------------
    # PROBE QC
    # --------------------------------------------------

    print()
    print("RECALCULATING PROBE MISSINGNESS")
    print("=" * 34)

    n_probes = matrix.shape[0]

    probe_missing = np.zeros(
        n_probes,
        dtype=np.float32,
    )

    for start in range(
        0,
        n_probes,
        BLOCK_SIZE
    ):

        end = min(
            start + BLOCK_SIZE,
            n_probes
        )

        block = np.asarray(
            matrix[
                start:end,
                keep_indices
            ],
            dtype=np.float32,
        )

        probe_missing[
            start:end
        ] = (
            np.isnan(block).mean(axis=1)
            * 100
        )

        print(
            f"Probe QC: "
            f"{end:,}/{n_probes:,}"
        )

    probe_pass = (
        probe_missing
        <= PROBE_MISSING_MAX
    )

    probe_qc = pd.DataFrame(
        {
            "cpg_id":
                probes["cpg_id"],
            "missing_percent":
                probe_missing,
            "qc_pass":
                probe_pass,
        }
    )

    probe_qc.to_csv(
        PROBE_FILTER_OUT,
        sep="\t",
        index=False,
    )

    n_pass = int(
        probe_pass.sum()
    )

    print()
    print(
        f"CpGs retained at <= "
        f"{PROBE_MISSING_MAX:.0f}% missing:"
    )

    print(
        f"{n_pass:,} / "
        f"{n_probes:,}"
    )

    # --------------------------------------------------
    # Differential methylation
    # --------------------------------------------------

    print()
    print("DIFFERENTIAL METHYLATION")
    print("=" * 28)

    passed_indices = np.where(
        probe_pass
    )[0]

    result_chunks = []

    with warnings.catch_warnings():

        warnings.simplefilter(
            "ignore"
        )

        for start in range(
            0,
            len(passed_indices),
            BLOCK_SIZE
        ):

            selected = passed_indices[
                start:
                start + BLOCK_SIZE
            ]

            beta = np.asarray(
                matrix[
                    selected,
                    :
                ][:, keep_indices],
                dtype=np.float32,
            )

            tumour_beta = beta[
                :,
                tumour_mask
            ]

            normal_beta = beta[
                :,
                normal_mask
            ]

            # Biological effect size on beta scale
            mean_tumour = np.nanmean(
                tumour_beta,
                axis=1,
            )

            mean_normal = np.nanmean(
                normal_beta,
                axis=1,
            )

            delta_beta = (
                mean_tumour
                - mean_normal
            )

            # Statistical testing on M-value scale
            tumour_m = beta_to_m(
                tumour_beta
            )

            normal_m = beta_to_m(
                normal_beta
            )

            test = stats.ttest_ind(
                tumour_m,
                normal_m,
                axis=1,
                equal_var=False,
                nan_policy="omit",
            )

            chunk = pd.DataFrame(
                {
                    "cpg_id":
                        probes.iloc[
                            selected
                        ]["cpg_id"]
                        .to_numpy(),

                    "mean_beta_tumour":
                        mean_tumour,

                    "mean_beta_normal":
                        mean_normal,

                    "delta_beta":
                        delta_beta,

                    "t_statistic":
                        test.statistic,

                    "p_value":
                        test.pvalue,

                    "tumour_n":
                        np.sum(
                            ~np.isnan(
                                tumour_beta
                            ),
                            axis=1,
                        ),

                    "normal_n":
                        np.sum(
                            ~np.isnan(
                                normal_beta
                            ),
                            axis=1,
                        ),
                }
            )

            result_chunks.append(
                chunk
            )

            done = min(
                start + BLOCK_SIZE,
                len(passed_indices)
            )

            print(
                f"Tested "
                f"{done:,}/"
                f"{len(passed_indices):,} CpGs"
            )

    results = pd.concat(
        result_chunks,
        ignore_index=True,
    )

    # --------------------------------------------------
    # Multiple testing correction
    # --------------------------------------------------

    print()
    print(
        "Applying Benjamini-Hochberg "
        "FDR correction..."
    )

    results["fdr"] = bh_adjust(
        results["p_value"].to_numpy()
    )

    results["significant"] = (
        (results["fdr"] < FDR_THRESHOLD)
        &
        (
            results["delta_beta"].abs()
            >= DELTA_BETA_THRESHOLD
        )
    )

    results["direction"] = "Not significant"

    results.loc[
        results["significant"]
        &
        (results["delta_beta"] > 0),
        "direction"
    ] = "Hypermethylated"

    results.loc[
        results["significant"]
        &
        (results["delta_beta"] < 0),
        "direction"
    ] = "Hypomethylated"

    results = results.sort_values(
        [
            "significant",
            "fdr",
        ],
        ascending=[
            False,
            True,
        ],
    )

    results.to_csv(
        RESULTS_OUT,
        sep="\t",
        index=False,
    )

    significant = results[
        results["significant"]
    ].copy()

    significant.to_csv(
        SIGNIFICANT_OUT,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # Final summary
    # --------------------------------------------------

    hyper = int(
        (
            significant["direction"]
            == "Hypermethylated"
        ).sum()
    )

    hypo = int(
        (
            significant["direction"]
            == "Hypomethylated"
        ).sum()
    )

    print()
    print("=" * 45)
    print("DIFFERENTIAL METHYLATION SUMMARY")
    print("=" * 45)

    print(
        f"Samples analysed: "
        f"{len(keep_indices)}"
    )

    print(
        f"  Tumour: {n_tumour}"
    )

    print(
        f"  Normal: {n_normal}"
    )

    print(
        f"CpGs tested: "
        f"{len(results):,}"
    )

    print()
    print(
        f"Significance criteria:"
    )

    print(
        f"  FDR < {FDR_THRESHOLD}"
    )

    print(
        f"  |delta-beta| >= "
        f"{DELTA_BETA_THRESHOLD}"
    )

    print()
    print(
        f"Significant DMPs: "
        f"{len(significant):,}"
    )

    print(
        f"Hypermethylated: "
        f"{hyper:,}"
    )

    print(
        f"Hypomethylated: "
        f"{hypo:,}"
    )

    print()
    print("Saved:")
    print(RESULTS_OUT)
    print(SIGNIFICANT_OUT)
    print(SAMPLE_FILTER_OUT)
    print(PROBE_FILTER_OUT)

    print()
    print(
        "DIFFERENTIAL METHYLATION PASSED"
    )


if __name__ == "__main__":
    main()
