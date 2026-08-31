"""
TCGA-KIRC ATAC-seq analysis.

Performs:
- normalized accessibility matrix validation
- peak-coordinate QC
- sample accessibility summaries
- variable-peak discovery
- PCA
- sample correlation
- unsupervised accessibility subgroup discovery
- exploratory subgroup-associated accessibility testing
- genomic peak annotation summaries
- reproducible result/figure export

The subgroup comparison is exploratory and unsupervised; it should not be
interpreted as a predefined clinical comparison.
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]

MATRIX_FILE = ROOT / "atac_seq/data/kirc/KIRC_Log2norm.txt"
PEAK_FILE = ROOT / "atac_seq/data/kirc/KIRC_peakCalls.txt"

RESULTS = ROOT / "atac_seq/results"
FIGURES = ROOT / "atac_seq/figures"

RESULTS.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)


def bh_fdr(pvalues):
    """Benjamini-Hochberg multiple-testing correction."""
    p = np.asarray(pvalues, dtype=float)

    valid = np.isfinite(p)
    adjusted = np.full(len(p), np.nan)

    pv = p[valid]

    if len(pv) == 0:
        return adjusted

    order = np.argsort(pv)
    ranked = pv[order]

    n = len(ranked)

    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    corrected = np.empty(n)
    corrected[order] = q

    adjusted[valid] = corrected

    return adjusted


def load_matrix(path):
    """
    Load TCGA cancer-type ATAC matrix.

    The first five columns contain genomic/peak metadata.
    Columns from position 6 onward contain normalized sample accessibility.
    """

    print("\nLoading normalized ATAC matrix...")

    df = pd.read_csv(
        path,
        sep="\t",
        low_memory=False,
    )

    if df.shape[1] <= 5:
        raise ValueError(
            f"Expected >5 columns in ATAC matrix; found {df.shape[1]}"
        )

    print(f"Raw matrix: {df.shape[0]:,} rows x {df.shape[1]:,} columns")

    # First five columns = genomic metadata
    metadata = df.iloc[:, :5].copy()

    # Remaining columns = biological ATAC samples
    X = df.iloc[:, 5:].apply(
        pd.to_numeric,
        errors="coerce",
    ).astype("float32")

    # Construct UNIQUE peak IDs from genomic coordinates.
    # Never use chromosome alone as an index because chr1/chr2 etc.
    # occur thousands of times.
    if metadata.shape[1] >= 3:
        peak_ids = (
            metadata.iloc[:, 0].astype(str)
            + ":"
            + metadata.iloc[:, 1].astype(str)
            + "-"
            + metadata.iloc[:, 2].astype(str)
        )
    else:
        peak_ids = pd.Series(
            [f"ATAC_peak_{i+1}" for i in range(len(metadata))]
        )

    # Guarantee uniqueness even if duplicate genomic coordinates occur
    duplicated = peak_ids.duplicated(keep=False)

    if duplicated.any():
        peak_ids = pd.Series(
            [
                f"{pid}|peak_{i+1}"
                if duplicated.iloc[i]
                else pid
                for i, pid in enumerate(peak_ids)
            ]
        )

    X.index = peak_ids.values

    X = X.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    X = X.dropna(
        axis=0,
        how="all",
    )

    X = X.dropna(
        axis=1,
        how="all",
    )

    # Missing values should be rare; use sample median if present
    if X.isna().values.any():
        X = X.fillna(
            X.median(axis=0)
        )

    print(
        f"Accessibility matrix: "
        f"{X.shape[0]:,} peaks x {X.shape[1]:,} samples"
    )

    print(
        f"Unique peak IDs: "
        f"{X.index.nunique():,}/{len(X):,}"
    )

    return X


def load_peaks(path):
    print("\nLoading KIRC peak annotations...")

    preview = pd.read_csv(
        path,
        sep="\t",
        header=None,
        nrows=2,
        dtype=str,
    )

    first = str(preview.iloc[0, 0]).strip().lower()

    # File appears to contain a real header
    if first in {
        "sequence",
        "chrom",
        "chr",
        "chromosome",
    }:

        peaks = pd.read_csv(
            path,
            sep="\t",
            low_memory=False,
        )

    else:

        peaks = pd.read_csv(
            path,
            sep="\t",
            header=None,
            low_memory=False,
        )

        default_columns = [
            "sequence",
            "start",
            "end",
            "name",
            "score",
            "annotation",
            "percentGC",
            "percentAT",
        ]

        peaks.columns = (
            default_columns[:peaks.shape[1]]
            + [
                f"extra_{i}"
                for i in range(
                    peaks.shape[1] - len(default_columns)
                )
            ]
        )

    # Normalize likely column names
    rename = {}

    for col in peaks.columns:

        c = str(col).lower()

        if c in {"chr", "chrom", "chromosome"}:
            rename[col] = "sequence"

        elif c == "start":
            rename[col] = "start"

        elif c == "end":
            rename[col] = "end"

        elif c in {"peak", "peak_name"}:
            rename[col] = "name"

    peaks = peaks.rename(columns=rename)

    for col in ["start", "end", "score", "percentGC", "percentAT"]:

        if col in peaks.columns:

            peaks[col] = pd.to_numeric(
                peaks[col],
                errors="coerce",
            )

    if {"start", "end"}.issubset(peaks.columns):

        peaks["peak_width"] = (
            peaks["end"] - peaks["start"]
        )

    print(
        f"Peak annotations: {len(peaks):,}"
    )

    return peaks


def annotation_class(value):
    x = str(value).lower()

    if (
        "promoter" in x
        or "tss" in x
    ):
        return "Promoter/TSS"

    if "5' utr" in x or "5utr" in x:
        return "5' UTR"

    if "3' utr" in x or "3utr" in x:
        return "3' UTR"

    if "exon" in x:
        return "Exonic"

    if "intron" in x:
        return "Intronic"

    if (
        "intergenic" in x
        or "distal" in x
    ):
        return "Intergenic/Distal"

    return "Other"


def main():

    print("=" * 72)
    print("TCGA-KIRC ATAC-SEQ ANALYSIS")
    print("=" * 72)

    if not MATRIX_FILE.exists():
        raise FileNotFoundError(MATRIX_FILE)

    if not PEAK_FILE.exists():
        raise FileNotFoundError(PEAK_FILE)

    X = load_matrix(MATRIX_FILE)
    peaks = load_peaks(PEAK_FILE)

    # ------------------------------------------------------
    # Basic matrix QC
    # ------------------------------------------------------

    missing_percent = (
        X.isna()
        .mean()
        .mean()
        * 100
    )

    sample_qc = pd.DataFrame(
        {
            "sample": X.columns,
            "mean_accessibility": X.mean(axis=0).values,
            "median_accessibility": X.median(axis=0).values,
            "sd_accessibility": X.std(axis=0).values,
        }
    )

    sample_qc.to_csv(
        RESULTS / "atac_sample_qc.tsv",
        sep="\t",
        index=False,
    )

    print("\n===== MATRIX QC =====")

    print(
        f"Peaks: {X.shape[0]:,}"
    )

    print(
        f"Samples: {X.shape[1]:,}"
    )

    print(
        f"Missing values: {missing_percent:.4f}%"
    )

    # ------------------------------------------------------
    # Peak-coordinate QC
    # ------------------------------------------------------

    if "peak_width" in peaks.columns:

        valid_width = peaks[
            peaks["peak_width"] > 0
        ]

        width_median = (
            valid_width["peak_width"].median()
        )

        width_mean = (
            valid_width["peak_width"].mean()
        )

        print("\n===== PEAK QC =====")

        print(
            f"Median peak width: "
            f"{width_median:,.1f} bp"
        )

        print(
            f"Mean peak width: "
            f"{width_mean:,.1f} bp"
        )

    # ------------------------------------------------------
    # Genomic annotations
    # ------------------------------------------------------

    if "annotation" in peaks.columns:

        peaks["annotation_class"] = (
            peaks["annotation"]
            .map(annotation_class)
        )

        annotation_summary = (
            peaks["annotation_class"]
            .value_counts()
            .rename_axis("annotation")
            .reset_index(name="peaks")
        )

        annotation_summary["percent"] = (
            annotation_summary["peaks"]
            / annotation_summary["peaks"].sum()
            * 100
        )

        annotation_summary.to_csv(
            RESULTS / "peak_annotation_summary.tsv",
            sep="\t",
            index=False,
        )

        print("\n===== GENOMIC ANNOTATION =====")

        print(
            annotation_summary
            .to_string(index=False)
        )

    # ------------------------------------------------------
    # Variable accessibility peaks
    # ------------------------------------------------------

    variance = X.var(axis=1)

    variance_table = pd.DataFrame(
        {
            "peak_id": variance.index,
            "variance": variance.values,
            "mean_accessibility": (
                X.mean(axis=1).values
            ),
        }
    ).sort_values(
        "variance",
        ascending=False,
    )

    variance_table.head(1000).to_csv(
        RESULTS / "top_variable_atac_peaks.tsv",
        sep="\t",
        index=False,
    )

    n_variable = min(
        5000,
        len(variance_table),
    )

    variable_ids = (
        variance_table
        .head(n_variable)["peak_id"]
        .tolist()
    )

    variable_matrix = (
        X.loc[variable_ids]
        .T
    )

    print("\n===== VARIABLE ACCESSIBILITY =====")

    print(
        f"Using top {n_variable:,} variable peaks"
    )

    # ------------------------------------------------------
    # PCA
    # ------------------------------------------------------

    scaled = StandardScaler().fit_transform(
        variable_matrix
    )

    n_components = min(
        5,
        scaled.shape[0] - 1,
        scaled.shape[1],
    )

    pca = PCA(
        n_components=n_components,
        random_state=42,
    )

    pcs = pca.fit_transform(scaled)

    pca_columns = [
        f"PC{i+1}"
        for i in range(n_components)
    ]

    pca_df = pd.DataFrame(
        pcs,
        index=variable_matrix.index,
        columns=pca_columns,
    )

    pca_df.index.name = "sample"

    pca_df.to_csv(
        RESULTS / "atac_sample_pca.tsv",
        sep="\t",
    )

    print("\n===== PCA =====")

    for i, ratio in enumerate(
        pca.explained_variance_ratio_,
        start=1,
    ):
        print(
            f"PC{i}: {ratio*100:.2f}%"
        )

    # ------------------------------------------------------
    # Exploratory ATAC subgroups
    # ------------------------------------------------------

    cluster_pcs = pcs[
        :,
        :min(5, pcs.shape[1])
    ]

    km = KMeans(
        n_clusters=2,
        random_state=42,
        n_init=50,
    )

    groups = km.fit_predict(
        cluster_pcs
    )

    cluster_df = pd.DataFrame(
        {
            "sample": variable_matrix.index,
            "ATAC_cluster": groups,
        }
    )

    cluster_df.to_csv(
        RESULTS / "atac_exploratory_clusters.tsv",
        sep="\t",
        index=False,
    )

    group_counts = (
        cluster_df["ATAC_cluster"]
        .value_counts()
        .sort_index()
    )

    print("\n===== EXPLORATORY ACCESSIBILITY CLUSTERS =====")

    print(group_counts.to_string())

    # ------------------------------------------------------
    # Cluster-associated accessibility testing
    # ------------------------------------------------------

    sample_by_peak = X.T

    group0 = sample_by_peak.iloc[
        np.where(groups == 0)[0]
    ]

    group1 = sample_by_peak.iloc[
        np.where(groups == 1)[0]
    ]

    if (
        len(group0) >= 2
        and len(group1) >= 2
    ):

        statistic, pvalue = ttest_ind(
            group1.to_numpy(),
            group0.to_numpy(),
            axis=0,
            equal_var=False,
            nan_policy="omit",
        )

        delta = (
            group1.mean(axis=0)
            - group0.mean(axis=0)
        )

        fdr = bh_fdr(pvalue)

        diff = pd.DataFrame(
            {
                "peak_id": X.index,
                "mean_cluster0": (
                    group0.mean(axis=0).values
                ),
                "mean_cluster1": (
                    group1.mean(axis=0).values
                ),
                "delta_accessibility": delta.values,
                "t_statistic": statistic,
                "pvalue": pvalue,
                "FDR": fdr,
            }
        )

        diff["significant"] = (
            (diff["FDR"] < 0.05)
            &
            (
                diff["delta_accessibility"]
                .abs()
                >= 1
            )
        )

        diff = diff.sort_values(
            ["FDR", "pvalue"],
            ascending=True,
        )

        significant = diff[
            diff["significant"]
        ].copy()

        # Attach genomic annotation where peak IDs match
        if "name" in peaks.columns:

            annotation_cols = [
                c
                for c in [
                    "name",
                    "sequence",
                    "start",
                    "end",
                    "annotation",
                    "annotation_class",
                ]
                if c in peaks.columns
            ]

            annotation_lookup = (
                peaks[annotation_cols]
                .drop_duplicates("name")
            )

            diff = diff.merge(
                annotation_lookup,
                left_on="peak_id",
                right_on="name",
                how="left",
            )

            significant = significant.merge(
                annotation_lookup,
                left_on="peak_id",
                right_on="name",
                how="left",
            )

        diff.head(1000).to_csv(
            RESULTS
            / "top_cluster_associated_atac_peaks.tsv",
            sep="\t",
            index=False,
        )

        significant.head(500).to_csv(
            RESULTS
            / "significant_cluster_associated_peaks.tsv",
            sep="\t",
            index=False,
        )

        print(
            "\nExploratory cluster-associated peaks:"
        )

        print(
            f"FDR < 0.05 and |delta| >= 1: "
            f"{len(significant):,}"
        )

    else:

        diff = None

        print(
            "\nCluster sizes too small for "
            "accessibility testing."
        )

    # ------------------------------------------------------
    # FIGURE 1: PCA
    # ------------------------------------------------------

    plt.figure(figsize=(8, 6))

    plt.scatter(
        pca_df["PC1"],
        pca_df["PC2"],
        s=70,
    )

    for sample in pca_df.index:

        plt.annotate(
            str(sample),
            (
                pca_df.loc[sample, "PC1"],
                pca_df.loc[sample, "PC2"],
            ),
            fontsize=6,
            xytext=(3, 3),
            textcoords="offset points",
        )

    plt.xlabel(
        f"PC1 "
        f"({pca.explained_variance_ratio_[0]*100:.1f}%)"
    )

    plt.ylabel(
        f"PC2 "
        f"({pca.explained_variance_ratio_[1]*100:.1f}%)"
    )

    plt.title(
        "TCGA-KIRC ATAC-seq accessibility PCA"
    )

    plt.tight_layout()

    plt.savefig(
        FIGURES / "kirc_atac_pca.png",
        dpi=200,
    )

    plt.close()

    # ------------------------------------------------------
    # FIGURE 2: sample correlation
    # ------------------------------------------------------

    corr = variable_matrix.T.corr()

    plt.figure(figsize=(8, 7))

    im = plt.imshow(
        corr,
        aspect="auto",
    )

    plt.colorbar(
        im,
        label="Pearson correlation",
    )

    plt.xticks(
        range(len(corr.columns)),
        corr.columns,
        rotation=90,
        fontsize=6,
    )

    plt.yticks(
        range(len(corr.index)),
        corr.index,
        fontsize=6,
    )

    plt.title(
        "KIRC ATAC-seq sample correlation"
    )

    plt.tight_layout()

    plt.savefig(
        FIGURES / "kirc_atac_sample_correlation.png",
        dpi=200,
    )

    plt.close()

    # ------------------------------------------------------
    # FIGURE 3: variable peak heatmap
    # ------------------------------------------------------

    heat_ids = (
        variance_table
        .head(min(50, len(variance_table)))
        ["peak_id"]
    )

    heat = X.loc[heat_ids].copy()

    heat = heat.sub(
        heat.mean(axis=1),
        axis=0,
    )

    sd = heat.std(axis=1).replace(
        0,
        1,
    )

    heat = heat.div(
        sd,
        axis=0,
    )

    plt.figure(figsize=(10, 8))

    im = plt.imshow(
        heat,
        aspect="auto",
        interpolation="nearest",
    )

    plt.colorbar(
        im,
        label="Peak accessibility z-score",
    )

    plt.xticks(
        range(len(heat.columns)),
        heat.columns,
        rotation=90,
        fontsize=6,
    )

    plt.yticks([])

    plt.xlabel("Samples")
    plt.ylabel("Top variable ATAC peaks")

    plt.title(
        "Top variable chromatin-accessibility peaks"
    )

    plt.tight_layout()

    plt.savefig(
        FIGURES / "kirc_atac_variable_peak_heatmap.png",
        dpi=200,
    )

    plt.close()

    # ------------------------------------------------------
    # FIGURE 4: peak annotation
    # ------------------------------------------------------

    if "annotation_class" in peaks.columns:

        counts = (
            peaks["annotation_class"]
            .value_counts()
        )

        plt.figure(figsize=(9, 5))

        counts.plot(
            kind="bar",
        )

        plt.ylabel("Number of peaks")
        plt.xlabel("Genomic annotation")

        plt.title(
            "KIRC ATAC-seq peak genomic distribution"
        )

        plt.xticks(
            rotation=35,
            ha="right",
        )

        plt.tight_layout()

        plt.savefig(
            FIGURES
            / "kirc_atac_peak_annotation.png",
            dpi=200,
        )

        plt.close()

    # ------------------------------------------------------
    # FIGURE 5: peak width
    # ------------------------------------------------------

    if "peak_width" in peaks.columns:

        width = peaks[
            (
                peaks["peak_width"] > 0
            )
            &
            (
                peaks["peak_width"] <
                peaks["peak_width"].quantile(0.99)
            )
        ]["peak_width"]

        plt.figure(figsize=(8, 5))

        plt.hist(
            width,
            bins=50,
        )

        plt.xlabel("Peak width (bp)")
        plt.ylabel("Number of peaks")

        plt.title(
            "KIRC ATAC-seq peak-width distribution"
        )

        plt.tight_layout()

        plt.savefig(
            FIGURES
            / "kirc_atac_peak_width_distribution.png",
            dpi=200,
        )

        plt.close()

    # ------------------------------------------------------
    # Analysis summary
    # ------------------------------------------------------

    summary = pd.DataFrame(
        [
            ["matrix_peaks", X.shape[0]],
            ["samples", X.shape[1]],
            ["annotated_peaks", len(peaks)],
            ["variable_peaks_used", n_variable],
            [
                "PC1_variance_percent",
                round(
                    pca.explained_variance_ratio_[0]
                    * 100,
                    3,
                ),
            ],
            [
                "PC2_variance_percent",
                round(
                    pca.explained_variance_ratio_[1]
                    * 100,
                    3,
                ),
            ],
            [
                "ATAC_cluster_0_samples",
                int((groups == 0).sum()),
            ],
            [
                "ATAC_cluster_1_samples",
                int((groups == 1).sum()),
            ],
        ],
        columns=[
            "metric",
            "value",
        ],
    )

    summary.to_csv(
        RESULTS / "atac_analysis_summary.tsv",
        sep="\t",
        index=False,
    )

    print("\n===== SAVED RESULTS =====")

    for path in sorted(
        RESULTS.glob("*")
    ):
        print(path)

    print("\n===== SAVED FIGURES =====")

    for path in sorted(
        FIGURES.glob("*")
    ):
        print(path)

    print(
        "\n=========================================="
    )

    print(
        "KIRC ATAC-SEQ ANALYSIS PASSED"
    )

    print(
        "=========================================="
    )


if __name__ == "__main__":
    main()
