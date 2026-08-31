from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "atac_seq" / "01_kirc_atac_analysis.py"


def test_atac_script_exists():
    assert SCRIPT.exists()


def test_atac_core_methods_present():
    text = SCRIPT.read_text()

    for item in [
        "PCA",
        "KMeans",
        "ttest_ind",
        "bh_fdr",
        "peak_width",
        "annotation_class",
    ]:
        assert item in text


def test_atac_summary_exists():
    assert (
        ROOT / "atac_seq" / "results" / "atac_analysis_summary.tsv"
    ).exists()


def test_atac_pca_exists():
    assert (
        ROOT / "atac_seq" / "results" / "atac_sample_pca.tsv"
    ).exists()


def test_atac_figure_exists():
    assert (
        ROOT / "atac_seq" / "figures" / "kirc_atac_pca.png"
    ).exists()
