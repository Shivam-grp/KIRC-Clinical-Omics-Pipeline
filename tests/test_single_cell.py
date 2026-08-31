from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "single_cell" / "01_scrna_pipeline.py"


def test_scrna_pipeline_exists():
    assert SCRIPT.exists()


def test_scrna_pipeline_contains_core_steps():
    text = SCRIPT.read_text()

    required = [
        "calculate_qc_metrics",
        "filter_cells",
        "normalize_total",
        "highly_variable_genes",
        "tl.pca",
        "pp.neighbors",
        "tl.umap",
        "tl.leiden",
        "rank_genes_groups",
    ]

    for item in required:
        assert item in text


def test_scrna_results_exist():
    results = ROOT / "single_cell" / "results"

    assert (results / "scrna_analysis_summary.tsv").exists()
    assert (results / "leiden_marker_genes.tsv").exists()
