from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "single_cell" / "02_annotate_celltypes.py"


def test_annotation_script_exists():
    assert SCRIPT.exists()


def test_annotation_uses_marker_scoring():
    text = SCRIPT.read_text()
    assert "score_genes" in text
    assert "marker_panels" in text


def test_annotation_outputs_exist():
    results = ROOT / "single_cell" / "results"

    required = [
        "cluster_celltype_annotations.tsv",
        "cluster_celltype_scores.tsv",
        "celltype_counts.tsv",
    ]

    for filename in required:
        assert (results / filename).exists()
