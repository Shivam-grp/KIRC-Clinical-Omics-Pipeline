from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_snakefile_exists():
    assert (ROOT / "workflow" / "Snakefile").exists()


def test_workflow_config_exists():
    assert (ROOT / "config" / "workflow.yaml").exists()


def test_snakefile_contains_core_rules():
    text = (
        ROOT /
        "workflow" /
        "Snakefile"
    ).read_text()

    expected_rules = [
        "rule methylation_matrix",
        "rule differential_methylation",
        "rule annotate_methylation",
        "rule integrate_multiomics",
        "rule multiomics_visualisation",
    ]

    for rule in expected_rules:
        assert rule in text
