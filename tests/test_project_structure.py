from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_required_directories_exist():
    required = [
        "scripts",
        "workflow",
        "config",
        "docs",
        "tests",
    ]

    for directory in required:
        assert (ROOT / directory).is_dir(), f"Missing directory: {directory}"


def test_required_project_files_exist():
    required = [
        "README.md",
        "pyproject.toml",
        "uv.lock",
        "workflow/Snakefile",
        "config/workflow.yaml",
    ]

    for filename in required:
        assert (ROOT / filename).exists(), f"Missing file: {filename}"


def test_analysis_scripts_exist():
    expected = [
        "01_discover_rnaseq.py",
        "02_prepare_rnaseq_cohort.py",
        "03_inspect_duplicate_aliquots.py",
        "04_download_rnaseq.py",
        "05_build_count_matrix.py",
        "06_rnaseq_qc.py",
        "07_differential_expression.py",
        "08_de_visualisation.py",
        "09_pathway_enrichment.py",
        "10_discover_methylation.py",
        "11_compare_methylation_platforms.py",
        "12_prepare_methylation_cohort.py",
        "13_download_methylation.py",
        "14_build_methylation_matrix.py",
        "15_differential_methylation.py",
        "16_annotate_methylation.py",
        "17_integrate_rna_methylation.py",
        "18_multiomics_visualisation.py",
    ]

    for script in expected:
        assert (ROOT / "scripts" / script).exists(), f"Missing script: {script}"
