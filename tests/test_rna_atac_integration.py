from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "19_integrate_rna_atac.py"
RESULTS = ROOT / "regulatory_integration" / "results"
FIGURES = ROOT / "regulatory_integration" / "figures"


def test_rna_atac_script_exists():
    assert SCRIPT.exists()


def test_rna_atac_core_methods_present():
    text = SCRIPT.read_text()

    for item in [
        "load_genes",
        "load_atac",
        "link_peaks_to_genes",
        "discover_rna_de_file",
        "combined_evidence_score",
    ]:
        assert item in text


def test_rna_atac_summary_exists():
    assert (RESULTS / "rna_atac_integration_summary.tsv").exists()


def test_shared_candidates_exist():
    assert (RESULTS / "shared_rna_atac_candidates.tsv").exists()


def test_rna_atac_figure_exists():
    assert (FIGURES / "rna_atac_evidence_counts.png").exists()
