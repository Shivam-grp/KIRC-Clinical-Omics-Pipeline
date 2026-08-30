from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEF = ROOT / "containers" / "kirc_pipeline.def"


def test_apptainer_definition_exists():
    assert DEF.exists()


def test_apptainer_definition_has_required_sections():
    text = DEF.read_text()

    required = [
        "Bootstrap: docker",
        "%files",
        "%post",
        "%environment",
        "%runscript",
        "%test",
    ]

    for section in required:
        assert section in text


def test_container_uses_locked_environment():
    text = DEF.read_text()

    assert "uv.lock" in text
    assert "uv sync --frozen" in text
