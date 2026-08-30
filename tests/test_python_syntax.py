from pathlib import Path
import py_compile


ROOT = Path(__file__).resolve().parents[1]


def test_all_pipeline_scripts_compile():
    scripts = sorted(
        (ROOT / "scripts").glob("*.py")
    )

    assert scripts, "No Python scripts found"

    for script in scripts:
        py_compile.compile(
            str(script),
            doraise=True
        )
