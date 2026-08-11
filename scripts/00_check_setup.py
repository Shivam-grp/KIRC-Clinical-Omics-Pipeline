"""Check the KIRC functional genomics project environment."""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pydeseq2
import pyarrow
import plotly
import streamlit


DATA_ROOT = Path("/mnt/e/KIRC_data")
REQUIRED = ("downloads", "raw", "interim", "processed", "results")


def main() -> None:
    print(f"Python: {sys.version.split()[0]}")
    print(f"Pandas: {pd.__version__}")
    print(f"NumPy: {np.__version__}")
    print(f"Data location: {DATA_ROOT}")

    missing = [
        name for name in REQUIRED
        if not (DATA_ROOT / name).is_dir()
    ]

    if missing:
        raise SystemExit(
            "Missing data directories: " + ", ".join(missing)
        )

    print("All data directories found.")
    print("All core packages imported.")
    print("SETUP CHECK PASSED")


if __name__ == "__main__":
    main()
