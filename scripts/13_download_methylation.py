"""Download TCGA-KIRC 450K methylation files with checksum validation."""

from pathlib import Path
import hashlib
import time

import pandas as pd
import requests


MANIFEST = Path(
    "/mnt/e/KIRC_data/downloads/"
    "tcga_kirc_methylation_450_rna_overlap.tsv"
)

OUT_DIR = Path("/mnt/e/KIRC_data/raw/methylation_450")
OUT_DIR.mkdir(parents=True, exist_ok=True)

LOG_FILE = Path(
    "/mnt/e/KIRC_data/downloads/"
    "tcga_kirc_methylation_download_log.tsv"
)

GDC_DATA = "https://api.gdc.cancer.gov/data"


def md5sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()

    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)

    return digest.hexdigest()


def download_file(file_id, file_name, expected_md5):

    folder = OUT_DIR / file_id
    folder.mkdir(parents=True, exist_ok=True)

    destination = folder / file_name

    if destination.exists():

        if md5sum(destination) == expected_md5:
            return "already_valid"

        destination.unlink()

    url = f"{GDC_DATA}/{file_id}"

    for attempt in range(1, 4):

        try:

            with requests.get(
                url,
                stream=True,
                timeout=(30, 300),
            ) as response:

                response.raise_for_status()

                with destination.open("wb") as handle:

                    for chunk in response.iter_content(
                        chunk_size=1024 * 1024
                    ):

                        if chunk:
                            handle.write(chunk)

            if md5sum(destination) != expected_md5:

                destination.unlink(missing_ok=True)

                raise RuntimeError(
                    "MD5 checksum mismatch"
                )

            return "downloaded"

        except Exception as exc:

            print(
                f"Attempt {attempt}/3 failed: {exc}"
            )

            if attempt == 3:
                return "failed"

            time.sleep(5)


def main():

    df = pd.read_csv(
        MANIFEST,
        sep="\t",
    )

    records = []

    print(f"Files to process: {len(df)}")
    print(f"Destination: {OUT_DIR}")
    print()

    for number, row in enumerate(
        df.itertuples(index=False),
        start=1,
    ):

        print(
            f"[{number}/{len(df)}] "
            f"{row.sample_submitter_id}"
        )

        status = download_file(
            row.file_id,
            row.file_name,
            row.md5sum,
        )

        print(f"    {status}")

        records.append(
            {
                "file_id": row.file_id,
                "sample_id": row.sample_submitter_id,
                "sample_type": row.sample_type,
                "status": status,
            }
        )

        pd.DataFrame(records).to_csv(
            LOG_FILE,
            sep="\t",
            index=False,
        )

    result = pd.DataFrame(records)

    print()
    print("Download summary")
    print("----------------")
    print(
        result["status"]
        .value_counts()
        .to_string()
    )

    print()
    print(f"Log: {LOG_FILE}")


if __name__ == "__main__":
    main()
