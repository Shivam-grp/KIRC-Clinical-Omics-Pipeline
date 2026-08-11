"""Discover harmonized TCGA-KIRC DNA methylation beta-value files."""

from datetime import datetime, timezone
from pathlib import Path
import json

import pandas as pd
import requests


API_URL = "https://api.gdc.cancer.gov/files"

OUT_DIR = Path("/mnt/e/KIRC_data/downloads")
OUT_DIR.mkdir(parents=True, exist_ok=True)

MANIFEST = OUT_DIR / "tcga_kirc_methylation_manifest.tsv"
PROVENANCE = OUT_DIR / "tcga_kirc_methylation_provenance.json"


filters = {
    "op": "and",
    "content": [
        {
            "op": "=",
            "content": {
                "field": "cases.project.project_id",
                "value": ["TCGA-KIRC"],
            },
        },
        {
            "op": "=",
            "content": {
                "field": "data_type",
                "value": ["Methylation Beta Value"],
            },
        },
        {
            "op": "=",
            "content": {
                "field": "access",
                "value": ["open"],
            },
        },
    ],
}


fields = [
    "file_id",
    "file_name",
    "file_size",
    "md5sum",
    "data_type",
    "data_format",
    "experimental_strategy",
    "platform",
    "cases.case_id",
    "cases.submitter_id",
    "cases.samples.sample_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_type",
]


params = {
    "filters": json.dumps(filters),
    "fields": ",".join(fields),
    "format": "JSON",
    "size": "2000",
}


print("Querying GDC for TCGA-KIRC methylation files...")

response = requests.get(
    API_URL,
    params=params,
    timeout=120,
)

response.raise_for_status()

hits = response.json()["data"]["hits"]

rows = []

for hit in hits:

    cases = hit.get("cases", [])

    if not cases:
        continue

    case = cases[0]

    samples = case.get("samples", [])

    sample = samples[0] if samples else {}

    rows.append(
        {
            "file_id": hit.get("file_id"),
            "file_name": hit.get("file_name"),
            "file_size": hit.get("file_size"),
            "md5sum": hit.get("md5sum"),
            "data_type": hit.get("data_type"),
            "data_format": hit.get("data_format"),
            "experimental_strategy": hit.get(
                "experimental_strategy"
            ),
            "platform": hit.get("platform"),
            "case_submitter_id": case.get(
                "submitter_id"
            ),
            "sample_submitter_id": sample.get(
                "submitter_id"
            ),
            "sample_type": sample.get(
                "sample_type"
            ),
        }
    )


df = pd.DataFrame(rows)

if df.empty:
    raise SystemExit(
        "ERROR: no TCGA-KIRC methylation files found."
    )


df = df.sort_values(
    [
        "platform",
        "sample_type",
        "case_submitter_id",
    ],
    na_position="last",
).reset_index(drop=True)


df.to_csv(
    MANIFEST,
    sep="\t",
    index=False,
)


provenance = {
    "project": "TCGA-KIRC",
    "data_type": "Methylation Beta Value",
    "source": "NCI Genomic Data Commons",
    "retrieved_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "number_of_files": int(len(df)),
    "manifest": str(MANIFEST),
    "filters": filters,
}


with PROVENANCE.open(
    "w",
    encoding="utf-8",
) as handle:
    json.dump(
        provenance,
        handle,
        indent=2,
    )


print()
print("METHYLATION DISCOVERY COMPLETE")
print("------------------------------")

print(f"Files discovered: {len(df)}")

print("\nPlatforms:")
print(
    df["platform"]
    .value_counts(dropna=False)
    .to_string()
)

print("\nSample types:")
print(
    df["sample_type"]
    .value_counts(dropna=False)
    .to_string()
)

total_gib = df["file_size"].sum() / (1024 ** 3)

print(
    f"\nEstimated total download: "
    f"{total_gib:.2f} GiB"
)

print()
print(f"Manifest: {MANIFEST}")
print(f"Provenance: {PROVENANCE}")
