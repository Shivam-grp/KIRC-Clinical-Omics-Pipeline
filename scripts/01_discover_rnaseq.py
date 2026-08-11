"""Discover TCGA-KIRC RNA-seq STAR-count files from the GDC API."""

from datetime import datetime, timezone
from pathlib import Path
import json

import pandas as pd
import requests


API_URL = "https://api.gdc.cancer.gov/files"
OUTPUT_DIR = Path("/mnt/e/KIRC_data/downloads")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MANIFEST_PATH = OUTPUT_DIR / "tcga_kirc_rnaseq_manifest.tsv"
PROVENANCE_PATH = OUTPUT_DIR / "tcga_kirc_rnaseq_provenance.json"


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
                "value": ["Gene Expression Quantification"],
            },
        },
        {
            "op": "=",
            "content": {
                "field": "analysis.workflow_type",
                "value": ["STAR - Counts"],
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
    "access",
    "data_type",
    "analysis.workflow_type",
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

print("Querying GDC API...")

response = requests.get(API_URL, params=params, timeout=120)
response.raise_for_status()

payload = response.json()
hits = payload["data"]["hits"]

rows = []

for hit in hits:
    case = hit.get("cases", [{}])[0]
    sample = case.get("samples", [{}])[0]

    rows.append(
        {
            "file_id": hit.get("file_id"),
            "file_name": hit.get("file_name"),
            "file_size": hit.get("file_size"),
            "md5sum": hit.get("md5sum"),
            "access": hit.get("access"),
            "workflow": hit.get("analysis", {}).get("workflow_type"),
            "case_id": case.get("case_id"),
            "case_submitter_id": case.get("submitter_id"),
            "sample_id": sample.get("sample_id"),
            "sample_submitter_id": sample.get("submitter_id"),
            "sample_type": sample.get("sample_type"),
        }
    )

df = pd.DataFrame(rows)

if df.empty:
    raise SystemExit("ERROR: GDC returned no RNA-seq files.")

df = df.sort_values(
    ["sample_type", "case_submitter_id"]
).reset_index(drop=True)

df.to_csv(MANIFEST_PATH, sep="\t", index=False)

provenance = {
    "project": "TCGA-KIRC",
    "data_type": "Gene Expression Quantification",
    "workflow": "STAR - Counts",
    "source": "NCI Genomic Data Commons",
    "api_endpoint": API_URL,
    "retrieved_utc": datetime.now(timezone.utc).isoformat(),
    "number_of_files": int(len(df)),
    "manifest": str(MANIFEST_PATH),
    "filters": filters,
}

with open(PROVENANCE_PATH, "w", encoding="utf-8") as handle:
    json.dump(provenance, handle, indent=2)

print()
print("RNA-seq discovery complete")
print("--------------------------")
print(f"Files discovered: {len(df)}")
print()
print("Sample types:")
print(df["sample_type"].value_counts(dropna=False).to_string())
print()
print(f"Manifest:   {MANIFEST_PATH}")
print(f"Provenance: {PROVENANCE_PATH}")
