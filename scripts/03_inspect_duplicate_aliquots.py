"""Inspect GDC aliquot metadata for duplicated TCGA-KIRC RNA-seq samples."""

from pathlib import Path
import json

import pandas as pd
import requests


MANIFEST = Path(
    "/mnt/e/KIRC_data/downloads/tcga_kirc_rnaseq_cohort.tsv"
)

API_URL = "https://api.gdc.cancer.gov/files"


df = pd.read_csv(MANIFEST, sep="\t")

dup = df[
    df.duplicated("sample_submitter_id", keep=False)
].copy()

duplicate_ids = dup["file_id"].tolist()

filters = {
    "op": "in",
    "content": {
        "field": "file_id",
        "value": duplicate_ids,
    },
}

fields = [
    "file_id",
    "file_name",
    "file_size",
    "created_datetime",
    "updated_datetime",
    "cases.submitter_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_type",
    "cases.samples.portions.submitter_id",
    "cases.samples.portions.analytes.submitter_id",
    "cases.samples.portions.analytes.analyte_type",
    "cases.samples.portions.analytes.aliquots.aliquot_id",
    "cases.samples.portions.analytes.aliquots.submitter_id",
]

params = {
    "filters": json.dumps(filters),
    "fields": ",".join(fields),
    "format": "JSON",
    "size": "100",
}

print("Querying GDC aliquot metadata...")

response = requests.get(API_URL, params=params, timeout=120)
response.raise_for_status()

hits = response.json()["data"]["hits"]

rows = []

for hit in hits:
    for case in hit.get("cases", []):
        for sample in case.get("samples", []):
            for portion in sample.get("portions", []):
                for analyte in portion.get("analytes", []):
                    for aliquot in analyte.get("aliquots", []):
                        rows.append(
                            {
                                "sample_id": sample.get("submitter_id"),
                                "sample_type": sample.get("sample_type"),
                                "portion": portion.get("submitter_id"),
                                "analyte": analyte.get("submitter_id"),
                                "analyte_type": analyte.get("analyte_type"),
                                "aliquot": aliquot.get("submitter_id"),
                                "aliquot_uuid": aliquot.get("aliquot_id"),
                                "file_id": hit.get("file_id"),
                                "file_size": hit.get("file_size"),
                                "created": hit.get("created_datetime"),
                                "updated": hit.get("updated_datetime"),
                            }
                        )

result = pd.DataFrame(rows)

result = result.sort_values(
    ["sample_id", "aliquot", "file_id"]
)

print()
print(result.to_string(index=False))

output = Path(
    "/mnt/e/KIRC_data/downloads/"
    "tcga_kirc_duplicate_aliquots.tsv"
)

result.to_csv(output, sep="\t", index=False)

print()
print("Saved:")
print(output)
