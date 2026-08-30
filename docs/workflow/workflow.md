# KIRC Clinical Multi-Omics Workflow

## Workflow overview

mermaid
flowchart TD
    A[TCGA-KIRC RNA-seq] --> B[RNA-seq QC]
    B --> C[Differential Expression]
    C --> D[Pathway and Biological Interpretation]

    E[TCGA-KIRC 450K Methylation] --> F[Methylation Matrix and QC]
    F --> G[Differential Methylation]
    G --> H[CpG and Gene Annotation]

    C --> I[RNA x Methylation Integration]
    H --> I

    I --> J[Concordant Epigenetic Candidates]
    J --> K[Multi-Omics Visualisation]


## Pipeline stages

1. TCGA-KIRC RNA-seq cohort discovery
2. RNA-seq cohort preparation and duplicate handling
3. Checksum-validated RNA-seq acquisition
4. Count-matrix construction
5. RNA-seq quality control and filtering
6. Tumour-versus-normal differential expression
7. Differential-expression visualisation
8. Pathway enrichment and biological interpretation
9. Illumina 450K methylation cohort discovery
10. Methylation platform comparison
11. RNA-overlapping methylation cohort preparation
12. Checksum-validated methylation acquisition
13. Methylation matrix construction
14. Methylation quality control
15. Differential methylation
16. CpG-to-gene and genomic annotation
17. RNA-seq x promoter-methylation integration
18. Multi-omics candidate prioritisation and visualisation

## Reproducibility

The downstream workflow is orchestrated using Snakemake.

Configuration:

config/workflow.yaml

Workflow:

workflow/Snakefile

Dependency management:

pyproject.toml and uv.lock

Large TCGA datasets are intentionally excluded from GitHub.
The repository contains analysis code, configuration, provenance,
workflow definitions and selected derived figures/results.
