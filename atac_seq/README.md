# TCGA-KIRC ATAC-seq Analysis

This module demonstrates reproducible chromatin-accessibility analysis
using public TCGA-KIRC ATAC-seq data.

## Analysis

The workflow performs:

1. ATAC accessibility-matrix validation
2. Peak-coordinate and width QC
3. Sample-level accessibility QC
4. Genomic peak annotation
5. Highly variable accessible-region discovery
6. PCA
7. Sample correlation analysis
8. Exploratory chromatin-accessibility clustering
9. Cluster-associated accessibility testing
10. Benjamini-Hochberg FDR correction
11. Publication-style visualisation

## Outputs

Results are written to:

atac_seq/results/

including:

- atac_analysis_summary.tsv
- atac_sample_qc.tsv
- atac_sample_pca.tsv
- atac_exploratory_clusters.tsv
- peak_annotation_summary.tsv
- top_variable_atac_peaks.tsv
- top_cluster_associated_atac_peaks.tsv
- significant_cluster_associated_peaks.tsv

Figures are written to:

atac_seq/figures/

including:

- ATAC PCA
- sample correlation heatmap
- variable accessibility heatmap
- genomic annotation distribution
- peak-width distribution

The discovered ATAC groups are unsupervised exploratory accessibility
subgroups and should not be interpreted as predefined clinical classes.

## Reproducibility

```bash
python atac_seq/01_kirc_atac_analysis.py
pytest -q
