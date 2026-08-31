# TCGA-KIRC RNA–ATAC Regulatory Integration

This module integrates RNA-seq differential-expression results with
ATAC-seq chromatin-accessibility evidence in TCGA-KIRC.

## Analysis strategy

The workflow:

1. Loads significant KIRC ATAC-seq peaks
2. Uses GENCODE v22 gene annotations
3. Maps accessible regions to nearby transcription start sites
4. Identifies promoter and proximal regulatory associations
5. Loads KIRC RNA differential-expression results
6. Identifies genes independently supported by RNA and ATAC evidence
7. Ranks shared regulatory candidates
8. Generates summary tables and figures

## Scientific interpretation

The RNA and ATAC analyses represent different comparisons.

RNA-seq compares tumour versus normal tissue.

The ATAC-seq groups were obtained through unsupervised clustering of
KIRC chromatin-accessibility profiles.

Therefore RNA and ATAC effect directions are not interpreted as direct
concordance. Instead, this module identifies genes supported by
independent transcriptional and regulatory-accessibility evidence.

## Results

Outputs are written to:

regulatory_integration/results/

including:

- atac_nearest_gene_links.tsv
- rna_atac_integration_summary.tsv
- shared_rna_atac_candidates.tsv
- significant_rna_genes.tsv
- top_50_rna_atac_candidates.tsv

## Figures

Figures are written to:

regulatory_integration/figures/

including:

- RNA/ATAC evidence counts
- top integrated regulatory candidates
- ATAC peak-to-TSS distance distribution

## Reproducibility

Run:

```bash
python scripts/19_integrate_rna_atac.py
pytest -q
