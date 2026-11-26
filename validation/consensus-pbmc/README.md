# Consensus vs flow cytometry (GEO GSE107572)

Data: blood-derived immune-cell mixtures (PBMC admixed with a small fraction of PMN) of 9 healthy
donors, from the quanTIseq study (Finotello et al. 2019, Genome Medicine). TPM from
`GSE107572_tpm_PBMC_RNAseq.txt.gz` (19,423 genes x 9 samples); flow cytometry fractions from the
sample annotation of `GSE107572_series_matrix.txt.gz`. Downloaded 2026-11-26 from the GEO FTP site.

The data files are not in the repository; download the two files named above from the GEO FTP site
(`ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107572/suppl/` and `.../matrix/`) into `data/`.

Run from the package root: `Rscript validation/consensus-pbmc/validate.R`
(10 built-in human methods; seed 20261126; 2000 donor bootstrap resamples).
Outputs: `results/consensus-vs-flow.csv` (every method and the consensus, per lineage, Spearman rho with
bootstrap 95% interval) and `results/consensus-vs-flow-summary.csv`.

## Result (Spearman rho against the flow fraction, n = 9 donors)

| lineage | methods | consensus | median method | best method | worst method | consensus 95% interval |
| --- | --- | --- | --- | --- | --- | --- |
| B cell | 9 | 0.57 | 0.56 | 0.78 | 0.42 | -0.16 to 0.93 |
| dendritic | 7 | 0.00 | -0.02 | 0.65 | -0.55 | -0.73 to 0.74 |
| monocyte | 7 | 0.65 | 0.49 | 0.83 | 0.02 | 0.08 to 0.89 |
| neutrophil | 10 | 0.38 | 0.39 | 0.68 | 0.08 | -0.33 to 0.80 |
| NK | 10 | 0.83 | 0.68 | 0.80 | -0.67 | 0.34 to 1.00 |
| CD4 T | 5 | 0.72 | 0.55 | 0.71 | 0.12 | 0.09 to 0.95 |
| CD8 T | 8 | 0.77 | 0.45 | 0.92 | 0.27 | 0.26 to 0.94 |

Facts: the consensus rho was at least the median method's in 6 of 7 lineages (neutrophil is 0.02 below), and at
least the best single method's in 2 of 7 (NK, CD4 T). Across the 7 lineages the Spearman correlation between the
methods' `agreement` and the consensus rho was -0.11.

## Reading it

- Exploratory, one data set, 9 samples: the bootstrap intervals are wide (several contain 0), so the data do not
  show that the consensus is better than a good single method. They are consistent with it being about as good as
  the median method and safer than the worst one (the worst methods have negative rho in dendritic, NK).
- In this test `agreement` did not predict accuracy (rho -0.11 across 7 lineages), so it should not be read as a
  quality score.
- Not tested: tumor samples; fibroblast, endothelial, macrophage lineages (no flow values); the tumor-content
  lineage; low-abundance lineages beyond the ones here. Blood mixtures contain no tumor cells, and the methods'
  tumor signatures were not designed for them.
- One xCell output (monocyte) is constant across the 9 samples and has no correlation; it is left out of the
  method summaries.
