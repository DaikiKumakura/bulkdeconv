# Consensus (mouse) vs flow cytometry: Petitprez et al. 2020 data

Data: `dataset_petitprez.rda` (1,132,974 bytes, SHA-256 `cc8e482b869e59f20913956aad8c7584a1e9b7db538c8217a9cdaa0afaf078ca`)
from `github.com/omnideconv/immunedeconv`, folder `data/`, downloaded 2026-11-26. It holds bulk RNA-seq (53,819 genes x 6
samples: blood and spleen of 3 TC-1 tumour-bearing mice; non-integer values, treated as TPM) and flow cytometry fractions
for T cell, CD8+ T cell, B cell, NK cell, macrophage, macrophage/monocyte, neutrophil, monocyte, eosinophil and granulocyte
(Petitprez et al. 2020, Genome Medicine; raw reads are E-MTAB-9271 in ArrayExpress). The file is not in this repository and
its licence is not stated on GitHub; it is used for this check only. Download it into `data/` to rerun.

Run from the package root: `Rscript validation/consensus-mouse/validate.R` (seed 20261126, 2000 bootstrap resamples).
Outputs: `results/consensus-vs-flow.csv`, `results/consensus-vs-flow-summary.csv`.

Methods run: `mmcp_counter`, `immucellai_mouse`, `base`. `seqimmucc` was skipped because it accepts counts only and the file
holds TPM-like values. Five lineages have flow values and at least two methods (B cell, CD8 T, NK, monocyte, neutrophil).

## Result (Spearman rho against flow, n = 6 samples)

| lineage | methods | consensus | median method | best method | worst method |
| --- | --- | --- | --- | --- | --- |
| B cell | 3 | 0.49 | 0.49 | 0.83 | 0.43 |
| CD8 T | 3 | 0.09 | -0.09 | 0.46 | -0.83 |
| NK | 3 | 0.49 | 0.77 | 0.83 | -0.31 |
| monocyte | 2 | 0.32 | 0.31 | 0.37 | 0.26 |
| neutrophil | 2 | 0.84 | 0.83 | 0.89 | 0.77 |

Facts: the consensus rho was at least the median method's in 4 of 5 lineages (three of them by 0.01 or less) and at least the
best single method's in 0 of 5. The bootstrap 95% intervals span almost the whole range (for example -1 to 1 for B cell,
monocyte and CD8 T).

## Reading it

- This check cannot tell the consensus from a single method. Six samples are 3 mice sampled in two tissues, so the samples
  are not independent, and a blood-versus-spleen difference can drive the correlations. Intervals are almost -1 to 1.
- The mouse table is therefore still not validated. A usable check needs more independent mouse samples with flow data. The
  Chen et al. 2017 set (ArrayExpress E-MTAB-6458, 12 samples of spleen, bone marrow, blood and lymph node, flow fractions for
  four cell types) is public only as raw FASTQ files; using it would need alignment and quantification, which was not done.
- Not tested: CD4 T, dendritic, macrophage (flow values near zero), fibroblast, endothelial, and any tumour tissue sample.
