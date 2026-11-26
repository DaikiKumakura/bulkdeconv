# bulkdeconv

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24210481.svg)](https://doi.org/10.5281/zenodo.24210481)
<!-- badges: end -->

**Cell-type deconvolution of human and mouse tumor bulk RNA-seq: 25
published methods through one function.**

## Installation

```r
install.packages(c("remotes", "BiocManager"))
remotes::install_github("DaikiKumakura/bulkdeconv", repos = BiocManager::repositories())
```

## Quick start

bulkdeconv takes a table of gene symbols by samples (TPM or counts) and returns, for each sample, the cell-type
scores of the method you choose. This example builds a small synthetic TPM table from the MCP-counter marker genes,
so it runs as it is, with no files to download. Sample `S3` is given eight times the expression of the T cell
markers:

```r
library(bulkdeconv)

markers <- readRDS(system.file("extdata", "mcp_counter", "markers.rds", package = "bulkdeconv"))
set.seed(1)
genes <- unique(c(markers[["HUGO symbols"]], sprintf("GENE%03d", 1:300)))
tpm <- matrix(rlnorm(length(genes) * 4, meanlog = 2, sdlog = 1), length(genes), 4,
              dimnames = list(genes, paste0("S", 1:4)))
t_genes <- markers[["HUGO symbols"]][markers[["Cell population"]] == "T cells"]
tpm[t_genes, "S3"] <- tpm[t_genes, "S3"] * 8

res <- bulkdeconv(tpm, expr = "tpm", species = "human", method = "mcp_counter")
round(res$score[, 1:5], 2)             # samples x cell types
res$metadata$output_class              # what the numbers are
res$metadata$comparability             # what may be compared
res$qc                                 # input and signature-coverage checks
```

```text
   T cells CD8 T cells Cytotoxic lymphocytes B lineage NK cells
S1   11.94        7.27                 13.26     10.09     8.42
S2   13.78        3.39                  6.14     13.57     5.90
S3  154.77        2.61                 33.51      8.32     7.17
S4   12.06        2.04                 21.17      5.06    21.56
```

`S3` stands out in the `T cells` column, as built. `output_class` is `"abundance_score"` and
`comparability` is "arbitrary units; comparable across samples for one cell type, not across cell types". The
`qc` list reports 411 genes, 4 samples and all 111 signature genes found.

With RSEM output, [rsemverse](https://daikikumakura.github.io/rsemverse/) reads the table in the right layout:

```r
tpm <- rsemverse::rsem_read("data", "metadata.csv", species = "human", expr = "tpm")
res <- bulkdeconv::bulkdeconv(tpm, "tpm", "human", "quantiseq")   # 25 methods: see bulkdeconv_methods()
```

## Documentation

[Package website](https://daikikumakura.github.io/bulkdeconv/) · [Function reference](https://daikikumakura.github.io/bulkdeconv/reference/) · [Changelog](https://daikikumakura.github.io/bulkdeconv/news/) · [Releases](https://github.com/DaikiKumakura/bulkdeconv/releases) · [R-universe](https://daikikumakura.r-universe.dev/bulkdeconv) · [Usage guide](https://github.com/DaikiKumakura/bulkdeconv/blob/main/README-GUIDE.md)

## Validation and limitations

- **Validation:** Compared built-in methods with reference implementations and tested representative external methods and synthetic mixtures.
- **Results:** The first 15 built-in methods reproduced the recorded reference outputs in tests.
- **Limitations:** Method scores have different units and interpretations; the consensus score was checked only on 9 blood mixtures with flow cytometry (seven lineages; small and not tumor tissue), see `validation/consensus-pbmc/README.md`.
- **Reproduction:** [inst/COPYRIGHTS](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS).

## Citation

Kumakura, D. (2026). *bulkdeconv: Cell-Type Deconvolution of Tumor Bulk RNA-Seq Through One Function* [Software]. [DOI: 10.5281/zenodo.24210481](https://doi.org/10.5281/zenodo.24210481).

Citation metadata: [CITATION.cff](CITATION.cff).

Also cite the methods used. Sources are listed by `bulkdeconv_methods()` and in [inst/COPYRIGHTS](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS).

## License

This project is licensed under GPL-3.0-only. See [LICENSE.md](LICENSE.md) for the license terms.

Third-party code and data are documented in [inst/COPYRIGHTS](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS). Separately installed methods retain their own licenses.
