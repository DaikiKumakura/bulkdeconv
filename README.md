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

With RSEM output, [rsemverse](https://daikikumakura.github.io/rsemverse/)
reads the table in the right layout:

```r
tpm <- rsemverse::rsem_read("data", "metadata.csv", species = "human", expr = "tpm")

res <- bulkdeconv::bulkdeconv(tpm, "tpm", "human", "mcp_counter")
res$score                    # samples x MCP-counter abundance scores
res$metadata$output_class    # what the numbers are
res$metadata$comparability   # what may be compared
res$qc                       # input and signature-coverage checks
```

## Documentation

[Package website](https://daikikumakura.github.io/bulkdeconv/) · [Function reference](https://daikikumakura.github.io/bulkdeconv/reference/) · [Changelog](https://daikikumakura.github.io/bulkdeconv/news/) · [Releases](https://github.com/DaikiKumakura/bulkdeconv/releases) · [R-universe](https://daikikumakura.r-universe.dev/bulkdeconv) · [Usage guide](https://github.com/DaikiKumakura/bulkdeconv/blob/main/README-GUIDE.md)

## Validation and limitations

- **Validation:** Compared built-in methods with reference implementations and tested representative external methods and synthetic mixtures.
- **Results:** The first 15 built-in methods reproduced the recorded reference outputs in tests.
- **Limitations:** Method scores have different units and interpretations; the consensus score has not been validated against known mixtures.
- **Reproduction:** [inst/COPYRIGHTS](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS).

## Citation

Kumakura, D. (2026). *bulkdeconv: Cell-Type Deconvolution of Tumor Bulk RNA-Seq Through One Function* [Software]. [DOI: 10.5281/zenodo.24210481](https://doi.org/10.5281/zenodo.24210481).

Citation metadata: [CITATION.cff](CITATION.cff).

Also cite the methods used. Sources are listed by `bulkdeconv_methods()` and in [inst/COPYRIGHTS](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS).

## License

This project is licensed under GPL-3.0-only. See [LICENSE.md](LICENSE.md) for the license terms.

Third-party code and data are documented in [inst/COPYRIGHTS](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS). Separately installed methods retain their own licenses.
