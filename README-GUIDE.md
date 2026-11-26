# bulkdeconv guide

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

This is all for the 19 built-in methods. For the three that run from their
original packages, whose licences do not allow shipping them here, install
once:

```r
bulkdeconv::bulkdeconv_install(c("epic", "dcq", "ssmd"))
```

They go into a private library (`tools::R_user_dir("bulkdeconv", "data")`,
changeable with `options(bulkdeconv.lib = "path")`) at the versions recorded
in the registry, and each method runs in its own R process, so your normal
library is not touched and methods cannot conflict.

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

[Package website](https://daikikumakura.github.io/bulkdeconv/) · [Function reference](https://daikikumakura.github.io/bulkdeconv/reference/) · [Changelog](https://daikikumakura.github.io/bulkdeconv/news/) · [Releases](https://github.com/DaikiKumakura/bulkdeconv/releases) · [R-universe](https://daikikumakura.r-universe.dev/bulkdeconv)

## Workflow




```r
res <- bulkdeconv(tpm, expr = "tpm", species = "human", method = "quantiseq")
res$score   # samples x cell types
```

- **One interface for 25 methods.** quanTIseq, MCP-counter, xCell, xCell2,
  ESTIMATE, EPIC, TIMER, GEDIT, PUREE (tumor purity) and others for human; mMCP-counter, seq-ImmuCC, BASE,
  DCQ and others for mouse. The call and the result look the same for all.
- **Installing the package is enough for 19 of them.** They are implemented
  inside bulkdeconv and reproduce the original output; their reference data
  come with the package. EPIC, DCQ, SSMD, PUREE, mySORT and DualSimplex need one extra install step.
- **The result says how to read it.** Every result states whether the numbers
  are fractions, abundance scores or enrichment scores, and what can be
  compared, with input and signature-coverage checks.
- **No registration-only methods.** CIBERSORT, CIBERSORTx, EcoTyper and ImmuCC
  are not included.

**Documentation:** <https://daikikumakura.github.io/bulkdeconv/> ·
[Releases](https://github.com/DaikiKumakura/bulkdeconv/releases)

```text
  gene symbol × sample table (TPM; counts for seq-ImmuCC)
                 │
   bulkdeconv_methods()      which methods fit the species and the data?
                 │
   bulkdeconv(data, expr, species, method, ...)
                 │
   res$score      samples × cell types
   res$metadata   what the numbers are and what may be compared
   res$qc         input checks and signature coverage
```

| Function | What it does |
|---|---|
| `bulkdeconv()` | Runs one method and returns `score`, `metadata` and `qc` |
| `bulkdeconv_methods()` | Lists the methods with species, input scale, minimum samples, output kind, licence and whether they are ready to run |
| `bulkdeconv_install()` | Installs the packages of EPIC, DCQ or SSMD into a private library |


## Input




One gene symbol per row and one sample per column, with the symbols either in
the first column or as row names (the table `rsemverse::rsem_read()`
returns); a numeric matrix with row names also works.

| gene  | S1   | S2   |
|-------|------|------|
| TP53  | 12.4 | 10.9 |
| ...   | ...  | ...  |

- Gene symbols and sample names must be unique.
- Values must be numeric, non-negative and **not log-transformed**; each
  sample needs at least one nonzero value.
- Every method takes TPM (`expr = "tpm"`) except `seqimmucc`, which takes
  counts. A wrong scale is refused with a message that says so.
- Pass the values as they are: do not filter or normalize them first, since
  each method needs its signature genes.


## Methods




| Method | Species | Built in | Output | Min. samples |
|---|---|---|---|---|
| `quantiseq` | human | yes | fraction | 1 |
| `mcp_counter` | human | yes | abundance score | 1 |
| `xcell` | human | yes | enrichment score | 2 |
| `xcell2` | human, mouse | yes (9 reference presets) | enrichment score | 1 |
| `estimate` | human | yes | enrichment score | 1 |
| `consensus_tme` | human | yes | enrichment score | 1 |
| `immucellai2` | human | yes | fraction | 1 |
| `imsig` | human | yes | abundance score | 2 |
| `spatialdecon_safetme` | human | yes | abundance estimate | 1 |
| `abis` | human | yes (blood-derived signatures) | abundance estimate | 1 |
| `timer` | human | yes | abundance score | 2 |
| `gedit` | human | yes (2 reference presets) | fraction | 1 |
| `danaher` | human | yes | abundance score | 1 |
| `charoentong` | human | yes | enrichment score | 1 |
| `bindea` | human | yes | enrichment score | 1 |
| `puree` | human | no: `bulkdeconv_install("puree")` (downloads the model) | fraction (tumor purity) | 1 |
| `mysort` | human | no: `bulkdeconv_install("mysort")` | fraction (percent) | 1 |
| `dualsimplex` | human, mouse | no: `bulkdeconv_install("dualsimplex")` | latent component (share per sample) | 5 |
| `epic` | human | no: `bulkdeconv_install("epic")` | fraction | 1 |
| `mmcp_counter` | mouse | yes | abundance score | 1 |
| `immucellai_mouse` | mouse | yes | fraction | 2 |
| `seqimmucc` | mouse | yes (LLSR; counts) | fraction | 1 |
| `base` | mouse | yes | enrichment score | 2 |
| `dcq` | mouse | no: `bulkdeconv_install("dcq")` | relative quantity | 1 |
| `ssmd` | mouse | no: `bulkdeconv_install("ssmd")` | marker modules | 6 |

`bulkdeconv_methods()` gives the full table. EPIC and PUREE are free for non-commercial
use only, and SSMD and mySORT declare no licence; this is reported when run.
PUREE's model is not bundled: `bulkdeconv_install("puree")` downloads it
from a recorded commit and checks it.


## Method options




Options are passed by name through `...`:

| Option | Methods |
|---|---|
| `preset` | `xcell2` (reference set, for example `"ImmGenData"` for mouse); `gedit` (`"HPCA"`, default, or `"ImmunoStates"`) |
| `cancer_type` | `consensus_tme`; required for `timer` |
| `tissue` | `immucellai2`, `ssmd` |
| `n_permutations`, `seed` | `base` |

```r
bulkdeconv(tpm, "tpm", "mouse", "xcell2", preset = "ImmGenData")
```


## Reading the result




- `res$metadata$output_class` says what the numbers are: fractions, abundance
  scores, enrichment scores or marker modules. Only fractions sum to a
  composition; scores compare one cell type across samples, not cell types
  with each other.
- Methods whose output depends on the other samples in the run (`xcell`,
  `xcell2`, `consensus_tme`, `imsig`, `timer`, `charoentong`, `bindea`,
  `mysort`, `immucellai_mouse`, `ssmd`, `base`, `dualsimplex`) warn when there are fewer than 10 samples. Use samples from one
  cohort for relative comparisons.
- A single sample is accepted only where the method can handle it. For xCell2,
  a one-sample run returns zero for every cell type, because the method
  subtracts the minimum score within the cohort.


## Consensus across methods

Methods disagree, and their numbers have different units. `bulkdeconv_consensus()`
combines several results for the same samples into one score per lineage (B cell,
CD4 T, CD8 T, NK, monocyte, macrophage, neutrophil, dendritic, fibroblast,
endothelial, and tumor content from PUREE and the inverted ESTIMATE score):

```r
res  <- lapply(c("quantiseq", "mcp_counter", "xcell"), function(m)
  bulkdeconv(tpm, expr = "tpm", species = "human", method = m))
cons <- bulkdeconv_consensus(res)
cons$score       # mean percentile rank of each sample, per lineage
cons$agreement   # how well the methods agree, per lineage (Spearman)
```

The value is a rank within the run, so it orders the samples for one lineage. It is
not a fraction, it cannot be compared between lineages, and it changes when samples
are added. It needs at least 5 samples, and a lineage appears only when two methods
cover it. Low or negative agreement means the methods order the samples differently:
look at the single results before using the consensus. Which output columns enter
which lineage is in `bulkdeconv_lineages()`; human only. The consensus has not been
checked on a small set of samples with known composition: 9 blood mixtures with flow cytometry (GEO GSE107572), where the consensus correlated with flow at least as well as the median method in 6 of 7 lineages, with wide uncertainty and no tumor samples.

## Related tools




- [rsemverse](https://daikikumakura.github.io/rsemverse/) reads RSEM files
  into the gene symbol by sample table that bulkdeconv takes.
- [rsemflow](https://daikikumakura.github.io/rsemflow/) runs DESeq2, GSEA and
  GSVA from the same RSEM files.


<a id="validation"></a>

## Validation and limitations

The 15 built-in methods of the first release are checked in the tests against
output recorded from the original implementations; on a clean installation
from GitHub, all 15 ran without any extra step, and EPIC, DCQ and SSMD ran after
`bulkdeconv_install()`. GEDIT is checked by recovering the proportions of
noise-free mixtures of its reference profiles, and the Danaher and Charoentong
scores against their published definitions (mean log2 expression of marker
genes; ssGSEA). PUREE reproduces the purities of its original test data
exactly; mySORT matches the authors' example output for 21 of 24 samples (the
others differ where two of its model settings fit almost equally well). R CMD check passes on Linux, macOS and Windows in CI.
Sources, licences and the few deliberate differences from the originals are
listed in [`inst/COPYRIGHTS`](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS) and the registry
([`inst/registry/`](https://github.com/DaikiKumakura/bulkdeconv/tree/main/inst/registry)).

## Citation

Kumakura, D. (2026). *bulkdeconv: Cell-Type Deconvolution of Tumor Bulk RNA-Seq Through One Function* [Software]. [DOI: 10.5281/zenodo.24210481](https://doi.org/10.5281/zenodo.24210481).

Citation metadata: [CITATION.cff](CITATION.cff).

Also cite the methods used. Sources are listed by `bulkdeconv_methods()` and in [inst/COPYRIGHTS](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS).

## License

This project is licensed under GPL-3.0-only. See [LICENSE.md](LICENSE.md) for the license terms.

Third-party code and data are documented in [inst/COPYRIGHTS](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS). Separately installed methods retain their own licenses.
