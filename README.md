# bulkdeconv

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24210482.svg)](https://doi.org/10.5281/zenodo.24210482)
<!-- badges: end -->

**Cell-type deconvolution of human and mouse tumor bulk RNA-seq: 18
published methods through one function.**

```r
res <- bulkdeconv(tpm, expr = "tpm", species = "human", method = "quantiseq")
res$score   # samples x cell types
```

- **One interface for 19 methods.** quanTIseq, MCP-counter, xCell, xCell2,
  ESTIMATE, EPIC, TIMER and others for human; mMCP-counter, seq-ImmuCC, BASE,
  BayesPrism, DCQ and others for mouse. The call and the result look the same for all.
- **Installing the package is enough for 15 of them.** They are implemented
  inside bulkdeconv and reproduce the original output; their reference data
  come with the package. Only EPIC, BayesPrism, DCQ and SSMD need one extra
  install step.
- **The result says how to read it.** Every result states whether the numbers
  are fractions, abundance scores or enrichment scores, and what can be
  compared, with input and signature-coverage checks.
- **No registration-only methods.** CIBERSORT, CIBERSORTx, EcoTyper and ImmuCC
  are not included.

**Documentation:** <https://daikikumakura.github.io/bulkdeconv/> ·
[Releases](https://github.com/DaikiKumakura/bulkdeconv/releases)

## Workflow

```text
  gene symbol × sample table (TPM; counts for seq-ImmuCC and BayesPrism)
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

## Installation

```r
install.packages(c("remotes", "BiocManager"))
remotes::install_github("DaikiKumakura/bulkdeconv", repos = BiocManager::repositories())
```

This is all for the 15 built-in methods. For the three that run from their
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
- Every method takes TPM (`expr = "tpm"`) except `seqimmucc` and `bayesprism`,
  which take counts. A wrong scale is refused with a message that says so.
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
| `epic` | human | no: `bulkdeconv_install("epic")` | fraction | 1 |
| `mmcp_counter` | mouse | yes | abundance score | 1 |
| `immucellai_mouse` | mouse | yes | fraction | 2 |
| `seqimmucc` | mouse | yes (LLSR; counts) | fraction | 1 |
| `base` | mouse | yes | enrichment score | 2 |
| `dcq` | mouse | no: `bulkdeconv_install("dcq")` | relative quantity | 1 |
| `ssmd` | mouse | no: `bulkdeconv_install("ssmd")` | marker modules | 6 |
| `bayesprism` | mouse | no: `bulkdeconv_install("bayesprism")`; references for 23 tumor models | fraction | 1 |

`bulkdeconv_methods()` gives the full table. EPIC is free for non-commercial
use only and SSMD declares no licence; both are reported when run.

## Method options

Options are passed by name through `...`:

| Option | Methods |
|---|---|
| `preset` | `xcell2` (reference set, for example `"ImmGenData"` for mouse) |
| `cancer_type` | `consensus_tme`; required for `timer` |
| `tissue` | `immucellai2`, `ssmd` |
| `n_permutations`, `seed` | `base` |
| `model` (required), `n_cores`, `seed` | `bayesprism` |

```r
bulkdeconv(tpm, "tpm", "mouse", "xcell2", preset = "ImmGenData")
```

## BayesPrism with tumor-model references

`bayesprism` runs BayesPrism (Chu et al. 2022) with a reference built for the
mouse tumor model: summed UMI counts per cell state from public scRNA-seq of
that model (GEO GSE310560 and other GEO series, listed in the table). Give raw counts and the model:

```r
bulkdeconv_install("bayesprism")
res <- bulkdeconv(counts, "count", "mouse", "bayesprism", model = "MC38", n_cores = 4)
```

| `model` | Status | Cell states | Tumor cells | Malignant RMSE, independent | Malignant RMSE, own dataset | Source |
|---|---|---|---|---|---|---|
| `B16F10` | validated | 34 | 3141 | - | 0.007 | GSE310560 |
| `MC38` | validated | 29 | 1846 | 0.448 | 0.019 | GSE310560 |
| `CT26` | validated | 27 | 6554 | 0.41 | 0.026 | GSE310560 |
| `LLC` | validated | 33 | 644 | - | 0.039 | GSE310560 |
| `RENCA` | validated | 33 | 2528 | 0.051 | 0.028 | GSE310560 |
| `4T1` | validated | 31 | 1488 | 0.097 | 0.207 | GSE310560 |
| `KPCY-F` | candidate | 26 | 95 | - | 0.072 | GSE310560 |
| `E0771` | validated | 11 | 2701 | 0.06 | 0.013 | GSE266617 |
| `F420` | validated | 9 | 697 | 0.033 | 0.012 | GSE179681 |
| `MB49` | candidate | 5 | 8776 | 0.137 | 0.011 | GSE259337 |
| `K7M2` | candidate | 13 | 409 | 0.393 | 0.027 | GSE179681 |
| `MB49_C5` | candidate | 8 | 358 | - | 0.017 | GSE290112 |
| `ID8` | candidate | 10 | 179 | - | 0.023 | GSE297132 |
| `MOC1` | candidate | 9 | 1333 | - | 0.01 | GSE273507 |
| `SB28` | candidate | 8 | 1553 | - | 0.045 | GSE226292 |
| `MycCaP` | candidate | 7 | 32942 | - | 0.011 | GSE240382 |
| `KPC1` | candidate | 11 | 3146 | - | 0.004 | GSE337208 |
| `GL261_HRE_UnaG` | candidate | 9 | 4542 | - | 0.036 | GSE179077 |
| `TC1` | candidate | 7 | 3779 | - | 0.026 | GSE312357 |
| `CT2A` | candidate | 6 | 1487 | - | 0.05 | GSE309632 |
| `AB1` | candidate | 9 | 3790 | - | 0.104 | GSE153942 |
| `A20` | candidate | 6 | 1323 | - | 0.111 | GSE173117 |
| `B16OVA` | candidate | 10 | 2438 | - | 0.002 | GSE261553 |

RMSE is the error of the malignant fraction on pseudo-bulk samples made from
scRNA-seq cells: "independent" from another dataset of the same model,
"own dataset" from a held-out mouse (or half of the cells) of the reference
dataset. `validated` models passed an independent check (or, for B16F10 and
LLC, a held-out-mouse check in GSE310560); `candidate` models have no passing
independent check (most were checked only within their own dataset) and give a
warning. For MC38 and CT26 the second
independent dataset (GSE222002, GSE311795) disagreed (RMSE 0.45, 0.41), likely
in part from its cell annotation; an earlier check on GSE121861 gave 0.077 and
0.098. The references of the GSE310560 models hold tumor cells, T cells (CD4,
regulatory, CD8), monocytes, macrophages, dendritic cells and fibroblasts; the
other references hold the cell types found in their dataset (column
`cell_types` of `inst/registry/bayesprism_models.csv`). Cells missing from a
reference (endothelial cells, B and NK cells,
neutrophils) are counted as other types. The values are fractions of reads
over the genes used, not of cells. A run takes several minutes; `n_cores`
splits it across cores.

## Reading the result

- `res$metadata$output_class` says what the numbers are: fractions, abundance
  scores, enrichment scores or marker modules. Only fractions sum to a
  composition; scores compare one cell type across samples, not cell types
  with each other.
- Methods whose output depends on the other samples in the run (`xcell`,
  `xcell2`, `consensus_tme`, `imsig`, `timer`, `immucellai_mouse`, `ssmd`,
  `base`) warn when there are fewer than 10 samples. Use samples from one
  cohort for relative comparisons.
- A single sample is accepted only where the method can handle it. For xCell2,
  a one-sample run returns zero for every cell type, because the method
  subtracts the minimum score within the cohort.

## Validation

The 15 built-in methods are checked in the tests against output recorded from
the original implementations. On a clean installation from GitHub, all 15 ran
without any extra step, and EPIC, DCQ and SSMD ran after
`bulkdeconv_install()`. R CMD check passes on Linux, macOS and Windows in CI.
Sources, licences and the few deliberate differences from the originals are
listed in [`inst/COPYRIGHTS`](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS) and the registry
([`inst/registry/`](https://github.com/DaikiKumakura/bulkdeconv/tree/main/inst/registry)).

## Related tools

- [rsemverse](https://daikikumakura.github.io/rsemverse/) reads RSEM files
  into the gene symbol by sample table that bulkdeconv takes.
- [rsemflow](https://daikikumakura.github.io/rsemflow/) runs DESeq2, GSEA and
  GSVA from the same RSEM files.

## Citation

Cite the archived version: <https://doi.org/10.5281/zenodo.24210482>, and the
method you used: its source is listed in `bulkdeconv_methods()` and
[`inst/COPYRIGHTS`](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS).

## Licence

GPL-3. The built-in methods are ports of published code (GPL-2 or later,
GPL-3, MIT); sources are in [`inst/COPYRIGHTS`](https://github.com/DaikiKumakura/bulkdeconv/blob/main/inst/COPYRIGHTS). The three
separately installed methods keep their own licences, which the registry lists
and which are reported when a method runs.
