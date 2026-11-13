# bulkdeconv

Cell-type deconvolution of **tumor bulk RNA-seq** (human or mouse) through one
function.

[Documentation and function reference](https://daikikumakura.github.io/bulkdeconv/) · [Releases](https://github.com/DaikiKumakura/bulkdeconv/releases)

```r
res <- bulkdeconv(df, expr = "tpm", species = "human", method = "quantiseq")
res$score    # samples x cell types
```

## Input

`df` is a table with one gene symbol per row and one column per sample (count
or TPM): the **first column is the gene symbol**, or the symbols are the row
names, as in the table `rsemverse::rsem_read()` returns.

| gene   | S1   | S2   |
|--------|------|------|
| TP53   | 12.4 | 10.9 |
| ...    | ...  | ...  |

Gene symbols must be unique; sample names must be unique; values numeric,
non-negative, not log-transformed. A numeric matrix with gene symbols as row
names also works. Each sample must have at least one nonzero value. Each method
accepts one scale: all need `expr = "tpm"`
except `seqimmucc`, which needs `"count"`. A wrong scale is refused with a
message that says so.

```r
tpm <- rsemverse::rsem_read(
  data = "data", metadata = "metadata.csv", species = "human", expr = "tpm"
)
res <- bulkdeconv::bulkdeconv(tpm, "tpm", "human", "mcp_counter")
res$score                    # samples x MCP-counter abundance scores
res$metadata$output_class   # interpretation of the score units
res$qc                       # input and signature-coverage checks
```

## Methods

`bulkdeconv_methods()` lists them with species, output kind, licence and
whether they are installed.

| species | methods |
|---|---|
| human | quantiseq, mcp_counter, xcell, estimate, consensus_tme, immucellai2, imsig, spatialdecon_safetme, epic, abis, timer |
| mouse | mmcp_counter, dcq, immucellai_mouse, ssmd, seqimmucc, base |
| both  | xcell2 (9 reference presets) |

Methods that need registration (CIBERSORT, CIBERSORTx, EcoTyper, ImmuCC) are
not included. EPIC is free for non-commercial use only; SSMD has no declared
licence. Both are reported when run. Fifteen methods are implemented inside bulkdeconv, so no method package has to be installed for them: ports of the published algorithms that reproduce the original output (`quantiseq`, `mcp_counter`, `xcell`, `estimate`, `consensus_tme`, `immucellai2`, `imsig`, `spatialdecon_safetme`, `abis`, `timer`, `mmcp_counter`, `immucellai_mouse`, `seqimmucc` (LLSR only), `base`, `xcell2`). Everything they need (reference data, and the packages `sva`, `MASS` and `quadprog`) is installed together with bulkdeconv. Sources, licences and the few deliberate differences are in `inst/COPYRIGHTS` and the registry. `abis` uses blood-derived signatures. The other three (`epic`, `dcq`, `ssmd`) run from their original packages, each in its own R process.

## Install and run

```r
install.packages(c("remotes", "BiocManager"))
remotes::install_github("DaikiKumakura/bulkdeconv", repos = BiocManager::repositories())

# tpm is a gene-by-sample table with gene symbols in the first column.
res <- bulkdeconv::bulkdeconv(tpm, "tpm", "human", "quantiseq")
res$score
res$metadata[c("output_class", "comparability")]
res$qc
```

That is all for the 15 methods implemented in the package: installing bulkdeconv
installs their reference data and the packages they use. Only `epic`, `dcq` and
`ssmd` run from their original packages, which their licences do not allow to be
shipped here; install those once with `bulkdeconv::bulkdeconv_install("epic")`
(or `"dcq"`, `"ssmd"`).

The packages of these three methods go into a **private library**
(`tools::R_user_dir("bulkdeconv", "data")`; change with
`options(bulkdeconv.lib = "path")`), at the versions recorded in the registry
(GitHub packages at a fixed commit). Each method runs in its own R process, so
dependency conflicts between methods do not occur and your normal library is
untouched.

Use `bulkdeconv::bulkdeconv_methods()` to check species, input scale,
minimum sample count, output class and package dependencies. The `installed`
column reports R-package availability. All reference data are installed with
the package.
If a method's dependencies are missing, run `bulkdeconv_install()` for that ID
and check any installation warnings for the names of the missing packages.

## Method options

Passed by name through `...`: `preset` (xcell2), `cancer_type`
(consensus_tme; required for timer), `tissue` (immucellai2, ssmd),
`n_permutations` and `seed` (base).

```r
bulkdeconv(df, "tpm", "mouse", "xcell2", preset = "ImmGenData")
```

## What the result says about comparability

`res$metadata$output_class` and `$comparability` state what the numbers are
(fractions, enrichment scores, marker modules) and what may be compared.
Methods whose output depends on the other samples in the run
(xcell, xcell2, consensus_tme, imsig, timer, immucellai_mouse, ssmd, base) warn when there
are fewer than 10 samples. A single sample is accepted only where the method
can handle it (some are run with an internal copy that is dropped again).
For xCell2, a one-sample run returns zero for every cell type because the
method subtracts the minimum score within the cohort. Use multiple samples
from the same cohort for relative comparisons.

## Files

* `R/`, `tests/`: the package.
* `inst/registry/`: the method registry (`methods.csv`, `presets.csv`,
  `sources.csv`); field definitions in `inst/registry/README.md`.
* `data-raw/`: script that builds the registry from the smoke-test results.
* `validation/`: scripts and report of the 2026-11-11 installation and smoke test.

## Licence

GPL-3. The natively implemented methods are ports of published code (GPL-2 or
later, GPL-3, MIT). Reference and code sources are in `inst/COPYRIGHTS`.
The methods that are installed separately keep their
own licences, which the registry lists and which are reported when a method runs.
