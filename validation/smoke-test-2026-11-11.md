# Smoke test of the 16 fixed methods (2026-11-11)

Purpose: before any package code is written, install each method in an isolated
library and find out what it really needs: input unit, behaviour with one
sample, dependence on the other samples, install success, and run time. This is
not an accuracy evaluation.

## Environment

* Windows 11, R 4.5.3, Bioconductor 3.22, Rtools 45. R and the main R library
  were not modified: everything was installed into an isolated library.
* Install size of the isolated library: **764 MB in 43 packages**. 539 MB of it is
  `macrophage`, a data package that `quantiseqr` lists under Suggests and that
  was pulled in because `dependencies = TRUE` was used (it can be avoided).
* Sources: CRAN (tidyestimate, imsig, ComICS), Bioconductor (quantiseqr, xCell2,
  SpatialDecon), GitHub (the rest). Versions and commit SHAs: `inst/registry/sources.csv`.

## Test data (bundled with the packages, no download)

| | Source | Size | Unit |
| --- | --- | --- | --- |
| Human tumour | `quantiseqr::dataset_racle` (melanoma RNA-seq) | 32,467 genes x 4 samples | TPM, HGNC symbols |
| Mouse | `SSMD::example_bulk` | 24,647 unique genes x 100 samples (first 6 used; SSMD used 2-100) | not stated (median 12.1), MGI symbols |

## Scenarios

* **multi**: all samples of the test matrix (4 human, 6 mouse).
* **one**: a one-column matrix.
* **dup**: the first sample duplicated into two columns (a workaround for methods
  that drop a one-column matrix to a vector).
* **independence**: maximum absolute difference between the values of sample 1 in the
  multi-sample run and when run alone (or duplicated). 0 means a sample is scored
  independently of the others.

## Results

| method | species | multi | one sample | independence | notes |
| --- | --- | --- | --- | --- | --- |
| quantiseq | human | ok | ok | 0 | |
| mcp_counter | human | ok | ok | 0 | downloads its signature from GitHub at call time |
| xcell | human | ok | **error** | **0.566** | depends on the other samples |
| xcell2 (6 human refs) | human | ok | **error** | 0.0002-0.053 | duplicate works; default `minSharedGenes = 0.9` rejects TabulaSapiensBlood (0.75 shared) |
| estimate | human | ok | ok | 0 | RNA-seq: stromal, immune, combined only |
| consensus_tme | human | ok | ok | **0.254** | ssGSEA is cohort-relative; needs `cancerType` ("Unfiltered" works) |
| immucellai2 | human | ok | ok | 0.068 | 15-19 s |
| imsig | human | ok | **error** | n/a | works with 2 distinct samples, fails with identical ones |
| spatialdecon_safetme | human | ok | **error** | 0 (duplicate) | convergence warning; needs a background matrix |
| epic | human | ok | ok | 0 | non-commercial licence |
| immucellai (initial) | human | **error** | **error** | n/a | calls the removed GSVA API; unusable with GSVA 2.4.9 |
| mmcp_counter | mouse | ok | **error** | 0 (duplicate) | warning: no genes for Basophils |
| dcq | mouse | ok | **error** | 0 (duplicate) | 61-gene bundled reference; 207 outputs |
| immucellai_mouse | mouse | ok | **error** | **2.62** | strongly cohort-dependent |
| ssmd (original repo) | mouse | ok for n >= 6 | error | n/a | 3 and 2 samples fail; 90-112 s per run |
| xcell2 (3 mouse refs) | mouse | ok | **error** | 0.16-0.22 | shared genes 0.83-0.88, so `minSharedGenes` was lowered to 0.5 |

("one sample: error" means the one-column matrix failed; for xcell2, spatialdecon,
mmcp_counter and dcq the duplicate-column workaround succeeded, so an adapter can
support a single sample.)

## Findings that change the design

1. **Single-sample support is an adapter job.** Five methods fail only because
   R drops a one-column matrix; duplicating the column makes them work and, for
   spatialdecon, mmcp_counter and dcq, gives identical values.
2. **Cohort dependence must be recorded.** xCell, ConsensusTME, ImmuCellAI-mouse,
   imsig and SSMD score a sample relative to the others, so a single-sample
   result is not comparable with a cohort result. The registry has a
   `sample_independent` column and `result$metadata` should carry it.
3. **Install problems were real.** SSMD declares only `bcv` but also needs `nnls`;
   the `zy26/SSMD` fork lacks `R2_two_mat` and fails, while the original
   `xiaoyulu95/SSMD` works. ImmuCellAI (initial) needs the old GSVA API.
   Isolated installs, pinned commits and per-method availability are necessary.
4. **MCPcounter reaches the internet at run time** (default `genes` and `probesets`
   arguments read `raw.githubusercontent.com`): an adapter must pass a pinned copy.
5. **Some "built-in signatures" are small.** DCQ ships 61 genes. xCell2 presets
   need at least 90% of their genes in the input by default.
6. **ImmuCellAI-mouse is MIT**, not ambiguous: the DESCRIPTION says
   `MIT + file LICENSE` and the LICENSE file is the standard CRAN stub. An earlier
   statement that its licence was unclear was wrong.
7. **ImmuCC is excluded** under the registration rule: its algorithm is
   CIBERSORT's SVR and its script tells users to request CIBERSORT from
   cibersort.stanford.edu.

## Not done

* Accuracy was not evaluated (the human data has FACS truth for only 4 samples).
* Only Windows and 4-6 samples (SSMD 2-100) were tested; Linux and macOS are
  untested. No method was tested with raw counts.
* The units of mcp_counter, xcell, consensus_tme, spatialdecon, mmcp_counter, dcq
  and ssmd are not stated in their documentation; they ran with TPM, which is not
  proof that TPM is the intended input.

## Reproduce

```bash
Rscript validation/01_install_cran_bioc.R  <library-dir>
Rscript validation/02_install_github.R     <library-dir>
Rscript validation/03_smoke_test.R         <library-dir>   # writes smoke_results2.csv next to the library
Rscript validation/04_ssmd_original.R      <library-dir>   # installs the original SSMD in a second library
Rscript validation/05_collect_pins.R       <library-dir>
python  data-raw/build_registry.py
```

`02_install_github.R` installs the `zy26/SSMD` fork that turned out to be broken;
`04_ssmd_original.R` installs and tests the original repository.
