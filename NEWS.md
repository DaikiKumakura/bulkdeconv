# bulkdeconv (development)

* `bulkdeconv_consensus()` combines results of several methods into one rank-based
  score per lineage (10 immune and stromal lineages, plus tumor content), with the
  agreement between the methods; `bulkdeconv_lineages()` shows the mapping from
  method outputs to lineages (human). The consensus orders samples within a lineage
  only and needs at least 5 samples.
* DualSimplex: the paper's year is 2026 (Genome Biology 27:25), not 2025.

# bulkdeconv 0.3.0

* `dualsimplex`: reference-free deconvolution with DualSimplex (Kleverov et al.
  2026). It needs no signature or single-cell reference: it factorizes the
  expression matrix of the run into `n_components` components (required, 2 to
  20) and returns each sample's share, summing to 1 (`output_class =
  latent_component`). The components are not cell types; assign them
  afterwards, for example with marker genes. Needs more samples than
  components; needs `bulkdeconv_install("dualsimplex")` (compiles C++ code).

# bulkdeconv 0.2.0

* Four signature-based human methods, implemented in the package with their
  data included: `gedit` (GEDIT, Nadel et al. 2021; presets `"HPCA"` with
  stromal cell types and `"ImmunoStates"`), `danaher` (marker-gene scores of
  14 immune populations, Danaher et al. 2017), `charoentong` (ssGSEA on the
  28 immune cell gene sets of Charoentong et al. 2017) and `bindea` (ssGSEA on
  the 24 immune cell gene sets of Bindea et al. 2013).
* `puree`: tumor purity with PUREE (Revkov et al. 2023), computed in the
  package; the model (CC BY-NC-SA 4.0) is downloaded by
  `bulkdeconv_install("puree")` from a recorded commit and checked.
* `mysort`: mySORT (Chen et al. 2018) with its signature matrix; needs
  `bulkdeconv_install("mysort")` (e1071, preprocessCore).
* BayesPrism was removed.

# bulkdeconv 0.1.0

First version.

* `bulkdeconv(data, expr, species, method, ...)` runs one method on a
  gene-by-sample table (first column gene symbols, other columns samples), the
  layout `rsemverse::rsem_read()` returns.
* 18 methods that need no registration: human (quantiseq, mcp_counter, xcell,
  estimate, consensus_tme, immucellai2, imsig, spatialdecon_safetme, epic,
  abis, timer), mouse (mmcp_counter, dcq, immucellai_mouse, ssmd, seqimmucc,
  base), both (xcell2). seqimmucc takes counts (LLSR only); the other 17 take TPM.
* 15 methods are implemented in the package (verified against the output of the
  original code): quantiseq, mcp_counter, xcell, estimate, consensus_tme,
  immucellai2, imsig, spatialdecon_safetme, abis, timer, mmcp_counter,
  immucellai_mouse, seqimmucc, base and xcell2. xCell2 inference includes the
  nine pre-trained presets. EPIC, dcq and ssmd use their original packages.
* immucellai2 takes 1.7 s instead of 18 s (the original reads a 41 MB text file
  on every call). estimate, immucellai2, spatialdecon_safetme and
  immucellai_mouse no longer need ggplot2, dplyr, GeomxTools, SeuratObject, GSVA
  or logNormReg.
* base is about 3x faster than the original code (51 s vs 175 s for 6 samples
  at 100 permutations).
* Licence: GPL-3.
* Each method runs in a separate R process from a private library at a pinned
  version (`bulkdeconv_install()`).
* Input validation, a registry (`bulkdeconv_methods()`), single-sample
  handling, and a standard result object.
* Installing bulkdeconv is enough for the 15 methods implemented in the package:
  their reference data ship with it and the packages they use (quadprog, MASS,
  sva) are installed as dependencies. Only epic, dcq and ssmd need
  `bulkdeconv_install()`.
