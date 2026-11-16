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
