# Extract the inference-only fields from the pinned Bioconductor reference objects.
.libPaths(c(commandArgs(trailingOnly = TRUE), .libPaths()))
stopifnot(as.character(getNamespaceVersion("xCell2")) == "1.2.3")

presets <- utils::read.csv("inst/registry/presets.csv", stringsAsFactors = FALSE)
out <- file.path("inst", "extdata", "xcell2")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (preset in presets$preset[presets$method == "xcell2"]) {
  name <- paste0(preset, ".xCell2Ref")
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "xCell2", envir = env)
  ref <- get(name, env, inherits = FALSE)
  saveRDS(list(
    signatures = xCell2::getSignatures(ref),
    params = xCell2::getParams(ref),
    spill = xCell2::getSpillMat(ref),
    genes = xCell2::getGenesUsed(ref)
  ), file.path(out, paste0(preset, ".rds")), compress = "xz")
}
