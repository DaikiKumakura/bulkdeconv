# From the package root: Rscript validation/06_xcell2_parity.R <xCell2 1.2.3 library>
ref_lib <- commandArgs(trailingOnly = TRUE)
stopifnot(length(ref_lib) == 1L, dir.exists(ref_lib))
pkgload::load_all(quiet = TRUE)

golden <- readRDS("tests/testthat/fixtures/native-golden.rds")
cases <- list(
  TMECompendium = list(X = golden$xcell$X, species = "human"),
  MouseRNAseqData = list(X = golden$immucellai_mouse$X, species = "mouse")
)
for (preset in names(cases)) {
  X <- cases[[preset]]$X
  species <- cases[[preset]]$species
  native_time <- system.time({
    native <- suppressWarnings(bulkdeconv(X, "tpm", species, "xcell2", preset = preset))
  })[["elapsed"]]
  original_time <- system.time({
    original <- callr::r(function(preset, X) {
      stopifnot(as.character(getNamespaceVersion("xCell2")) == "1.2.3")
      env <- new.env(parent = emptyenv())
      name <- paste0(preset, ".xCell2Ref")
      utils::data(list = name, package = "xCell2", envir = env)
      t(suppressWarnings(suppressMessages(
        xCell2::xCell2Analysis(X, get(name, env, inherits = FALSE), minSharedGenes = 0)
      )))
    }, args = list(preset = preset, X = X), libpath = c(ref_lib, .libPaths()))
  })[["elapsed"]]
  error <- max(abs(native$score - original))
  stopifnot(identical(dimnames(native$score), dimnames(original)), error < 1e-6)
  cat(sprintf("%s: %d genes x %d samples, %.3fs original -> %.3fs native, max abs diff %.3g\n",
              preset, nrow(X), ncol(X), original_time, native_time, error))
}
