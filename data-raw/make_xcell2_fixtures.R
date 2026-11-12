# Capture xCell2 1.2.3 output for each bundled preset with tied ranks and partial overlap.
.libPaths(c(commandArgs(trailingOnly = TRUE), .libPaths()))
stopifnot(as.character(getNamespaceVersion("xCell2")) == "1.2.3")

presets <- utils::read.csv("inst/registry/presets.csv", stringsAsFactors = FALSE)
fixtures <- setNames(vector("list", nrow(presets)), presets$preset)
for (preset in names(fixtures)) {
  name <- paste0(preset, ".xCell2Ref")
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "xCell2", envir = env)
  ref <- get(name, env, inherits = FALSE)
  genes <- xCell2::getGenesUsed(ref)
  genes <- genes[seq_along(genes) %% 5L != 0L]
  i <- seq_along(genes)
  X <- cbind(S3 = i %% 113L, S1 = (i * 19L) %% 127L, S2 = (i * 29L) %% 131L)
  rownames(X) <- genes
  expected <- suppressWarnings(suppressMessages(
    t(xCell2::xCell2Analysis(X, ref, minSharedGenes = 0))
  ))
  fixtures[[preset]] <- list(X = X, expected = expected)
}
saveRDS(fixtures["LM22"], "tests/testthat/fixtures/xcell2-LM22-golden.rds", compress = "xz")
fixtures$LM22 <- NULL
saveRDS(fixtures, "tests/testthat/fixtures/xcell2-golden.rds", compress = "xz")
