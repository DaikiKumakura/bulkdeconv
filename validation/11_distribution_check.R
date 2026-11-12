# Run from the package root with the dependencies from 07_clean_install.R.
# Builds the source archive, installs it into a new library and runs every method from the
# installed package alone (no reference directory, no development tree).
stage <- tempfile("bulkdeconv-distribution-")
dir.create(stage)
archive <- pkgbuild::build(".", dest_path = stage, vignettes = FALSE, binary = FALSE)
files <- utils::untar(archive, list = TRUE)
excluded <- "^bulkdeconv/(docs(/|$)|AGENTS[.]md$|data-raw(/|$)|validation(/|$))"
required <- paste0("bulkdeconv/", c(
  "inst/COPYRIGHTS", "inst/extdata/abis/sigmatrixRNAseq.txt", "inst/extdata/seqimmucc/sig_matr_seqImmuCC.txt",
  "inst/extdata/base/compendium.rds", "inst/extdata/base/cell_correspondance.txt", "inst/extdata/xcell2/LM22.rds"))
stopifnot(!any(grepl(excluded, files)), all(required %in% files))
cat("Verified that the source archive holds all reference files and COPYRIGHTS.\n")

lib <- file.path(stage, "lib")
dir.create(lib)
install.packages(archive, lib = lib, repos = NULL, type = "source")
.libPaths(c(lib, .libPaths()))
library(bulkdeconv, lib.loc = lib)
native <- readRDS("tests/testthat/fixtures/native-golden.rds")
lm22 <- readRDS("tests/testthat/fixtures/xcell2-LM22-golden.rds")$LM22
cases <- list(
  abis = list(X = native$abis$X, expr = "tpm", species = "human"),
  seqimmucc = list(X = native$seqimmucc$X, expr = "count", species = "mouse"),
  base = list(X = native$base$X, expr = "tpm", species = "mouse", n_permutations = 5L),
  xcell2 = list(X = lm22$X, expr = "tpm", species = "human", preset = "LM22")
)
for (method in names(cases)) {
  args <- c(unname(cases[[method]][1L]), cases[[method]][-1L], list(method = method))
  result <- suppressWarnings(do.call(bulkdeconv, args))
  stopifnot(identical(rownames(result$score), colnames(cases[[method]]$X)),
            ncol(result$score) > 0L, all(is.finite(result$score)))
  if (method == "xcell2") stopifnot(max(abs(result$score - lm22$expected)) < 1e-6)
  cat(method, ": ran from the installed package without any extra reference.\n", sep = "")
}

Sys.setenv(BULKDECONV_TEST_LIB = paste(.libPaths(), collapse = ";"),
           BULKDECONV_TEST_SLOW = "1", NOT_CRAN = "true")
testthat::test_dir("tests/testthat", package = "bulkdeconv", load_package = "installed",
                   filter = "integration|xcell2", stop_on_failure = TRUE)
