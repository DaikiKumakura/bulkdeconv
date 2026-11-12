# Run from the package root, passing the extracted SpatialDecon 1.20.1 directory.
# Only the four upstream R files used for this call and logNormReg are needed.
upstream_dir <- commandArgs(trailingOnly = TRUE)[1L]
stopifnot(read.dcf(file.path(upstream_dir, "DESCRIPTION"), fields = "Version")[1L] == "1.20.1")
upstream <- new.env()
for (file in c("spatialdecon.R", "algorithm2.R", "deconLNR.R", "flagOutliers.R")) {
  sys.source(file.path(upstream_dir, "R", file), envir = upstream)
}
source("R/native.R")
g <- readRDS("tests/testthat/fixtures/native-golden.rds")$spatialdecon
signature <- readRDS("inst/extdata/spatialdecon/safeTME.rds")
background <- matrix(0.1, nrow(g$X), ncol(g$X), dimnames = dimnames(g$X))
expected <- t(upstream$spatialdecon(g$X, bg = background, X = signature)$beta)
actual <- bd_native$spatialdecon_fit(g$X, signature)
error <- max(abs(actual - expected))
cat("Maximum absolute difference from same-platform SpatialDecon:", error, "\n")
stopifnot(identical(dimnames(actual), dimnames(expected)), error < 1e-5)
