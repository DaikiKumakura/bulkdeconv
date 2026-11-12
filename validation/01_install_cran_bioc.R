# Usage: Rscript 01_install_cran_bioc.R <isolated-library-directory>
# Run the scripts in numeric order; each takes the same library directory.
lib <- normalizePath(commandArgs(TRUE)[1], winslash = "/", mustWork = FALSE)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(timeout = 1200, repos = c(CRAN = "https://cloud.r-project.org"))
cat("install target (isolated):", lib, "\n")
cat("main library untouched:", .libPaths()[2], "\n\n")

before <- rownames(installed.packages(lib.loc = lib))
install_one <- function(p, bioc) {
  t0 <- Sys.time()
  ok <- tryCatch({
    if (bioc) {
      BiocManager::install(p, lib = lib, update = FALSE, ask = FALSE, dependencies = TRUE, quiet = TRUE)
    } else {
      install.packages(p, lib = lib, dependencies = TRUE, quiet = TRUE)
    }
    requireNamespace(p, quietly = TRUE)
  }, error = function(e) { message("  error: ", conditionMessage(e)); FALSE })
  v <- if (isTRUE(ok)) as.character(packageVersion(p)) else "-"
  cat(sprintf("%-14s %-5s %5.0fs  version=%s\n", p, if (isTRUE(ok)) "OK" else "FAIL", as.numeric(difftime(Sys.time(), t0, units = "secs")), v))
}
for (p in c("tidyestimate", "imsig", "ComICS")) install_one(p, FALSE)
for (p in c("quantiseqr", "xCell2", "SpatialDecon")) install_one(p, TRUE)

after <- rownames(installed.packages(lib.loc = lib))
cat("\nnewly installed in the isolated library:", length(setdiff(after, before)), "packages\n")
cat(paste(sort(setdiff(after, before)), collapse = ", "), "\n")
