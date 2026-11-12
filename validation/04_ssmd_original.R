# Usage: Rscript 04_ssmd_original.R <isolated-library-directory>
# Run the scripts in numeric order; each takes the same library directory.
base <- dirname(normalizePath(commandArgs(TRUE)[1], winslash = "/", mustWork = FALSE))
lib <- file.path(base, "lib"); lib2 <- file.path(base, "lib_ssmd_original"); dir.create(lib2, showWarnings = FALSE)
.libPaths(c(lib2, lib, .libPaths()))
options(warn = 1, timeout = 600, repos = BiocManager::repositories())
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
cat("install xiaoyulu95/SSMD into:", lib2, "\n")
ok <- tryCatch({ remotes::install_github("xiaoyulu95/SSMD", lib = lib2, upgrade = "never", dependencies = TRUE, quiet = TRUE); TRUE }, error = function(e) { cat("install error:", conditionMessage(e), "\n"); FALSE })
cat("installed:", ok, " loaded from:", find.package("SSMD"), "\n")
suppressWarnings(suppressMessages(library(SSMD)))
cat("exports:", paste(sort(getNamespaceExports("SSMD")), collapse = ", "), "\n")
print(args(SSMD::SSMD))
load_obj <- function(name, pkg) { e <- new.env(); utils::data(list = name, package = pkg, envir = e); get(ls(e)[1], e) }
M <- tryCatch(load_obj("example_bulk", "SSMD"), error = function(e) NULL)
if (is.null(M)) { M <- load_obj("example_bulk", "SSMD") }
cat("example_bulk dim:", dim(M), "\n")
M <- M[!duplicated(rownames(M)), , drop = FALSE]
for (n in c(100, 50, 20, 10, 6, 3, 2)) {
  t0 <- Sys.time()
  r <- tryCatch(suppressMessages(SSMD::SSMD(M[, seq_len(n), drop = FALSE], "Inflammatory")), error = function(e) e)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "error")) cat(sprintf("n=%-3d error (%.0fs): %s\n", n, secs, substr(gsub("\\s+", " ", conditionMessage(r)), 1, 150)))
  else {
    p <- if (is.list(r)) r[[1]] else r
    cat(sprintf("n=%-3d ok (%.0fs) names=%s | first element class=%s dim=%s\n", n, secs, paste(names(r), collapse = ","), class(p)[1], paste(dim(p), collapse = "x")))
    if (n == 20) { saveRDS(r, file.path(base, "ssmd_n20.rds")); if (!is.null(dim(p))) { cat("   cell types:", paste(utils::head(rownames(p), 12), collapse = ", "), "\n"); cat("   columns:", paste(utils::head(colnames(p), 4), collapse = ", "), "\n") } }
  }
}
