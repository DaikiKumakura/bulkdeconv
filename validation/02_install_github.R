# Usage: Rscript 02_install_github.R <isolated-library-directory>
# Run the scripts in numeric order; each takes the same library directory.
lib <- normalizePath(commandArgs(TRUE)[1], winslash = "/", mustWork = FALSE)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(timeout = 1200, repos = BiocManager::repositories())
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
cat("install target (isolated):", lib, "\n\n")
before <- rownames(installed.packages(lib.loc = lib))

targets <- list(
  list(name = "MCPcounter",       repo = "ebecht/MCPcounter",                subdir = "Source"),
  list(name = "xCell",            repo = "dviraran/xCell",                   subdir = NULL),
  list(name = "ConsensusTME",     repo = "cansysbio/ConsensusTME",           subdir = NULL),
  list(name = "ImmuCellAI2.0",    repo = "GuoBioinfoLab/ImmuCellAI-2.0",     subdir = "R-package/ImmuCellAI2.0"),
  list(name = "EPIC",             repo = "GfellerLab/EPIC",                  subdir = NULL),
  list(name = "ImmuCellAI",       repo = "lydiaMyr/ImmuCellAI",              subdir = NULL),
  list(name = "mMCPcounter",      repo = "cit-bioinfo/mMCP-counter",         subdir = NULL),
  list(name = "ImmuCellAImouse",  repo = "BioInfoCloud/ImmuCellAI_mouse",    subdir = NULL),
  list(name = "SSMD",             repo = "zy26/SSMD",                        subdir = NULL)
)
for (t in targets) {
  t0 <- Sys.time()
  ok <- tryCatch({
    args <- list(t$repo, lib = lib, upgrade = "never", dependencies = TRUE, build_vignettes = FALSE, quiet = TRUE)
    if (!is.null(t$subdir)) args$subdir <- t$subdir
    do.call(remotes::install_github, args)
    requireNamespace(t$name, quietly = TRUE)
  }, error = function(e) { message("  error (", t$name, "): ", gsub("\\s+", " ", substr(conditionMessage(e), 1, 400))); FALSE })
  v <- if (isTRUE(ok)) as.character(packageVersion(t$name)) else "-"
  cat(sprintf("%-16s %-5s %5.0fs  version=%s   [%s]\n", t$name, if (isTRUE(ok)) "OK" else "FAIL", as.numeric(difftime(Sys.time(), t0, units = "secs")), v, t$repo))
}
after <- rownames(installed.packages(lib.loc = lib))
cat("\nnewly installed in this step:", length(setdiff(after, before)), "packages\n")
cat(paste(sort(setdiff(after, before)), collapse = ", "), "\n")
