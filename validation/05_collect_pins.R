# Usage: Rscript 05_collect_pins.R <isolated-library-directory>
# Run the scripts in numeric order; each takes the same library directory.
base <- dirname(normalizePath(commandArgs(TRUE)[1], winslash = "/", mustWork = FALSE))
lib <- file.path(base, "lib")
pk <- c("quantiseqr", "MCPcounter", "xCell", "xCell2", "tidyestimate", "ConsensusTME", "ImmuCellAI2.0", "imsig", "SpatialDecon",
        "EPIC", "ImmuCellAI", "mMCPcounter", "ComICS", "ImmuCellAImouse", "SSMD", "nnls", "bcv", "GSVA", "limSolve", "e1071", "preprocessCore")
rows <- lapply(pk, function(p) {
  d <- tryCatch(packageDescription(p, lib.loc = c(lib, .libPaths())), error = function(e) NULL)
  if (is.null(d) || identical(d, NA)) return(data.frame(package = p, version = NA, source = "not installed", sha = "", license = ""))
  src <- if (!is.null(d$RemoteType)) paste0(d$RemoteType, ":", d$RemoteUsername, "/", d$RemoteRepo, if (!is.null(d$RemoteSubdir)) paste0("/", d$RemoteSubdir) else "")
         else if (!is.null(d$biocViews) || !is.null(d$git_url)) "bioconductor" else "CRAN"
  data.frame(package = p, version = d$Version, source = src, sha = substr(if (!is.null(d$RemoteSha)) d$RemoteSha else if (!is.null(d$git_last_commit)) d$git_last_commit else "", 1, 12), license = gsub("\\s+", " ", d$License))
})
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE, right = FALSE)
write.csv(tab, file.path(base, "pins.csv"), row.names = FALSE)

est <- readRDS(file.path(base, "smoke_estimates2.rds"))
cat("\n== output features per method (multi-sample run)\n")
for (n in grep(":multi$", names(est), value = TRUE)) {
  m <- est[[n]]; cat(sprintf("%-32s %2d features: %s\n", sub(":multi$", "", n), ncol(m), paste(colnames(m), collapse = ", ")))
}
sz <- function(p) { d <- file.path(lib, p); if (!dir.exists(d)) return(NA_real_); f <- list.files(d, recursive = TRUE, full.names = TRUE); sum(file.info(f)$size) / 1e6 }
cat("\n== installed size (MB) in the isolated library\n")
for (p in c("xCell2", "SpatialDecon", "quantiseqr", "macrophage", "GeomxTools", "ImmuCellAI2.0", "ImmuCellAI", "ImmuCellAImouse", "EPIC", "xCell", "mMCPcounter", "MCPcounter", "tidyestimate", "imsig", "ComICS", "ConsensusTME")) cat(sprintf("  %-16s %7.1f\n", p, sz(p)))
cat(sprintf("  TOTAL %.0f MB in %d packages\n", sum(file.info(list.files(lib, recursive = TRUE, full.names = TRUE))$size) / 1e6, length(list.files(lib))))
