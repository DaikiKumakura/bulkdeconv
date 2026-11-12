# "Install bulkdeconv and it works": installs the source archive into a NEW,
# EMPTY library exactly as the README says (one remotes::install_local() call),
# then runs every natively implemented method through bulkdeconv() without
# calling bulkdeconv_install() and without the development tree on the path.
#
# Run from the package root in a session that sees only the new library:
#   R_LIBS=<empty dir> R_LIBS_USER=<empty dir> Rscript validation/12_install_only.R <empty dir> [archive]
args <- commandArgs(trailingOnly = TRUE)
lib <- normalizePath(args[1L], winslash = "/")
archive <- if (length(args) > 1L) args[2L] else "bulkdeconv_0.1.0.tar.gz"
stopifnot(dir.exists(lib), !length(list.files(lib)), file.exists(archive))
.libPaths(c(lib, .Library), include.site = FALSE)
stopifnot(identical(.libPaths(), c(lib, normalizePath(.Library, winslash = "/"))))
options(timeout = 600, repos = c(CRAN = "https://cloud.r-project.org"), warn = 1,
        bulkdeconv.lib = file.path(lib, "no-private-library"))

t0 <- Sys.time()
install.packages(c("remotes", "BiocManager"), lib = lib)
remotes::install_local(archive, lib = lib, repos = BiocManager::repositories(), upgrade = "never")
cat("install minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")
cat("packages in the new library:", length(list.files(lib)), "\n")

library(bulkdeconv, lib.loc = lib)
reg <- bulkdeconv_methods()
native <- reg$method[reg$engine == "native"]
stopifnot(length(native) == 15L, all(reg$installed[reg$engine == "native"]))

g <- readRDS("tests/testthat/fixtures/native-golden.rds")
x2 <- readRDS("tests/testthat/fixtures/xcell2-golden.rds")
lm22 <- readRDS("tests/testthat/fixtures/xcell2-LM22-golden.rds")$LM22
cases <- list(
  quantiseq = list(g$quantiseq$X, "tpm", "human"), mcp_counter = list(g$mcp_counter$X, "tpm", "human"),
  xcell = list(g$xcell$X, "tpm", "human"), estimate = list(g$estimate$X, "tpm", "human"),
  consensus_tme = list(g$consensus_tme$X, "tpm", "human"), immucellai2 = list(g$immucellai2$X, "tpm", "human"),
  imsig = list(g$imsig$X, "tpm", "human"), spatialdecon_safetme = list(g$spatialdecon$X, "tpm", "human"),
  abis = list(g$abis$X, "tpm", "human"), timer = list(g$timer$X, "tpm", "human", cancer_type = "LUAD"),
  mmcp_counter = list(g$mmcp_counter$X, "tpm", "mouse"), immucellai_mouse = list(g$immucellai_mouse$X, "tpm", "mouse"),
  seqimmucc = list(g$seqimmucc$X, "count", "mouse"), base = list(g$base$X, "tpm", "mouse", n_permutations = 5L),
  xcell2 = list(x2$TMECompendium$X, "tpm", "human"),
  `xcell2 (LM22)` = list(lm22$X, "tpm", "human", preset = "LM22"))
stopifnot(setequal(sub(" .*", "", names(cases)), native))
for (nm in names(cases)) {
  a <- cases[[nm]]
  names(a)[1:3] <- c("data", "expr", "species")
  res <- suppressWarnings(do.call(bulkdeconv, c(a, list(method = sub(" .*", "", nm)))))
  stopifnot(identical(rownames(res$score), colnames(a$data)), ncol(res$score) > 0L, all(is.finite(res$score)))
  cat(sprintf("%-22s ok  %s  %.1f s\n", nm, paste(dim(res$score), collapse = "x"), res$metadata$seconds))
}
cat("All 15 native methods ran after installing bulkdeconv alone.\n")
