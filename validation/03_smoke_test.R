# Usage: Rscript 03_smoke_test.R <isolated-library-directory>
# Run the scripts in numeric order; each takes the same library directory.
lib <- normalizePath(commandArgs(TRUE)[1], winslash = "/", mustWork = FALSE)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(warn = 1, stringsAsFactors = FALSE)
suppressWarnings(suppressMessages({
  library(quantiseqr); library(MCPcounter); library(xCell); library(xCell2); library(tidyestimate)
  library(ConsensusTME); library(imsig); library(SpatialDecon); library(EPIC); library(ImmuCellAI)
  library(mMCPcounter); library(ComICS); library(ImmuCellAImouse); library(SSMD)
}))
load_obj <- function(name, pkg) { e <- new.env(); utils::data(list = name, package = pkg, envir = e); get(ls(e)[1], e) }
racle <- new.env(); utils::data("dataset_racle", package = "quantiseqr", envir = racle)
H <- get("dataset_racle", racle)$expr_mat; H <- H[!duplicated(rownames(H)), , drop = FALSE]
M <- load_obj("example_bulk", "SSMD"); colnames(M) <- paste0("S", seq_len(ncol(M))); M <- M[!duplicated(rownames(M)), , drop = FALSE]
M6 <- M[, 1:6, drop = FALSE]
cm <- new.env(); utils::data("commons", package = "ComICS", envir = cm)

num <- function(x) {                                       # -> numeric matrix
  if (inherits(x, "data.frame")) x <- as.matrix(x[, vapply(x, is.numeric, logical(1)), drop = FALSE])
  x
}
# each adapter: function(X) -> sample x feature numeric matrix (X = genes x samples matrix, symbols as rownames)
ad <- list(
  quantiseq = list(sp = "human", f = function(X) num(run_quantiseq(X, "TIL10", FALSE, TRUE, TRUE)), pick = "rows_are_samples"),
  mcp_counter = list(sp = "human", f = function(X) t(MCPcounter.estimate(X, featuresType = "HUGO_symbols"))),
  xcell = list(sp = "human", f = function(X) t(xCellAnalysis(X, rnaseq = TRUE, parallel.sz = 1))),
  estimate = list(sp = "human", f = function(X) { d <- tidyestimate::filter_common_genes(as.data.frame(X), id = "hgnc_symbol", tidy = FALSE, tell_missing = FALSE, find_alias = TRUE); num(tidyestimate::estimate_score(d, is_affymetrix = FALSE)) }),
  consensus_tme = list(sp = "human", f = function(X) t(ConsensusTME::consensusTMEAnalysis(X, cancerType = "Unfiltered", statMethod = "ssgsea"))),
  immucellai2 = list(sp = "human", f = function(X) num(ImmuCellAI2.0::run_immucellai2(X, tissue.type = "tumor"))),
  imsig = list(sp = "human", f = function(X) num(imsig::imsig(as.data.frame(X), r = 0.6))),
  spatialdecon_safetme = list(sp = "human", f = function(X) t(SpatialDecon::spatialdecon(as.matrix(X), matrix(0.1, nrow(X), ncol(X), dimnames = dimnames(X)), X = SpatialDecon::safeTME, align_genes = TRUE)$beta)),
  epic = list(sp = "human", f = function(X) num(EPIC::EPIC(bulk = as.matrix(X))$cellFractions)),
  immucellai = list(sp = "human", f = function(X) { r <- ImmuCellAI::ImmuCellAI_new(as.data.frame(X), "rnaseq", 0, 0); num(r[[1]]) }),
  mmcp_counter = list(sp = "mouse", f = function(X) t(mMCPcounter::mMCPcounter.estimate(X, features = "Gene.Symbol", genomeVersion = "GCRm39"))),
  dcq = list(sp = "mouse", f = function(X) num(ComICS::dcq(get("immgen_dat", cm), as.data.frame(X), get("DCQ_mar", cm))$average)),
  immucellai_mouse = list(sp = "mouse", f = function(X) { r <- ImmuCellAImouse::ImmuCellAI_mouse(as.data.frame(X), "rnaseq", 0, 0); num(r[[1]]) }),
  ssmd = list(sp = "mouse", f = function(X) { r <- SSMD::SSMD(X, "Inflammatory"); num(r) })
)
refs_h <- c("TMECompendium", "PanCancer", "ImmuneCompendium", "BlueprintEncode", "LM22", "TabulaSapiensBlood")
refs_m <- c("MouseRNAseqData", "ImmGenData", "TabulaMurisBlood")
for (r in refs_h) local({ rr <- r; ad[[paste0("xcell2[", rr, "]")]] <<- list(sp = "human", f = function(X) t(xCell2::xCell2Analysis(as.matrix(X), load_obj(paste0(rr, ".xCell2Ref"), "xCell2"), minSharedGenes = 0.5))) })
for (r in refs_m) local({ rr <- r; ad[[paste0("xcell2[", rr, "]")]] <<- list(sp = "mouse", f = function(X) t(xCell2::xCell2Analysis(as.matrix(X), load_obj(paste0(rr, ".xCell2Ref"), "xCell2"), minSharedGenes = 0.5))) })

fix <- function(x) {                                       # -> numeric matrix (sample x feature) or error
  if (is.list(x) && !is.data.frame(x)) {
    cand <- Filter(function(z) is.matrix(z) || is.data.frame(z), x)
    if (!length(cand)) stop("non-matrix output: ", paste(names(x), collapse = ","))
    x <- cand[[1]]
  }
  x <- num(x)
  if (is.null(dim(x))) x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  x
}
try_run <- function(f, X) {
  warns <- character(); t0 <- Sys.time()
  out <- tryCatch(withCallingHandlers(suppressMessages(fix(f(X))), warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }), error = function(e) e)
  list(out = out, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")), warns = warns)
}
sink_null <- function(expr) { zz <- file(tempfile(), open = "wt"); sink(zz); on.exit({ sink(); close(zz) }); force(expr) }
rows <- list(); ests <- list()
note_of <- function(r) if (inherits(r$out, "error")) substr(gsub("\\s+", " ", conditionMessage(r$out)), 1, 140) else if (length(r$warns)) substr(gsub("\\s+", " ", r$warns[1]), 1, 90) else ""
for (id in names(ad)) {
  a <- ad[[id]]; Xfull <- if (a$sp == "human") H else if (id == "ssmd") M else M6
  X1 <- Xfull[, 1, drop = FALSE]
  Xdup <- cbind(Xfull[, 1], Xfull[, 1]); colnames(Xdup) <- c("A", "B"); rownames(Xdup) <- rownames(Xfull)
  Xfull <- as.matrix(Xfull)
  res <- list(multi = sink_null(try_run(a$f, Xfull)), one = sink_null(try_run(a$f, X1)), dup = sink_null(try_run(a$f, Xdup)))
  st <- vapply(res, function(r) if (inherits(r$out, "error")) "error" else "ok", "")
  for (sc in names(res)) {
    r <- res[[sc]]
    rows[[length(rows) + 1L]] <- data.frame(method = id, species = a$sp, scenario = sc, status = st[[sc]], seconds = round(r$secs, 1),
      out_dim = if (st[[sc]] == "ok") paste(dim(r$out), collapse = "x") else "", features = if (st[[sc]] == "ok") paste(utils::head(colnames(r$out), 3), collapse = "|") else "", note = note_of(r))
    if (st[[sc]] == "ok") ests[[paste(id, sc, sep = ":")]] <- r$out
  }
  # independence: sample 1 in the multi-sample run vs when run alone (or duplicated)
  indep <- NA_character_
  if (st[["multi"]] == "ok") for (alt in c("one", "dup")) if (is.na(indep) && st[[alt]] == "ok") {
    a1 <- res$multi$out[1, ]; b1 <- res[[alt]]$out[1, ]
    nm <- intersect(names(a1), names(b1)); d <- suppressWarnings(max(abs(a1[nm] - b1[nm]), na.rm = TRUE))
    indep <- sprintf("%s: max|diff|=%.3g", alt, d)
  }
  cat(sprintf("%-26s %-5s multi=%-5s one=%-5s dup=%-5s  independence[%s]  %s\n", id, a$sp, st[["multi"]], st[["one"]], st[["dup"]], ifelse(is.na(indep), "n/a", indep), note_of(res$multi)))
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(dirname(lib), "smoke_results2.csv"), row.names = FALSE)
saveRDS(ests, file.path(dirname(lib), "smoke_estimates2.rds"))

# ---- SSMD / imsig: minimum number of samples ---------------------------------------------
cat("\n--- SSMD minimum samples (mouse test matrix) ---\n")
for (n in c(100, 50, 30, 20, 10, 6, 3, 2)) {
  r <- sink_null(try_run(function(X) num(SSMD::SSMD(X, "Inflammatory")), M[, seq_len(n), drop = FALSE]))
  cat(sprintf("  n=%-3d %s %s\n", n, if (inherits(r$out, "error")) "error" else "ok   ", if (inherits(r$out, "error")) substr(gsub("\\s+", " ", conditionMessage(r$out)), 1, 110) else paste0("dim=", paste(dim(r$out), collapse = "x"), sprintf(" (%.0fs)", r$secs))))
}
cat("\n--- imsig minimum samples (human test matrix, 4 samples) ---\n")
for (n in c(4, 3, 2)) {
  r <- sink_null(try_run(function(X) num(imsig::imsig(as.data.frame(X), r = 0.6)), H[, seq_len(n), drop = FALSE]))
  cat(sprintf("  n=%-3d %s\n", n, if (inherits(r$out, "error")) paste("error:", substr(gsub("\\s+", " ", conditionMessage(r$out)), 1, 110)) else paste0("ok dim=", paste(dim(r$out), collapse = "x"))))
}
