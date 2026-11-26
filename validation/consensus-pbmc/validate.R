# Validate bulkdeconv_consensus() against flow cytometry (quanTIseq PBMC/PMN mixtures,
# GEO GSE107572: 9 healthy donors, TPM + flow fractions in the sample annotation).
# Run from the package root:  Rscript validation/consensus-pbmc/validate.R
suppressMessages(devtools::load_all(".", quiet = TRUE))

# ---- parameters ---------------------------------------------------------------
SEED <- 20261126
N_BOOT <- 2000L
METHODS <- c("quantiseq", "mcp_counter", "xcell", "consensus_tme", "abis",
             "spatialdecon_safetme", "danaher", "imsig", "bindea", "charoentong")
# flow annotation name -> consensus lineage
FLOW <- c("b cells" = "B_cell", "cd4+ t cells" = "T_CD4", "cd8+ t cells" = "T_CD8",
          "natural killer cells" = "NK", "monocytes" = "monocyte",
          "neutrophils" = "neutrophil", "myeloid dendritic cells" = "dendritic")
dir <- file.path("validation", "consensus-pbmc")

# ---- data ---------------------------------------------------------------------
tpm <- utils::read.delim(gzfile(file.path(dir, "data", "GSE107572_tpm_PBMC_RNAseq.txt.gz")),
                         check.names = FALSE, stringsAsFactors = FALSE)
cat("TPM:", nrow(tpm), "genes x", ncol(tpm) - 1L, "samples\n")
donor <- sub("^pbmc_([0-9]+).*$", "donor \\1", colnames(tpm)[-1])
stopifnot(!anyDuplicated(donor), !anyDuplicated(tpm$GENE))
colnames(tpm)[-1] <- donor

sm <- readLines(gzfile(file.path(dir, "data", "GSE107572_series_matrix.txt.gz")))
ch <- sm[grepl("^!Sample_characteristics_ch1", sm)]
ch <- lapply(strsplit(ch, "\t"), function(x) gsub('"', "", x[-1]))
key <- vapply(ch, function(x) sub(":.*$", "", x[1]), "")
flow <- do.call(cbind, lapply(ch[key %in% names(FLOW)], function(x) as.numeric(sub("^[^:]*: ", "", x))))
colnames(flow) <- key[key %in% names(FLOW)]
sid <- sub("^subject id: ", "", ch[[which(key == "subject id")]])
rownames(flow) <- sid
stopifnot(setequal(rownames(flow), donor), nrow(flow) == 9L)
flow <- flow[donor, , drop = FALSE]
colnames(flow) <- unname(FLOW[colnames(flow)])
cat("flow:", paste(colnames(flow), collapse = ", "), "\n")

# ---- run the methods ----------------------------------------------------------
res <- list()
for (m in METHODS) {
  r <- tryCatch(suppressWarnings(suppressMessages(bulkdeconv(tpm, expr = "tpm", species = "human", method = m))),
                error = function(e) { cat("  skipped", m, ":", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) res[[m]] <- r
}
cat("methods run:", paste(names(res), collapse = ", "), "\n")
cons <- bulkdeconv_consensus(unname(res))
print(cons)

# ---- compare with flow --------------------------------------------------------
lin_of <- bulkdeconv_lineages()
rank01 <- function(x) (rank(x, ties.method = "average") - 0.5) / length(x)
method_lineage <- function(r, lin) {          # one method's rank score for a lineage
  rows <- lin_of[lin_of$method == r$metadata$method & lin_of$lineage == lin, ]
  cols <- intersect(rows$cell_type, colnames(r$score))
  if (!length(cols)) return(NULL)
  rowMeans(sapply(cols, function(cn) rank01(r$score[donor, cn])))
}
boot_ci <- function(x, y, B = N_BOOT) {
  n <- length(x)
  b <- replicate(B, { i <- sample.int(n, n, TRUE); if (length(unique(x[i])) < 2L || length(unique(y[i])) < 2L) NA else stats::cor(x[i], y[i], method = "spearman") })
  stats::quantile(b, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}
set.seed(SEED)
rows <- list()
for (lin in intersect(colnames(cons$score), colnames(flow))) {
  y <- flow[, lin]
  one <- function(x, who) {
    ci <- boot_ci(x, y)
    data.frame(lineage = lin, source = who, rho = stats::cor(x, y, method = "spearman"),
               ci_low = ci[1], ci_high = ci[2], n = length(y))
  }
  rows[[length(rows) + 1L]] <- one(cons$score[donor, lin], "consensus")
  for (m in names(res)) {
    x <- method_lineage(res[[m]], lin)
    if (!is.null(x)) rows[[length(rows) + 1L]] <- one(x, m)
  }
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
tab$agreement <- cons$agreement[tab$lineage]
utils::write.csv(tab, file.path(dir, "results", "consensus-vs-flow.csv"), row.names = FALSE)

summ <- do.call(rbind, lapply(split(tab, tab$lineage), function(d) {
  ind <- d[d$source != "consensus", ]
  ind <- ind[!is.na(ind$rho), ]   # a method with a constant output has no correlation
  data.frame(lineage = d$lineage[1], n_methods = nrow(ind),
             consensus_rho = d$rho[d$source == "consensus"],
             median_method_rho = stats::median(ind$rho), best_method_rho = max(ind$rho),
             worst_method_rho = min(ind$rho), agreement = d$agreement[1])
}))
summ$consensus_minus_median <- summ$consensus_rho - summ$median_method_rho
utils::write.csv(summ, file.path(dir, "results", "consensus-vs-flow-summary.csv"), row.names = FALSE)
print(summ, digits = 2, row.names = FALSE)
cat("\nlineages where consensus rho >= median method rho:", sum(summ$consensus_minus_median >= 0), "of", nrow(summ), "\n")
cat("lineages where consensus rho >= best method rho:", sum(summ$consensus_rho >= summ$best_method_rho), "of", nrow(summ), "\n")
cat("Spearman(agreement, consensus rho) across lineages:", round(stats::cor(summ$agreement, summ$consensus_rho, method = "spearman"), 2), "\n")
