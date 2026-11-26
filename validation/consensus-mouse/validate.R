# Validate bulkdeconv_consensus() (mouse) against flow cytometry, using the Petitprez et al.
# 2020 data set as packaged in immunedeconv (dataset_petitprez: bulk RNA-seq of blood and
# spleen of 3 TC-1 tumour-bearing mice, 6 samples, with flow cytometry fractions).
# Source of the file: github.com/omnideconv/immunedeconv, data/dataset_petitprez.rda
# (not redistributed here). Run from the package root:
#   Rscript validation/consensus-mouse/validate.R
suppressMessages(devtools::load_all(".", quiet = TRUE))

# ---- parameters ---------------------------------------------------------------
SEED <- 20261126
N_BOOT <- 2000L
METHODS <- c("mmcp_counter", "immucellai_mouse", "seqimmucc", "base")
# flow column -> consensus lineage
FLOW <- c("B cell" = "B_cell", "T cell CD8+" = "T_CD8", "NK cell" = "NK",
          "Monocyte" = "monocyte", "Neutrophil" = "neutrophil")
dir <- file.path("validation", "consensus-mouse")

# ---- data ---------------------------------------------------------------------
e <- new.env()
load(file.path(dir, "data", "dataset_petitprez.rda"), envir = e)
expr <- e$dataset_petitprez$expr_mat
ref <- as.data.frame(e$dataset_petitprez$ref, check.names = FALSE)
cat("expression:", nrow(expr), "genes x", ncol(expr), "samples\n")
stopifnot(!anyDuplicated(rownames(expr)), setequal(ref[["Sample Name"]], colnames(expr)),
          all(names(FLOW) %in% colnames(ref)))
ref <- ref[match(colnames(expr), ref[["Sample Name"]]), ]
flow <- as.matrix(ref[, names(FLOW)]); colnames(flow) <- unname(FLOW)
rownames(flow) <- colnames(expr)
data <- data.frame(gene = rownames(expr), expr, check.names = FALSE)
samples <- colnames(expr)

# ---- run the methods ----------------------------------------------------------
res <- list()
for (m in METHODS) {
  r <- tryCatch(suppressWarnings(suppressMessages(bulkdeconv(data, expr = "tpm", species = "mouse", method = m))),
                error = function(e) { cat("  skipped", m, ":", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) res[[m]] <- r
}
cat("methods run:", paste(names(res), collapse = ", "), "\n")
cons <- bulkdeconv_consensus(unname(res))
print(cons)

# ---- compare with flow --------------------------------------------------------
lin_of <- bulkdeconv_lineages("mouse")
rank01 <- function(x) (rank(x, ties.method = "average") - 0.5) / length(x)
method_lineage <- function(r, lin) {
  rows <- lin_of[lin_of$method == r$metadata$method & lin_of$lineage == lin, ]
  cols <- intersect(rows$cell_type, colnames(r$score))
  if (!length(cols)) return(NULL)
  rowMeans(sapply(cols, function(cn) rank01(r$score[samples, cn])))
}
boot_ci <- function(x, y, B = N_BOOT) {
  n <- length(x)
  b <- replicate(B, { i <- sample.int(n, n, TRUE)
    if (length(unique(x[i])) < 2L || length(unique(y[i])) < 2L) NA else stats::cor(x[i], y[i], method = "spearman") })
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
  rows[[length(rows) + 1L]] <- one(cons$score[samples, lin], "consensus")
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
  ind <- ind[!is.na(ind$rho), ]
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
