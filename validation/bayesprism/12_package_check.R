# Check bulkdeconv(method = "bayesprism") end to end with the bundled MC38
# reference on pseudo-bulk samples from GSE121861 (cells not used to build the
# reference): agreement with the truth, reproducibility with the same seed,
# dependence on the other samples in the run, and one sample.
#
#   Rscript validation/bayesprism/12_package_check.R <work_dir> <private_lib> [n_cores]

suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
work <- args[[1]]
options(bulkdeconv.lib = args[[2]])
n_cores <- if (length(args) >= 3) as.integer(args[[3]]) else 10L
# Source the package code instead of load_all() so that imports needed only by
# other methods (sva) need not be installed; system.file() points at inst/.
env <- new.env()
env$system.file <- function(..., package = "base") {
  if (identical(package, "bulkdeconv")) normalizePath(file.path("inst", ...), "/", mustWork = FALSE)
  else base::system.file(..., package = package)
}
for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = env)
env$bd_pkg_version <- function(...) "dev"
bulkdeconv <- env$bulkdeconv

SEED <- 20261115
N_PB <- 6
N_CELLS_PB <- 2000

counts <- readRDS(file.path(work, "gse121861", "counts.rds"))
lab <- utils::read.csv(file.path(work, "gse121861", "annot", "labels.csv"), stringsAsFactors = FALSE)
lab <- lab[lab$model == "MC38" & lab$cell_type != "low_quality", ]
set.seed(SEED)
pick <- lapply(seq_len(N_PB), function(i) sample(seq_len(nrow(lab)), N_CELLS_PB, replace = TRUE))
X <- sapply(pick, function(ix) Matrix::rowSums(counts[, lab$cell[ix]]))
colnames(X) <- sprintf("pb%d", seq_len(N_PB))
truth <- sapply(pick, function(ix) mean(lab$cell_type[ix] == "malignant"))  # cell share, for orientation

run <- function(x, seed = 1L) {
  t0 <- Sys.time()
  r <- bulkdeconv(x, expr = "count", species = "mouse", method = "bayesprism", model = "MC38",
                  n_cores = n_cores, seed = seed)
  list(score = r$score, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")), res = r)
}
a <- run(X)
b <- run(X)
c2 <- run(X[, 1:2])
one <- run(X[, 1, drop = FALSE])

cat(sprintf("6 samples: %.0f s; 2 samples: %.0f s; 1 sample: %.0f s\n", a$secs, c2$secs, one$secs))
cat("output columns:", paste(colnames(a$score), collapse = ", "), "\n")
cat(sprintf("same seed, max difference: %.3g\n", max(abs(a$score - b$score))))
cat(sprintf("2 samples alone vs in 6, max difference: %.3g\n", max(abs(c2$score - a$score[1:2, ]))))
cat(sprintf("1 sample alone vs in 6, max difference: %.3g\n", max(abs(one$score - a$score[1, , drop = FALSE]))))
cat(sprintf("Tumor fraction (reads) vs malignant cell share: r = %.3f\n", cor(a$score[, "Tumor"], truth)))
print(round(cbind(Tumor = a$score[, "Tumor"], malignant_cell_share = truth), 3))
print(a$res$qc)
