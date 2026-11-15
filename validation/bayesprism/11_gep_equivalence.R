# Does a GEP reference (counts summed per cell state) give the same BayesPrism
# result as the per-cell count matrix it was summed from? new.prism() collapses
# a count matrix by cell state before normalising (R/new_prism.R), so the two
# should agree; this checks it numerically on one model.
#
#   Rscript validation/bayesprism/11_gep_equivalence.R <work_dir> [model] [n_cores]

suppressPackageStartupMessages({
  library(Matrix)
  library(BayesPrism)
})

SEED <- 20261115
N_PB <- 3
N_CELLS_PB <- 2000
MIN_CELLS_STATE <- 20

args <- commandArgs(trailingOnly = TRUE)
work <- args[[1]]
model <- if (length(args) >= 2) args[[2]] else "MC38"
n_cores <- if (length(args) >= 3) as.integer(args[[3]]) else 10L

m <- readRDS(file.path(work, "gse310560", model, "counts.rds"))
meta <- utils::read.csv(file.path(work, "gse310560", model, "meta.csv"), stringsAsFactors = FALSE)
meta$cell_state <- paste(meta$cell_type, meta$cell_state, sep = "|")
mice <- sort(unique(meta$mouse))

ref_idx <- which(meta$mouse == mice[1])
n_state <- table(meta$cell_state[ref_idx])
ref_idx <- ref_idx[meta$cell_state[ref_idx] %in% names(n_state)[n_state >= MIN_CELLS_STATE]]
test_idx <- which(meta$mouse == mice[2])

set.seed(SEED)
mixture <- t(sapply(seq_len(N_PB), function(i) Matrix::rowSums(m[, sample(test_idx, N_CELLS_PB, replace = TRUE)])))
rownames(mixture) <- sprintf("pb%02d", seq_len(N_PB))
genes <- rownames(m)[Matrix::rowSums(m[, ref_idx]) > 0 & colSums(mixture) > 0]
genes <- genes[!grepl("^(Rp[sl]|Mrp[sl]|mt-)", genes)]

ref_cells <- Matrix::t(m[genes, ref_idx])                       # cells x genes (dgCMatrix)
state <- factor(meta$cell_state[ref_idx])
gep <- as.matrix(Matrix::fac2sparse(state) %*% ref_cells)      # states x genes
type_cells <- meta$cell_type[ref_idx]
type_state <- tapply(type_cells, state, unique)[rownames(gep)]

run <- function(reference, input.type, types, states) {
  t0 <- Sys.time()
  utils::capture.output({
    p <- new.prism(reference = reference, mixture = mixture[, genes], input.type = input.type,
                   cell.type.labels = types, cell.state.labels = states, key = "Tumor",
                   outlier.cut = 0.01, outlier.fraction = 0.1)
    bp <- run.prism(p, n.cores = n_cores)
  })
  list(theta = get.fraction(bp, "final", "type"), seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

a <- run(ref_cells, "count.matrix", type_cells, as.character(state))
b <- run(gep, "GEP", unname(type_state), rownames(gep))
types <- sort(colnames(a$theta))
d <- abs(a$theta[, types] - b$theta[, types])
cat(sprintf("%s: %d reference cells in %d states, %d genes, %d pseudo-bulk samples\n",
            model, length(ref_idx), nrow(gep), length(genes), N_PB))
cat(sprintf("count.matrix input: %.0f s; GEP input: %.0f s\n", a$seconds, b$seconds))
cat(sprintf("largest absolute difference in any cell-type fraction: %.3g\n", max(d)))
print(round(rbind(count_matrix = colMeans(a$theta[, types]), GEP = colMeans(b$theta[, types])), 4))
cat(sprintf("GEP size in memory: %.1f MB (%d x %d doubles)\n", object.size(gep) / 1e6, nrow(gep), ncol(gep)))
