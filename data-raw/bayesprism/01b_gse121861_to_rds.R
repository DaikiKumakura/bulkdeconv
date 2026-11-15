# Turn the recovered GSE121861 counts (01_gse121861_counts.py) into the input of
# 02_annotate.R: a dgCMatrix genes x cells with gene symbols as row names (rows
# that share a symbol are summed) and a cell table with `cell` and `model`.
#
#   Rscript data-raw/bayesprism/01b_gse121861_to_rds.R <work_dir>

suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
stopifnot("usage: 01b_gse121861_to_rds.R <work_dir>" = length(args) == 1)
dir <- args[[1]]

m <- as(Matrix::readMM(file.path(dir, "counts.mtx")), "CsparseMatrix")
genes <- utils::read.csv(file.path(dir, "genes.csv"), stringsAsFactors = FALSE)
cells <- utils::read.csv(file.path(dir, "cells.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(m) == nrow(genes), ncol(m) == nrow(cells))
cat("matrix:", nrow(m), "genes x", ncol(m), "cells; total UMI", sum(m), "\n")

symbol <- factor(genes$GeneSymbol, levels = unique(genes$GeneSymbol))
cat("gene symbols used by more than one row:", sum(duplicated(genes$GeneSymbol)), "\n")
out <- as(Matrix::fac2sparse(symbol) %*% m, "CsparseMatrix")
rownames(out) <- levels(symbol)
colnames(out) <- cells$cell
stopifnot(sum(out) == sum(m), all(Matrix::colSums(out) == Matrix::colSums(m)))
cat("after summing shared symbols:", nrow(out), "genes\n")

cells$model <- cells$Model
saveRDS(out, file.path(dir, "counts.rds"))
utils::write.csv(cells, file.path(dir, "cells_model.csv"), row.names = FALSE)
cat("written counts.rds and cells_model.csv to", dir, "\n")
