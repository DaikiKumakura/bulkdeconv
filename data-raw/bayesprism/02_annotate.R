# Coarse cell-type annotation of a mouse tumor scRNA-seq count matrix by fixed
# marker rules, applied per model. The same rules are used for every dataset,
# so that references built from different datasets share one label set and can
# be validated against each other without transferring labels between them.
#
#   Rscript data-raw/bayesprism/02_annotate.R <counts.rds> <cells.csv> <out_dir>
#
# <counts.rds>: dgCMatrix genes x cells (gene symbols as row names, cell IDs as
# column names). <cells.csv>: a table with columns `cell` and `model`.
# Output: <out_dir>/labels.csv (cell, model, cluster, cell_type) and
# <out_dir>/clusters.csv (per cluster: size, marker detection rates, the call),
# which is the evidence to review by eye.

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(scuttle)
  library(scran)
  library(scater)
  library(bluster)
})

SEED <- 20261115
N_HVG <- 2000
N_PC <- 30
K_SNN <- 20
MIN_CELLS_MODEL <- 200
# A cluster belongs to a compartment when the mean detection rate of its
# markers (share of cells with a nonzero count) reaches this value.
DETECTION <- 0.4

PANELS <- list(
  immune = "Ptprc",
  # Endothelial cells and fibroblasts are called by genes that tumor cells of
  # mesenchymal-like lines rarely express: CT26, EMT6 and LL2 cells express
  # collagens (Col1a2, Col3a1) and Dcn, but not Lum and Pdgfra (GSE121861).
  endothelial = c("Pecam1", "Cdh5"),
  fibroblast = c("Lum", "Pdgfra"),
  T = c("Cd3e", "Cd3d", "Cd3g"),
  NK = c("Nkg7", "Klrb1c", "Ncr1"),
  B = c("Cd79a", "Ms4a1", "Cd19"),
  plasma = c("Jchain", "Mzb1", "Xbp1"),
  monocyte_macrophage = c("Lyz2", "C1qa", "Adgre1", "Csf1r"),
  neutrophil = c("S100a8", "S100a9", "Retnlg"),
  DC = c("Flt3", "Xcr1", "Cd209a", "Siglech"),
  mast_basophil = c("Cpa3", "Mcpt8", "Cma1")
)
IMMUNE_TYPES <- c("T", "NK", "B", "plasma", "monocyte_macrophage", "neutrophil", "DC", "mast_basophil")

args <- commandArgs(trailingOnly = TRUE)
stopifnot("usage: 02_annotate.R <counts.rds> <cells.csv> <out_dir>" = length(args) == 3)
counts <- readRDS(args[[1]])
cells <- utils::read.csv(args[[2]], stringsAsFactors = FALSE)
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(inherits(counts, "dgCMatrix"), identical(colnames(counts), cells$cell))
missing_markers <- setdiff(unlist(PANELS), rownames(counts))
cat("markers not in the data:", if (length(missing_markers)) paste(missing_markers, collapse = ", ") else "none", "\n")

detect_rate <- function(m, genes) {
  genes <- intersect(genes, rownames(m))
  if (!length(genes)) return(NA_real_)
  mean(Matrix::rowMeans(m[genes, , drop = FALSE] > 0))
}

# A non-immune cluster in which almost nothing is detected (Vim, expressed by
# nearly every tumor and stromal cell here, in under 30% of cells, and a median
# UMI below half of the model's median) is treated as low quality.
call_cluster <- function(rates, vim, umi, model_umi) {
  if (rates[["immune"]] < 0.1 && vim < 0.3 && umi < model_umi / 2) return("low_quality")
  if (rates[["immune"]] >= DETECTION) {
    imm <- rates[IMMUNE_TYPES]
    if (all(is.na(imm)) || max(imm, na.rm = TRUE) < DETECTION / 2) return("immune_unassigned")
    return(names(which.max(imm)))
  }
  if (rates[["endothelial"]] >= DETECTION) return("endothelial")
  if (rates[["fibroblast"]] >= DETECTION) return("fibroblast")
  "malignant"
}

labels <- list()
summary <- list()
for (model in sort(unique(cells$model))) {
  idx <- which(cells$model == model)
  if (length(idx) < MIN_CELLS_MODEL) {
    cat(sprintf("%s: %d cells, below %d; skipped\n", model, length(idx), MIN_CELLS_MODEL))
    next
  }
  m <- counts[, idx]
  m <- m[Matrix::rowSums(m) > 0, ]
  sce <- SingleCellExperiment(list(counts = m))
  set.seed(SEED)
  sce <- logNormCounts(sce)
  hvg <- getTopHVGs(modelGeneVar(sce), n = N_HVG)
  sce <- runPCA(sce, subset_row = hvg, ncomponents = N_PC, BSPARAM = BiocSingular::IrlbaParam())
  cl <- clusterRows(reducedDim(sce, "PCA"), NNGraphParam(k = K_SNN, cluster.fun = "louvain"))
  cl <- as.integer(cl)
  model_umi <- stats::median(Matrix::colSums(m))
  rows <- lapply(sort(unique(cl)), function(k) {
    sub <- m[, cl == k, drop = FALSE]
    rates <- vapply(PANELS, function(g) detect_rate(sub, g), numeric(1))
    data.frame(model = model, cluster = k, n = ncol(sub), t(round(rates, 2)),
               Cd8 = round(detect_rate(sub, c("Cd8a", "Cd8b1")), 2),
               Cd4 = round(detect_rate(sub, "Cd4"), 2),
               Mki67 = round(detect_rate(sub, "Mki67"), 2),
               Vim = round(detect_rate(sub, "Vim"), 2),
               median_umi = stats::median(Matrix::colSums(sub)),
               call = call_cluster(rates, detect_rate(sub, "Vim"), stats::median(Matrix::colSums(sub)), model_umi),
               check.names = FALSE)
  })
  tab <- do.call(rbind, rows)
  # Split T cells into CD8 and CD4 by the cluster's detection rates.
  is_t <- tab$call == "T"
  tab$call[is_t] <- ifelse(tab$Cd8[is_t] >= tab$Cd4[is_t], "T_CD8", "T_CD4")
  summary[[model]] <- tab
  labels[[model]] <- data.frame(cell = colnames(m), model = model, cluster = cl,
                                cell_type = tab$call[match(cl, tab$cluster)])
  cat(sprintf("%s: %d cells, %d clusters; cells per type: %s\n", model, ncol(m), nrow(tab),
              paste(sprintf("%s %d", names(table(labels[[model]]$cell_type)), table(labels[[model]]$cell_type)), collapse = ", ")))
}
labels <- do.call(rbind, labels)
summary <- do.call(rbind, summary)
utils::write.csv(labels, file.path(out_dir, "labels.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "clusters.csv"), row.names = FALSE)
cat("written:", nrow(labels), "cells,", nrow(summary), "clusters to", out_dir, "\n")
