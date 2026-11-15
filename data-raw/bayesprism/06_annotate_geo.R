# Read the GEO samples of one model (sources.csv), run QC, cluster, and label
# clusters by marker detection. Clusters that match no normal cell type are
# labelled Tumor; a per-model override handles models whose malignant cells
# look like a normal type (A20: B-cell lymphoma). Labels use the same coarse
# vocabulary as GSE310560 where possible.
#
#   Rscript data-raw/bayesprism/06_annotate_geo.R <geo_dir> <out_dir> <model>
#
# Writes <out_dir>/<model>/counts.rds (genes x cells, integer dgCMatrix),
# meta.csv (cell, sample, role, cluster, cell_type, cell_state) and
# clusters.csv (per cluster: size, label, marker detection).

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
})

SEED <- 20261115
MIN_GENES <- 200
MIN_UMI <- 500
MAX_MITO <- 0.20
DETECTION <- 0.30   # share of cells in a cluster expressing a marker
N_HVG <- 2000
N_PC <- 30
K_NN <- 20

# Panels: a cluster gets the label whose markers have the highest mean
# detection, if that mean is >= DETECTION. Immune labels need Ptprc.
PANELS <- list(
  CD8_T_cells = c("Cd3e", "Cd8a", "Cd8b1"), CD4_Treg = c("Cd3e", "Cd4", "Foxp3"),
  CD4_Tconv = c("Cd3e", "Cd4", "Cd40lg"), gdT_other_T_cells = c("Cd3e", "Trdc", "Tcrg-C1"),
  NK = c("Nkg7", "Klrb1c", "Ncr1"), B_cells = c("Cd79a", "Ms4a1", "Cd19"), Plasma = c("Jchain", "Mzb1", "Xbp1"),
  Neutrophils = c("S100a8", "S100a9", "Retnlg"), Mast_Basophils = c("Cpa3", "Cyp11a1", "Ms4a2"),
  pDC = c("Siglech", "Bst2", "Ccr9"), cDC = c("Flt3", "Cd209a", "Xcr1", "Ccr7"),
  Monocytes = c("Ly6c2", "Plac8", "Ccr2"), Macrophages = c("C1qa", "C1qb", "Adgre1"),
  Microglia = c("P2ry12", "Tmem119", "Sall1"),
  Endothelial = c("Pecam1", "Cdh5", "Kdr"), Fibroblasts = c("Lum", "Pdgfra", "Dcn"),
  Pericytes = c("Rgs5", "Kcnj8", "Higd1b"), Astrocytes = c("Aqp4", "Aldh1l1", "Gja1"),
  Oligodendrocytes = c("Mog", "Mbp", "Plp1"), Neurons = c("Snap25", "Syt1", "Rbfox3"),
  Erythroid = c("Hba-a1", "Hbb-bs", "Alas2"),
  Lung_epithelial = c("Sftpc", "Sftpb", "Scgb1a1"), Osteoblasts = c("Bglap", "Sp7", "Ibsp")
)
IMMUNE <- c("CD8_T_cells", "CD4_Treg", "CD4_Tconv", "gdT_other_T_cells", "NK", "B_cells", "Plasma", "Neutrophils",
            "Mast_Basophils", "pDC", "cDC", "Monocytes", "Macrophages", "Microglia")
DROP <- c("Erythroid", "Low_quality")
TUMOR_AS <- list(A20 = c("B_cells"))  # malignant cells of these models carry a normal label
# Tissue-resident panels are used only where that tissue is sampled (tumor cells
# can express single markers such as Gja1 or Sp7, so the panels are not global).
BRAIN <- c("Microglia", "Astrocytes", "Oligodendrocytes", "Neurons")
TISSUE_PANELS <- list(SB28 = BRAIN, GL261_HRE_UnaG = BRAIN, CT2A = BRAIN,
                      F420 = "Lung_epithelial", K7M2 = "Lung_epithelial", TC1 = "Lung_epithelial")

args <- commandArgs(trailingOnly = TRUE)
geo <- args[[1]]; out <- args[[2]]; model <- args[[3]]
PANELS <- PANELS[setdiff(names(PANELS), setdiff(c(BRAIN, "Lung_epithelial", "Osteoblasts"), TISSUE_PANELS[[model]]))]
IMMUNE <- intersect(IMMUNE, names(PANELS))
src <- utils::read.csv(file.path("data-raw", "bayesprism", "sources.csv"), stringsAsFactors = FALSE)
src <- src[src$model == model, ]
stopifnot(nrow(src) > 0)

# ---- readers ----------------------------------------------------------------
read_mtx_dir <- function(mtx, feat, bc) {
  m <- Matrix::readMM(mtx)
  f <- utils::read.delim(feat, header = FALSE, stringsAsFactors = FALSE)
  b <- utils::read.delim(bc, header = FALSE, stringsAsFactors = FALSE)[[1]]
  if (ncol(f) >= 3) f <- f[f[[3]] %in% c("Gene Expression", f[[3]][1]), , drop = FALSE]
  m <- m[seq_len(nrow(f)), , drop = FALSE]
  rownames(m) <- sub("^[A-Za-z0-9]+___", "", f[[min(2, ncol(f))]]); colnames(m) <- b  # strip genome prefixes (mm10___Xkr4)
  as(m, "CsparseMatrix")
}
read_h5 <- function(path) {
  g <- rhdf5::h5ls(path); grp <- if ("/matrix" %in% paste0(g$group, g$name)) "matrix" else g$name[g$group == "/"][1]
  r <- function(x) rhdf5::h5read(path, paste0(grp, "/", x))
  feat <- if (paste0("/", grp, "/features") %in% paste0(g$group, "/", g$name) | "features" %in% g$name) r("features") else NULL
  sym <- if (!is.null(feat)) feat$name else r("gene_names")
  sym <- as.vector(sym)
  keep <- if (!is.null(feat) && !is.null(feat$feature_type)) as.vector(feat$feature_type) == "Gene Expression" else rep(TRUE, length(sym))
  m <- Matrix::sparseMatrix(i = as.vector(r("indices")) + 1, p = as.vector(r("indptr")), x = as.numeric(r("data")),
                            dims = as.integer(r("shape")), dimnames = list(sym, as.vector(r("barcodes"))))
  m[keep, , drop = FALSE]
}
read_dense <- function(path) {
  d <- data.table::fread(path, data.table = FALSE)
  g <- sub("^ENSMUSG[0-9]+_", "", d[[1]]); m <- as.matrix(d[, -1]); rownames(m) <- g  # ENSMUSG..._Symbol
  if (mean(m == round(m)) < 1) stop("non-integer dense matrix: ", path)
  Matrix::Matrix(m, sparse = TRUE)
}
read_sample <- function(gse, files) {
  p <- file.path(geo, gse, files)
  if (length(p) == 3 && any(grepl("matrix.mtx", p))) {
    return(read_mtx_dir(p[grepl("matrix.mtx", p)], p[grepl("features|genes", p)], p[grepl("barcodes", p)]))
  }
  if (length(p) > 1) p <- p[grepl("filtered", p)]  # filtered and raw h5 offered: use the filtered one
  stopifnot(length(p) == 1)
  if (grepl("\\.h5$", p)) return(read_h5(p))
  if (grepl("\\.tar(\\.gz)?$", p)) {
    d <- file.path(tempdir(), basename(p)); dir.create(d, showWarnings = FALSE)
    utils::untar(p, exdir = d)
    fs <- list.files(d, recursive = TRUE, full.names = TRUE)
    h5 <- fs[grepl("\\.h5$", fs)]
    if (length(h5)) return(read_h5(h5[1]))
    mtx <- fs[grepl("matrix.mtx", fs)]
    if (length(mtx) != 1) mtx <- mtx[grepl("filtered", mtx)][1]
    dd <- dirname(mtx)
    return(read_mtx_dir(mtx, list.files(dd, "features|genes", full.names = TRUE)[1],
                        list.files(dd, "barcodes", full.names = TRUE)[1]))
  }
  read_dense(p)
}

# ---- read, QC ------------------------------------------------------------------
mats <- list(); meta <- list()
for (gsm in unique(src$gsm)) {
  s <- src[src$gsm == gsm, ]
  m <- read_sample(s$gse[1], s$file)
  m <- m[!is.na(rownames(m)) & nzchar(rownames(m)), , drop = FALSE]
  if (anyDuplicated(rownames(m))) m <- Matrix::fac2sparse(factor(rownames(m))) %*% m
  stopifnot(all(m@x == round(m@x)))
  umi <- Matrix::colSums(m); ng <- Matrix::colSums(m > 0)
  mito <- Matrix::colSums(m[grepl("^mt-", rownames(m)), , drop = FALSE]) / pmax(umi, 1)
  keep <- umi >= MIN_UMI & ng >= MIN_GENES & mito <= MAX_MITO
  cat(sprintf("%s %s: %d barcodes, %d pass QC (UMI>=%d, genes>=%d, mito<=%.2f)\n", model, gsm, ncol(m), sum(keep),
              MIN_UMI, MIN_GENES, MAX_MITO))
  m <- m[, keep, drop = FALSE]
  colnames(m) <- paste(gsm, colnames(m), sep = "_")
  mats[[gsm]] <- m
  meta[[gsm]] <- data.frame(cell = colnames(m), sample = gsm, role = s$role[1], umi = umi[keep])
}
genes <- Reduce(union, lapply(mats, rownames))
mats <- lapply(mats, function(m) {
  miss <- setdiff(genes, rownames(m))
  if (length(miss)) m <- rbind(m, Matrix::Matrix(0, length(miss), ncol(m), dimnames = list(miss, colnames(m)), sparse = TRUE))
  m[genes, , drop = FALSE]
})
counts <- do.call(cbind, mats)
meta <- do.call(rbind, meta)

# ---- doublets, clustering --------------------------------------------------------
set.seed(SEED)
sce <- SingleCellExperiment(list(counts = counts))
sce <- scDblFinder::scDblFinder(sce, samples = meta$sample, verbose = FALSE)
single <- sce$scDblFinder.class == "singlet"
cat(sprintf("%s: %d cells after QC, %d doublets removed\n", model, ncol(sce), sum(!single)))
sce <- sce[, single]; meta <- meta[single, ]
sce <- scuttle::logNormCounts(sce)
hvg <- scran::getTopHVGs(scran::modelGeneVar(sce, block = meta$sample), n = N_HVG)
sce <- scater::runPCA(sce, ncomponents = N_PC, subset_row = hvg, BSPARAM = BiocSingular::IrlbaParam())
cl <- scran::clusterCells(sce, use.dimred = "PCA", BLUSPARAM = bluster::NNGraphParam(k = K_NN, cluster.fun = "louvain"))
meta$cluster <- as.integer(cl)

# ---- labels --------------------------------------------------------------------
det <- function(g, idx) if (g %in% rownames(sce)) mean(counts(sce)[g, idx] > 0) else NA_real_
rows <- list()
for (k in sort(unique(meta$cluster))) {
  idx <- which(meta$cluster == k)
  ptprc <- det("Ptprc", idx)
  score <- vapply(PANELS, function(gs) mean(vapply(gs, det, numeric(1), idx = idx), na.rm = TRUE), numeric(1))
  if (is.na(ptprc) || ptprc < 0.3) score[IMMUNE] <- NA            # immune labels need Ptprc
  if (!is.na(ptprc) && ptprc >= 0.5) score[setdiff(names(score), IMMUNE)] <- NA
  best <- if (all(is.na(score))) "none" else names(which.max(score))
  label <- if (best != "none" && score[[best]] >= DETECTION) best else "Tumor"
  if (label %in% TUMOR_AS[[model]]) label <- "Tumor"
  # ambient hemoglobin reaches every cell; call a cluster erythroid only when nearly all cells carry it
  if (label == "Erythroid" && score[["Erythroid"]] < 0.8) {
    score[["Erythroid"]] <- NA
    best <- if (all(is.na(score))) "none" else names(which.max(score))
    label <- if (best != "none" && score[[best]] >= DETECTION) best else "Tumor"
  }
  med_umi <- stats::median(meta$umi[idx])
  if (label == "Tumor" && med_umi < 0.3 * stats::median(meta$umi)) label <- "Low_quality"
  rows[[k]] <- data.frame(cluster = k, n = length(idx), label = label, ptprc = round(ptprc, 2),
                          best = best, best_score = if (best == "none") NA else round(max(score, na.rm = TRUE), 2), median_umi = med_umi,
                          samples = paste(names(table(meta$sample[idx])), table(meta$sample[idx]), sep = ":", collapse = " "))
}
cls <- do.call(rbind, rows)
meta$cell_type <- cls$label[match(meta$cluster, cls$cluster)]
keep <- !meta$cell_type %in% DROP
meta <- meta[keep, ]
meta$cell_state <- meta$cell_type
meta$cell_state[meta$cell_type == "Tumor"] <- paste0("Tumor_c", meta$cluster[meta$cell_type == "Tumor"])
m <- counts(sce)[, meta$cell]
m <- m[Matrix::rowSums(m) > 0, ]
dir <- file.path(out, model); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(m, file.path(dir, "counts.rds"))
utils::write.csv(meta, file.path(dir, "meta.csv"), row.names = FALSE)
utils::write.csv(cls, file.path(dir, "clusters.csv"), row.names = FALSE)
print(cls[, c("cluster", "n", "label", "ptprc", "best", "best_score", "median_umi")], row.names = FALSE)
cat(model, "cells per type:", paste(names(table(meta$cell_type)), table(meta$cell_type), collapse = ", "), "\n")
