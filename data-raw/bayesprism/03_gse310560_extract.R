# Extract per-model count matrices and author annotations from the GSE310560
# (HuMu) Seurat object: singlets only, annotated doublets removed, young
# untreated mice only (the two control groups of the dataset).
#
#   Rscript data-raw/bayesprism/03_gse310560_extract.R <seurat.rds> <out_dir> [models]
#
# Writes <out_dir>/<model>/counts.rds (dgCMatrix genes x cells) and meta.csv
# (cell, mouse, cell_type = Coarse_Annot, cell_state = Fine_Annot, ...).

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
})

CONTROL_CONDITIONS <- c("Young/H2O_LFD/CTRL", "Young/Regular/Regular")

args <- commandArgs(trailingOnly = TRUE)
stopifnot("usage: 03_gse310560_extract.R <seurat.rds> <out_dir> [models]" = length(args) >= 2)
models <- if (length(args) >= 3) strsplit(args[[3]], ",")[[1]] else c("B16F10", "MC38", "CT26", "LLC")

x <- readRDS(args[[1]])
md <- x@meta.data
md$cell <- rownames(md)
md$condition <- paste(md$Age, md$Drink_Diet, md$Microbiota, sep = "/")
md$mouse <- sub("--.*", "", md$exp_keys)
counts <- LayerData(x, assay = "RNA", layer = "counts")
stopifnot(identical(colnames(counts), md$cell), all(counts@x == round(counts@x)))
cat("all cells:", nrow(md), "\n")

for (model in models) {
  keep <- md$Tumor_Line == model & md$doubletFinder == "Singlet" &
    !grepl("^Doublets", md$Fine_Annot) & md$condition %in% CONTROL_CONDITIONS
  n_line <- sum(md$Tumor_Line == model)
  sub <- md[keep, c("cell", "mouse", "tissue_cell", "condition", "Coarse_Annot", "Fine_Annot", "nCount_RNA")]
  names(sub)[names(sub) == "Coarse_Annot"] <- "cell_type"
  names(sub)[names(sub) == "Fine_Annot"] <- "cell_state"
  m <- counts[, sub$cell]
  m <- m[Matrix::rowSums(m) > 0, ]
  stopifnot(all(Matrix::colSums(m) == sub$nCount_RNA))
  dir <- file.path(args[[2]], model)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(m, file.path(dir, "counts.rds"))
  utils::write.csv(sub, file.path(dir, "meta.csv"), row.names = FALSE)
  cat(sprintf("%s: %d of %d cells kept (%d genes); mice %s; cells per type: %s\n", model, nrow(sub), n_line,
              nrow(m), paste(sort(unique(sub$mouse)), collapse = ", "),
              paste(sprintf("%s %d", names(table(sub$cell_type)), table(sub$cell_type)), collapse = ", ")))
}
