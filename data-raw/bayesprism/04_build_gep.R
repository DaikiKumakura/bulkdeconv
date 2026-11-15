# Build the bundled BayesPrism references: one GEP (summed UMI counts per cell
# state, states x genes) per mouse tumor model, from GSE310560 (both mice of
# the young control groups, author annotation). Gene filtering is done here so
# that the run time needs neither org.Mm.eg.db nor BayesPrism::cleanup.genes:
# protein-coding genes (org.Mm.eg.db GENETYPE), without the ribosomal,
# mitochondrial and chrY groups BayesPrism recommends removing, without Xist.
#
#   Rscript data-raw/bayesprism/04_build_gep.R <work_dir>
#
# Writes inst/extdata/bayesprism/<model>.rds and inst/registry/bayesprism_models.csv.
# Models are listed only after validation (validation/bayesprism/10_cross_dataset.R).

suppressPackageStartupMessages(library(Matrix))

MIN_CELLS_STATE <- 20
GENE_GROUPS <- c("Rb", "Mrp", "other_Rb", "chrM", "chrY")
# Malignant-cell RMSE of the cross-dataset check (HuMu reference -> GSE121861
# pseudo-bulk; validation/bayesprism, 2026-11-15)
MODELS <- c(B16F10 = "0.007", MC38 = "0.077", CT26 = "0.098", LLC = "0.039")

args <- commandArgs(trailingOnly = TRUE)
stopifnot("usage: 04_build_gep.R <work_dir>" = length(args) == 1)
work <- args[[1]]
out_dir <- file.path("inst", "extdata", "bayesprism")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

a <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
  keys = AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, "ENTREZID"), keytype = "ENTREZID",
  columns = c("SYMBOL", "GENETYPE")))
protein_coding <- unique(a$SYMBOL[a$GENETYPE == "protein-coding"])

rows <- list()
for (model in names(MODELS)) {
  m <- readRDS(file.path(work, "gse310560", model, "counts.rds"))
  meta <- utils::read.csv(file.path(work, "gse310560", model, "meta.csv"), stringsAsFactors = FALSE)
  stopifnot(identical(colnames(m), meta$cell))
  meta$cell_state <- paste(meta$cell_type, meta$cell_state, sep = "|")
  n_state <- table(meta$cell_state)
  keep <- meta$cell_state %in% names(n_state)[n_state >= MIN_CELLS_STATE]
  state <- factor(meta$cell_state[keep])
  gep <- as.matrix(Matrix::fac2sparse(state) %*% Matrix::t(m[, keep]))
  stopifnot(all(gep == round(gep)))
  type_of_state <- tapply(meta$cell_type[keep], state, unique)
  stopifnot(all(lengths(type_of_state) == 1), "Tumor" %in% unlist(type_of_state))

  gep <- gep[, colSums(gep) > 0]
  utils::capture.output(g <- BayesPrism::cleanup.genes(input = gep, input.type = "GEP", species = "mm",
                                                         gene.group = GENE_GROUPS, exp.cells = 1))
  genes <- setdiff(intersect(colnames(g), protein_coding), "Xist")
  gep <- gep[, genes]
  storage.mode(gep) <- "integer"

  ref <- list(gep = gep, state_type = unname(unlist(type_of_state)[rownames(gep)]), key = "Tumor",
              model = model, source = "GEO GSE310560 (young control groups, author annotation Coarse_Annot/Fine_Annot)",
              n_cells = sum(keep), n_mice = length(unique(meta$mouse[keep])))
  saveRDS(ref, file.path(out_dir, paste0(model, ".rds")), compress = "xz")
  rows[[model]] <- data.frame(model = model, status = "validated", states = nrow(gep), genes = ncol(gep),
                              cells = sum(keep), cell_types = paste(sort(unique(ref$state_type)), collapse = "|"),
                              source = "GSE310560", cross_rmse_malignant = MODELS[[model]])
  cat(sprintf("%s: %d cells, %d states, %d genes, %.2f MB\n", model, sum(keep), nrow(gep), ncol(gep),
              file.size(file.path(out_dir, paste0(model, ".rds"))) / 1e6))
}
utils::write.csv(do.call(rbind, rows), file.path("inst", "registry", "bayesprism_models.csv"), row.names = FALSE)
