# Build the bundled BayesPrism references: one GEP (summed UMI counts per cell
# state, states x genes) per mouse tumor model. Gene filtering is done here so
# that the run time needs neither org.Mm.eg.db nor BayesPrism::cleanup.genes:
# protein-coding genes (org.Mm.eg.db GENETYPE), without the ribosomal,
# mitochondrial and chrY groups BayesPrism recommends removing, without Xist.
#
#   Rscript data-raw/bayesprism/04_build_gep.R <work_dir>
#
# The models, their reference data and validation results are listed in
# data-raw/bayesprism/model_table.csv (results from validation/bayesprism/).
# Writes inst/extdata/bayesprism/<model>.rds and inst/registry/bayesprism_models.csv.

suppressPackageStartupMessages(library(Matrix))

MIN_CELLS_STATE <- 20
GENE_GROUPS <- c("Rb", "Mrp", "other_Rb", "chrM", "chrY")

args <- commandArgs(trailingOnly = TRUE)
stopifnot("usage: 04_build_gep.R <work_dir>" = length(args) == 1)
work <- args[[1]]
out_dir <- file.path("inst", "extdata", "bayesprism")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tab <- utils::read.csv(file.path("data-raw", "bayesprism", "model_table.csv"), stringsAsFactors = FALSE,
                       check.names = FALSE)

a <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
  keys = AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, "ENTREZID"), keytype = "ENTREZID",
  columns = c("SYMBOL", "GENETYPE")))
protein_coding <- unique(a$SYMBOL[a$GENETYPE == "protein-coding"])

load_ref <- function(model, kind) {
  if (kind == "humu") {
    m <- readRDS(file.path(work, "gse310560", model, "counts.rds"))
    meta <- utils::read.csv(file.path(work, "gse310560", model, "meta.csv"), stringsAsFactors = FALSE)
    meta$cell_state <- paste(meta$cell_type, meta$cell_state, sep = "|")
  } else if (kind == "kumar") {
    m <- readRDS(file.path(work, "gse121861", "counts.rds"))
    meta <- utils::read.csv(file.path(work, "gse121861", "annot", "labels.csv"), stringsAsFactors = FALSE)
    meta <- meta[meta$model == model & meta$cell_type != "low_quality", ]
    meta$cell_type[meta$cell_type == "malignant"] <- "Tumor"
    meta$cell_state <- ifelse(meta$cell_type == "Tumor", paste0("Tumor|", meta$cluster), meta$cell_type)
    m <- m[, meta$cell]
  } else {
    m <- readRDS(file.path(work, "geo_annot", model, "counts.rds")); m@x <- as.numeric(m@x)
    meta <- utils::read.csv(file.path(work, "geo_annot", model, "meta.csv"), stringsAsFactors = FALSE)
    meta <- meta[meta$role == "ref", ]
    m <- m[, meta$cell]
  }
  stopifnot(identical(colnames(m), meta$cell))
  list(counts = m, meta = meta)
}

rows <- list()
for (i in seq_len(nrow(tab))) {
  model <- tab$model[i]
  d <- load_ref(model, tab$ref[i])
  n_state <- table(d$meta$cell_state)
  small <- d$meta$cell_state %in% names(n_state)[n_state < MIN_CELLS_STATE]
  d$meta$cell_state[small] <- paste0(d$meta$cell_type[small], "|other")  # merge small states within a type
  n_state <- table(d$meta$cell_state)
  keep <- d$meta$cell_state %in% names(n_state)[n_state >= MIN_CELLS_STATE]
  state <- factor(d$meta$cell_state[keep])
  gep <- as.matrix(Matrix::fac2sparse(state) %*% Matrix::t(d$counts[, keep]))
  stopifnot(all(gep == round(gep)))
  type_of_state <- tapply(d$meta$cell_type[keep], state, unique)
  stopifnot(all(lengths(type_of_state) == 1), "Tumor" %in% unlist(type_of_state))

  gep <- gep[, colSums(gep) > 0]
  up <- toupper(colnames(gep)); gep <- gep[, !(duplicated(up) | duplicated(up, fromLast = TRUE))]
  utils::capture.output(g <- BayesPrism::cleanup.genes(input = gep, input.type = "GEP", species = "mm",
                                                         gene.group = GENE_GROUPS, exp.cells = 1))
  genes <- setdiff(intersect(colnames(g), protein_coding), "Xist")
  gep <- gep[, genes]
  storage.mode(gep) <- "integer"
  types <- unname(unlist(type_of_state)[rownames(gep)])

  ref <- list(gep = gep, state_type = types, key = "Tumor", model = model, source = tab$source[i],
              status = tab$status[i], n_cells = sum(keep))
  f <- file.path(out_dir, paste0(model, ".rds"))
  saveRDS(ref, f, compress = "xz")
  rows[[model]] <- data.frame(model = model, status = tab$status[i], states = nrow(gep), genes = ncol(gep),
                              cells = sum(keep), cell_types = paste(sort(unique(types)), collapse = "|"),
                              source = tab$source[i], tumor_cells = sum(d$meta$cell_type[keep] == "Tumor"),
                              rmse_malignant_independent = tab$rmse_malignant_independent[i],
                              rmse_malignant_within = tab$rmse_malignant_within[i], note = tab$note[i])
  cat(sprintf("%s: %d cells, %d states, %d genes, %.2f MB\n", model, sum(keep), nrow(gep), ncol(gep), file.size(f) / 1e6))
}
utils::write.csv(do.call(rbind, rows), file.path("inst", "registry", "bayesprism_models.csv"), row.names = FALSE)
