# Builds inst/extdata/ (reference data of the four natively implemented methods:
# abis, seqimmucc, base, timer) from the files of the omnideconv/immunedeconv
# repository at a fixed commit. immunedeconv itself is neither installed nor
# called; only these data files are read. Provenance: inst/COPYRIGHTS.
#
#   Rscript data-raw/build_extdata.R
#
# Changes made to the data (all documented in inst/COPYRIGHTS):
#   * base: the compendium weights are rounded to 6 decimals (max change 5e-7)
#     and stored xz-compressed, 38 MB -> 11 MB.
#   * timer: the 32 per-cancer marker-gene vectors are stored as one named list.

sha <- "e625e6c28ed14a30f9f40f159925cc9f0df4fa49"
raw <- function(path) paste0("https://raw.githubusercontent.com/omnideconv/immunedeconv/", sha, "/inst/extdata/", path)
root <- if (file.exists("DESCRIPTION")) "." else stop("run from the package root")
out <- file.path(root, "inst", "extdata")
tmp <- tempfile(); dir.create(tmp)
fetch <- function(path) {
  f <- file.path(tmp, basename(path))
  utils::download.file(raw(path), f, mode = "wb", quiet = TRUE)
  f
}
for (d in c("seqimmucc", "abis", "base", "timer")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

file.copy(fetch("mouse_deconvolution/sig_matr_seqImmuCC.txt"), file.path(out, "seqimmucc", "sig_matr_seqImmuCC.txt"), overwrite = TRUE)
file.copy(fetch("abis/sigmatrixRNAseq.txt"), file.path(out, "abis", "sigmatrixRNAseq.txt"), overwrite = TRUE)
file.copy(fetch("mouse_deconvolution/cell_correspondance.txt"), file.path(out, "base", "cell_correspondance.txt"), overwrite = TRUE)

comp <- as.matrix(readRDS(fetch("mouse_deconvolution/BASE_cell_compendium.rds")))
stopifnot(min(comp) >= 0, max(comp) <= 1)
saveRDS(round(comp, 6), file.path(out, "base", "compendium.rds"), compress = "xz")

cur <- base::get(load(fetch("timer/precalculated/immune.expression.curated.RData")))
stopifnot(is.list(cur), all(c("genes", "celltypes") %in% names(cur)))
saveRDS(list(genes = as.matrix(cur$genes), celltypes = cur$celltypes), file.path(out, "timer", "immune_reference.rds"), compress = "xz")

types <- c("kich", "blca", "brca", "cesc", "gbm", "hnsc", "kirp", "lgg", "lihc", "luad", "lusc", "prad", "sarc", "pcpg",
           "paad", "tgct", "ucec", "ov", "skcm", "dlbc", "kirc", "acc", "meso", "thca", "uvm", "ucs", "thym", "esca",
           "stad", "read", "coad", "chol")
markers <- lapply(setNames(types, toupper(types)), function(t) base::get(load(fetch(paste0("timer/precalculated/genes_", t, ".RData")))))
saveRDS(markers, file.path(out, "timer", "markers_by_cancer.rds"), compress = "xz")

# ---- immucellai2: reference restricted to the 5510 marker genes (the original reads
# a 41 MB text file and keeps these rows). Source: GuoBioinfoLab/ImmuCellAI-2.0 (MIT).
icai <- "https://raw.githubusercontent.com/GuoBioinfoLab/ImmuCellAI-2.0/0135016ee966fbd974f9f694f0f468e627a09ab9/R-package/ImmuCellAI2.0/inst/extdata/"
dir.create(file.path(out, "immucellai2"), showWarnings = FALSE)
mk_file <- tempfile(); utils::download.file(paste0(icai, "MarkerUsedDeconvolution_5510.txt"), mk_file, mode = "wb", quiet = TRUE)
markers <- trimws(readLines(mk_file, warn = FALSE)); markers <- markers[nzchar(markers)]
if (tolower(markers[1]) %in% c("gene", "genes", "symbol")) markers <- markers[-1]
markers <- unique(markers)
ref_file <- tempfile(); utils::download.file(paste0(icai, "reference_53celltypesTPM20260518.txt"), ref_file, mode = "wb", quiet = TRUE)
x <- utils::read.table(ref_file, header = TRUE, sep = "	", check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE, row.names = NULL)
genes <- trimws(as.character(x[[1L]]))
m <- as.matrix(x[-1L]); storage.mode(m) <- "double"
stopifnot(!anyNA(genes), all(is.finite(m)), all(m >= 0))
if (anyDuplicated(genes)) {  # same rule as the original reader: duplicates averaged
  m <- rowsum(m, group = genes, reorder = FALSE) / as.numeric(table(factor(genes, levels = unique(genes))))
  genes <- unique(genes)
}
rownames(m) <- genes
m <- m[markers[markers %in% rownames(m)], , drop = FALSE]
saveRDS(list(reference = m, markers = markers), file.path(out, "immucellai2", "reference_markers.rds"), compress = "xz")
cat("immucellai2 reference:", paste(dim(m), collapse = " x "), "of", length(markers), "markers
")

for (f in list.files(out, recursive = TRUE, full.names = TRUE)) cat(sprintf("%-60s %8.2f MB\n", sub(paste0(out, "/"), "", f, fixed = TRUE), file.size(f) / 1e6))
