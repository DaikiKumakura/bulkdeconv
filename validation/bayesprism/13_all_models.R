# Pseudo-bulk check of the BayesPrism reference of each model.
#   independent: reference from the reference dataset, pseudo-bulk from another
#                dataset of the same model (different lab), when one exists
#   within     : reference from all but one sample (mouse), pseudo-bulk from the
#                held-out sample; with a single sample, a random half of the
#                cells (a weak check: same mouse and run)
# Truth: share of UMIs per category over the genes used (BayesPrism's theta).
#
#   Rscript validation/bayesprism/13_all_models.R <work_dir> <out_dir> <models> [n_cores]

suppressPackageStartupMessages({
  library(Matrix)
  library(BayesPrism)
})

SEED <- 20261115
N_PB <- 10
N_CELLS_PB <- 2000
MALIGNANT_RANGE <- c(0.2, 0.9)
MIN_CELLS_STATE <- 20
GENE_GROUPS <- c("Rb", "Mrp", "other_Rb", "chrM", "chrY")
LYMPHOID <- c("CD4_Tconv", "CD4_Treg", "CD8_T_cells", "gdT_other_T_cells", "NK", "B_cells", "Plasma",
              "T_CD8", "T_CD4", "B", "plasma", "immune_unassigned")
MYELOID <- c("cDC", "pDC", "Mono_Mac", "Monocytes", "TAMs", "Macrophages", "Neutrophils", "Mast_Basophils",
             "Microglia", "monocyte_macrophage", "neutrophil", "DC", "mast_basophil")
category_of <- function(type) {
  ifelse(type %in% c("Tumor", "malignant"), "malignant",
         ifelse(type %in% LYMPHOID, "lymphoid", ifelse(type %in% MYELOID, "myeloid", "stroma")))
}

# reference dataset and independent check dataset per model
PLAN <- list(
  B16F10 = list(ref = "humu"), MC38 = list(ref = "humu", check = "geo"), CT26 = list(ref = "humu", check = "geo"),
  LLC = list(ref = "humu"), RENCA = list(ref = "humu", check = "geo"), `4T1` = list(ref = "humu", check = "geo"),
  `KPCY-F` = list(ref = "humu"), `KPCY-C2` = list(ref = "humu"),
  EMT6 = list(ref = "kumar"), Sa1N = list(ref = "kumar"),
  E0771 = list(ref = "geo", check = "geo"), MB49 = list(ref = "geo", check = "geo"),
  F420 = list(ref = "geo", check = "geo"), K7M2 = list(ref = "geo", check = "geo")
)
plan_of <- function(m) PLAN[[m]] %||% list(ref = "geo")
`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
work <- args[[1]]; out_dir <- args[[2]]; models <- strsplit(args[[3]], ",")[[1]]
n_cores <- if (length(args) >= 4) as.integer(args[[4]]) else 10L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_set <- function(model, kind, role = NULL) {
  if (kind == "humu") {
    m <- readRDS(file.path(work, "gse310560", model, "counts.rds"))
    meta <- utils::read.csv(file.path(work, "gse310560", model, "meta.csv"), stringsAsFactors = FALSE)
    meta$sample <- meta$mouse
    meta$cell_state <- paste(meta$cell_type, meta$cell_state, sep = "|")
  } else if (kind == "kumar") {
    m <- readRDS(file.path(work, "gse121861", "counts.rds"))
    lab <- utils::read.csv(file.path(work, "gse121861", "annot", "labels.csv"), stringsAsFactors = FALSE)
    meta <- lab[lab$model == model & lab$cell_type != "low_quality", ]
    cells <- utils::read.csv(file.path(work, "gse121861", "cells_model.csv"), stringsAsFactors = FALSE)
    meta$sample <- cells$Sample_ID[match(meta$cell, cells$cell)]
    meta$cell_type[meta$cell_type == "malignant"] <- "Tumor"
    meta$cell_state <- ifelse(meta$cell_type == "Tumor", paste0("Tumor|", meta$cluster), meta$cell_type)
    m <- m[, meta$cell]
  } else {
    m <- readRDS(file.path(work, "geo_annot", model, "counts.rds")); m@x <- as.numeric(m@x)
    meta <- utils::read.csv(file.path(work, "geo_annot", model, "meta.csv"), stringsAsFactors = FALSE)
    if (!is.null(role)) { keep <- meta$role == role; meta <- meta[keep, ]; m <- m[, meta$cell] }
  }
  stopifnot(identical(colnames(m), meta$cell))
  meta$category <- category_of(meta$cell_type)
  list(counts = m[Matrix::rowSums(m) > 0, ], meta = meta)
}
subset_cells <- function(d, keep) list(counts = d$counts[, keep], meta = d$meta[keep, ])

merge_small_states <- function(meta) {
  n_state <- table(meta$cell_state)
  small <- meta$cell_state %in% names(n_state)[n_state < MIN_CELLS_STATE]
  meta$cell_state[small] <- paste0(meta$cell_type[small], "|other")
  meta
}
make_gep <- function(d) {
  d$meta <- merge_small_states(d$meta)
  n_state <- table(d$meta$cell_state)
  keep <- d$meta$cell_state %in% names(n_state)[n_state >= MIN_CELLS_STATE]
  state <- factor(d$meta$cell_state[keep])
  gep <- as.matrix(Matrix::fac2sparse(state) %*% Matrix::t(d$counts[, keep]))
  type <- tapply(d$meta$cell_type[keep], state, unique)
  list(gep = gep, state_type = unlist(type)[rownames(gep)])
}
protein_coding <- local({
  pc <- NULL
  function(g) {
    if (is.null(pc)) {
      a <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
        keys = AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, "ENTREZID"), keytype = "ENTREZID", columns = c("SYMBOL", "GENETYPE")))
      pc <<- unique(a$SYMBOL[a$GENETYPE == "protein-coding"])
    }
    g[g %in% pc]
  }
})
select_genes <- function(gep, mixture_genes) {
  shared <- intersect(colnames(gep), mixture_genes)
  shared <- shared[!duplicated(toupper(shared)) & !duplicated(toupper(shared), fromLast = TRUE)]
  utils::capture.output(g <- BayesPrism::cleanup.genes(input = gep[, shared, drop = FALSE], input.type = "GEP",
                                                         species = "mm", gene.group = GENE_GROUPS, exp.cells = 1))
  setdiff(protein_coding(colnames(g)), "Xist")
}
make_pseudobulk <- function(d, seed) {
  set.seed(seed)
  cats <- names(which(table(d$meta$category) >= MIN_CELLS_STATE))
  others <- setdiff(cats, "malignant")
  frac <- t(vapply(seq_len(N_PB), function(i) {
    p <- stats::setNames(numeric(length(cats)), cats)
    mal <- if ("malignant" %in% cats) stats::runif(1, MALIGNANT_RANGE[1], MALIGNANT_RANGE[2]) else 0
    w <- stats::rgamma(length(others), 1); p[others] <- (1 - mal) * w / sum(w)
    if ("malignant" %in% cats) p[["malignant"]] <- mal
    p
  }, numeric(length(cats))))
  idx <- split(seq_len(nrow(d$meta)), d$meta$category)
  lapply(stats::setNames(cats, cats), function(k) {
    x <- sapply(seq_len(N_PB), function(i) {
      n <- round(frac[i, k] * N_CELLS_PB)
      if (n == 0) return(numeric(nrow(d$counts)))
      Matrix::rowSums(d$counts[, sample(idx[[k]], n, replace = TRUE), drop = FALSE])
    })
    rownames(x) <- rownames(d$counts); x
  })
}
evaluate <- function(ref_d, test_d, label, seed) {
  ref <- make_gep(ref_d)
  if (!"Tumor" %in% ref$state_type) return(NULL)
  pb <- make_pseudobulk(test_d, seed)
  genes <- select_genes(ref$gep, rownames(pb[[1]]))
  mixture <- t(Reduce(`+`, lapply(pb, function(x) x[genes, , drop = FALSE])))
  rownames(mixture) <- sprintf("pb%02d", seq_len(nrow(mixture)))
  truth <- sapply(pb, function(x) colSums(x[genes, , drop = FALSE])); truth <- truth / rowSums(truth)
  t0 <- Sys.time()
  utils::capture.output({
    p <- new.prism(reference = ref$gep[, genes], mixture = mixture, input.type = "GEP", cell.type.labels = ref$state_type,
                   cell.state.labels = rownames(ref$gep), key = "Tumor", outlier.cut = 0.01, outlier.fraction = 0.1)
    bp <- run.prism(p, n.cores = n_cores)
  })
  th <- get.fraction(bp, "final", "type")
  est <- t(rowsum(t(th), category_of(colnames(th))))
  cats <- union(colnames(est), colnames(truth))
  fill <- function(m) { o <- matrix(0, nrow(m), length(cats), dimnames = list(NULL, cats)); o[, colnames(m)] <- m; o }
  est <- fill(est); truth <- fill(truth)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  %-40s %d states, %d genes, %.0f s\n", label, nrow(ref$gep), length(genes), secs))
  do.call(rbind, lapply(cats, function(k) data.frame(
    comparison = label, category = k, truth_mean = mean(truth[, k]), est_mean = mean(est[, k]),
    bias = mean(est[, k] - truth[, k]), rmse = sqrt(mean((est[, k] - truth[, k])^2)),
    r = if (stats::sd(truth[, k]) > 0 && stats::sd(est[, k]) > 0) stats::cor(est[, k], truth[, k]) else NA_real_,
    seconds = secs)))
}

for (model in models) {
  cat(sprintf("\n== %s\n", model))
  pl <- plan_of(model)
  ref <- load_set(model, pl$ref, if (pl$ref == "geo") "ref" else NULL)
  res <- list()
  samples <- sort(unique(ref$meta$sample))
  if (length(samples) >= 2) {
    held <- samples[length(samples)]
    res$within <- evaluate(subset_cells(ref, ref$meta$sample != held), subset_cells(ref, ref$meta$sample == held),
                           sprintf("%s within (held-out sample %s)", model, held), SEED)
  } else {
    set.seed(SEED); half <- sample(c(TRUE, FALSE), nrow(ref$meta), replace = TRUE)
    res$within <- evaluate(subset_cells(ref, half), subset_cells(ref, !half), sprintf("%s within (cell halves)", model), SEED)
  }
  if (!is.null(pl$check)) {
    chk <- load_set(model, pl$check, "check")
    res$independent <- evaluate(ref, chk, sprintf("%s independent (%s -> other dataset)", model, pl$ref), SEED + 1)
  }
  met <- do.call(rbind, res)
  if (is.null(met)) { cat("  no Tumor state with enough cells; skipped\n"); next }
  met <- cbind(model = model, met)
  utils::write.csv(met, file.path(out_dir, sprintf("metrics_%s.csv", model)), row.names = FALSE)
  print(met[, c("comparison", "category", "truth_mean", "est_mean", "rmse", "r")], digits = 3, row.names = FALSE)
}
