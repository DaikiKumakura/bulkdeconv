# Stage-1 validation of BayesPrism with GEP references built from public mouse
# tumor scRNA-seq (GSE310560 "HuMu" and GSE121861 "Kumar").
#
# For each model, pseudo-bulk samples with known composition are made from one
# set of cells and deconvolved with a GEP reference (summed counts per cell
# state) built from another set:
#   within  : HuMu mouse 1 -> HuMu mouse 2 (same lab, different animal)
#   cross   : HuMu (both mice) -> Kumar (different lab and experiment)
#   reverse : Kumar -> HuMu (both mice)
# The truth is the share of UMIs of each category over the genes used, which is
# what BayesPrism estimates (README: "fraction of reads over the given gene
# set"). Cell types are compared after mapping to four categories shared by
# both datasets: malignant, myeloid, lymphoid, fibroblast.
#
#   Rscript validation/bayesprism/10_cross_dataset.R <work_dir> <out_dir> [models] [n_cores]
#
# <work_dir> holds gse310560/<model>/{counts.rds,meta.csv} (03_gse310560_extract.R)
# and gse121861/{counts.rds,annot/labels.csv} (01, 01b, 02 in data-raw/bayesprism).

suppressPackageStartupMessages({
  library(Matrix)
  library(BayesPrism)
})

SEED <- 20261115
N_PB <- 20            # pseudo-bulk samples per comparison
N_CELLS_PB <- 2000    # cells per pseudo-bulk sample
MALIGNANT_RANGE <- c(0.2, 0.9)
MIN_CELLS_STATE <- 20 # BayesPrism recommends at least 20 cells per state
GENE_GROUPS <- c("Rb", "Mrp", "other_Rb", "chrM", "chrY")  # mouse has no MALAT1 group; Malat1 is non-coding and dropped below
KUMAR_NAME <- c(B16F10 = "B16F10", MC38 = "MC38", CT26 = "CT26", LLC = "LL2")

CATEGORY <- c(
  # HuMu (Coarse_Annot)
  Tumor = "malignant", CD4_Tconv = "lymphoid", CD4_Treg = "lymphoid", CD8_T_cells = "lymphoid",
  gdT_other_T_cells = "lymphoid", cDC = "myeloid", pDC = "myeloid", Mono_Mac = "myeloid",
  Monocytes = "myeloid", TAMs = "myeloid", Fibroblasts = "fibroblast",
  # Kumar (02_annotate.R)
  malignant = "malignant", monocyte_macrophage = "myeloid", neutrophil = "myeloid", DC = "myeloid",
  mast_basophil = "myeloid", T_CD8 = "lymphoid", T_CD4 = "lymphoid", NK = "lymphoid", B = "lymphoid",
  plasma = "lymphoid", immune_unassigned = "lymphoid", fibroblast = "fibroblast", endothelial = "endothelial"
)

args <- commandArgs(trailingOnly = TRUE)
stopifnot("usage: 10_cross_dataset.R <work_dir> <out_dir> [models] [n_cores]" = length(args) >= 2)
work <- args[[1]]
out_dir <- args[[2]]
models <- if (length(args) >= 3) strsplit(args[[3]], ",")[[1]] else names(KUMAR_NAME)
n_cores <- if (length(args) >= 4) as.integer(args[[4]]) else 10L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- data -------------------------------------------------------------------

load_humu <- function(model) {
  m <- readRDS(file.path(work, "gse310560", model, "counts.rds"))
  meta <- utils::read.csv(file.path(work, "gse310560", model, "meta.csv"), stringsAsFactors = FALSE)
  stopifnot(identical(colnames(m), meta$cell))
  meta$cell_state <- paste(meta$cell_type, meta$cell_state, sep = "|")
  meta$category <- unname(CATEGORY[meta$cell_type])
  stopifnot(!anyNA(meta$category))
  list(counts = m, meta = meta, malignant = "Tumor")
}

kumar_all <- NULL
load_kumar <- function(model) {
  if (is.null(kumar_all)) {
    kumar_all <<- list(
      counts = readRDS(file.path(work, "gse121861", "counts.rds")),
      labels = utils::read.csv(file.path(work, "gse121861", "annot", "labels.csv"), stringsAsFactors = FALSE)
    )
  }
  lab <- kumar_all$labels[kumar_all$labels$model == KUMAR_NAME[[model]] & kumar_all$labels$cell_type != "low_quality", ]
  m <- kumar_all$counts[, lab$cell]
  m <- m[Matrix::rowSums(m) > 0, ]
  lab$cell_state <- ifelse(lab$cell_type == "malignant", paste0("malignant|", lab$cluster), lab$cell_type)
  lab$category <- unname(CATEGORY[lab$cell_type])
  stopifnot(!anyNA(lab$category))
  list(counts = m, meta = lab, malignant = "malignant")
}

subset_cells <- function(d, keep) {
  list(counts = d$counts[, keep], meta = d$meta[keep, ], malignant = d$malignant)
}

# ---- reference ----------------------------------------------------------------

# GEP: summed counts per cell state (states x genes), with the type of each state.
make_gep <- function(d) {
  n_state <- table(d$meta$cell_state)
  keep_state <- names(n_state)[n_state >= MIN_CELLS_STATE]
  keep <- d$meta$cell_state %in% keep_state
  meta <- d$meta[keep, ]
  state <- factor(meta$cell_state)
  gep <- as.matrix(Matrix::fac2sparse(state) %*% Matrix::t(d$counts[, keep]))
  type_of_state <- tapply(meta$cell_type, state, function(v) unique(v))
  stopifnot(all(lengths(type_of_state) == 1))
  list(gep = gep, state_type = unlist(type_of_state)[rownames(gep)],
       dropped_states = setdiff(names(n_state), keep_state), malignant = d$malignant,
       n_cells = sum(keep))
}

protein_coding <- local({
  pc <- NULL
  function(genes) {
    if (is.null(pc)) {
      a <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
        keys = AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, "ENTREZID"), keytype = "ENTREZID",
        columns = c("SYMBOL", "GENETYPE")))
      pc <<- unique(a$SYMBOL[a$GENETYPE == "protein-coding"])
    }
    genes[genes %in% pc]
  }
})

# Genes used: shared by reference and mixture, protein-coding, without the
# groups BayesPrism recommends removing, without Xist.
select_genes <- function(gep, mixture_genes) {
  shared <- intersect(colnames(gep), mixture_genes)
  cleaned <- utils::capture.output(g <- BayesPrism::cleanup.genes(input = gep[, shared, drop = FALSE],
    input.type = "GEP", species = "mm", gene.group = GENE_GROUPS, exp.cells = 1))
  genes <- setdiff(protein_coding(colnames(g)), "Xist")
  genes
}

# ---- pseudo-bulk ----------------------------------------------------------------

make_pseudobulk <- function(d, seed) {
  set.seed(seed)
  cats <- names(which(table(d$meta$category) >= MIN_CELLS_STATE))
  others <- setdiff(cats, "malignant")
  frac <- t(vapply(seq_len(N_PB), function(i) {
    p <- stats::setNames(numeric(length(cats)), cats)
    mal <- if ("malignant" %in% cats) stats::runif(1, MALIGNANT_RANGE[1], MALIGNANT_RANGE[2]) else 0
    w <- stats::rgamma(length(others), 1)
    p[others] <- (1 - mal) * w / sum(w)
    if ("malignant" %in% cats) p[["malignant"]] <- mal
    p
  }, numeric(length(cats))))
  idx <- split(seq_len(nrow(d$meta)), d$meta$category)
  by_cat <- lapply(stats::setNames(cats, cats), function(k) {
    sapply(seq_len(N_PB), function(i) {
      n <- round(frac[i, k] * N_CELLS_PB)
      if (n == 0) return(numeric(nrow(d$counts)))
      Matrix::rowSums(d$counts[, sample(idx[[k]], n, replace = TRUE), drop = FALSE])
    })
  })
  for (k in cats) rownames(by_cat[[k]]) <- rownames(d$counts)
  list(by_cat = by_cat, cell_frac = frac)
}

# ---- BayesPrism -------------------------------------------------------------------

run_bayesprism <- function(ref, mixture, input.type = "GEP", reference = ref$gep,
                           types = ref$state_type, states = rownames(ref$gep)) {
  t0 <- Sys.time()
  log <- utils::capture.output(type = "output", {
    prism <- BayesPrism::new.prism(reference = reference, mixture = mixture, input.type = input.type,
                                   cell.type.labels = types, cell.state.labels = states,
                                   key = ref$malignant, outlier.cut = 0.01, outlier.fraction = 0.1)
    bp <- BayesPrism::run.prism(prism = prism, n.cores = n_cores)
  })
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(final = BayesPrism::get.fraction(bp, which.theta = "final", state.or.type = "type"),
       first = BayesPrism::get.fraction(bp, which.theta = "first", state.or.type = "type"),
       seconds = secs, n_genes = ncol(prism@mixture))
}

to_category <- function(theta) {
  cat_of <- unname(CATEGORY[colnames(theta)])
  stopifnot(!anyNA(cat_of))
  t(rowsum(t(theta), cat_of))
}

score <- function(est, truth, label) {
  cats <- union(colnames(est), colnames(truth))
  fill <- function(m) { out <- matrix(0, nrow(m), length(cats), dimnames = list(rownames(m), cats)); out[, colnames(m)] <- m; out }
  est <- fill(est); truth <- fill(truth)
  do.call(rbind, lapply(cats, function(k) {
    e <- est[, k]; t <- truth[, k]
    data.frame(comparison = label, category = k, truth_mean = mean(t), est_mean = mean(e),
               bias = mean(e - t), rmse = sqrt(mean((e - t)^2)),
               r = if (stats::sd(t) > 0 && stats::sd(e) > 0) stats::cor(e, t) else NA_real_)
  }))
}

evaluate <- function(ref_data, test_data, label, seed) {
  ref <- make_gep(ref_data)
  pb <- make_pseudobulk(test_data, seed)
  mixture_genes <- rownames(pb$by_cat[[1]])
  genes <- select_genes(ref$gep, mixture_genes)
  mixture <- t(Reduce(`+`, lapply(pb$by_cat, function(x) x[genes, , drop = FALSE])))
  rownames(mixture) <- sprintf("pb%02d", seq_len(nrow(mixture)))
  truth <- sapply(pb$by_cat, function(x) colSums(x[genes, , drop = FALSE]))
  truth <- truth / rowSums(truth)
  rownames(truth) <- rownames(mixture)
  res <- run_bayesprism(ref, mixture[, genes])
  est <- to_category(res$final)
  est_first <- to_category(res$first)
  cat(sprintf("  %-34s ref %d states (%d cells), %d genes used, %.0f s\n", label, nrow(ref$gep),
              ref$n_cells, res$n_genes, res$seconds))
  list(metrics = rbind(score(est, truth, paste(label, "final")), score(est_first, truth, paste(label, "first"))),
       timing = data.frame(comparison = label, states = nrow(ref$gep), genes = res$n_genes,
                           samples = N_PB, cores = n_cores, seconds = res$seconds),
       est = est, truth = truth, ref = ref)
}

# ---- run ----------------------------------------------------------------------------

metrics <- list(); timing <- list(); equivalence <- list()
for (model in models) {
  cat(sprintf("\n== %s\n", model))
  humu <- load_humu(model)
  kumar <- load_kumar(model)
  mice <- sort(unique(humu$meta$mouse))
  stopifnot(length(mice) == 2)
  runs <- list(
    within = evaluate(subset_cells(humu, humu$meta$mouse == mice[1]), subset_cells(humu, humu$meta$mouse == mice[2]),
                      sprintf("%s within (HuMu %s -> %s)", model, mice[1], mice[2]), SEED),
    cross = evaluate(humu, kumar, sprintf("%s cross (HuMu -> Kumar)", model), SEED + 1),
    reverse = evaluate(kumar, humu, sprintf("%s reverse (Kumar -> HuMu)", model), SEED + 2)
  )
  for (r in runs) { metrics[[length(metrics) + 1]] <- cbind(model = model, r$metrics); timing[[length(timing) + 1]] <- cbind(model = model, r$timing) }
  saveRDS(lapply(runs, function(r) list(est = r$est, truth = r$truth, states = r$ref$state_type,
                                         dropped_states = r$ref$dropped_states)),
          file.path(out_dir, sprintf("runs_%s.rds", model)))
}

metrics <- do.call(rbind, metrics)
timing <- do.call(rbind, timing)
utils::write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)
utils::write.csv(timing, file.path(out_dir, "timing.csv"), row.names = FALSE)
cat("\n== metrics (final theta)\n")
print(metrics[grepl("final$", metrics$comparison), c("comparison", "category", "truth_mean", "est_mean", "bias", "rmse", "r")],
      digits = 3, row.names = FALSE)
cat("\n== timing\n"); print(timing, digits = 3, row.names = FALSE)
