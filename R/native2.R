# More natively implemented methods (see native.R): estimate and
# immucellai_mouse, with the helpers they need (ssGSEA, non-negative least
# squares). Same conventions: X is genes x samples, results samples x features,
# each function self-contained.

# ---- ssGSEA ----
# As in GSVA (Hanzelmann et al. 2013; Barbie et al. 2009) with GSVA's defaults
# (alpha = 0.25, scores divided by the range over all scores). Checked against
# GSVA 2.4.9: agrees to 1e-15 (tests/testthat/test-native.R).
bd_native$ssgsea <- function(X, sets, alpha = 0.25, normalize = TRUE) {
  X <- as.matrix(X)
  sets <- lapply(sets, function(g) intersect(g, rownames(X)))
  sets <- sets[lengths(sets) >= 1L]
  R <- apply(X, 2L, function(x) as.integer(rank(x)))
  if (!is.matrix(R)) R <- matrix(R, ncol = ncol(X))
  es <- vapply(seq_len(ncol(X)), function(j) {
    r <- R[, j]
    o <- order(r, decreasing = TRUE)
    rs <- r[o]
    g <- rownames(X)[o]
    vapply(sets, function(gs) {
      hit <- g %in% gs
      w <- (rs * hit)^alpha * hit
      sum(cumsum(w) / sum(w) - cumsum(!hit) / sum(!hit))
    }, 0)
  }, numeric(length(sets)))
  es <- matrix(es, nrow = length(sets), dimnames = list(names(sets), colnames(X)))
  if (normalize) es <- es / (max(es) - min(es))
  es
}

# ---- non-negative least squares ----
# min ||C x - d||^2 subject to x >= 0 (Lawson & Hanson 1974, active set).
bd_native$nnls <- function(C, d, tol = 1e-10, max_iter = 500L) {
  C <- as.matrix(C)
  n <- ncol(C)
  x <- numeric(n)
  passive <- rep(FALSE, n)
  w <- drop(crossprod(C, d - C %*% x))
  iter <- 0L
  while (any(!passive) && max(w[!passive]) > tol && iter < max_iter) {
    iter <- iter + 1L
    t <- which(!passive)[which.max(w[!passive])]
    passive[t] <- TRUE
    repeat {
      z <- numeric(n)
      z[passive] <- qr.solve(C[, passive, drop = FALSE], d)
      if (all(z[passive] > 0)) break
      neg <- passive & z <= 0
      alpha <- min(x[neg] / (x[neg] - z[neg]))
      x <- x + alpha * (z - x)
      passive <- passive & x > tol
      x[!passive] <- 0
    }
    x <- z
    w <- drop(crossprod(C, d - C %*% x))
  }
  x
}

# ---- ESTIMATE (GPL >= 2; Yoshihara K, et al. Nat Commun 2013, via tidyestimate) ----
bd_native$estimate_fit <- function(X, common_genes, gene_sets) {
  # Keep the genes of the ESTIMATE common-gene list (plus genes found through
  # unambiguous aliases), rank the genes within each sample, and compute the
  # stromal and immune enrichment scores and their sum.
  sym <- rownames(X)
  cg <- common_genes
  direct <- which(sym %in% cg$hgnc_symbol)
  missing <- cg$hgnc_symbol[!cg$hgnc_symbol %in% sym]
  cand <- cg[cg$hgnc_symbol %in% missing, , drop = FALSE]
  m <- match(cand$external_synonym, sym)
  ok <- !is.na(m)
  j <- data.frame(hgnc_symbol = cand$hgnc_symbol[ok], row_id = m[ok], stringsAsFactors = FALSE)
  if (nrow(j)) {
    one_gene <- as.vector(table(j$hgnc_symbol)[j$hgnc_symbol]) == 1L
    one_row <- as.vector(table(j$row_id)[as.character(j$row_id)]) == 1L
    j <- j[one_gene & one_row, , drop = FALSE]
    j <- j[order(j$row_id), , drop = FALSE]  # the original's right_join returns them in data order
  }
  rows <- c(direct, j$row_id)
  mat <- X[rows, , drop = FALSE]
  rownames(mat) <- c(sym[direct], j$hgnc_symbol)
  if (nrow(mat) < 1000L) stop("only ", nrow(mat), " of the 10,412 ESTIMATE genes are in the data", call. = FALSE)
  rk <- apply(mat, 2L, rank, ties.method = "average")
  rk <- 10000 * rk / nrow(rk)
  rownames(rk) <- rownames(mat)
  scores <- matrix(NA_real_, ncol(rk), 2L, dimnames = list(colnames(rk), c("stromal", "immune")))
  for (i in 1:2) {
    gs <- intersect(gene_sets[[i]], rownames(rk))
    if (length(gs) == 0L) next
    for (smp in seq_len(ncol(rk))) {
      o <- order(rk[, smp], decreasing = TRUE)
      v <- rk[o, smp]
      hit <- names(v) %in% gs
      miss <- 1 - hit
      v <- v^0.25
      scores[smp, i] <- sum(cumsum(hit * v / sum(v[hit])) - cumsum(miss / sum(miss)))
    }
  }
  cbind(scores, estimate = scores[, "stromal"] + scores[, "immune"])
}

# ---- ImmuCellAI-mouse (MIT; Miao Y-R, et al. Adv Sci 2023; github.com/BioInfoCloud/ImmuCellAI_mouse) ----
bd_native$immucellai_mouse_fit <- function(X, ref) {
  # Three layers (7 lineages, 20 sub-lineages, 9 T-cell states): marker genes
  # are weighted by their expression in the reference, scored per sample with
  # ssGSEA, compensated for spill-over between related cell types by
  # non-negative least squares, and the layers are combined by multiplying the
  # fractions. Also returns an overall infiltration score. Port of
  # ImmuCellAI_mouse(sample, "rnaseq", group_tag = 0, ...)$abundance.
  # Deviation from the original: the per-sample gene sets of the infiltration
  # score are built by exact sample-name matching; the original uses grep(), so
  # a sample name that is part of another one (e.g. S1 and S10) mixes their
  # genes. The two agree whenever no sample name occurs inside another.
  marker_exp <- ref$marker_exp
  sample_exp <- as.matrix(X)
  mode(sample_exp) <- "numeric"
  sample_ratio <- function(data, mexp, tag) {
    acc <- 0
    for (cell in colnames(mexp)) acc <- acc + log2(data + 1) / mexp[, cell] * tag[, cell]
    acc
  }
  tag_matrix <- function(genes, markers) {
    out <- vapply(names(markers), function(cell) as.numeric(genes %in% markers[[cell]]), numeric(length(genes)))
    matrix(out, nrow = length(genes), dimnames = list(genes, names(markers)))
  }
  rescale <- function(x) round(x - min(x) + min(abs(x)), 3)
  compensate <- function(raw, comp) {
    diag(comp) <- 1
    rows <- rownames(raw)[rownames(raw) %in% rownames(comp)]
    cm <- comp[rows, rows, drop = FALSE]
    if (ncol(raw) == 1L) {
      sc <- as.matrix(nnls(cm, raw[rows, 1L]))
    } else {
      sc <- apply(raw[rows, , drop = FALSE], 2L, function(x) nnls(cm, x))
    }
    sc <- apply(sc, 1L, rescale)
    if (!is.matrix(sc)) sc <- matrix(sc, nrow = 1L)
    colnames(sc) <- rows
    t(sc)
  }
  layer_scores <- function(markers, comp) {
    genes <- intersect(intersect(rownames(sample_exp), unlist(markers)), rownames(marker_exp))
    sam <- sample_exp[genes, , drop = FALSE]
    mexp <- marker_exp[genes, names(markers), drop = FALSE]
    tags <- tag_matrix(genes, markers)
    ratio <- apply(sam, 2L, function(x) sample_ratio(x, mexp, tags))
    ratio <- matrix(ratio, nrow = length(genes), dimnames = list(genes, colnames(sam)))
    es <- ssgsea(ratio, markers)
    res <- apply(es, 1L, rescale)
    res <- matrix(res, nrow = ncol(ratio), dimnames = list(colnames(ratio), rownames(es)))
    out <- compensate(t(res), comp)
    if (ncol(out) > 1L) apply(out, 1L, function(x) round(x, 3)) else t(round(out, 3))
  }
  l1 <- ref$l1_marker
  l2 <- ref$l2_marker
  l3 <- ref$l3_marker
  # infiltration score
  genes <- intersect(intersect(rownames(sample_exp), unlist(l1)), rownames(marker_exp))
  sam <- sample_exp[genes, , drop = FALSE]
  tags1 <- tag_matrix(genes, l1)
  dev <- apply(sam, 2L, function(x) sample_ratio(x, marker_exp[genes, names(l1), drop = FALSE], tags1))
  dev <- matrix(dev, nrow = length(genes), dimnames = list(genes, colnames(sam)))
  infil_mt <- NULL
  for (cell in names(l1)[-c(5, 6)]) {
    gn <- intersect(l1[[cell]], rownames(dev))
    nm <- as.vector(unlist(lapply(colnames(dev), function(x) paste(x, gn, sep = "_"))))
    ex <- matrix(as.vector(dev[gn, , drop = FALSE]), ncol = 1L, dimnames = list(nm, "Ratio"))
    sets <- lapply(colnames(dev), function(smp) paste(smp, gn, sep = "_"))
    names(sets) <- colnames(dev)
    infil_mt <- cbind(infil_mt, ssgsea(ex, sets))
  }
  infil_mt[infil_mt < 0] <- 0
  infil <- round(rowSums(infil_mt), 3)
  # layer 1
  layer1 <- round(t(apply(layer_scores(l1, ref$l1_correction), 1L, function(x) x / sum(x))), 4)
  # layer 2
  t2 <- t(layer_scores(l2, ref$l2_correction))
  groups2 <- list(T_cell = c("CD4_T_cell", "CD8_T_cell", "NKT", "Tgd"),
                  B_cell = c("B1_cell", "Follicular_B", "Germinal_center_B", "Marginal_Zone_B", "Memory_B", "Plasma_cell"),
                  Dendritic_cells = c("cDC1", "cDC2", "MoDC", "pDC"),
                  Granulocytes = c("Basophil", "Eosinophil", "mast_cell", "Neutrophils"),
                  Macrophage = c("M1_macrophage", "M2_macrophage"))
  part <- function(t_mat, cells, parent) {
    sub <- t(t_mat[rownames(t_mat) %in% cells, , drop = FALSE])
    round(t(apply(sub, 1L, function(x) x / sum(x))), 4) * parent
  }
  l2_parts <- lapply(names(groups2), function(g) part(t2, groups2[[g]], layer1[, g]))
  names(l2_parts) <- names(groups2)
  layer2 <- do.call(cbind, l2_parts)
  # layer 3
  t3 <- t(layer_scores(l3, ref$l3_correction))
  groups3 <- list(CD4_T_cell = c("CD4_Tm", "Naive_CD4_T", "T_helper_cell", "Treg"),
                  CD8_T_cell = c("CD8_Tc", "CD8_Tcm", "CD8_Tem", "CD8_Tex", "Naive_CD8_T"))
  l3_parts <- lapply(names(groups3), function(g) part(t3, groups3[[g]], l2_parts$T_cell[, g]))
  final <- round(cbind(layer1, layer2, do.call(cbind, l3_parts)), 4)
  final[is.na(final)] <- 0
  cbind(final, Infiltration_score = infil[rownames(final)])
}
