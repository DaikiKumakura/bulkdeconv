# Methods implemented inside bulkdeconv (engine "native"): abis, seqimmucc
# (LLSR), base and timer. They are ports of the published algorithms, written
# to give the same numbers as the original code (see tests/testthat/test-native.R
# and inst/COPYRIGHTS for sources and licences), and they call no other
# deconvolution package. Each function is self-contained: it runs in a
# separate R process together with the other helpers (see bd_isolate()).
#
# Conventions: X is genes x samples (symbols as row names); results are
# samples x features.

bd_native <- list(

  # seq-ImmuCC, linear least-squares regression (LLSR) variant
  # (Chen et al. 2018, doi:10.3389/fimmu.2018.01286). Values are percent and
  # sum to 100 per sample. The original code computes a quantile
  # normalisation but never uses it; it is therefore not applied here either.
  llsr = function(X, signature) {
    signature <- data.matrix(signature)
    X <- data.matrix(X)
    if (max(signature) > 50) signature <- log2(signature + 1)
    if (max(X) > 50) X <- log2(X + 1)
    signature <- (signature - mean(signature)) / stats::sd(as.vector(signature))
    common <- intersect(rownames(signature), rownames(X))
    if (length(common) < 10L) stop("fewer than 10 genes of the signature are in the data", call. = FALSE)
    signature <- signature[common, , drop = FALSE]
    X <- X[common, , drop = FALSE]
    out <- matrix(NA_real_, ncol(X), ncol(signature), dimnames = list(colnames(X), colnames(signature)))
    for (i in seq_len(ncol(X))) {
      y <- X[, i]
      y <- (y - mean(y)) / stats::sd(y)
      b <- stats::lsfit(signature, y, intercept = FALSE)$coefficients
      b[b < 0] <- 0
      out[i, ] <- b / sum(b) * 100
    }
    out
  },

  # ABIS, robust regression on a signature matrix built from sorted blood
  # cells (Monaco et al. 2019, doi:10.1016/j.celrep.2019.01.041), RNA-seq
  # version. Values are regression coefficients times 100, 3 significant digits.
  abis_fit = function(X, signature) {
    genes <- intersect(rownames(X), rownames(signature))
    if (length(genes) < 10L) stop("fewer than 10 genes of the signature are in the data", call. = FALSE)
    sig <- as.matrix(signature[genes, , drop = FALSE])
    Xg <- X[genes, , drop = FALSE]
    co <- vapply(seq_len(ncol(Xg)), function(j) stats::coef(MASS::rlm(sig, Xg[, j], maxit = 100)),
                 numeric(ncol(sig)))
    co <- matrix(co, nrow = ncol(sig), dimnames = list(colnames(sig), colnames(Xg)))
    signif(t(co) * 100, 3)
  },

  # BASE (Varn et al. 2017, doi:10.1158/0008-5472.CAN-16-2490) with mouse
  # cell compendium: quantile normalisation, per-gene median centring,
  # enrichment of up/down weight sets, normalisation by permutations, then
  # the sub-types are combined by the median of their z-scored values.
  # Same arithmetic and same random-number calls as the original; the
  # cumulative sums are vectorised, which makes it faster: 51 s instead of 175 s
  # for 6 samples and 100 permutations (measured 2026-11-11, same machine).
  base_scores = function(X, compendium, mapping, perm = 100L, seed = 1L) {
    set.seed(seed)
    data <- as.matrix(X)
    rk <- apply(data, 2L, rank)
    srt <- apply(data, 2L, sort)
    med <- apply(srt, 1L, stats::median, na.rm = TRUE)
    data[] <- med[rk]
    reg <- as.matrix(compendium)
    common <- intersect(rownames(data), rownames(reg))
    if (length(common) < 100L) stop("fewer than 100 genes of the compendium are in the data", call. = FALSE)
    data <- data[common, , drop = FALSE]
    reg <- reg[common, , drop = FALSE]
    data <- data - apply(data, 1L, stats::median)
    cnum <- ncol(data)
    rnum <- nrow(data)
    # column-wise cumulative sums from one cumsum over the whole matrix
    # (agrees with per-column cumsum to about 1e-11 relative)
    col_cumsum <- function(m) {
      cs <- cumsum(as.vector(m))
      off <- c(0, cs[seq_len(ncol(m) - 1L) * rnum])
      matrix(cs - rep(off, each = rnum), rnum)
    }
    one_es <- function(expv, w) {
      ae <- abs(expv)
      fg <- col_cumsum(w * ae)
      bg <- col_cumsum((1 - w) * ae)
      d <- fg / rep(fg[rnum, ], each = rnum) - bg / rep(bg[rnum, ], each = rnum)
      list(pos = pmax(apply(d, 2L, max), 0), neg = pmin(apply(d, 2L, min), 0))
    }
    es <- matrix(0, cnum, ncol(reg))
    for (k in seq_len(cnum)) {
      o <- order(data[, k], decreasing = TRUE)
      r <- one_es(data[o, k], reg[o, , drop = FALSE])
      es[k, ] <- ifelse(r$pos > abs(r$neg), r$pos, r$neg)
    }
    pos_es <- neg_es <- matrix(0, ncol(reg), perm)
    for (k in seq_len(perm)) {
      se <- sample(1:cnum, 1)
      cur <- sample(data[, se])
      r <- one_es(cur, reg)
      pos_es[, k] <- r$pos
      neg_es[, k] <- r$neg
    }
    pavg <- apply(pos_es, 1L, mean)
    navg <- abs(apply(neg_es, 1L, mean))
    for (k in seq_len(nrow(es))) {
      tmp <- es[k, ]
      es[k, ] <- ifelse(tmp > 0, tmp / pavg, tmp / navg)
    }
    colnames(es) <- paste0(colnames(reg), ".ES")
    rownames(es) <- colnames(data)
    half <- ncol(es) / 2
    up <- es[, seq_len(half), drop = FALSE]
    dn <- es[, half + seq_len(half), drop = FALSE]
    scores <- up - dn
    colnames(scores) <- make.names(gsub("_UP\\.ES", "", colnames(scores)), allow_ = FALSE)
    scores <- scores[, is.finite(colSums(scores)), drop = FALSE]
    z <- scale(scores)
    subtypes <- unique(mapping$Subtype)
    out <- list()
    for (s in subtypes) {
      cols <- colnames(z) %in% mapping[mapping$Subtype == s, 1L]
      if (any(cols)) out[[s]] <- apply(z[, cols, drop = FALSE], 1L, stats::median)
    }
    res <- do.call(cbind, out)
    rownames(res) <- rownames(scores)
    res
  },

  # TIMER (Li et al. 2016, doi:10.1186/s13059-016-1028-7): remove the batch
  # effect between the samples and a curated immune-cell reference (ComBat),
  # take the median reference profile per cell type, and fit cancer-type
  # specific marker genes by constrained regression (Abbas et al. 2009).
  timer_fit = function(X, reference, markers) {
    expr <- X[nchar(rownames(X)) > 1L, , drop = FALSE]
    outlier <- sort(unique(unlist(lapply(seq_len(ncol(expr)), function(j) {
      rownames(expr)[utils::tail(order(expr[, j]), 5L)]
    }))))
    expr <- expr[!rownames(expr) %in% outlier, , drop = FALSE]
    imm <- reference$genes
    cell_types <- reference$celltypes
    shared <- intersect(rownames(expr), rownames(imm))
    n1 <- ncol(expr)
    n2 <- ncol(imm)
    dd <- cbind(expr[shared, , drop = FALSE], imm[shared, , drop = FALSE])
    mode(dd) <- "numeric"
    batch <- c(rep(1, n1), rep(2, n2))
    dd0 <- sva::ComBat(dd, batch, NULL, BPPARAM = BiocParallel::bpparam("SerialParam"))
    expr_n <- dd0[, seq_len(n1), drop = FALSE]
    imm_n <- dd0[, n1 + seq_len(n2), drop = FALSE]
    types <- unique(names(cell_types))
    agg <- vapply(types, function(k) apply(imm_n[, names(cell_types) == k, drop = FALSE], 1L,
                                           stats::median, na.rm = TRUE), numeric(nrow(imm_n)))
    rownames(agg) <- rownames(imm_n)
    sel <- intersect(markers, rownames(expr_n))
    XX <- agg[sel, -4L, drop = FALSE]  # the fourth reference type (NK) is not reported
    YY <- expr_n[sel, , drop = FALSE]
    abbas <- function(XX, y) {
      remove <- integer()
      names_kept <- colnames(XX)
      cur <- XX
      repeat {
        if (length(remove)) {
          if (is.null(ncol(cur))) return(structure(rep(0, ncol(XX)), names = colnames(XX)))
          cur <- cur[, -remove, drop = FALSE]
          names_kept <- names_kept[-remove]
          if (length(names_kept) == 0L) return(structure(rep(0, ncol(XX)), names = colnames(XX)))
        }
        fit <- stats::lsfit(cur, y, intercept = FALSE)
        beta <- fit$coefficients[seq_len(ncol(cur))]
        if (min(beta > 0)) break
        remove <- which.min(beta)
      }
      out <- structure(rep(0, ncol(XX)), names = colnames(XX))
      out[names_kept] <- beta
      out
    }
    res <- vapply(seq_len(ncol(YY)), function(j) abbas(XX, YY[, j]), numeric(ncol(XX)))
    res <- matrix(res, nrow = ncol(XX), dimnames = list(colnames(XX), colnames(YY)))
    t(res)
  }
)

# ---- ImmuCellAI 2.0 (MIT; Yang J, github.com/GuoBioinfoLab/ImmuCellAI-2.0) ----
bd_native$immucellai2_fit <- function(X, reference, markers, tissue = "tumor",
                                      pseudo_depth = 1e5, n_iter = 50L, tol = 1e-4,
                                      alpha_major = 10, alpha_sub = 5, alpha_state = 1) {
  # Variational Bayes deconvolution on a 53-state reference with a T-cell
  # hierarchy (lineage > subtype > state). Port of run_immucellai2() with its
  # default settings (hierarchy.mode = "tcell", no UNKNOWN state). The
  # reference matrix passed in holds only the marker genes (5510); the
  # original reads a 41 MB text file and then keeps the same rows.
  normalize_prob <- function(x) {
    x[is.na(x)] <- 0
    x[x < 0] <- 0
    s <- sum(x)
    if (s <= 0) return(rep(1 / length(x), length(x)))
    x / s
  }
  if (tissue == "non-tumor") reference <- reference[, !colnames(reference) %in% c("TAM", "MDSC"), drop = FALSE]
  bulk <- as.matrix(X)
  mode(bulk) <- "numeric"
  bulk[is.na(bulk)] <- 0
  bulk[bulk < 0] <- 0
  common <- intersect(markers, intersect(rownames(reference), rownames(bulk)))
  if (length(common) == 0L) stop("no marker genes overlap the data; check the gene symbols", call. = FALSE)
  if (length(common) < 1000L) warning("only ", length(common), " marker genes overlap the data", call. = FALSE)
  bulk <- bulk[common, , drop = FALSE]
  ref <- reference[common, , drop = FALSE]
  # hierarchy of the 53 states (tcell mode)
  st <- colnames(ref)
  hier <- data.frame(state_name = st, subtype = st, major_lineage = st, stringsAsFactors = FALSE)
  cd4_memory <- c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm")
  cd4_naive <- c("CD4Tn", "CD4Tnaive")
  cd4_helper <- c("Th1", "Th1/17", "Th1/Th17", "Th1.Th17", "Th1_Th17", "Th17", "Th2")
  cd4_special <- c("Tfh", "CD4Tfh", "Treg", "Tr1")
  cd8_memory <- c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm")
  cd8_naive <- c("CD8Tn", "CD8Tnaive")
  cd8_effector <- c("Tc", "Tex", "exhausted_T")
  hier$major_lineage[st %in% c(cd4_memory, cd4_naive, cd4_helper, cd4_special)] <- "CD4T"
  hier$subtype[st %in% cd4_memory] <- "CD4 Memory T"
  hier$subtype[st %in% cd4_naive] <- "CD4 Naive T"
  hier$subtype[st %in% cd4_helper] <- "CD4 Helper T"
  hier$subtype[st %in% c("Tfh", "CD4Tfh")] <- "Tfh"
  hier$major_lineage[st %in% c(cd8_memory, cd8_naive, cd8_effector)] <- "CD8T"
  hier$subtype[st %in% cd8_memory] <- "CD8 Memory T"
  hier$subtype[st %in% cd8_naive] <- "CD8 Naive T"
  hier$subtype[st %in% "Tc"] <- "Tc"
  hier$subtype[st %in% c("Tex", "exhausted_T")] <- "Tex"
  state_labels <- hier$state_name
  subtype_labels <- hier$subtype
  lineage_labels <- hier$major_lineage
  keep <- rowSums(bulk) > 0 & rowSums(ref) > 0
  bulk <- bulk[keep, , drop = FALSE]
  ref <- ref[keep, , drop = FALSE]
  phi <- ref + 1e-8
  phi <- sweep(phi, 2L, colSums(phi), "/")
  lineages <- unique(lineage_labels)
  subtypes <- unique(subtype_labels)
  # prior over states: uniform over lineages, then subtypes, then states
  prior <- rep(0, length(state_labels))
  names(prior) <- state_labels
  major_prior <- rep(1 / length(lineages), length(lineages))
  names(major_prior) <- lineages
  for (lin in lineages) {
    idx_lin <- which(lineage_labels == lin)
    sbts <- unique(subtype_labels[idx_lin])
    sub_prior <- rep(1 / length(sbts), length(sbts))
    names(sub_prior) <- sbts
    for (sbt in sbts) {
      idx_sbt <- idx_lin[subtype_labels[idx_lin] == sbt]
      prior[idx_sbt] <- major_prior[lin] * sub_prior[sbt] / length(idx_sbt)
    }
  }
  prior <- normalize_prob(prior)
  fit_one <- function(bulk_vec) {
    counts_in <- round(bulk_vec / sum(bulk_vec) * pseudo_depth)
    theta <- prior
    for (iter in seq_len(n_iter)) {
      weighted <- sweep(phi, 2L, normalize_prob(theta[colnames(phi)]), "*")
      denom <- rowSums(weighted)
      denom[is.na(denom) | denom <= 0] <- 1e-300
      sc <- colSums((weighted / denom) * counts_in)
      names(sc) <- colnames(phi)
      sc <- sc[state_labels]
      sc[is.na(sc)] <- 0
      lc <- vapply(lineages, function(l) sum(sc[lineage_labels == l]), 0)
      pm <- normalize_prob(vapply(lineages, function(l) sum(prior[lineage_labels == l]), 0))
      th_major <- normalize_prob(alpha_major * pm + lc)
      names(th_major) <- lineages
      th_sub <- list()
      for (lin in lineages) {
        sbts <- unique(subtype_labels[lineage_labels == lin])
        sc_sub <- vapply(sbts, function(b) sum(sc[subtype_labels == b]), 0)
        ps <- normalize_prob(vapply(sbts, function(b) sum(prior[subtype_labels == b]), 0))
        th_sub[[lin]] <- normalize_prob(alpha_sub * ps + sc_sub)
        names(th_sub[[lin]]) <- sbts
      }
      th_state <- list()
      for (sbt in subtypes) {
        idx <- which(subtype_labels == sbt)
        th_state[[sbt]] <- normalize_prob(alpha_state * normalize_prob(prior[idx]) + sc[idx])
        names(th_state[[sbt]]) <- state_labels[idx]
      }
      new <- rep(0, length(state_labels))
      names(new) <- state_labels
      for (i in seq_along(state_labels)) {
        new[state_labels[i]] <- th_major[lineage_labels[i]] * th_sub[[lineage_labels[i]]][subtype_labels[i]] *
          th_state[[subtype_labels[i]]][state_labels[i]]
      }
      new <- normalize_prob(new)
      delta <- max(abs(new[state_labels] - theta[state_labels]))
      theta <- new[state_labels]
      if (delta < tol) break
    }
    theta
  }
  out <- list()
  for (j in colnames(bulk)) {
    v <- bulk[, j]
    if (sum(v) <= 0) next
    out[[j]] <- fit_one(v)
  }
  if (!length(out)) stop("no sample has expression in the marker genes", call. = FALSE)
  do.call(rbind, out)
}

# ---- SpatialDecon (MIT; Danaher P, et al. Nat Commun 2022, doi:10.1038/s41467-022-28020-5) ----
bd_native$spatialdecon_fit <- function(X, signature, background = 0.1, resid_thresh = 3,
                                       lower_thresh = 0.5, maxit = 1000L) {
  # Log-normal regression of each sample on the cell-type profiles (here the
  # SAFE-TME matrix), run twice: the second time without the data points whose
  # log2 residual in the first fit exceeds `resid_thresh`. Port of
  # spatialdecon() with its defaults and a constant background; the original
  # fits with logNormReg::lognlm(), whose objective and gradient are written
  # out below using the same L-BFGS-B configuration.
  # The standard errors of the original are not needed and not computed.
  norm <- as.matrix(X)
  sig <- as.matrix(signature)
  sig <- sig * 2 / stats::quantile(sig, 0.99)
  # The upstream refit orders genes by the reference, independently of input order.
  shared <- intersect(rownames(sig), rownames(norm))
  if (length(shared) < 100L) stop("only ", length(shared), " genes are shared with the signature; at least 100 are needed", call. = FALSE)
  Y <- norm[shared, , drop = FALSE]
  Xs <- sig[shared, , drop = FALSE]
  bg <- matrix(background, nrow(Y), ncol(Y), dimnames = dimnames(Y))
  epsilon <- min(Y[(Y > 0) & !is.na(Y)])
  K <- ncol(Xs)
  fit_all <- function(Y, eps, stage) {
    fits <- lapply(seq_len(ncol(Y)), function(j) {
      y <- Y[, j]
      b <- bg[, j]
      use <- !is.na(y)
      y0 <- y[use]
      y <- pmax(y0, eps)
      Xd <- cbind(b[use], Xs[use, , drop = FALSE])
      ly <- log(y)
      dmin <- function(par) {
        mu <- pmax(drop(Xd %*% par), 1e-04)
        sum((ly - log(mu))^2)
      }
      grad <- function(par) {
        mu <- pmax(drop(Xd %*% par), 1e-04)
        -drop(crossprod(Xd, (ly - log(mu)) / mu))
      }
      init <- rep(mean(y0) / (mean(Xs) * K), K)
      stats::optim(c(1, init), dmin, grad, method = "L-BFGS-B", lower = c(1, rep(0, K)),
                   upper = c(1, rep(Inf, K)), control = list(maxit = maxit))
    })
    failed <- vapply(fits, function(fit) fit$convergence != 0L, logical(1L))
    if (any(failed)) {
      warning("SpatialDecon ", stage, " did not converge for sample(s): ",
              paste(colnames(Y)[failed], collapse = ", "), ".", call. = FALSE)
    }
    beta <- vapply(fits, function(fit) fit$par[-1L], numeric(K))
    beta <- matrix(beta, nrow = K, dimnames = list(colnames(Xs), colnames(Y)))
    pmax(beta, 0)
  }
  b0 <- fit_all(Y, epsilon, "initial fit")
  yhat <- Xs %*% b0 + bg
  resids <- log2(pmax(Y, lower_thresh)) - log2(pmax(yhat, lower_thresh))
  Y2 <- replace(Y, abs(resids) > resid_thresh, NA)
  t(fit_all(Y2, epsilon, "outlier refit"))
}
