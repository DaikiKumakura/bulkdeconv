# Natively implemented signature-based methods added after the first release:
# GEDIT and the Danaher cell-type scores (the Charoentong gene sets are scored
# with ssgsea() in native2.R). Same conventions as native.R: X is genes x
# samples, results are samples x features, each function self-contained.

# ---- GEDIT ----
# Port of GEDIT 3.0 (Nadel et al. 2021; GEDIT3.py and scripts/ at BNadel/GEDIT
# be9b33f) with its defaults: 50 signature genes per cell type chosen by
# negative entropy, no row scaling power (each gene rescaled to [0, 1] over
# reference and mixture together), non-negative least squares per sample and
# the coefficients divided by their sum. As in the original, the mixture is
# quantile normalised to the overall value distribution of the reference, the
# reference itself is not, and gene names are matched case-insensitively.
bd_native$gedit_fit <- function(X, ref, sigs_per_type = 50L) {
  X <- as.matrix(X)
  ref <- ref[rowSums(ref) > 0, , drop = FALSE]
  target <- sort(as.vector(ref))
  qnorm_col <- function(x) {
    # value at the same quantile of the target (linear interpolation), ranks
    # broken by position as numpy's argsort does for distinct values
    o <- order(x)
    pos <- (seq_along(x) - 1) / (length(x) - 1) * (length(target) - 1)
    lo <- floor(pos)
    val <- target[lo + 1] + (pos - lo) * (target[pmin(lo + 2, length(target))] - target[lo + 1])
    out <- numeric(length(x)); out[o] <- val
    out
  }
  mix <- apply(X, 2L, qnorm_col)
  if (!is.matrix(mix)) mix <- matrix(mix, ncol = ncol(X))
  dimnames(mix) <- dimnames(X)
  up <- function(m) { u <- toupper(rownames(m)); m <- m[!duplicated(u, fromLast = TRUE), , drop = FALSE]; rownames(m) <- toupper(rownames(m)); m }
  mix <- up(mix); ref <- up(ref)
  shared <- sort(intersect(rownames(mix), rownames(ref)), method = "radix")
  if (length(shared) < ncol(ref)) stop("only ", length(shared), " genes of the GEDIT reference are in the data", call. = FALSE)
  r <- ref[shared, , drop = FALSE]
  # signature genes: genes with at most half (and fewer than all) zero entries;
  # zeros replaced by the smallest non-zero value; score = negative entropy
  n_zero <- rowSums(r == 0)
  r <- r[n_zero <= min(0.5 * ncol(r), ncol(r) - 1), , drop = FALSE]
  r[r == 0] <- min(r[r > 0])
  p <- r / rowSums(r)
  score <- rowSums(p * log(p))
  top <- colnames(r)[max.col(r, ties.method = "first")]
  pick <- unlist(lapply(split(seq_along(score), top), function(i) {
    o <- i[order(score[i], rownames(r)[i], method = "radix")]  # ascending (score, name) as Python sorts tuples; the best are last
    utils::tail(o, sigs_per_type)
  }), use.names = FALSE)
  total <- sigs_per_type * ncol(ref)
  if (length(pick) < total) {
    rest <- setdiff(seq_along(score), pick)
    rest <- rest[order(score[rest], rownames(r)[rest], method = "radix")]
    pick <- c(pick, utils::tail(rest, total - length(pick)))
  }
  genes <- sort(rownames(r)[pick], method = "radix")
  sref <- ref[genes, , drop = FALSE]; smix <- mix[genes, , drop = FALSE]
  both <- cbind(sref, smix)
  lo <- apply(both, 1L, min); hi <- apply(both, 1L, max)
  scaled <- (both - lo) / (hi - lo)
  scaled[!is.finite(scaled)] <- 0
  A <- scaled[, seq_len(ncol(sref)), drop = FALSE]
  B <- scaled[, -seq_len(ncol(sref)), drop = FALSE]
  est <- vapply(seq_len(ncol(B)), function(j) {
    f <- nnls(A, B[, j])
    f / sum(f)
  }, numeric(ncol(A)))
  est <- matrix(est, nrow = ncol(A), dimnames = list(colnames(ref), colnames(X)))
  t(est)
}

# ---- Danaher cell-type scores ----
# Danaher et al. 2017: each population's score is the mean log2 expression of
# its marker genes (the paper's TCGA analysis used log2 RSEM values; here
# log2(TPM + 1)). Genes absent from the data are skipped; a population with no
# marker present gives NA.
bd_native$danaher_fit <- function(X, markers) {
  L <- log2(as.matrix(X) + 1)
  rownames(L) <- toupper(rownames(L))
  s <- vapply(markers, function(g) {
    g <- intersect(toupper(g), rownames(L))
    if (!length(g)) return(rep(NA_real_, ncol(L)))
    colMeans(L[g, , drop = FALSE])
  }, numeric(ncol(L)))
  s <- matrix(s, nrow = ncol(L), dimnames = list(colnames(L), names(markers)))
  s
}

# ---- PUREE tumour purity ----
# Port of PUREE 1.0 (Revkov et al. 2023; predict_purity.py at skandlab/PUREE
# 5a0a702) with HGNC symbols as input: symbols are mapped to the Ensembl IDs of
# the ranking universe, each sample is ranked over the universe genes present
# (percentile rank, ties to the minimum, as pandas rank(method = "min",
# pct = TRUE)), the 158 model genes are taken, missing ones get the training
# medians, and the linear model gives the purity. `d` is the RDS made by
# bulkdeconv_install("puree").
bd_native$puree_fit <- function(X, d) {
  X <- as.matrix(X)
  U <- matrix(NA_real_, length(d$universe), ncol(X), dimnames = list(d$universe, colnames(X)))
  if (mean(grepl("^ENSG[0-9]+", rownames(X))) > 0.5) {
    # Ensembl gene IDs (version suffixes dropped), as PUREE's default input
    ids <- sub("\\.[0-9]+$", "", rownames(X))
    keep <- ids %in% d$universe & !duplicated(ids)
    U[ids[keep], ] <- X[keep, , drop = FALSE]
  } else {
    conv <- d$conversion[d$conversion$ENSEMBL_ID %in% d$universe, , drop = FALSE]
    conv <- conv[!duplicated(conv), , drop = FALSE]
    conv <- conv[conv$HGNC %in% rownames(X), , drop = FALSE]
    conv <- conv[!duplicated(conv$ENSEMBL_ID), , drop = FALSE]
    U[conv$ENSEMBL_ID, ] <- X[conv$HGNC, , drop = FALSE]
  }
  R <- apply(U, 2L, function(v) {
    ok <- !is.na(v)
    r <- rep(NA_real_, length(v))
    r[ok] <- rank(v[ok], ties.method = "min") / sum(ok)
    r
  })
  if (!is.matrix(R)) R <- matrix(R, ncol = ncol(X))
  rownames(R) <- d$universe
  F <- R[d$genes, , drop = FALSE]
  miss <- is.na(F)
  F[miss] <- matrix(d$medians, nrow(F), ncol(F))[miss]
  purity <- d$intercept + drop(crossprod(F, d$coef))
  matrix(purity, ncol = 1L, dimnames = list(colnames(X), "tumor_purity"))
}
