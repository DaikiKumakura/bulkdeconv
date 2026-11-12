# xCell2 1.2.3 inference with pre-trained references. See inst/COPYRIGHTS.
bd_native$xcell2_fit <- function(X, ref) {
  shared <- intersect(rownames(X), ref$genes)
  fraction <- round(length(shared) / length(ref$genes), 2)
  if (fraction < 0.85) {
    warning("This xCell2 reference shares only ", fraction,
            " of its genes with the mixture; at least 85% is recommended.", call. = FALSE)
  }

  ranks <- vapply(seq_len(ncol(X)), function(j) {
    rank(X[shared, j], ties.method = "average")
  }, numeric(length(shared)))
  dimnames(ranks) <- list(shared, colnames(X))
  signatures <- ref$signatures
  matched <- utils::relist(
    match(unlist(signatures, use.names = FALSE), shared), signatures
  )
  signature_types <- sub("#.*$", "", names(signatures))
  types <- unique(signature_types)
  raw <- vapply(types, function(type) {
    is_type <- signature_types == type
    sets <- lapply(matched[is_type], function(genes) unique(genes[!is.na(genes)]))
    found <- lengths(sets)
    if (sum(found > 0L) < 3L) {
      warning("Cannot calculate enrichment for '", type,
              "': fewer than three signatures overlap the mixture.", call. = FALSE)
      return(numeric(ncol(X)))
    }
    if (any(found != lengths(signatures[is_type]))) {
      warning("Some signature genes are missing for '", type, "'.", call. = FALSE)
    }
    scores <- vapply(sets[found > 0L], function(genes) {
      k <- length(genes)
      (colMeans(ranks[genes, , drop = FALSE]) - (k + 1) / 2) / (nrow(ranks) - k)
    }, numeric(ncol(X)))
    rowMeans(scores)
  }, numeric(ncol(X)))
  rownames(raw) <- colnames(X)
  if (any(!is.finite(raw) | raw < 0)) {
    stop("The mixture produced non-finite or negative xCell2 enrichment scores.", call. = FALSE)
  }

  for (type in types) {
    p <- ref$params[match(type, ref$params$celltype), , drop = FALSE]
    score <- (raw[, type]^(1 / p$b)) / p$a * p$m
    raw[, type] <- round(score - min(score), 5)
  }
  spill <- ref$spill * 0.5
  diag(spill) <- 1
  cells <- intersect(types, rownames(spill))
  spill <- spill[cells, cells, drop = FALSE]
  score <- t(vapply(seq_len(nrow(raw)), function(i) {
    nnls(spill, raw[i, cells])
  }, numeric(length(cells))))
  dimnames(score) <- list(rownames(raw), cells)
  score
}
