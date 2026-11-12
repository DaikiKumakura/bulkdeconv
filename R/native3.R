# More natively implemented methods (see native.R, native2.R): consensus_tme,
# xcell, mcp_counter, mmcp_counter. Same conventions: X is genes x samples,
# results samples x features, each function self-contained (they call the
# helpers ssgsea() and nnls() of native2.R by name).

# ---- ConsensusTME (GPL-3; Jimenez-Sanchez A, et al. Cancer Res 2019, doi:10.1158/0008-5472.CAN-18-3711) ----
bd_native$consensus_tme_fit <- function(X, gene_sets, cancer_type = "Unfiltered") {
  # ssGSEA of the consensus gene sets of one cancer type (the package ships
  # them precomputed for all TCGA types and "Unfiltered"). Port of
  # consensusTMEAnalysis(statMethod = "ssgsea").
  sets <- gene_sets[[cancer_type]]
  if (is.null(sets)) stop("unknown cancer type '", cancer_type, "'", call. = FALSE)
  # Immune_Score: all genes of the immune signatures (all signatures when unfiltered)
  immune <- if (cancer_type == "Unfiltered") sets else sets[!names(sets) %in% c("Endothelial", "Fibroblasts")]
  sets$Immune_Score <- sort(unique(unlist(immune)))
  X <- as.matrix(X)
  if (nrow(X) < 1000L) warning("low number of genes; the enrichment may not be reliable", call. = FALSE)
  sets <- sets[lengths(sets) > 0L]  # signatures left empty by the purity filter are dropped
  n_in <- vapply(sets, function(g) length(intersect(rownames(X), g)), 0L)
  if (all(n_in == 0L)) stop("no gene of the gene sets matches the gene symbols of the data", call. = FALSE)
  if (any(n_in == 0L)) {
    warning("these signatures were not run (no gene matched): ", paste(names(sets)[n_in == 0L], collapse = ", "), call. = FALSE)
    sets <- sets[n_in > 0L]
  }
  t(ssgsea(X, sets))
}

# ---- xCell (GPL-3; Aran D, Hu Z, Butte AJ. Genome Biol 2017, doi:10.1186/s13059-017-1349-1) ----
bd_native$xcell_fit <- function(X, signatures, genes, fit_values, spill) {
  # 489 gene signatures scored by ssGSEA on gene ranks, averaged per cell type,
  # shifted and transformed with fitted parameters, then corrected for
  # spill-over between cell types by non-negative least squares (alpha = 0.5);
  # the immune, stroma and microenvironment scores are added. Port of
  # xCellAnalysis(rnaseq = TRUE) with its defaults.
  X <- as.matrix(X)
  shared <- intersect(rownames(X), genes)
  if (length(shared) < 5000L) stop("only ", length(shared), " genes are shared with xCell's gene list; at least 5000 are needed", call. = FALSE)
  ranked <- apply(X[shared, , drop = FALSE], 2L, rank)
  scores <- ssgsea(ranked, signatures, normalize = FALSE)
  scores <- scores - apply(scores, 1L, min)
  cell_type <- vapply(strsplit(rownames(scores), "%", fixed = TRUE), `[`, "", 1L)
  types <- sort(unique(cell_type))
  agg <- t(matrix(vapply(types, function(ct) colMeans(scores[cell_type == ct, , drop = FALSE]), numeric(ncol(scores))),
                  nrow = ncol(scores)))
  rownames(agg) <- types
  colnames(agg) <- colnames(scores)
  rows <- rownames(agg)[rownames(agg) %in% rownames(fit_values)]
  ts <- agg[rows, , drop = FALSE]
  ts <- (ts - apply(ts, 1L, min)) / 5000
  ts[ts < 0] <- 0
  ts <- (ts^fit_values[rows, 2]) / (fit_values[rows, 3] * 2)
  K <- spill * 0.5
  diag(K) <- 1
  rows2 <- rownames(ts)[rownames(ts) %in% rownames(K)]
  adj <- apply(ts[rows2, , drop = FALSE], 2L, function(x) nnls(K[rows2, rows2, drop = FALSE], x))
  adj[adj < 0] <- 0
  rownames(adj) <- rows2
  immune <- colSums(adj[c("B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Eosinophils", "Macrophages",
                          "Monocytes", "Mast cells", "Neutrophils", "NK cells"), , drop = FALSE]) / 1.5
  stroma <- colSums(adj[c("Adipocytes", "Endothelial cells", "Fibroblasts"), , drop = FALSE]) / 2
  t(rbind(adj, ImmuneScore = immune, StromaScore = stroma, MicroenvironmentScore = immune + stroma))
}

# ---- MCP-counter (GPL-3; Becht E, et al. Genome Biol 2016, doi:10.1186/s13059-016-1070-5) ----
bd_native$mcp_counter_fit <- function(X, markers) {
  # Score of a cell population = mean expression of its marker genes (HUGO
  # symbols) in the sample. `markers` has columns "HUGO symbols" and
  # "Cell population". Port of MCPcounter.estimate(featuresType = "HUGO_symbols").
  X <- as.matrix(X)
  present <- markers[markers[["HUGO symbols"]] %in% rownames(X), , drop = FALSE]
  pops <- unique(present[["Cell population"]])
  missing <- setdiff(unique(markers[["Cell population"]]), pops)
  if (length(missing)) warning("no markers found for: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!length(pops)) stop("no marker gene is in the data; check the gene symbols", call. = FALSE)
  res <- vapply(pops, function(p) {
    g <- intersect(rownames(X), present[["HUGO symbols"]][present[["Cell population"]] == p])
    colMeans(X[g, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(X)))
  matrix(res, nrow = ncol(X), dimnames = list(colnames(X), pops))
}

# ---- mMCP-counter (GPL-3; Petitprez F, et al. Genome Med 2020, doi:10.1186/s13073-020-00783-w) ----
bd_native$mmcp_counter_fit <- function(X, signatures) {
  # Score of a cell population = median expression of its marker genes (mouse
  # symbols, genome GCRm39). Port of mMCPcounter.estimate(features = "Gene.Symbol").
  X <- as.matrix(X)
  found <- intersect(signatures[["Gene.Symbol"]], rownames(X))
  if (!length(found)) stop("no signature gene is in the data; check the gene symbols", call. = FALSE)
  absent <- setdiff(unique(signatures$Denomination),
                    unique(signatures$Denomination[signatures[["Gene.Symbol"]] %in% rownames(X)]))
  if (length(absent)) warning("no genes found for: ", paste(absent, collapse = ", "), call. = FALSE)
  loc <- signatures[signatures[["Gene.Symbol"]] %in% found, , drop = FALSE]
  groups <- split(seq_len(nrow(loc)), loc$Denomination)
  agg <- t(matrix(vapply(groups, function(i) {
    m <- X[loc[["Gene.Symbol"]][i], , drop = FALSE]
    apply(m, 2L, stats::median, na.rm = TRUE)
  }, numeric(ncol(X))), nrow = ncol(X)))
  rownames(agg) <- names(groups)
  colnames(agg) <- colnames(X)
  want <- c("T cells", "CD8 T cells", "NK cells", "B derived", "Memory B cells", "Monocytes / macrophages",
            "Monocytes", "Granulocytes", "Mast cells", "Eosinophils", "Neutrophils", "Basophils", "Vessels",
            "Lymphatics", "Endothelial cells", "Fibroblasts")
  agg <- agg[match(want, rownames(agg)), , drop = FALSE]
  rownames(agg) <- want
  agg <- agg[rowSums(is.na(agg)) < ncol(agg), , drop = FALSE]
  t(agg)
}

# ---- ImSig (GPL-3; Nirmal AJ, et al. Cancer Immunol Res 2018, doi:10.1158/2326-6066.CIR-17-0568) ----
bd_native$imsig_fit <- function(X, sig, r = 0.6) {
  # Keep the signature genes of each cell type whose expression correlates
  # (Pearson, above r) with at least one other gene of the same signature, and
  # score a cell type as the mean expression of its retained genes.
  # Port of imsig(exp, r) without the sorting of the samples.
  X <- as.matrix(X)
  g <- intersect(rownames(X), as.character(sig$gene))
  if (length(g) < 10L) stop("only ", length(g), " ImSig genes are in the data", call. = FALSE)
  e <- X[g, , drop = FALSE]
  sg <- sig[as.character(sig$gene) %in% g, , drop = FALSE]
  keep <- character()
  for (ct in levels(sg$cell)) {
    gs <- as.character(sg$gene[sg$cell %in% ct])
    if (length(gs) < 2L) next
    cm <- suppressWarnings(stats::cor(t(e[gs, , drop = FALSE])))
    cm[cm <= r] <- 0
    diag(cm) <- 0
    keep <- c(keep, gs[apply(cm, 1L, function(x) !all(x == 0))])
  }
  sg_keep <- sg[as.character(sg$gene) %in% keep, , drop = FALSE]
  avg <- e[rownames(e) %in% keep, , drop = FALSE]
  out <- vapply(levels(sg_keep$cell), function(ct) {
    gs <- as.character(sg_keep$gene[sg_keep$cell %in% ct])
    colMeans(avg[gs, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(X)))
  matrix(out, nrow = ncol(X), dimnames = list(colnames(X), levels(sg_keep$cell)))
}

# ---- quanTIseq (GPL-3; Finotello F, et al. Genome Med 2019, doi:10.1186/s13073-019-0638-6; quantiseqr, Federico A, Ponzoni L, Marini F) ----
bd_native$quantiseq_fit <- function(X, ref, tumor = TRUE, scale_mrna = TRUE) {
  # 10 immune cell types from the TIL10 signature by constrained least squares
  # (fractions >= 0, summing to at most 1; the rest is "Other"), with
  # re-annotation of gene symbols, TPM normalisation, removal of noisy genes and
  # of genes high in tumours, mRNA-content scaling and the Treg/CD4 correction.
  # Port of run_quantiseq(signature_matrix = "TIL10", is_arraydata = FALSE).
  mix <- as.matrix(X)
  hg <- ref$hgnc
  cur <- rownames(mix)
  ind <- match(cur, hg$ApprovedSymbol)
  newg <- rep(NA_character_, length(cur))
  newg2 <- rep(NA_character_, length(cur))
  ok <- which(!is.na(ind))
  st <- as.character(hg$Status[ind[ok]])
  approved <- ok[st == "Approved"]
  newg[approved] <- cur[approved]
  renamed <- ok[!st %in% c("Approved", "EntryWithdrawn")]
  newg2[renamed] <- gsub("symbolwithdrawn,see", "", as.character(hg$ApprovedName[ind[renamed]]), fixed = TRUE)
  # genes not found as approved symbols: look them up among the previous symbols.
  # (The original also looks at a "Synonyms" column, which the HGNC table it
  # ships does not have, so that step never finds anything; it is left out.)
  tok <- function(col) {
    parts <- strsplit(as.character(hg[[col]]), ",", fixed = TRUE)
    data.frame(token = unlist(parts), approved = rep(as.character(hg$ApprovedSymbol), lengths(parts)), stringsAsFactors = FALSE)
  }
  prev <- tok("PreviousSymbols")
  last_match <- function(tab, genes) {
    idx <- match(genes, rev(tab$token))
    tab_rev <- rev(tab$approved)
    ifelse(is.na(idx), NA_character_, tab_rev[idx])
  }
  miss <- which(is.na(ind))
  newg2[miss] <- last_match(prev, cur[miss])
  newg2[newg2 %in% setdiff(newg, NA)] <- NA
  fill <- intersect(which(is.na(newg)), which(!is.na(newg2)))
  newg[fill] <- newg2[fill]
  has <- which(!is.na(newg))
  m <- mix[has, , drop = FALSE]
  newg <- newg[has]
  if (anyDuplicated(newg)) {
    grp <- factor(newg, levels = sort(unique(newg)))
    m <- t(vapply(levels(grp), function(gn) apply(m[grp == gn, , drop = FALSE], 2L, stats::median), numeric(ncol(m))))
    if (ncol(mix) == 1L) m <- matrix(m, ncol = 1L)
    rownames(m) <- levels(grp)
  } else {
    rownames(m) <- newg
  }
  colnames(m) <- colnames(mix)
  if (max(m) < 50) m <- 2^m
  m <- t(t(m) * 1e+06 / apply(m, 2L, sum))
  sig <- as.matrix(ref$signature)
  sig <- sig[!rownames(sig) %in% ref$rmgenes, , drop = FALSE]
  if (tumor) sig <- sig[!rownames(sig) %in% ref$abgenes, , drop = FALSE]
  if (length(intersect(rownames(sig), rownames(m))) < 20L) stop("too few TIL10 signature genes are in the data", call. = FALSE)
  mrna <- if (scale_mrna) ref$scaling[match(colnames(sig), names(ref$scaling))] else rep(1, ncol(sig))
  solve_all <- function(sg, mx, scaling) {
    cg <- intersect(rownames(sg), rownames(mx))
    A <- sg[cg, , drop = FALSE]
    B <- mx[cg, , drop = FALSE]
    sc <- norm(A, "2")
    A <- A / sc
    k <- ncol(A)
    Dmat <- crossprod(A)
    Amat <- cbind(diag(k), rep(-1, k))
    res <- t(apply(B, 2L, function(b) {
      b <- b / sc
      est <- quadprog::solve.QP(Dmat, drop(crossprod(A, b)), Amat, c(rep(0, k), -1))$solution
      est.sum <- sum(est)
      est <- est / scaling
      est <- est / sum(est) * est.sum
      c(est, pmax(0, 1 - sum(est)))
    }))
    if (ncol(B) == 1L) res <- matrix(res, nrow = 1L)
    dimnames(res) <- list(colnames(B), c(colnames(sg), "Other"))
    res
  }
  r1 <- solve_all(sig, m, mrna)
  if ("Tregs" %in% colnames(sig) && "T.cells.CD4" %in% colnames(sig)) {
    i <- which(colnames(sig) == "T.cells.CD4")
    r2 <- solve_all(sig[, -i, drop = FALSE], m, mrna[-i])
    low <- which(r1[, "Tregs"] < 0.02)
    if (length(low)) {
      r1[low, "Tregs"] <- (r2[low, "Tregs"] + r1[low, "Tregs"]) / 2
      r1[low, "T.cells.CD4"] <- pmax(0, r1[low, "T.cells.CD4"] - (r2[low, "Tregs"] + r1[low, "Tregs"]) / 2)
    }
  }
  r1 / rowSums(r1)
}
