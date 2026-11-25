#' Consensus across methods
#'
#' Combines the results of several methods run on the same samples into one
#' score per cell lineage (B cell, CD4 T, CD8 T, NK, monocyte, macrophage,
#' neutrophil, dendritic, fibroblast, endothelial) and one for tumor content
#' (PUREE purity and the inverted ESTIMATE score).
#'
#' Methods report different quantities (fractions, abundance scores,
#' enrichment scores), so their numbers cannot be averaged. Each mapped output
#' column is therefore replaced by its percentile rank among the samples of the
#' run, `(rank - 0.5) / n`, inverted where a higher method score means less of
#' the lineage. The ranks of the columns that belong to one lineage are
#' averaged within a method, and the method values are averaged across
#' methods. The mapping from method outputs to lineages is in
#' [bulkdeconv_lineages()].
#'
#' The result orders the samples within one lineage. It is not a fraction,
#' and it cannot be compared between lineages or between runs: adding or
#' removing samples changes every value.
#'
#' @param results A list of results from [bulkdeconv()], for the same samples
#'   and species, one per method (each method at most once).
#' @param min_methods A lineage is reported only when at least this many
#'   methods cover it. Default 2.
#' @param min_samples The smallest number of samples accepted, because ranks
#'   need a cohort. Default 5.
#'
#' @return An object of class `bulkdeconv_consensus`, a list with `score`
#'   (samples x lineages, mean percentile rank in (0, 1)), `agreement` (per
#'   lineage, the mean pairwise Spearman correlation between the methods that
#'   cover it, `NA` with one method), `n_methods`, `sources` (which output
#'   columns of which method entered each lineage) and `dropped` (lineages,
#'   methods and columns left out, with the reason).
#'
#' @examples
#' \dontrun{
#' res <- lapply(c("quantiseq", "mcp_counter", "xcell"), function(m) {
#'   bulkdeconv(tpm, expr = "tpm", species = "human", method = m)
#' })
#' cons <- bulkdeconv_consensus(res)
#' cons$score
#' cons$agreement
#' }
#' @export
bulkdeconv_consensus <- function(results, min_methods = 2L, min_samples = 5L) {
  if (!is.list(results) || !length(results) ||
      !all(vapply(results, inherits, NA, what = "bulkdeconv_result"))) {
    bd_abort("`results` must be a list of results returned by bulkdeconv().")
  }
  for (v in list(min_methods, min_samples)) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 1) {
      bd_abort("`min_methods` and `min_samples` must be whole numbers of at least 1.")
    }
  }
  methods <- vapply(results, function(r) r$metadata$method, "")
  if (anyDuplicated(methods)) {
    bd_abort("Each method may appear once; repeated: ",
             paste(unique(methods[duplicated(methods)]), collapse = ", "), ".")
  }
  species <- unique(vapply(results, function(r) r$metadata$species, ""))
  if (length(species) != 1L) bd_abort("All results must be for the same species.")
  samples <- rownames(results[[1L]]$score)
  same <- vapply(results, function(r) setequal(rownames(r$score), samples), NA)
  if (!all(same)) bd_abort("All results must hold the same samples.")
  if (length(samples) < min_samples) {
    bd_abort("The consensus is rank based and needs at least ", min_samples,
             " samples; the results have ", length(samples), ".")
  }

  map <- bd_read_registry("lineages.csv")
  map <- map[map$species == species, , drop = FALSE]
  map$direction <- as.integer(map$direction)
  dropped <- data.frame(what = character(), name = character(), reason = character(),
                        stringsAsFactors = FALSE)
  drop_row <- function(what, name, reason) {
    dropped[nrow(dropped) + 1L, ] <<- list(what, name, reason)
  }

  per <- list()   # method x lineage rank vectors
  src <- list()
  for (i in seq_along(results)) {
    m <- methods[i]
    sc <- results[[i]]$score[samples, , drop = FALSE]
    mm <- map[map$method == m, , drop = FALSE]
    if (!nrow(mm)) {
      drop_row("method", m, "no lineage mapping")
      next
    }
    for (lin in unique(mm$lineage)) {
      rows <- mm[mm$lineage == lin, , drop = FALSE]
      ranks <- list()
      for (k in seq_len(nrow(rows))) {
        col <- rows$cell_type[k]
        if (!col %in% colnames(sc)) {
          drop_row("column", paste0(m, ": ", col), "not in the output")
          next
        }
        x <- sc[, col]
        if (anyNA(x)) { drop_row("column", paste0(m, ": ", col), "has missing values"); next }
        if (length(unique(x)) < 2L) { drop_row("column", paste0(m, ": ", col), "constant"); next }
        r <- (rank(x, ties.method = "average") - 0.5) / length(x)
        ranks[[col]] <- if (rows$direction[k] < 0L) 1 - r else r
      }
      if (!length(ranks)) next
      per[[paste(lin, m, sep = "\r")]] <- rowMeans(do.call(cbind, ranks))
      src[[paste(lin, m, sep = "\r")]] <- data.frame(
        lineage = lin, method = m, columns = paste(names(ranks), collapse = "; "),
        stringsAsFactors = FALSE)
    }
  }
  keys <- strsplit(names(per), "\r", fixed = TRUE)
  lin_of <- vapply(keys, `[`, "", 1L)
  lineages <- sort(unique(lin_of))

  score <- matrix(NA_real_, length(samples), 0L, dimnames = list(samples, NULL))
  agreement <- numeric(); n_methods <- integer()
  for (lin in lineages) {
    vs <- per[lin_of == lin]
    if (length(vs) < min_methods) {
      drop_row("lineage", lin, paste0("covered by ", length(vs), " method(s); needs ", min_methods))
      next
    }
    mat <- do.call(cbind, vs)
    score <- cbind(score, rowMeans(mat))
    colnames(score)[ncol(score)] <- lin
    n_methods[lin] <- ncol(mat)
    agreement[lin] <- if (ncol(mat) < 2L) NA_real_ else {
      cc <- suppressWarnings(stats::cor(mat, method = "spearman"))
      mean(cc[upper.tri(cc)])
    }
  }
  if (!ncol(score)) {
    bd_abort("No lineage is covered by at least ", min_methods, " of the methods given. ",
             "See bulkdeconv_lineages() for the mapped outputs.")
  }
  structure(list(
    score = score, agreement = agreement[colnames(score)], n_methods = n_methods[colnames(score)],
    sources = do.call(rbind, unname(src[names(src) %in% names(per)[lin_of %in% colnames(score)]])),
    dropped = dropped,
    metadata = list(methods = methods, species = species, n_samples = length(samples),
                    bulkdeconv_version = bd_pkg_version())
  ), class = "bulkdeconv_consensus")
}

#' @export
print.bulkdeconv_consensus <- function(x, ...) {
  cat("<bulkdeconv_consensus> ", length(x$metadata$methods), " methods (",
      paste(x$metadata$methods, collapse = ", "), "), ", x$metadata$n_samples, " samples\n", sep = "")
  cat("  mean percentile rank within this run; orders samples within a lineage only\n")
  tab <- data.frame(lineage = colnames(x$score), n_methods = unname(x$n_methods),
                    agreement = round(unname(x$agreement), 2), row.names = NULL)
  print(tab, row.names = FALSE)
  if (nrow(x$dropped)) cat("  left out: ", nrow(x$dropped), " item(s) (see $dropped)\n", sep = "")
  invisible(x)
}

#' Lineages used by the consensus
#'
#' The table that maps the output columns of each method to the lineages of
#' [bulkdeconv_consensus()].
#'
#' @param species `"human"` (the only species with a mapping so far).
#' @return A data frame with `method`, `species`, `cell_type` (the output
#'   column), `lineage`, and `direction` (1, or -1 when a higher method score
#'   means less of the lineage).
#' @examples
#' head(bulkdeconv_lineages())
#' @export
bulkdeconv_lineages <- function(species = "human") {
  species <- bd_match_species(species)
  out <- bd_read_registry("lineages.csv")
  out <- out[out$species == species, , drop = FALSE]
  out$direction <- as.integer(out$direction)
  rownames(out) <- NULL
  out
}
