#' Cell-type deconvolution of tumor bulk RNA-seq with one function
#'
#' Runs one published deconvolution or cell-abundance method on a gene-by-sample
#' table. The input layout is the one `rsemverse::rsem_read()` returns: the
#' first column holds gene symbols and every other column is one sample.
#' The method runs in a separate R process. External methods use the original
#' package at the version recorded in the registry (see [bulkdeconv_install()]).
#'
#' @param data A data frame whose first column is the gene symbol (character)
#'   and whose other columns are samples (numeric). A numeric matrix with gene
#'   symbols as row names is also accepted. Gene symbols must be unique and
#'   every sample must contain at least one nonzero value.
#' @param expr `"count"` or `"tpm"`: the scale of the values in `data`.
#'   Each method accepts one scale (see `expr_accepted` in
#'   [bulkdeconv_methods()]): `"count"` for seqimmucc, `"tpm"` for the others.
#'   The wrong one is refused with a message that says so.
#' @param species `"human"` or `"mouse"`.
#' @param method A method ID; see [bulkdeconv_methods()].
#' @param ... Method-specific options, named: `preset` (xcell2),
#'   `cancer_type` (consensus_tme; required for timer), `tissue` (immucellai2,
#'   ssmd), `n_permutations` and `seed` (base). Unknown or invalid options are errors.
#' @details
#' Use [bulkdeconv_methods()] to check the input scale, minimum sample count
#' and output class. Scores from different methods can have different units.
#' Cohort-dependent scores are interpreted within the run.
#'
#' All reference data the methods need are installed with the package.
#' @return An object of class `bulkdeconv_result`, a list with
#'   * `score`: numeric matrix, samples in rows, features in columns;
#'   * `metadata`: method, species, expr, options, what kind of value the
#'     method returns (`output_class`) and how it can be compared
#'     (`comparability`), licence, package versions, runtime;
#'   * `qc`: gene and sample counts, signature coverage when the method
#'     exposes its genes, and the notes and warnings raised.
#' @examples
#' \dontrun{
#' # df: first column gene symbols, other columns samples (TPM)
#' res <- bulkdeconv(df, expr = "tpm", species = "human", method = "quantiseq")
#' res$score
#' }
#' @export
bulkdeconv <- function(data, expr, species, method, ...) {
  if (missing(data) || missing(expr) || missing(species) || missing(method)) {
    bd_abort("`data`, `expr`, `species` and `method` are all required; see ?bulkdeconv.")
  }
  expr <- bd_match_expr(expr)
  species <- bd_match_species(species)
  row <- bd_method_row(method)

  input <- bd_prepare_input(data, expr, species)
  plan <- bd_check_method(row, species, expr, input$n_samples, list(...))

  libs <- bd_lib_paths()
  if (!bd_is_installed(method, libs)) {
    bd_abort("The packages for '", method, "' are not installed. Run bulkdeconv_install(\"",
             method, "\") first.")
  }

  warnings_out <- c(plan$warnings)
  X <- input$matrix
  if (plan$duplicate) {
    X <- cbind(X, X)
    colnames(X) <- make.unique(c(colnames(input$matrix), ".bulkdeconv_copy"))
  }
  for (w in warnings_out) warning(w, call. = FALSE)

  t0 <- Sys.time()
  out <- bd_run_isolated(bd_get_adapter(method), X, plan$opts, bd_split(row$packages), libs)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  est <- bd_finish_score(out$result$score, colnames(X))
  if (plan$duplicate) est <- est[colnames(input$matrix), , drop = FALSE]
  for (w in unique(out$warnings)) warning("[", method, "] ", w, call. = FALSE)

  cov <- out$result$coverage
  structure(list(
    score = est,
    metadata = list(
      method = method, species = species, expr = expr, options = plan$opts, engine = row$engine,
      output_class = row$output_class, comparability = row$comparability,
      cohort_dependent = row$cohort_dependent, license = row$license,
      registration = row$registration, packages = as.list(out$versions),
      r_version = out$r_version, bulkdeconv_version = bd_pkg_version(),
      seconds = round(secs, 1)),
    qc = list(
      n_genes = input$n_genes, n_samples = input$n_samples,
      signature_genes = if (is.null(cov)) NA_integer_ else unname(cov[["signature"]]),
      signature_genes_found = if (is.null(cov)) NA_integer_ else unname(cov[["found"]]),
      single_sample_workaround = plan$duplicate,
      notes = input$notes,
      warnings = c(plan$warnings, if (length(out$warnings)) paste0("[", method, "] ", unique(out$warnings)))
    )
  ), class = "bulkdeconv_result")
}

#' @export
print.bulkdeconv_result <- function(x, ...) {
  m <- x$metadata
  cat("<bulkdeconv_result> ", m$method, " (", m$species, ", ", m$expr, ")\n", sep = "")
  cat("  ", nrow(x$score), " samples x ", ncol(x$score), " features; ",
      x$qc$n_genes, " input genes; ", m$seconds, " s\n", sep = "")
  cat("  output: ", m$output_class, " -- ", m$comparability, "\n", sep = "")
  if (!is.na(x$qc$signature_genes)) {
    cat("  signature genes found: ", x$qc$signature_genes_found, "/", x$qc$signature_genes, "\n", sep = "")
  }
  if (length(m$options)) {
    cat("  options: ", paste(names(m$options), unlist(m$options), sep = " = ", collapse = ", "), "\n", sep = "")
  }
  if (length(x$qc$notes)) cat("  notes: ", paste(x$qc$notes, collapse = " "), "\n", sep = "")
  if (length(x$qc$warnings)) cat("  warnings: ", length(x$qc$warnings), " (see $qc$warnings)\n", sep = "")
  cat("  $score (first rows):\n")
  print(utils::head(round(x$score[, seq_len(min(5L, ncol(x$score))), drop = FALSE], 4L), 4L))
  invisible(x)
}
