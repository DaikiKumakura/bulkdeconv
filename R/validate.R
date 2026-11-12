# Input validation. The expected layout is the one rsemverse::rsem_read()
# returns: first column = gene symbols, remaining columns = one per sample
# (count or TPM). A matrix with gene symbols as row names is also accepted.

bd_prepare_input <- function(data, expr, species) {
  notes <- character()

  if (is.matrix(data)) {
    if (is.null(rownames(data))) bd_abort("A matrix `data` needs gene symbols as row names.")
    genes <- rownames(data)
    if (is.null(colnames(data))) bd_abort("A matrix `data` needs sample IDs as column names.")
    mat <- data
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(df) < 2L) bd_abort("`data` needs a gene column plus at least one sample column.")
    if (is.character(df[[1L]]) || is.factor(df[[1L]])) {
      genes <- as.character(df[[1L]])
      mat <- df[, -1L, drop = FALSE]
      names(mat) <- names(df)[-1L]  # `[.data.frame` silently renames duplicates
    } else if (!identical(rownames(df), as.character(seq_len(nrow(df))))) {
      genes <- rownames(df)
      mat <- df
      notes <- c(notes, "Using row names as gene symbols (the first column is numeric).")
    } else {
      bd_abort("The first column of `data` must hold gene symbols (character). ",
               "It is numeric and there are no row names to fall back on.")
    }
  } else {
    bd_abort("`data` must be a data frame (first column = gene symbols, other columns = samples) ",
             "or a numeric matrix with gene symbols as row names.")
  }

  genes <- trimws(genes)
  if (anyNA(genes) || any(!nzchar(genes))) {
    bd_abort(sum(is.na(genes) | !nzchar(genes)), " gene symbol(s) are missing or empty.")
  }
  if (anyDuplicated(genes)) {
    d <- unique(genes[duplicated(genes)])
    bd_abort(length(d), " gene symbol(s) occur more than once (e.g. ",
             paste(utils::head(d, 3L), collapse = ", "),
             "). Merge duplicates before calling bulkdeconv(); rsemverse::rsem_read() already does.")
  }

  samples <- colnames(mat)
  if (is.null(samples) || anyNA(samples) || any(!nzchar(samples))) {
    bd_abort("Every sample column needs a non-empty name.")
  }
  if (anyDuplicated(samples)) {
    bd_abort("Sample names must be unique; duplicated: ",
             paste(utils::head(unique(samples[duplicated(samples)]), 3L), collapse = ", "), ".")
  }
  if (!length(samples)) bd_abort("`data` needs at least one sample column.")

  is_num <- if (is.matrix(mat)) rep(is.numeric(mat), ncol(mat)) else vapply(mat, is.numeric, logical(1))
  if (!all(is_num)) {
    bd_abort("Non-numeric sample column(s): ",
             paste(utils::head(samples[!is_num], 3L), collapse = ", "),
             ". All columns after the gene column must be numeric.")
  }
  x <- as.matrix(mat)
  storage.mode(x) <- "double"
  dimnames(x) <- list(genes, samples)

  if (anyNA(x)) bd_abort(sum(is.na(x)), " value(s) are NA or NaN; remove or impute them first.")
  if (any(is.infinite(x))) bd_abort("`data` contains infinite values.")
  if (any(x < 0)) {
    bd_abort("`data` contains negative values (", sum(x < 0),
             "). Counts and TPM are non-negative; log-transformed or batch-corrected data is not supported.")
  }

  if (nrow(x) < 10L) bd_abort("Only ", nrow(x), " genes in `data`; at least 10 are required.")

  if (mean(grepl("^ENS[A-Z]*G[0-9]+", genes)) > 0.5) {
    bd_abort("Most gene IDs look like Ensembl IDs. bulkdeconv() needs gene symbols ",
             "(rsemverse::rsem_read() converts them).")
  }

  has_letter <- grepl("[A-Za-z]", genes)
  if (!any(has_letter)) {
    bd_abort("No gene IDs look like gene symbols. Supply symbols such as ACTB or Actb.")
  }
  frac_lower <- mean(grepl("[a-z]", genes[has_letter]))
  if (species == "human" && frac_lower > 0.5) {
    bd_abort("`species` is \"human\" but ", round(100 * frac_lower),
             "% of the gene symbols contain lower-case letters, which is typical of mouse symbols (e.g. Actb).")
  }
  if (species == "mouse" && frac_lower < 0.5) {
    bd_abort("`species` is \"mouse\" but only ", round(100 * frac_lower),
             "% of the gene symbols contain lower-case letters; human symbols are upper case (e.g. ACTB).")
  }

  zero_samples <- samples[colSums(x) == 0]
  if (length(zero_samples)) {
    bd_abort("Sample column(s) contain only zeros: ",
             paste(utils::head(zero_samples, 3L), collapse = ", "), ".")
  }

  non_integer <- any(abs(x - round(x)) > 1e-8)
  if (expr == "count" && non_integer) {
    notes <- c(notes, "Counts are not all integers (expected counts, e.g. from RSEM); they are used as they are.")
  }
  if (expr == "tpm" && max(x) < 30 && non_integer) {
    bd_warn("Maximum value is ", signif(max(x), 3),
            ": the data may be log-transformed. Methods expect linear-scale TPM.")
  }
  if (expr == "tpm") {
    sums <- colSums(x)
    if (any(sums < 1e5 | sums > 1.2e6)) {
      notes <- c(notes, paste0("Column sums of TPM range from ", format(min(sums), big.mark = ","),
                               " to ", format(max(sums), big.mark = ","), " (expected about 1,000,000)."))
    }
  }

  list(matrix = x, notes = notes, n_genes = nrow(x), n_samples = ncol(x))
}
