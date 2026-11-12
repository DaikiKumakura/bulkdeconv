# Runs an adapter in a separate R process whose library path starts with the
# private library. Warnings raised inside are collected and returned.

bd_cache_dir <- function() tools::R_user_dir("bulkdeconv", "cache")

bd_run_isolated <- function(adapter, X, opts, packages, libs = bd_lib_paths()) {
  ctx <- list(cache_dir = bd_cache_dir(), extdata = system.file("extdata", package = "bulkdeconv"))
  out <- tryCatch(
    callr::r(
      function(adapter, X, opts, ctx, packages) {
        warns <- character()
        res <- withCallingHandlers(
          suppressMessages(adapter(X, opts, ctx)),
          warning = function(w) {
            warns <<- c(warns, conditionMessage(w))
            invokeRestart("muffleWarning")
          })
        versions <- vapply(packages, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), "")
        list(result = res, warnings = warns, versions = versions, r_version = R.version.string)
      },
      args = list(adapter = adapter, X = X, opts = opts, ctx = ctx, packages = packages),
      libpath = c(libs, .libPaths()), show = FALSE, spinner = FALSE, error = "error"),
    error = function(e) {
      parent <- e$parent %||% e
      msg <- gsub("\\s+", " ", conditionMessage(parent))
      bd_abort("The method failed: ", msg)
    })
  out
}

# Turn the adapter output into a samples x features matrix with the right names.
bd_finish_score <- function(est, samples) {
  if (!is.matrix(est) || !is.numeric(est)) bd_abort("The method returned an unexpected result.")
  if (is.null(rownames(est)) || !all(samples %in% rownames(est))) {
    if (nrow(est) != length(samples)) bd_abort("The method returned ", nrow(est), " rows for ",
                                               length(samples), " samples.")
    rownames(est) <- samples
  }
  est <- est[samples, , drop = FALSE]
  est
}
