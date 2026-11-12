#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

bd_abort <- function(...) stop(paste0(...), call. = FALSE)

bd_warn <- function(...) warning(paste0(...), call. = FALSE, immediate. = FALSE)

bd_note <- function(...) message("bulkdeconv: ", paste0(...))

bd_split <- function(x, sep = "|") {
  if (is.na(x) || !nzchar(x)) character() else strsplit(x, sep, fixed = TRUE)[[1]]
}

bd_pkg_version <- function() as.character(utils::packageVersion("bulkdeconv"))

# The native methods call these packages inside functions that run in a
# separate R process (see R/native*.R). This never-called function states the
# dependency where R CMD check looks for it, without loading the packages when
# bulkdeconv itself is loaded.
bd_runtime_imports <- function() {
  list(MASS::rlm, quadprog::solve.QP, sva::ComBat, BiocParallel::bpparam)
}
