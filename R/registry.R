bd_registry_path <- function(file) {
  p <- system.file("registry", file, package = "bulkdeconv")
  if (!nzchar(p)) bd_abort("Registry file '", file, "' was not found in the installed package.")
  p
}

bd_read_registry <- function(file = "methods.csv") {
  utils::read.csv(bd_registry_path(file), stringsAsFactors = FALSE,
                  na.strings = character(), colClasses = "character",
                  encoding = "UTF-8", check.names = FALSE)
}

bd_registry <- function() {
  reg <- bd_read_registry("methods.csv")
  reg$min_samples_n <- suppressWarnings(as.integer(reg$min_samples_n))
  reg$min_samples_n[is.na(reg$min_samples_n)] <- 1L
  reg
}

# One registry row, or an informative error for an unknown method.
bd_method_row <- function(method) {
  reg <- bd_registry()
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    bd_abort("`method` must be a single method ID; see bulkdeconv_methods().")
  }
  row <- reg[reg$method == method, , drop = FALSE]
  if (nrow(row) == 0L) {
    close <- utils::head(reg$method[agrepl(method, reg$method, max.distance = 0.3)], 3L)
    bd_abort("Unknown method '", method, "'.",
             if (length(close)) paste0(" Did you mean: ", paste(close, collapse = ", "), "?") else "",
             " See bulkdeconv_methods() for the available IDs.")
  }
  row
}

#' List the methods bulkdeconv can run
#'
#' Returns the registry: one row per method, with the species it supports,
#' the expression scale it accepts, what it needs, how its output can be
#' compared, its licence, and whether it is installed in the private library.
#'
#' @param species Optional. `"human"` or `"mouse"`; keeps only methods that
#'   support it.
#' @param usable If `TRUE` (default), drops methods that cannot be run
#'   (status `unavailable` or `excluded`).
#' @return A data frame with one row per method. `installed` is `TRUE` when
#'   every package the method needs is found in the private library or in the
#'   normal library paths.
#' @examples
#' bulkdeconv_methods()
#' bulkdeconv_methods("mouse")[, c("method", "output_class", "installed")]
#' @export
bulkdeconv_methods <- function(species = NULL, usable = TRUE) {
  reg <- bd_registry()
  if (!is.null(species)) {
    species <- bd_match_species(species)
    reg <- reg[vapply(reg$species, function(s) species %in% bd_split(s), logical(1)), , drop = FALSE]
  }
  if (usable) reg <- reg[!reg$status %in% c("unavailable", "excluded"), , drop = FALSE]
  libs <- bd_lib_paths()
  reg$installed <- vapply(reg$method, function(m) bd_is_installed(m, libs), logical(1))
  rownames(reg) <- NULL
  reg
}

bd_match_species <- function(species) {
  if (!is.character(species) || length(species) != 1L || is.na(species) ||
      !species %in% c("human", "mouse")) {
    bd_abort("`species` must be \"human\" or \"mouse\".")
  }
  species
}

bd_match_expr <- function(expr) {
  if (!is.character(expr) || length(expr) != 1L || is.na(expr) ||
      !expr %in% c("count", "tpm")) {
    bd_abort("`expr` must be \"count\" or \"tpm\".")
  }
  expr
}
