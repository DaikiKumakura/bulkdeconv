# The private library: where the third-party method packages live, apart from
# the user's normal library, so versions pinned here never disturb it.

bd_default_lib <- function() {
  file.path(tools::R_user_dir("bulkdeconv", "data"), "lib",
            paste0("R-", R.version$major, ".", strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L]))
}

# Library paths used for method packages. Option `bulkdeconv.lib` (or the
# environment variable BULKDECONV_LIB, separated by ";") may name several
# directories; the first is where installs go.
bd_lib_paths <- function() {
  o <- getOption("bulkdeconv.lib")
  if (is.null(o)) {
    e <- Sys.getenv("BULKDECONV_LIB", "")
    if (nzchar(e)) o <- strsplit(e, ";", fixed = TRUE)[[1L]]
  }
  o <- o %||% bd_default_lib()
  normalizePath(o, winslash = "/", mustWork = FALSE)
}

bd_is_installed <- function(method, libs = bd_lib_paths()) {
  row <- bd_registry()
  row <- row[row$method == method, , drop = FALSE]
  pkgs <- bd_split(row$packages)
  if (!length(pkgs)) return(identical(row$engine, "native"))
  all(vapply(pkgs, function(p) nzchar(system.file(package = p, lib.loc = c(libs, .libPaths()))), logical(1)))
}

# "github::owner/repo@sha#subdir" -> list(type, name, ref, subdir)
bd_parse_spec <- function(spec) {
  type <- sub("::.*$", "", spec)
  rest <- sub("^[^:]+::", "", spec)
  subdir <- if (grepl("#", rest, fixed = TRUE)) sub("^.*#", "", rest) else NULL
  rest <- sub("#.*$", "", rest)
  ref <- if (grepl("@", rest, fixed = TRUE)) sub("^.*@", "", rest) else NULL
  name <- sub("@.*$", "", rest)
  list(type = type, name = name, ref = ref, subdir = subdir)
}

#' Install the packages a method needs into the private library
#'
#' Installs the original method packages at the versions recorded in the
#' registry (GitHub packages at a fixed commit) into a private library, so they
#' do not touch your normal library. Interactively it asks before installing.
#'
#' @param method One or more method IDs, or `"all"` for every usable method.
#' @param lib Library directory. Default: the first path of the option
#'   `bulkdeconv.lib`, else a per-user directory (see Details).
#' @param ask Ask for confirmation first. Default: `interactive()`.
#' @details The default library is `tools::R_user_dir("bulkdeconv", "data")`
#'   plus `lib/R-<major.minor>`. Set `options(bulkdeconv.lib = "path")` to use
#'   another place. The packages `remotes` and `BiocManager` are needed for the
#'   installation itself.
#' @return Invisibly, the library path used.
#' @examples
#' \dontrun{
#' bulkdeconv_install("quantiseq")
#' bulkdeconv_install("all")
#' }
#' @export
bulkdeconv_install <- function(method, lib = NULL, ask = interactive()) {
  reg <- bulkdeconv_methods(usable = TRUE)
  if (identical(method, "all")) method <- reg$method
  rows <- do.call(rbind, lapply(method, bd_method_row))
  rows <- rows[rows$method %in% reg$method, , drop = FALSE]
  if (nrow(rows) == 0L) bd_abort("None of the requested methods can be installed.")
  for (p in c("remotes", "BiocManager")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      bd_abort("Installing method packages needs the '", p, "' package: install.packages(\"", p, "\")")
    }
  }
  lib <- normalizePath(lib %||% bd_lib_paths()[1L], winslash = "/", mustWork = FALSE)
  specs <- unique(unlist(lapply(rows$install, bd_split, sep = ";")))
  if (ask) {
    cat("bulkdeconv will install into:\n  ", lib, "\n", sep = "")
    cat(paste0("  - ", specs, collapse = "\n"), "\n", sep = "")
    ans <- readline("Proceed? [y/N] ")
    if (!tolower(trimws(ans)) %in% c("y", "yes")) bd_abort("Installation cancelled.")
  }
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  old <- .libPaths()
  on.exit(.libPaths(old), add = TRUE)
  .libPaths(c(lib, old))
  repos <- BiocManager::repositories()
  deps <- c("Depends", "Imports", "LinkingTo")
  for (s in specs) {
    p <- bd_parse_spec(s)
    bd_note("installing ", s)
    switch(p$type,
      cran = ,
      bioc = BiocManager::install(p$name, lib = lib, update = FALSE, ask = FALSE,
                                  dependencies = deps),
      github = remotes::install_github(p$name, ref = p$ref, subdir = p$subdir, lib = lib,
                                       upgrade = "never", dependencies = deps, repos = repos,
                                       quiet = TRUE),
      bd_abort("Unknown install type '", p$type, "' in the registry."))
  }
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    pkgs <- bd_split(r$packages)
    missing <- pkgs[!vapply(pkgs, function(p) nzchar(system.file(package = p, lib.loc = lib)), logical(1))]
    if (length(missing)) {
      bd_warn("Method '", r$method, "': package(s) not found in ", lib, " after installation: ",
              paste(missing, collapse = ", "))
    }
  }
  invisible(lib)
}
