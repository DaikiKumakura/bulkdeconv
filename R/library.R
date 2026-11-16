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
  data <- bd_data_specs(row$install)
  if (length(data) && !all(vapply(data, function(d) nzchar(bd_find_data(d, libs)), logical(1)))) return(FALSE)
  pkgs <- bd_split(row$packages)
  if (!length(pkgs) && length(data)) return(TRUE)
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

# ---- method data downloaded at install time ----
# Some methods' data cannot be distributed with the package (licence); they are
# downloaded from a fixed commit, checked by size and MD5, converted to an RDS
# file and kept in <lib>/_bulkdeconv_data/<name>/. Install spec: "data::<name>".
bd_data_sources <- list(
  puree = list(
    url = "https://raw.githubusercontent.com/skandlab/PUREE/5a0a702535e79e37b071971063e72fa697540818/",
    files = data.frame(
      path = c("models/model.joblib", "models/median_imputer.joblib", "data/predictive_genes.csv",
               "data/ranking_universe.csv", "data/gene_id_conversion_table.csv"),
      md5 = c("b7058fc32dae966d41f165c40d4fde12", "cad65d68b8ce78a08fb680b3f30bf3b6", "db18c84300eb7346a4f28c626eb9393f",
              "9faa30744be80a02c49454531d0bd3ac", "403d119d13878559d4963fe81c1b80fb"),
      size = c(6799, 5468, 3062, 199536, 515657), stringsAsFactors = FALSE)
  )
)
bd_data_specs <- function(install) {
  s <- bd_split(install, sep = ";")
  sub("^data::", "", s[startsWith(s, "data::")])
}
bd_find_data <- function(name, libs = bd_lib_paths()) {
  f <- file.path(libs, "_bulkdeconv_data", name, paste0(name, ".rds"))
  f <- f[file.exists(f)]
  if (length(f)) f[1L] else ""
}
bd_install_data <- function(name, lib) {
  src <- bd_data_sources[[name]]
  if (is.null(src)) bd_abort("Unknown data source '", name, "' in the registry.")
  dir <- file.path(lib, "_bulkdeconv_data", name)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, basename(src$files$path))
  for (i in seq_len(nrow(src$files))) {
    utils::download.file(paste0(src$url, src$files$path[i]), paths[i], mode = "wb", quiet = TRUE)
    if (file.size(paths[i]) != src$files$size[i] || unname(tools::md5sum(paths[i])) != src$files$md5[i]) {
      unlink(paths[i])
      bd_abort("Downloaded file '", src$files$path[i], "' does not match the recorded size and checksum.")
    }
  }
  saveRDS(bd_build_data(name, paths), file.path(dir, paste0(name, ".rds")))
  invisible(dir)
}
# PUREE (Revkov et al. 2023): a linear model on percentile ranks of 158 genes.
# The coefficients, intercept and imputation medians are read as little-endian
# doubles at fixed offsets of the checked joblib files (no unpickling); the
# offsets belong to the files of the recorded commit.
bd_build_data <- function(name, paths) {
  stopifnot(identical(name, "puree"))
  rd <- function(f, offset, n) {
    con <- file(f, "rb"); on.exit(close(con))
    seek(con, offset)
    readBin(con, "double", n, size = 8L, endian = "little")
  }
  genes <- utils::read.csv(paths[3], stringsAsFactors = FALSE)$ENSEMBL_ID
  conv <- utils::read.csv(paths[5], stringsAsFactors = FALSE)
  list(genes = genes, coef = rd(paths[1], 4105L, 158L), intercept = rd(paths[1], 6756L, 1L),
       medians = rd(paths[2], 4165L, 158L),
       universe = utils::read.csv(paths[4], stringsAsFactors = FALSE)$ENSEMBL_ID,
       conversion = conv[, c("ENSEMBL_ID", "HGNC")])
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
  all_specs <- unique(unlist(lapply(rows$install, bd_split, sep = ";")))
  data_specs <- all_specs[startsWith(all_specs, "data::")]
  if (length(data_specs) == length(all_specs)) {
    lib <- normalizePath(lib %||% bd_lib_paths()[1L], winslash = "/", mustWork = FALSE)
    if (ask) {
      cat("bulkdeconv will download into:\n  ", file.path(lib, "_bulkdeconv_data"), "\n", sep = "")
      cat(paste0("  - ", data_specs, collapse = "\n"), "\n", sep = "")
      ans <- readline("Proceed? [y/N] ")
      if (!tolower(trimws(ans)) %in% c("y", "yes")) bd_abort("Installation cancelled.")
    }
    for (d in sub("^data::", "", data_specs)) { bd_note("downloading data for ", d); bd_install_data(d, lib) }
    return(invisible(lib))
  }
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
    if (p$type == "data") { bd_note("downloading data for ", p$name); bd_install_data(p$name, lib); next }
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
