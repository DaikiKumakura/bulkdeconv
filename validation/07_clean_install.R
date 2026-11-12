# Run in a clean R session with R_LIBS and R_LIBS_USER set to an empty directory.
lib <- normalizePath(commandArgs(trailingOnly = TRUE)[1L], winslash = "/")
stopifnot(dir.exists(lib), !length(list.files(lib)))
.libPaths(c(lib, .Library), include.site = FALSE)
stopifnot(identical(.libPaths(), c(lib, normalizePath(.Library, winslash = "/"))))
options(timeout = 600, repos = c(CRAN = "https://cloud.r-project.org"), warn = 1)
deps <- c("Depends", "Imports", "LinkingTo")
install.packages(c("remotes", "BiocManager", "testthat"), lib = lib, dependencies = deps)
# One step, as in the README: bulkdeconv with everything its native methods need.
remotes::install_local(".", lib = lib, dependencies = deps, upgrade = "never", repos = BiocManager::repositories())
stopifnot(with(bulkdeconv::bulkdeconv_methods(), all(installed[engine == "native"])))
stopifnot(as.character(BiocManager::version()) == "3.22")
bulkdeconv::bulkdeconv_install("all", lib = lib, ask = FALSE)

registry <- bulkdeconv::bulkdeconv_methods()
stopifnot(nrow(registry) == 18L, all(registry$installed))
packages <- unique(unlist(strsplit(registry$packages, "|", fixed = TRUE)))
packages <- packages[nzchar(packages)]
versions <- vapply(packages, function(pkg) {
  stopifnot(requireNamespace(pkg, quietly = TRUE))
  as.character(getNamespaceVersion(pkg))
}, character(1L))
print(versions)

# quantiseqr supplies the human integration fixture. Install it only after
# confirming that bulkdeconv_install("all") can install all runtime packages.
BiocManager::install("quantiseqr", lib = lib, dependencies = deps, ask = FALSE, update = FALSE)
Sys.setenv(BULKDECONV_TEST_LIB = lib, BULKDECONV_TEST_SLOW = "1", NOT_CRAN = "true")
testthat::test_local(filter = "integration", stop_on_failure = TRUE)
