#' bulkdeconv: cell-type deconvolution of tumor bulk RNA-seq through one function
#'
#' bulkdeconv runs 19 published cell-type deconvolution and cell-abundance
#' methods for human and mouse tumor bulk RNA-seq through one call, and
#' returns the same kind of result for all of them.
#'
#' @section Workflow:
#' 1. [bulkdeconv_methods()] lists the methods with species, input scale,
#'    minimum number of samples, output kind and licence.
#' 2. [bulkdeconv()] runs one method on a gene symbol by sample table (TPM;
#'    counts for `seqimmucc` and `bayesprism`) and returns `score` (samples x cell types),
#'    `metadata` (what the numbers are and what may be compared) and `qc`.
#'
#' Fifteen methods are implemented in the package and run after installing
#' it; their reference data come with the package. EPIC, BayesPrism, DCQ and SSMD run
#' from their original packages: install them once with
#' [bulkdeconv_install()] into a private library.
#'
#' @keywords internal
"_PACKAGE"
