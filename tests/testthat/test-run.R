sys.source(test_path("fixtures", "input-helpers.R"), envir = environment())

# End-to-end through bulkdeconv() with a stand-in adapter (no method package needed).

local_dummy <- function(adapter = dummy_adapter, env = parent.frame()) {
  testthat::local_mocked_bindings(
    bd_is_installed = function(...) TRUE,
    bd_get_adapter = function(method) bulkdeconv:::bd_isolate(adapter),
    .env = env)
}

test_that("bulkdeconv() returns a standard result object", {
  skip_on_cran()
  local_dummy()
  d <- make_input()
  res <- suppressWarnings(bulkdeconv(d, expr = "tpm", species = "human", method = "epic"))
  expect_s3_class(res, "bulkdeconv_result")
  expect_equal(rownames(res$score), c("S1", "S2", "S3"))
  expect_equal(colnames(res$score), c("total", "mean"))
  expect_equal(unname(res$score[, "total"]), unname(colSums(d[, -1])))
  expect_equal(res$metadata$method, "epic")
  expect_equal(res$metadata$output_class, "fraction")
  expect_equal(res$qc$n_genes, 40L)
  expect_equal(res$qc$signature_genes_found, 8)
  expect_false(res$qc$single_sample_workaround)
  expect_output(print(res), "bulkdeconv_result")
})

test_that("the licence notice for EPIC is raised as a warning", {
  skip_on_cran()
  local_dummy()
  expect_warning(bulkdeconv(make_input(), "tpm", "human", "epic"), "non-commercial")
})

test_that("a single sample is duplicated for runs that need two and the copy is dropped", {
  skip_on_cran()
  seen <- NULL
  adapter <- function(X, opts, ctx) {
    list(score = cbind(n = rep(ncol(X), ncol(X)), total = colSums(X)))
  }
  local_dummy(adapter)
  d <- make_input(n_samples = 1)
  res <- suppressWarnings(bulkdeconv(d, "tpm", "human", "xcell2"))
  expect_equal(nrow(res$score), 1L)
  expect_equal(rownames(res$score), "S1")
  expect_equal(unname(res$score[1, "n"]), 2)
  expect_true(res$qc$single_sample_workaround)
  w <- testthat::capture_warnings(bulkdeconv(d, "tpm", "human", "xcell2"))
  expect_true(any(grepl("duplicated internally", w)))
})

test_that("errors inside the method surface with the original message", {
  skip_on_cran()
  local_dummy(function(X, opts, ctx) stop("boom from the method", call. = FALSE))
  expect_error(suppressWarnings(bulkdeconv(make_input(), "tpm", "human", "epic")),
               "The method failed: .*boom from the method")
})

test_that("warnings inside the method are carried into the result", {
  skip_on_cran()
  local_dummy(function(X, opts, ctx) {
    warning("careful")
    list(score = cbind(a = colSums(X)))
  })
  res <- suppressWarnings(bulkdeconv(make_input(), "tpm", "human", "quantiseq"))
  expect_true(any(grepl("\\[quantiseq\\] careful", res$qc$warnings)))
})

test_that("arguments are validated before any method runs", {
  d <- make_input()
  expect_error(bulkdeconv(d, "tpm", "human"), "all required")
  expect_error(bulkdeconv(d, "rpkm", "human", "epic"), "count")
  expect_error(bulkdeconv(d, "tpm", "rat", "epic"), "human")
  expect_error(bulkdeconv(d, "tpm", "human", "nope"), "Unknown method")
  expect_error(bulkdeconv(d, "count", "human", "quantiseq"), "accepts expr")
})

test_that("a method that is not installed gives an install hint", {
  old <- options(bulkdeconv.lib = file.path(tempdir(), "empty-lib"))
  on.exit(options(old))
  testthat::local_mocked_bindings(bd_is_installed = function(...) FALSE)
  expect_error(bulkdeconv(make_input(), "tpm", "human", "quantiseq"), "bulkdeconv_install")
})
