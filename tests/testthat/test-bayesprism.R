check_bp <- function(...) {
  bulkdeconv:::bd_check_method(bulkdeconv:::bd_method_row("bayesprism"), "mouse", "count", 3L, list(...))
}

test_that("bayesprism needs counts and a listed model", {
  row <- bulkdeconv:::bd_method_row("bayesprism")
  expect_error(bulkdeconv:::bd_check_method(row, "mouse", "tpm", 3L, list(model = "MC38")), "accepts expr = \"count\"")
  expect_error(check_bp(), "`model` is required")
  expect_error(check_bp(model = "EMT6"), "must be one of")
  p <- check_bp(model = "MC38")
  expect_identical(p$opts, list(model = "MC38", n_cores = 1L, seed = 1L))
  expect_error(check_bp(model = "MC38", n_cores = 0), "whole number")
  expect_false(any(grepl("candidate", p$warnings)))
  expect_true(any(grepl("candidate", check_bp(model = "MOC1")$warnings)))
})

test_that("every listed model has a consistent bundled reference", {
  models <- bulkdeconv:::bd_read_registry("bayesprism_models.csv")
  expect_length(models$model, 23L)
  expect_true(all(c("B16F10", "MC38", "CT26", "LLC") %in% models$model))
  expect_true(all(models$status %in% c("validated", "candidate")))
  for (i in seq_len(nrow(models))) {
    ref <- readRDS(system.file("extdata", "bayesprism", paste0(models$model[i], ".rds"), package = "bulkdeconv"))
    expect_type(ref$gep, "integer")
    expect_identical(nrow(ref$gep), as.integer(models$states[i]))
    expect_identical(ncol(ref$gep), as.integer(models$genes[i]))
    expect_length(ref$state_type, nrow(ref$gep))
    expect_true(ref$key %in% ref$state_type)
    expect_false(anyDuplicated(colnames(ref$gep)) > 0)
    expect_false(any(grepl("^(Rp[sl][0-9]+[a-z]?|mt-.*)$", colnames(ref$gep))))
    expect_true(all(rowSums(ref$gep) > 0))
  }
})
