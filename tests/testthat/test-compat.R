row_of <- function(m) bulkdeconv:::bd_method_row(m)
check <- function(m, species = "human", expr = "tpm", n = 5L, ...) {
  bulkdeconv:::bd_check_method(row_of(m), species, expr, n, list(...))
}

test_that("excluded and unavailable methods are refused", {
  expect_error(check("immucc", "mouse"), "excluded")
  expect_error(check("immucellai"), "cannot be run")
})

test_that("species mismatch lists methods that fit", {
  expect_error(check("quantiseq", "mouse"), "supports human, not mouse")
  expect_error(check("mmcp_counter", "human"), "supports mouse")
  expect_no_error(check("xcell2", "mouse"))
})

test_that("count input is refused for tpm-only methods", {
  expect_error(check("quantiseq", expr = "count"), "accepts expr = \"tpm\"")
})

test_that("single-sample handling follows the registry", {
  p <- check("xcell2", n = 1L)
  expect_true(p$duplicate)
  expect_match(p$warnings[1], "duplicated internally")
  expect_false(check("quantiseq", n = 1L)$duplicate)
  expect_error(check("xcell", n = 1L), "at least 2")
  expect_error(check("ssmd", "mouse", n = 3L), "at least 6")
})

test_that("cohort-dependent methods warn on small n", {
  p <- check("xcell", n = 4L)
  expect_true(any(grepl("depends on the other samples", p$warnings)))
  expect_false(any(grepl("depends on the other samples", check("quantiseq", n = 4L)$warnings)))
})

test_that("licence notices are attached", {
  expect_true(any(grepl("non-commercial", check("epic")$warnings)))
  expect_true(any(grepl("no declared licence", check("ssmd", "mouse", n = 8L)$warnings)))
})

test_that("options are validated and defaulted", {
  expect_equal(check("xcell2")$opts$preset, "TMECompendium")
  expect_equal(check("xcell2", "mouse")$opts$preset, "MouseRNAseqData")
  expect_equal(check("xcell2", preset = "LM22")$opts$preset, "LM22")
  expect_error(check("xcell2", preset = "ImmGenData"), "must be one of")
  expect_equal(check("consensus_tme", cancer_type = "LUAD")$opts$cancer_type, "LUAD")
  expect_error(check("consensus_tme", cancer_type = "XYZ"), "must be one of")
  expect_equal(check("immucellai2")$opts$tissue, "tumor")
  expect_error(check("quantiseq", tissue = "tumor"), "takes no options")
  expect_error(check("immucellai2", foo = 1), "Unknown option")
  expect_error(bulkdeconv:::bd_check_options("immucellai2", list(1), "human"), "must be named")
})

test_that("seqimmucc is the method that takes counts", {
  expect_no_error(check("seqimmucc", "mouse", expr = "count"))
  expect_error(check("seqimmucc", "mouse", expr = "tpm"), "accepts expr = \"count\"")
  expect_false(check("seqimmucc", "mouse", expr = "count", n = 1L)$duplicate)  # native: works on one sample
})

test_that("timer needs a valid cancer type; base takes a permutation count", {
  expect_error(check("timer"), "required")
  expect_equal(check("timer", cancer_type = "LUAD")$opts$cancer_type, "LUAD")
  expect_error(check("timer", cancer_type = "LAML"), "must be one of")
  expect_error(check("timer", cancer_type = "LUAD", n = 1L), "at least 2")
  expect_equal(check("base", "mouse")$opts$n_permutations, 100L)
  expect_equal(check("base", "mouse")$opts$seed, 1L)
  expect_equal(check("base", "mouse", seed = 7)$opts$seed, 7L)
  expect_equal(check("base", "mouse", n_permutations = 10)$opts$n_permutations, 10L)
  expect_error(check("base", "mouse", n_permutations = 0), "whole number")
  expect_error(check("base", "mouse", n_permutations = 2.5), "whole number")
  expect_true(any(grepl("permutations", check("base", "mouse")$warnings)))
})

test_that("abis warns that its signature is from blood", {
  expect_true(any(grepl("blood", check("abis")$warnings)))
})
