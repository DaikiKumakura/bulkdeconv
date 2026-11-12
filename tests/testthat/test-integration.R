# Runs the real methods. Skipped unless BULKDECONV_TEST_LIB names a library
# (";"-separated) that holds the method packages, e.g. after
# bulkdeconv_install("all"). Test data come from the packages themselves.

libs <- strsplit(Sys.getenv("BULKDECONV_TEST_LIB"), ";", fixed = TRUE)[[1L]]
local_method_library <- function() {
  skip_if(!length(libs), "BULKDECONV_TEST_LIB is not set")
  skip_on_cran()
  withr::local_options(bulkdeconv.lib = libs, .local_envir = parent.frame())
}

get_data <- function(species) {
  out <- callr::r(function(species) {
    e <- new.env()
    if (species == "human") {
      pk <- "quantiseqr"; utils::data(list = "dataset_racle", package = pk, envir = e)
      m <- get("dataset_racle", e)$expr_mat
    } else {
      pk <- "SSMD"; utils::data(list = "example_bulk", package = pk, envir = e)
      m <- get("example_bulk", e)
      colnames(m) <- paste0("S", seq_len(ncol(m)))
    }
    m <- m[!duplicated(rownames(m)), , drop = FALSE]
    data.frame(gene = rownames(m), as.data.frame(m, check.names = FALSE),
               check.names = FALSE, stringsAsFactors = FALSE)
  }, args = list(species = species), libpath = c(libs, .libPaths()))
  rownames(out) <- NULL
  out
}

human <- c("quantiseq", "mcp_counter", "xcell", "estimate", "consensus_tme", "immucellai2",
           "imsig", "spatialdecon_safetme", "epic", "xcell2", "abis", "timer")
mouse <- c("mmcp_counter", "dcq", "immucellai_mouse", "xcell2", "base")

test_that("human methods run on real data and return samples x features", {
  local_method_library()
  d <- get_data("human")
  for (m in human) {
    if (!bulkdeconv:::bd_is_installed(m)) next
    res <- if (m == "timer") suppressWarnings(bulkdeconv(d, "tpm", "human", m, cancer_type = "LUAD"))
           else suppressWarnings(bulkdeconv(d, "tpm", "human", m))
    expect_s3_class(res, "bulkdeconv_result")
    expect_equal(rownames(res$score), colnames(d)[-1], info = m)
    expect_true(ncol(res$score) >= 3L, info = m)
    expect_false(anyNA(res$score), info = m)
  }
})

test_that("mouse methods run on real data", {
  local_method_library()
  d <- get_data("mouse")
  d6 <- d[, 1:7]
  for (m in mouse) {
    if (!bulkdeconv:::bd_is_installed(m)) next
    res <- if (m == "base") suppressWarnings(bulkdeconv(d6, "tpm", "mouse", m, n_permutations = 5))
           else suppressWarnings(bulkdeconv(d6, "tpm", "mouse", m))
    expect_equal(rownames(res$score), colnames(d6)[-1], info = m)
    expect_true(ncol(res$score) >= 5L, info = m)
  }
})

test_that("seqimmucc runs on counts, alone and in a cohort", {
  local_method_library()
  skip_if(!bulkdeconv:::bd_is_installed("seqimmucc"))
  d <- get_data("mouse")[, 1:21]
  d[, -1] <- round(d[, -1])
  res <- suppressWarnings(bulkdeconv(d, "count", "mouse", "seqimmucc"))
  expect_equal(ncol(res$score), 10L)
  expect_equal(unname(rowSums(res$score)), rep(100, 20), tolerance = 1e-6)
  one <- suppressWarnings(bulkdeconv(d[, 1:2], "count", "mouse", "seqimmucc"))
  expect_equal(one$score[1, ], res$score[1, ], tolerance = 1e-8)
  expect_false(one$qc$single_sample_workaround)
  expect_equal(one$metadata$engine, "native")
})

test_that("ssmd runs (slow: about two minutes; set BULKDECONV_TEST_SLOW=1)", {
  local_method_library()
  skip_if(!nzchar(Sys.getenv("BULKDECONV_TEST_SLOW")), "BULKDECONV_TEST_SLOW is not set")
  skip_if(!bulkdeconv:::bd_is_installed("ssmd"))
  d <- get_data("mouse")[, 1:9]
  res <- suppressWarnings(bulkdeconv(d, "tpm", "mouse", "ssmd"))
  expect_equal(rownames(res$score), colnames(d)[-1])
})

test_that("a single sample works where the registry says it does", {
  local_method_library()
  d <- get_data("human")[, 1:2]
  for (m in c("quantiseq", "epic", "xcell2", "spatialdecon_safetme", "mcp_counter", "consensus_tme", "estimate", "immucellai2")) {
    if (!bulkdeconv:::bd_is_installed(m)) next
    res <- if (m == "timer") suppressWarnings(bulkdeconv(d, "tpm", "human", m, cancer_type = "LUAD"))
           else suppressWarnings(bulkdeconv(d, "tpm", "human", m))
    expect_equal(nrow(res$score), 1L, info = m)
  }
  dm <- get_data("mouse")[, 1:2]
  for (m in c("mmcp_counter", "xcell2")) {
    res <- suppressWarnings(bulkdeconv(dm, "tpm", "mouse", m))
    expect_equal(nrow(res$score), 1L, info = m)
  }
})

test_that("rsemverse-style output (non-integer values) is accepted", {
  local_method_library()
  d <- get_data("human")
  d[, -1] <- d[, -1] + 0.37
  skip_if(!bulkdeconv:::bd_is_installed("quantiseq"))
  expect_no_error(suppressWarnings(bulkdeconv(d, "tpm", "human", "quantiseq")))
})
