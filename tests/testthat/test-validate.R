sys.source(test_path("fixtures", "input-helpers.R"), envir = environment())

test_that("a data frame with a gene column is accepted and converted", {
  d <- make_input()
  p <- bulkdeconv:::bd_prepare_input(d, "tpm", "human")
  expect_equal(dim(p$matrix), c(40L, 3L))
  expect_equal(rownames(p$matrix)[1], "GENE1")
  expect_equal(colnames(p$matrix), c("S1", "S2", "S3"))
})

test_that("a matrix with row names is accepted", {
  d <- make_input()
  m <- as.matrix(d[, -1]); rownames(m) <- d$gene
  p <- bulkdeconv:::bd_prepare_input(m, "tpm", "human")
  expect_equal(p$n_genes, 40L)
})

test_that("a tibble-like data frame with a numeric first column and row names is accepted", {
  d <- make_input()
  m <- d[, -1]; rownames(m) <- d$gene
  p <- bulkdeconv:::bd_prepare_input(m, "tpm", "human")
  expect_true(any(grepl("row names", p$notes)))
})

test_that("bad layouts fail with clear messages", {
  d <- make_input()
  expect_error(bulkdeconv:::bd_prepare_input(list(1), "tpm", "human"), "data frame")
  expect_error(bulkdeconv:::bd_prepare_input(d[, 1, drop = FALSE], "tpm", "human"), "at least one sample")
  expect_error(bulkdeconv:::bd_prepare_input(d[, -1], "tpm", "human"), "first column")
  m <- as.matrix(d[, -1])
  expect_error(bulkdeconv:::bd_prepare_input(m, "tpm", "human"), "row names")
})

test_that("gene symbol problems are errors", {
  d <- make_input(); d$gene[2] <- d$gene[1]
  expect_error(bulkdeconv:::bd_prepare_input(d, "tpm", "human"), "more than once")
  d <- make_input(); d$gene[2] <- ""
  expect_error(bulkdeconv:::bd_prepare_input(d, "tpm", "human"), "missing or empty")
  d <- make_input(); d$gene <- sprintf("ENSG%011d", seq_len(40))
  expect_error(bulkdeconv:::bd_prepare_input(d, "tpm", "human"), "Ensembl")
})

test_that("sample column problems are errors", {
  d <- make_input(); names(d)[3] <- "S1"
  expect_error(bulkdeconv:::bd_prepare_input(d, "tpm", "human"), "unique")
  d <- make_input(); d$S2 <- as.character(d$S2)
  expect_error(bulkdeconv:::bd_prepare_input(d, "tpm", "human"), "Non-numeric")
})

test_that("values are checked", {
  d <- make_input(); d$S1[3] <- NA
  expect_error(bulkdeconv:::bd_prepare_input(d, "tpm", "human"), "NA")
  d <- make_input(); d$S1[3] <- Inf
  expect_error(bulkdeconv:::bd_prepare_input(d, "tpm", "human"), "infinite")
  d <- make_input(); d$S1[3] <- -1
  expect_error(bulkdeconv:::bd_prepare_input(d, "tpm", "human"), "negative")
  expect_error(bulkdeconv:::bd_prepare_input(make_input(n_genes = 5), "tpm", "human"), "at least 10")
})

test_that("species is checked against the symbol case", {
  expect_error(bulkdeconv:::bd_prepare_input(make_input(species = "mouse"), "tpm", "human"), "mouse symbols")
  expect_error(bulkdeconv:::bd_prepare_input(make_input(species = "human"), "tpm", "mouse"), "upper case")
  expect_no_error(bulkdeconv:::bd_prepare_input(make_input(species = "mouse"), "tpm", "mouse"))
})

test_that("the public entry point rejects empty samples and numeric gene IDs", {
  d <- make_input()
  d$S2 <- 0
  expect_error(bulkdeconv(d, "tpm", "human", "xcell"), "only zeros: S2")

  d <- make_input()
  d$gene <- as.character(seq_len(nrow(d)))
  expect_error(bulkdeconv(d, "tpm", "human", "xcell"), "No gene IDs look like gene symbols")
})

test_that("non-integer counts are noted; log-like TPM warns", {
  d <- make_input(scale = 1000)
  p <- bulkdeconv:::bd_prepare_input(d, "count", "human")
  expect_match(p$notes[1], "not all integers")
  d <- make_input(scale = 3)
  d[, -1] <- pmin(d[, -1], 15)
  expect_warning(bulkdeconv:::bd_prepare_input(d, "tpm", "human"), "log-transformed")
})
