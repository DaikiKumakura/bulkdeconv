# GEDIT, Danaher and Charoentong (data from data-raw/build_signatures.R).

ext <- function(method, filename) system.file("extdata", method, filename, package = "bulkdeconv", mustWork = TRUE)
run_native <- function(fn, ...) {
  f <- bulkdeconv:::bd_native[[fn]]
  environment(f) <- list2env(bulkdeconv:::bd_native, parent = globalenv())
  f(...)
}

test_that("bundled signature data have the published sizes", {
  expect_identical(dim(readRDS(ext("gedit", "HPCA.rds"))), c(19715L, 25L))
  expect_identical(dim(readRDS(ext("gedit", "ImmunoStates.rds"))), c(317L, 20L))
  mk <- readRDS(ext("danaher", "markers.rds"))
  expect_length(mk, 14L)
  expect_identical(sum(lengths(mk)), 60L)
  sets <- readRDS(ext("charoentong", "gene_sets.rds"))
  expect_length(sets, 28L)
  expect_identical(length(unique(unlist(sets))), 782L)
})

test_that("gedit recovers the proportions of noise-free mixtures of its reference", {
  ref <- readRDS(ext("gedit", "HPCA.rds"))
  set.seed(1)
  P <- matrix(stats::runif(ncol(ref) * 6), ncol(ref))
  P <- t(t(P) / colSums(P)); rownames(P) <- colnames(ref)
  X <- ref %*% P; colnames(X) <- paste0("s", 1:6)
  est <- run_native("gedit_fit", X, ref)
  expect_identical(dimnames(est), list(colnames(X), colnames(ref)))
  expect_equal(unname(rowSums(est)), rep(1, 6), tolerance = 1e-12)
  expect_true(all(est >= 0))
  expect_gt(stats::cor(as.vector(est), as.vector(t(P))), 0.9)
  expect_lt(max(abs(est - t(P))), 0.05)
})

test_that("gedit matches gene names case-insensitively and needs enough shared genes", {
  ref <- readRDS(ext("gedit", "ImmunoStates.rds"))
  X <- ref %*% rep(1 / ncol(ref), ncol(ref)); colnames(X) <- "s1"
  rownames(X) <- tolower(rownames(X))
  expect_equal(run_native("gedit_fit", X, ref), run_native("gedit_fit", `rownames<-`(X, toupper(rownames(X))), ref))
  expect_error(run_native("gedit_fit", X[1:5, , drop = FALSE], ref), "genes of the GEDIT reference")
})

test_that("danaher scores are mean log2 expression of the markers and sample independent", {
  mk <- list(A = c("G1", "G2"), B = "G3", C = "MISSING")
  X <- matrix(c(1, 3, 7, 0, 15, 1), 3, 2, dimnames = list(c("G1", "G2", "G3"), c("s1", "s2")))
  s <- run_native("danaher_fit", X, mk)
  expect_equal(s[, "A"], c(s1 = mean(log2(c(2, 4))), s2 = mean(log2(c(1, 16)))))
  expect_equal(s[, "B"], c(s1 = 3, s2 = 1))
  expect_true(all(is.na(s[, "C"])))
  expect_equal(run_native("danaher_fit", X[, 1, drop = FALSE], mk), s[1, , drop = FALSE])
})

test_that("the three methods run through bulkdeconv()", {
  mk <- readRDS(ext("danaher", "markers.rds"))
  ref <- readRDS(ext("gedit", "ImmunoStates.rds"))
  genes <- unique(c(rownames(ref), unlist(mk), unlist(readRDS(ext("charoentong", "gene_sets.rds")))))
  set.seed(2)
  X <- matrix(stats::rexp(length(genes) * 3, 0.1), length(genes), 3, dimnames = list(genes, paste0("s", 1:3)))
  for (m in c("gedit", "danaher", "charoentong")) {
    r <- suppressWarnings(bulkdeconv(X, "tpm", "human", m))  # charoentong warns: cohort-dependent with few samples
    expect_s3_class(r, "bulkdeconv_result")
    expect_identical(rownames(r$score), colnames(X))
  }
  expect_identical(ncol(bulkdeconv(X, "tpm", "human", "gedit", preset = "ImmunoStates")$score), 20L)
  expect_error(bulkdeconv(X, "tpm", "human", "gedit", preset = "LM22"), "must be one of")
})
