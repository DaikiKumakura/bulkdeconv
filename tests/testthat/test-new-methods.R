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

test_that("puree ranks over the universe, imputes missing genes and applies the linear model", {
  d <- list(universe = paste0("ENSG", 1:6), genes = c("ENSG2", "ENSG5", "ENSG6"), coef = c(0.5, -0.2, 1),
            intercept = 0.1, medians = c(0.4, 0.6, 0.9),
            conversion = data.frame(ENSEMBL_ID = paste0("ENSG", 1:6), HGNC = paste0("G", 1:6)))
  X <- matrix(c(5, 1, 3, 3, 9, 0,  1, 2, 3, 4, 5, 6), 6, 2, dimnames = list(paste0("G", 1:6), c("a", "b")))
  X <- X[-6, , drop = FALSE]  # G6 missing: imputed with its median
  got <- run_native("puree_fit", X, d)
  rk_a <- c(1, 5) / 5
  expect_equal(got["a", "tumor_purity"], 0.1 + 0.5 * rk_a[1] - 0.2 * rk_a[2] + 1 * 0.9)
  expect_equal(got["b", "tumor_purity"], 0.1 + 0.5 * 2 / 5 - 0.2 * 5 / 5 + 1 * 0.9)
  ens <- X; rownames(ens) <- paste0("ENSG", 1:5, ".3")
  expect_equal(run_native("puree_fit", ens, d), got)
  expect_equal(run_native("puree_fit", X[, 1, drop = FALSE], d), got[1, , drop = FALSE])
})

test_that("bulkdeconv_install('puree') downloads, checks and reads the model", {
  skip_on_cran()
  skip_if_offline("raw.githubusercontent.com")
  lib <- file.path(tempdir(), "puree-lib")
  bulkdeconv_install("puree", lib = lib, ask = FALSE)
  d <- readRDS(file.path(lib, "_bulkdeconv_data", "puree", "puree.rds"))
  expect_length(d$genes, 158L)
  expect_length(d$coef, 158L)
  expect_length(d$medians, 158L)
  expect_equal(d$intercept, 0.780836832988957, tolerance = 1e-14)
  expect_true(all(d$genes %in% d$universe))
})

test_that("bindea gene sets are the 24 published immune cell types", {
  sets <- readRDS(ext("bindea", "gene_sets.rds"))
  expect_length(sets, 24L)
  expect_identical(length(unique(unlist(sets))), 501L)
  expect_setequal(sets$DC, c("CCL13", "CCL17", "CCL22", "CD209", "HSD11B1", "NPR1", "PPFIBP2"))
  expect_true("MARCH6" %in% sets$`NK CD56bright cells`)
  expect_false(any(c("Angiogenesis", "Antigen presenting machinery") %in% names(sets)))
})

test_that("mysort recovers the composition of mixtures of its signature", {
  S <- readRDS(ext("mysort", "signature_matrix.rds"))
  expect_identical(dim(S), c(603L, 21L))
  ad <- bulkdeconv:::bd_get_adapter("mysort")
  set.seed(3)
  P <- matrix(stats::rexp(ncol(S) * 4), ncol(S)); P <- t(t(P) / colSums(P))
  X <- S %*% P; colnames(X) <- paste0("s", 1:4)
  r <- ad(X, list(), list(extdata = dirname(ext("mysort", "signature_matrix.rds")) |> dirname()))$score
  expect_equal(unname(rowSums(r)), rep(100, 4), tolerance = 1e-8)
  expect_gt(stats::cor(as.vector(r), as.vector(t(P)) * 100), 0.8)
})

test_that("the native quantile normalisation equals preprocessCore, ties included", {
  skip_if_not_installed("preprocessCore")
  set.seed(1)
  M <- matrix(stats::rpois(1000 * 4, 3), 1000, 4)
  f <- bulkdeconv:::bd_isolate(function(M) quantile_normalize(M))
  expect_equal(f(M), preprocessCore::normalize.quantiles(M), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("the primal nu-SVR matches the exact quadratic-programming solution", {
  skip_if_not_installed("quadprog")
  set.seed(7)
  n <- 60; r <- 4
  S <- matrix(stats::rnorm(n * r), n, r)
  y <- as.vector(S %*% c(1, -0.5, 0.3, 0.8)) + stats::rnorm(n, 0, 0.6); y <- (y - mean(y)) / stats::sd(y)
  K <- tcrossprod(S)
  f <- bulkdeconv:::bd_isolate(function(S, y, nu) nusvr_linear_w(S, y, nu))
  for (nu in c(0.25, 0.75)) {
    D <- rbind(cbind(K, -K), cbind(-K, K)) + diag(1e-8, 2 * n)
    A <- cbind(c(rep(1, n), rep(-1, n)), rep(1, 2 * n), diag(2 * n), -diag(2 * n))
    sol <- quadprog::solve.QP(D, c(y, -y), A, c(0, n * nu, rep(0, 2 * n), rep(-1, 2 * n)), meq = 2)$solution
    w_qp <- as.vector(crossprod(S, sol[1:n] - sol[(n + 1):(2 * n)]))
    expect_lt(max(abs(f(S, y, nu) - w_qp)), 1e-5)
  }
})

test_that("mysort stays close to libsvm (e1071) on the same input", {
  skip_if_not_installed("e1071")
  set.seed(2)
  S <- readRDS(ext("mysort", "signature_matrix.rds"))
  P <- matrix(stats::rexp(ncol(S) * 3), ncol(S)); P <- t(t(P) / colSums(P))
  X <- S %*% P + matrix(stats::rnorm(nrow(S) * 3, 0, 0.05 * stats::sd(S)), nrow(S)); colnames(X) <- paste0("s", 1:3)
  f <- bulkdeconv:::bd_isolate(function(S, y, nu) nusvr_linear_w(S, y, nu))
  Z <- (S - mean(S)) / stats::sd(S)
  y <- X[, 1]; y <- (y - mean(y)) / stats::sd(y)
  m <- e1071::svm(Z, y, scale = FALSE, type = "nu-regression", kernel = "linear", nu = 0.5)
  expect_lt(max(abs(f(Z, y, 0.5) - as.vector(t(m$coefs) %*% m$SV))), 5e-3)
})
