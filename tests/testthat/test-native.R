# The natively implemented methods must reproduce what the original
# implementations give (fixtures made by data-raw/make_fixtures.R).

golden <- readRDS(test_path("fixtures", "native-golden.rds"))
ext <- function(method, filename) system.file("extdata", method, filename, package = "bulkdeconv", mustWork = TRUE)
run_native <- function(fn, ...) {
  f <- bulkdeconv:::bd_native[[fn]]
  environment(f) <- list2env(bulkdeconv:::bd_native, parent = globalenv())
  f(...)
}

test_that("seqimmucc (LLSR) reproduces the original output", {
  g <- golden$seqimmucc
  sig <- utils::read.table(ext("seqimmucc", "sig_matr_seqImmuCC.txt"), header = TRUE, sep = "\t",
                           row.names = 1, check.names = FALSE)
  got <- run_native("llsr", g$X, sig)
  expect_equal(dim(got), dim(g$expected))
  expect_equal(colnames(got), colnames(g$expected))
  expect_equal(got, g$expected, tolerance = 1e-8)
  expect_equal(unname(rowSums(got)), rep(100, nrow(got)), tolerance = 1e-8)
})

test_that("seqimmucc gives a sample the same values alone and in a cohort", {
  g <- golden$seqimmucc
  sig <- utils::read.table(ext("seqimmucc", "sig_matr_seqImmuCC.txt"), header = TRUE, sep = "\t",
                           row.names = 1, check.names = FALSE)
  one <- run_native("llsr", g$X[, 1, drop = FALSE], sig)
  expect_equal(one[1, ], run_native("llsr", g$X, sig)[1, ], tolerance = 1e-10)
})

test_that("abis reproduces the original output", {
  skip_if_not_installed("MASS")
  g <- golden$abis
  sig <- utils::read.delim(ext("abis", "sigmatrixRNAseq.txt"), check.names = FALSE)
  got <- run_native("abis_fit", g$X, sig)
  expect_equal(dim(got), dim(g$expected))
  expect_equal(got, g$expected, tolerance = 1e-8)
  one <- run_native("abis_fit", g$X[, 1, drop = FALSE], sig)
  expect_equal(one[1, ], got[1, ], tolerance = 1e-10)
})

test_that("base reproduces the original output (same permutations)", {
  g <- golden$base
  comp <- readRDS(ext("base", "compendium.rds"))
  map <- utils::read.table(ext("base", "cell_correspondance.txt"), sep = "\t", header = TRUE)
  got <- run_native("base_scores", g$X, comp, map, perm = g$perm, seed = g$seed)
  expect_equal(dim(got), dim(g$expected))
  expect_equal(colnames(got), colnames(g$expected))
  # the compendium is stored rounded to 6 decimals, hence the tolerance
  expect_equal(got, g$expected, tolerance = 1e-3)
  expect_gt(stats::cor(as.vector(got), as.vector(g$expected)), 0.9999)
})

test_that("base is reproducible with a fixed seed and changes with another", {
  g <- golden$base
  comp <- readRDS(ext("base", "compendium.rds"))
  map <- utils::read.table(ext("base", "cell_correspondance.txt"), sep = "\t", header = TRUE)
  a <- run_native("base_scores", g$X, comp, map, perm = 5L, seed = 3L)
  b <- run_native("base_scores", g$X, comp, map, perm = 5L, seed = 3L)
  c <- run_native("base_scores", g$X, comp, map, perm = 5L, seed = 4L)
  expect_identical(a, b)
  expect_false(isTRUE(all.equal(a, c)))
})

test_that("timer reproduces the original output", {
  skip_if_not_installed("sva")
  g <- golden$timer
  ref <- readRDS(ext("timer", "immune_reference.rds"))
  markers <- readRDS(ext("timer", "markers_by_cancer.rds"))[[g$cancer]]
  got <- suppressMessages(run_native("timer_fit", g$X, ref, markers))
  expect_equal(dim(got), dim(g$expected))
  expect_equal(colnames(got), colnames(g$expected))
  expect_equal(got, g$expected, tolerance = 1e-6)
})

test_that("the bundled reference data are what the registry says", {
  expect_length(readRDS(ext("timer", "markers_by_cancer.rds")), 32L)
  expect_true(nzchar(system.file("COPYRIGHTS", package = "bulkdeconv")))
})

test_that("estimate reproduces the original output", {
  g <- golden$estimate
  d <- readRDS(ext("estimate", "estimate_data.rds"))
  got <- run_native("estimate_fit", g$X, d$common_genes, d$gene_sets)
  expect_equal(colnames(got), colnames(g$expected))
  expect_equal(rownames(got), rownames(g$expected))
  expect_equal(got, g$expected, tolerance = 1e-8)
})

test_that("spatialdecon (SAFE-TME) reproduces converged reference scores", {
  g <- golden$spatialdecon
  safe <- readRDS(ext("spatialdecon", "safeTME.rds"))
  got <- suppressWarnings(run_native("spatialdecon_fit", g$X, safe))
  expect_equal(dim(got), dim(g$expected))
  expect_equal(colnames(got), colnames(g$expected))
  expect_true(all(is.finite(got) & got >= 0))
  # LAU125 converges in both stages on Windows and Linux. The full fixture
  # records the upstream Windows fit, including iteration-limited samples.
  # The reference was recorded on Windows. On other platforms the optimiser
  # (L-BFGS-B) stops at a point that differs in the 4th-5th significant digit
  # (largest difference seen on Ubuntu and macOS: 8e-5), hence the tolerance.
  expect_equal(got["LAU125", ], g$expected["LAU125", ],
               tolerance = if (.Platform$OS.type == "windows") 1e-5 else 1e-3)
  if (.Platform$OS.type == "windows") {
    expect_equal(got, g$expected, tolerance = 1e-5)
  }
  reversed <- g$X[rev(seq_len(nrow(g$X))), , drop = FALSE]
  expect_identical(suppressWarnings(run_native("spatialdecon_fit", reversed, safe)), got)
  expect_warning(
    expect_warning(run_native("spatialdecon_fit", g$X, safe, maxit = 1L),
                   "initial fit did not converge"),
    "outlier refit did not converge"
  )
})

test_that("immucellai2 reproduces the original output (tumor and non-tumor)", {
  g <- golden$immucellai2
  d <- readRDS(ext("immucellai2", "reference_markers.rds"))
  got <- run_native("immucellai2_fit", g$X, d$reference, d$markers, tissue = "tumor")
  expect_equal(colnames(got), colnames(g$expected))
  expect_equal(got, g$expected, tolerance = 1e-8)
  nt <- run_native("immucellai2_fit", g$X, d$reference, d$markers, tissue = "non-tumor")
  expect_equal(colnames(nt), colnames(g$expected_nontumor))
  expect_equal(nt, g$expected_nontumor, tolerance = 1e-8)
  expect_equal(unname(rowSums(got)), rep(1, nrow(got)), tolerance = 1e-8)
})

test_that("immucellai_mouse reproduces the original output", {
  g <- golden$immucellai_mouse
  ref <- readRDS(ext("immucellai_mouse", "reference.rds"))
  got <- run_native("immucellai_mouse_fit", g$X, ref)
  expect_equal(colnames(got), colnames(g$expected))
  expect_equal(rownames(got), rownames(g$expected))
  # values are rounded to 3-4 decimals in the method: allow one unit in the last place
  expect_equal(got, g$expected, tolerance = 2e-3)
  expect_gt(stats::cor(as.vector(got), as.vector(g$expected)), 0.9999)
})

test_that("the ssGSEA helper agrees with GSVA", {
  skip_if_not_installed("GSVA")
  set.seed(1)
  X <- matrix(rexp(800 * 5, 0.1), 800, 5, dimnames = list(paste0("g", 1:800), paste0("s", 1:5)))
  X[sample(length(X), 300)] <- 0
  sets <- list(a = sample(rownames(X), 40), b = sample(rownames(X), 15), c = sample(rownames(X), 100))
  ref <- GSVA::gsva(GSVA::ssgseaParam(X, sets))
  got <- run_native("ssgsea", X, sets)
  expect_equal(unname(got), unname(matrix(as.numeric(as.matrix(ref)), nrow(ref))), tolerance = 1e-10)
})

test_that("the NNLS helper agrees with pracma::lsqlincon and is non-negative", {
  skip_if_not_installed("pracma")
  skip_if_not_installed("quadprog")
  set.seed(2)
  C <- diag(6) + matrix(runif(36, 0, 0.3), 6)
  for (i in 1:5) {
    d <- runif(6, -0.5, 2)
    x <- run_native("nnls", C, d)
    expect_true(all(x >= 0))
    expect_equal(x, as.vector(pracma::lsqlincon(C, d, lb = 0)), tolerance = 1e-6)
  }
})

# GPL-3 originals ported to R/native3.R. Columns are compared by name (the
# original orders cell types with the locale's sort()).
by_name <- function(x, ref) x[, colnames(ref), drop = FALSE]

test_that("consensus_tme reproduces the original output (two cancer types)", {
  g <- golden$consensus_tme
  sets <- readRDS(ext("consensus_tme", "gene_sets.rds"))
  for (ct in names(g$expected)) {
    got <- suppressWarnings(run_native("consensus_tme_fit", g$X, sets, ct))
    expect_setequal(colnames(got), colnames(g$expected[[ct]]))
    expect_equal(unname(by_name(got, g$expected[[ct]])), unname(matrix(as.numeric(g$expected[[ct]]), nrow(got))), tolerance = 1e-8, info = ct)
  }
})

test_that("xcell reproduces the original output", {
  g <- golden$xcell
  d <- readRDS(ext("xcell", "xcell_data.rds"))
  got <- run_native("xcell_fit", g$X, d$signatures, d$genes, d$fit_values, d$spill)
  expect_setequal(colnames(got), colnames(g$expected))
  expect_equal(by_name(got, g$expected), g$expected, tolerance = 1e-6)
})

test_that("mcp_counter reproduces the original output", {
  g <- golden$mcp_counter
  got <- suppressWarnings(run_native("mcp_counter_fit", g$X, readRDS(ext("mcp_counter", "markers.rds"))))
  expect_setequal(colnames(got), colnames(g$expected))
  expect_equal(by_name(got, g$expected), g$expected, tolerance = 1e-10)
  one <- suppressWarnings(run_native("mcp_counter_fit", g$X[, 1, drop = FALSE], readRDS(ext("mcp_counter", "markers.rds"))))
  expect_equal(one[1, ], got[1, colnames(one)], tolerance = 1e-10)
})

test_that("mmcp_counter reproduces the original output", {
  g <- golden$mmcp_counter
  got <- suppressWarnings(run_native("mmcp_counter_fit", g$X, readRDS(ext("mmcp_counter", "signatures_GCRm39.rds"))))
  expect_equal(colnames(got), colnames(g$expected))
  expect_equal(got, g$expected, tolerance = 1e-10)
})

test_that("imsig reproduces the original output", {
  g <- golden$imsig
  got <- run_native("imsig_fit", g$X, readRDS(ext("imsig", "sig.rds")))
  expect_setequal(colnames(got), colnames(g$expected))
  exp <- g$expected[rownames(got), colnames(got)]
  expect_equal(got, exp, tolerance = 1e-8)
})

test_that("quantiseq reproduces the original output", {
  skip_if_not_installed("quadprog")
  g <- golden$quantiseq
  got <- run_native("quantiseq_fit", g$X, readRDS(ext("quantiseq", "til10.rds")))
  expect_setequal(colnames(got), colnames(g$expected))
  expect_equal(got[, colnames(g$expected)], g$expected, tolerance = 1e-5)
  expect_equal(unname(rowSums(got)), rep(1, nrow(got)), tolerance = 1e-8)
})
