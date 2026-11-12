test_that("all xCell2 presets reproduce the original scores through bulkdeconv", {
  fixtures <- readRDS(test_path("fixtures", "xcell2-golden.rds"))
  presets <- utils::read.csv(system.file("registry", "presets.csv", package = "bulkdeconv"))
  fixtures <- c(fixtures, readRDS(test_path("fixtures", "xcell2-LM22-golden.rds")))
  expect_setequal(names(fixtures), presets$preset)
  for (preset in names(fixtures)) {
    g <- fixtures[[preset]]
    species <- presets$species[match(preset, presets$preset)]
    got <- suppressWarnings(bulkdeconv(g$X, "tpm", species, "xcell2", preset = preset))
    expect_identical(dimnames(got$score), dimnames(g$expected), info = preset)
    expect_true(max(abs(got$score - g$expected)) < 1e-6, info = preset)
    expect_identical(got$metadata$engine, "native")
  }
})

test_that("xCell2 keeps the single-sample workaround and enforces reference coverage", {
  g <- readRDS(test_path("fixtures", "xcell2-golden.rds"))$TMECompendium
  one <- g$X[, 1L, drop = FALSE]
  got <- suppressWarnings(bulkdeconv(one, "tpm", "human", "xcell2"))
  expect_identical(rownames(got$score), colnames(one))
  expect_true(got$qc$single_sample_workaround)
  expect_true(all(got$score == 0))
  expect_error(
    suppressWarnings(bulkdeconv(one[1:100, , drop = FALSE], "tpm", "human", "xcell2")),
    "check the species and the gene symbols"
  )
})
