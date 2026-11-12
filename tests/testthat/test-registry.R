test_that("the registry loads and holds the 18 usable methods", {
  reg <- bulkdeconv_methods()
  expect_setequal(reg$method, c(
    "quantiseq", "mcp_counter", "xcell", "estimate", "consensus_tme", "immucellai2", "imsig",
    "spatialdecon_safetme", "epic", "xcell2", "abis", "timer", "mmcp_counter", "dcq",
    "immucellai_mouse", "ssmd", "seqimmucc", "base"))
  expect_type(reg$installed, "logical")
})

test_that("unusable methods are listed only on request", {
  all <- bulkdeconv_methods(usable = FALSE)
  expect_true(all(c("immucellai", "immucc") %in% all$method))
  expect_false(any(c("immucellai", "immucc") %in% bulkdeconv_methods()$method))
})

test_that("species filter works", {
  expect_true(all(grepl("mouse", bulkdeconv_methods("mouse")$species)))
  expect_true("xcell2" %in% bulkdeconv_methods("mouse")$method)
  expect_false("quantiseq" %in% bulkdeconv_methods("mouse")$method)
  expect_error(bulkdeconv_methods("rat"), "human")
})

test_that("every usable method has an adapter and a parsable install spec", {
  reg <- bulkdeconv_methods()
  for (m in reg$method) {
    expect_true(is.function(bulkdeconv:::bd_adapter_fns[[m]]), info = m)
    specs <- bulkdeconv:::bd_split(reg$install[reg$method == m], ";")
    native <- reg$engine[reg$method == m] == "native"
    if (!native) expect_gt(length(specs), 0)
    for (s in specs) {
      p <- bulkdeconv:::bd_parse_spec(s)
      expect_true(p$type %in% c("cran", "bioc", "github"), info = s)
      if (p$type == "github") expect_match(p$ref, "^[0-9a-f]{40}$", info = s)
    }
  }
})

test_that("install specs are parsed", {
  p <- bulkdeconv:::bd_parse_spec("github::GuoBioinfoLab/ImmuCellAI-2.0@abc#R-package/ImmuCellAI2.0")
  expect_equal(p$name, "GuoBioinfoLab/ImmuCellAI-2.0")
  expect_equal(p$ref, "abc")
  expect_equal(p$subdir, "R-package/ImmuCellAI2.0")
  p <- bulkdeconv:::bd_parse_spec("bioc::quantiseqr")
  expect_null(p$ref); expect_null(p$subdir)
})

test_that("unknown methods suggest close names", {
  expect_error(bulkdeconv:::bd_method_row("quantiseqq"), "Did you mean: quantiseq")
  expect_error(bulkdeconv:::bd_method_row("nope"), "Unknown method")
})

test_that("the private library can be set by option", {
  withr::local_options(bulkdeconv.lib = c(tempdir(), tempdir()))
  expect_length(bulkdeconv:::bd_lib_paths(), 2L)
})

test_that("native methods need nothing beyond the packages bulkdeconv itself imports", {
  reg <- bulkdeconv_methods()
  nat <- reg[reg$engine == "native", ]
  expect_setequal(nat$method, c("abis", "timer", "seqimmucc", "base", "estimate", "immucellai2",
                                "spatialdecon_safetme", "immucellai_mouse", "consensus_tme", "xcell",
                                "mcp_counter", "mmcp_counter", "imsig", "quantiseq", "xcell2"))
  expect_true(all(nat$install == ""))
  used <- unlist(strsplit(nat$packages, "|", fixed = TRUE))
  imports <- read.dcf(system.file("DESCRIPTION", package = "bulkdeconv"), "Imports")[1L]
  for (p in used) expect_match(imports, p, fixed = TRUE)
  expect_true(all(nat$installed))
  expect_setequal(reg$method[reg$engine == "external"], c("epic", "dcq", "ssmd"))
})
