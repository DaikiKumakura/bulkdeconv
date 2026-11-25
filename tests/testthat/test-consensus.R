fake_result <- function(method, score, species = "human") {
  structure(list(score = score, metadata = list(method = method, species = species)),
            class = "bulkdeconv_result")
}
mk <- function(n = 10) {
  set.seed(3)
  s <- paste0("S", seq_len(n))
  z <- seq_len(n) + rnorm(n, 0, 0.1)        # CD8 abundance
  p <- rev(seq_len(n)) + rnorm(n, 0, 0.1)   # tumor purity, opposite order
  list(
    quantiseq = fake_result("quantiseq", cbind(T.cells.CD8 = z / 100, B.cells = runif(n), Other = runif(n)) |> `rownames<-`(s)),
    mcp_counter = fake_result("mcp_counter", cbind(`CD8 T cells` = z^2, `B lineage` = runif(n)) |> `rownames<-`(s)),
    estimate = fake_result("estimate", cbind(estimate = -p) |> `rownames<-`(s)),
    puree = fake_result("puree", cbind(tumor_purity = p / 100) |> `rownames<-`(s))
  )
}

test_that("lineages are ranked, combined and checked for agreement", {
  r <- mk()
  cons <- bulkdeconv_consensus(unname(r))
  expect_s3_class(cons, "bulkdeconv_consensus")
  expect_setequal(colnames(cons$score), c("T_CD8", "B_cell", "tumor_content"))
  expect_equal(dim(cons$score), c(10, 3))
  expect_true(all(cons$score > 0 & cons$score < 1))
  # T_CD8: both methods rank the samples identically
  expect_equal(unname(cons$agreement["T_CD8"]), 1)
  expect_equal(unname(cons$n_methods["T_CD8"]), 2L)
  expect_equal(unname(rank(cons$score[, "T_CD8"])), as.numeric(1:10))
  # tumor content: ESTIMATE is inverted, so the two sources agree
  expect_equal(unname(cons$agreement["tumor_content"]), 1)
  expect_equal(unname(rank(cons$score[, "tumor_content"])), as.numeric(10:1))
  expect_output(print(cons), "bulkdeconv_consensus")
})

test_that("unused outputs and lineages with one method are reported", {
  cons <- bulkdeconv_consensus(unname(mk()), min_methods = 2L)
  # B.cells is in quantiseq and B lineage in mcp_counter -> two methods, kept
  expect_true("B_cell" %in% colnames(cons$score))
  r <- mk(); r$mcp_counter <- NULL
  cons1 <- bulkdeconv_consensus(unname(r))
  expect_false("T_CD8" %in% colnames(cons1$score))
  expect_true(any(cons1$dropped$what == "lineage" & cons1$dropped$name == "T_CD8"))
  cons2 <- bulkdeconv_consensus(unname(r), min_methods = 1L)
  expect_true(is.na(cons2$agreement["T_CD8"]))
})

test_that("bad input is refused with a clear message", {
  r <- mk()
  expect_error(bulkdeconv_consensus(list(1, 2)), "list of results")
  expect_error(bulkdeconv_consensus(list(r$quantiseq, r$quantiseq)), "once")
  few <- lapply(r, function(x) { x$score <- x$score[1:3, , drop = FALSE]; x })
  expect_error(bulkdeconv_consensus(unname(few)), "at least 5 samples")
  other <- r; other$mcp_counter$score <- other$mcp_counter$score[10:1, , drop = FALSE]
  rownames(other$mcp_counter$score) <- paste0("X", 1:10)
  expect_error(bulkdeconv_consensus(unname(other)), "same samples")
  mouse <- r; mouse$mcp_counter$metadata$species <- "mouse"
  expect_error(bulkdeconv_consensus(unname(mouse)), "same species")
  expect_error(bulkdeconv_consensus(list(r$quantiseq, fake_result("timer", r$quantiseq$score))), "No lineage")
})

test_that("the lineage table only names real outputs of known methods", {
  tab <- bulkdeconv_lineages()
  expect_true(all(tab$method %in% bulkdeconv_methods()$method))
  expect_true(all(tab$direction %in% c(1L, -1L)))
  expect_false(anyDuplicated(tab[c("method", "cell_type")]) > 0)
})
