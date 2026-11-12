# Builds tests/testthat/fixtures/native-golden.rds: small inputs and the output
# the ORIGINAL implementations (omnideconv/immunedeconv at the commit in
# data-raw/build_registry.py) give for them. test-native.R checks that the
# natively implemented methods in bulkdeconv reproduce these numbers.
# immunedeconv is used here only to create the fixtures; the package itself
# never calls it.
#
#   Rscript data-raw/make_fixtures.R <library-with-immunedeconv-quantiseqr-SSMD>

args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args, .libPaths()))
suppressMessages(library(immunedeconv))
root <- if (file.exists("DESCRIPTION")) "." else stop("run from the package root")
extdata <- file.path(root, "inst", "extdata")
set.seed(20261111)

e <- new.env(); utils::data(list = "dataset_racle", package = "quantiseqr", envir = e)
H <- as.matrix(get("dataset_racle", e)$expr_mat); H <- H[!duplicated(rownames(H)), ]
e2 <- new.env(); utils::data(list = "example_bulk", package = "SSMD", envir = e2)
M <- as.matrix(get("example_bulk", e2)); colnames(M) <- paste0("S", seq_len(ncol(M))); M <- M[!duplicated(rownames(M)), ]

quiet <- function(expr) { zz <- file(tempfile(), "wt"); sink(zz); sink(zz, type = "message"); on.exit({ sink(type = "message"); sink(); close(zz) }); tryCatch(force(expr), error = function(e) { sink(type = "message"); sink(); message("ERROR inside quiet(): ", conditionMessage(e)); stop(e) }) }
pick <- function(have, want, n_extra = 0L) {
  g <- intersect(have, want)
  extra <- setdiff(have, g)
  sort(c(g, if (n_extra > 0L) sample(extra, min(n_extra, length(extra)))))
}

# seqimmucc (LLSR): counts, 20 mouse samples, genes of the signature
seq_sig <- read.table(file.path(extdata, "seqimmucc", "sig_matr_seqImmuCC.txt"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
g <- pick(rownames(M), rownames(seq_sig), 200L)
Xl <- round(M[g, 1:20])
ref_l <- quiet(deconvolute_seqimmucc(Xl, "LLSR"))
ref_l <- t(ref_l[!rownames(ref_l) %in% c("Correlation", "RSEM"), ])

# abis: human TPM, 4 samples, genes of the signature
asig <- read.delim(file.path(extdata, "abis", "sigmatrixRNAseq.txt"), check.names = FALSE)
Xa <- H[pick(rownames(H), rownames(asig), 200L), ]
ref_a <- t(quiet(deconvolute_abis(Xa, arrays = FALSE)))

# base: mouse TPM, 6 samples, a subset of the compendium genes, 20 permutations
comp <- readRDS(file.path(extdata, "base", "compendium.rds"))
gb <- pick(rownames(M), sample(rownames(comp), 5000L))
Xb <- M[gb, 1:6]
set.seed(1)
ref_b <- t(quiet(deconvolute_base_algorithm(Xb, n_permutations = 20)))

# timer: human TPM, 4 samples, LUAD, marker genes plus a random subset of the reference genes
tref <- readRDS(file.path(extdata, "timer", "immune_reference.rds"))
mk <- readRDS(file.path(extdata, "timer", "markers_by_cancer.rds"))$LUAD
gt <- pick(rownames(H), c(mk, sample(rownames(tref$genes), 3000L)))
Xt <- H[gt, ]
ref_t <- t(quiet(deconvolute_timer(Xt, indications = rep("luad", ncol(Xt)))))

# ---- methods ported from MIT / GPL-2+ packages: the originals are run here ----
suppressMessages({ library(tidyestimate); library(SpatialDecon); library(ImmuCellAImouse); library(ImmuCellAI2.0) })
as_num <- function(x) { x <- as.data.frame(x); as.matrix(x[, vapply(x, is.numeric, logical(1)), drop = FALSE]) }

# estimate: human TPM, 4 samples, genes of the ESTIMATE list plus a random subset of the rest
ecd <- readRDS(file.path(extdata, "estimate", "estimate_data.rds"))
ge <- pick(rownames(H), c(ecd$common_genes$hgnc_symbol, ecd$common_genes$external_synonym), 500L)
Xe <- H[ge, ]
d <- quiet(tidyestimate::filter_common_genes(as.data.frame(Xe), id = "hgnc_symbol", tidy = FALSE, tell_missing = FALSE, find_alias = TRUE))
ref_e <- quiet(as_num(tidyestimate::estimate_score(d, is_affymetrix = FALSE)))
rownames(ref_e) <- colnames(Xe)

# spatialdecon (SAFE-TME): human TPM, 4 samples, genes of the signature plus others
safe <- readRDS(file.path(extdata, "spatialdecon", "safeTME.rds"))
gs <- pick(rownames(H), rownames(safe), 300L)
Xs <- H[gs, ]
ref_s <- t(quiet(SpatialDecon::spatialdecon(norm = Xs, bg = matrix(0.1, nrow(Xs), ncol(Xs), dimnames = dimnames(Xs)), X = safe, align_genes = TRUE)$beta))

# immucellai2: human TPM, 4 samples, marker genes of the method
ic2 <- readRDS(file.path(extdata, "immucellai2", "reference_markers.rds"))
Xi <- H[pick(rownames(H), ic2$markers), ]
ref_i <- quiet(ImmuCellAI2.0::run_immucellai2(Xi, tissue.type = "tumor")$subtype.fraction)
ref_i2 <- quiet(ImmuCellAI2.0::run_immucellai2(Xi, tissue.type = "non-tumor")$subtype.fraction)

# immucellai_mouse: mouse TPM, 8 samples (names without overlaps), marker genes plus others
icm <- readRDS(file.path(extdata, "immucellai_mouse", "reference.rds"))
gm <- pick(rownames(M), c(unlist(icm$l1_marker), unlist(icm$l2_marker), unlist(icm$l3_marker), rownames(icm$marker_exp)), 300L)
Xm <- M[gm, 1:8]; colnames(Xm) <- sprintf("S%02d", 1:8)
ref_m <- quiet(ImmuCellAImouse::ImmuCellAI_mouse(as.data.frame(Xm), "rnaseq", 0, 0)[[1]])

# ---- GPL-3 originals ported to R/native3.R ----
suppressMessages({ library(ConsensusTME); library(xCell); library(MCPcounter); library(mMCPcounter) })
# consensus_tme: human TPM, 4 samples; two cancer types
cons <- readRDS(file.path(extdata, "consensus_tme", "gene_sets.rds"))
Xc <- H[pick(rownames(H), unlist(cons), 2000L), ]
ref_c <- lapply(c(LUAD = "LUAD", Unfiltered = "Unfiltered"), function(ct)
  t(quiet(suppressWarnings(ConsensusTME::consensusTMEAnalysis(Xc, cancerType = ct, statMethod = "ssgsea")))))
# xcell: human TPM, 4 samples, genes of xCell's list
xd <- readRDS(file.path(extdata, "xcell", "xcell_data.rds"))
Xx <- H[pick(rownames(H), xd$genes), ]
ref_x <- t(quiet(xCell::xCellAnalysis(Xx, rnaseq = TRUE, parallel.sz = 1)))
# mcp_counter: human TPM, 4 samples, marker genes plus others
mcp <- readRDS(file.path(extdata, "mcp_counter", "markers.rds"))
Xp <- H[pick(rownames(H), mcp[["HUGO symbols"]], 300L), ]
ref_p <- t(quiet(MCPcounter::MCPcounter.estimate(Xp, featuresType = "HUGO_symbols", genes = mcp)))
# mmcp_counter: mouse TPM, 6 samples
msig <- readRDS(file.path(extdata, "mmcp_counter", "signatures_GCRm39.rds"))
Xmm <- M[pick(rownames(M), msig$Gene.Symbol, 300L), 1:6]
ref_mm <- t(quiet(suppressWarnings(mMCPcounter::mMCPcounter.estimate(Xmm, features = "Gene.Symbol", genomeVersion = "GCRm39"))))

# imsig: human, 12 samples (the 4 real ones plus 8 noisy copies), genes of the signature
suppressMessages(library(imsig))
isig <- readRDS(file.path(extdata, "imsig", "sig.rds"))
# (imsig looks up an object called `sig` along the search path: never define one here)
gi <- intersect(rownames(H), as.character(isig$gene))
Xi0 <- H[gi, ]
Xim <- cbind(Xi0, vapply(1:8, function(i) Xi0[, (i %% 4) + 1] * exp(rnorm(nrow(Xi0), 0, 0.3)), numeric(nrow(Xi0))))
colnames(Xim) <- c(colnames(Xi0), paste0("N", 1:8))
ref_im <- as.matrix(quiet(imsig::imsig(as.data.frame(Xim), r = 0.6, sort = FALSE)))
# quantiseq: human TPM, 4 samples, signature / removed genes plus random others
suppressMessages(library(quantiseqr))
til <- readRDS(file.path(extdata, "quantiseq", "til10.rds"))
Xq <- H[pick(rownames(H), c(rownames(til$signature), til$rmgenes, til$abgenes), 3000L), ]
ref_q <- as.matrix(quiet(as.data.frame(quantiseqr::run_quantiseq(Xq, "TIL10", FALSE, TRUE, TRUE)))[, -1])
rownames(ref_q) <- colnames(Xq)

out <- list(
  imsig = list(X = Xim, expected = ref_im),
  quantiseq = list(X = Xq, expected = ref_q),
  consensus_tme = list(X = Xc, expected = ref_c),
  xcell = list(X = Xx, expected = ref_x),
  mcp_counter = list(X = Xp, expected = ref_p),
  mmcp_counter = list(X = Xmm, expected = ref_mm),
  estimate = list(X = Xe, expected = ref_e),
  spatialdecon = list(X = Xs, expected = ref_s),
  immucellai2 = list(X = Xi, expected = ref_i, expected_nontumor = ref_i2),
  immucellai_mouse = list(X = Xm, expected = as.matrix(ref_m)),
  seqimmucc = list(X = Xl, expected = ref_l),
  abis = list(X = Xa, expected = ref_a),
  base = list(X = Xb, expected = ref_b, perm = 20L, seed = 1L),
  timer = list(X = Xt, expected = ref_t, cancer = "LUAD"),
  immunedeconv_version = as.character(utils::packageVersion("immunedeconv")))
dir.create(file.path(root, "tests", "testthat", "fixtures"), showWarnings = FALSE)
f <- file.path(root, "tests", "testthat", "fixtures", "native-golden.rds")
saveRDS(out, f, compress = "xz")
cat(sprintf("wrote %s (%.2f MB)\n", f, file.size(f) / 1e6))
for (n in c("seqimmucc", "abis", "base", "timer", "estimate", "spatialdecon", "immucellai2", "immucellai_mouse", "xcell", "mcp_counter", "mmcp_counter", "imsig", "quantiseq")) cat(sprintf("  %-10s input %s  expected %s\n", n, paste(dim(out[[n]]$X), collapse = "x"), paste(dim(out[[n]]$expected), collapse = "x")))
