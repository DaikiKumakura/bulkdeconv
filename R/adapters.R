# Adapters: one function per method. Each runs in a separate R process, in the
# private library, and sees only its arguments plus the helpers below (it must
# not refer to bulkdeconv's namespace). Contract:
#   adapter(X, opts, ctx) -> list(score = samples x features numeric matrix,
#                                 coverage = NULL or c(signature = n, found = n))
# X is genes x samples (symbols as row names); ctx$cache_dir is a writable dir.

bd_helpers <- list(
  reference_file = function(ctx, method, filename) {
    path <- file.path(ctx$extdata, method, filename)
    if (!utils::file_test("-f", path)) {
      stop("Reference file '", file.path(method, filename), "' is missing from the installed ",
           "package; reinstall bulkdeconv.", call. = FALSE)
    }
    path
  },
  as_score = function(x) {
    if (is.list(x) && !is.data.frame(x)) {
      cand <- Filter(function(z) is.matrix(z) || is.data.frame(z), x)
      if (!length(cand)) stop("method returned no matrix", call. = FALSE)
      x <- cand[[1L]]
    }
    if (inherits(x, "data.frame")) {
      x <- as.matrix(x[, vapply(x, is.numeric, logical(1)), drop = FALSE])
    }
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L, dimnames = list(NULL, names(x)))
    storage.mode(x) <- "double"
    x
  },
  coverage = function(sig_genes, X) {
    sig_genes <- unique(sig_genes)
    c(signature = length(sig_genes), found = sum(sig_genes %in% rownames(X)))
  }
)

bd_adapter_fns <- list(
  quantiseq = function(X, opts, ctx) {
    ref <- readRDS(file.path(ctx$extdata, "quantiseq", "til10.rds"))
    list(score = quantiseq_fit(X, ref), coverage = coverage(rownames(ref$signature), X))
  },
  mcp_counter = function(X, opts, ctx) {
    mk <- readRDS(file.path(ctx$extdata, "mcp_counter", "markers.rds"))
    list(score = mcp_counter_fit(X, mk), coverage = coverage(mk[["HUGO symbols"]], X))
  },
  xcell = function(X, opts, ctx) {
    d <- readRDS(file.path(ctx$extdata, "xcell", "xcell_data.rds"))
    list(score = xcell_fit(X, d$signatures, d$genes, d$fit_values, d$spill), coverage = coverage(d$genes, X))
  },
  xcell2 = function(X, opts, ctx) {
    ref <- readRDS(reference_file(ctx, "xcell2", paste0(opts$preset, ".rds")))
    cov <- coverage(ref$genes, X)
    if (cov[["found"]] / cov[["signature"]] < 0.5) {
      stop("only ", round(100 * cov[["found"]] / cov[["signature"]]), "% of the genes of preset '",
           opts$preset, "' are in the data; check the species and the gene symbols", call. = FALSE)
    }
    list(score = xcell2_fit(X, ref), coverage = cov)
  },
  estimate = function(X, opts, ctx) {
    d <- readRDS(file.path(ctx$extdata, "estimate", "estimate_data.rds"))
    list(score = estimate_fit(X, d$common_genes, d$gene_sets),
         coverage = coverage(unique(d$common_genes$hgnc_symbol), X))
  },
  consensus_tme = function(X, opts, ctx) {
    sets <- readRDS(file.path(ctx$extdata, "consensus_tme", "gene_sets.rds"))
    list(score = consensus_tme_fit(X, sets, opts$cancer_type), coverage = coverage(unlist(sets[[opts$cancer_type]]), X))
  },
  immucellai2 = function(X, opts, ctx) {
    d <- readRDS(file.path(ctx$extdata, "immucellai2", "reference_markers.rds"))
    list(score = immucellai2_fit(X, d$reference, d$markers, tissue = opts$tissue),
         coverage = coverage(d$markers, X))
  },
  imsig = function(X, opts, ctx) {
    sig <- readRDS(file.path(ctx$extdata, "imsig", "sig.rds"))
    list(score = imsig_fit(X, sig), coverage = coverage(as.character(sig$gene), X))
  },
  spatialdecon_safetme = function(X, opts, ctx) {
    safe <- readRDS(file.path(ctx$extdata, "spatialdecon", "safeTME.rds"))
    list(score = spatialdecon_fit(X, safe), coverage = coverage(rownames(safe), X))
  },
  epic = function(X, opts, ctx) {
    list(score = as_score(EPIC::EPIC(bulk = as.matrix(X))$cellFractions))
  },
  mmcp_counter = function(X, opts, ctx) {
    sig <- readRDS(file.path(ctx$extdata, "mmcp_counter", "signatures_GCRm39.rds"))
    list(score = mmcp_counter_fit(X, sig), coverage = coverage(sig$Gene.Symbol, X))
  },
  dcq = function(X, opts, ctx) {
    e <- new.env()
    utils::data("commons", package = "ComICS", envir = e)
    r <- ComICS::dcq(get("immgen_dat", e), as.data.frame(X), get("DCQ_mar", e))
    list(score = as_score(r$average))
  },
  immucellai_mouse = function(X, opts, ctx) {
    ref <- readRDS(file.path(ctx$extdata, "immucellai_mouse", "reference.rds"))
    list(score = immucellai_mouse_fit(X, ref))
  },
  abis = function(X, opts, ctx) {
    sig <- utils::read.delim(reference_file(ctx, "abis", "sigmatrixRNAseq.txt"), check.names = FALSE)
    list(score = abis_fit(X, sig), coverage = coverage(rownames(sig), X))
  },
  timer = function(X, opts, ctx) {
    ref <- readRDS(file.path(ctx$extdata, "timer", "immune_reference.rds"))
    markers <- readRDS(file.path(ctx$extdata, "timer", "markers_by_cancer.rds"))[[opts$cancer_type]]
    list(score = timer_fit(X, ref, markers), coverage = coverage(markers, X))
  },
  seqimmucc = function(X, opts, ctx) {
    sig <- utils::read.table(reference_file(ctx, "seqimmucc", "sig_matr_seqImmuCC.txt"), header = TRUE,
                             sep = "	", row.names = 1, check.names = FALSE)
    list(score = llsr(X, sig), coverage = coverage(rownames(sig), X))
  },
  base = function(X, opts, ctx) {
    comp_file <- reference_file(ctx, "base", "compendium.rds")
    map_file <- reference_file(ctx, "base", "cell_correspondance.txt")
    comp <- readRDS(comp_file)
    map <- utils::read.table(map_file, sep = "	", header = TRUE)
    list(score = base_scores(X, comp, map, perm = opts$n_permutations, seed = opts$seed),
         coverage = coverage(rownames(comp), X))
  },
  ssmd = function(X, opts, ctx) {
    r <- SSMD::SSMD(X, opts$tissue)
    list(score = t(as_score(r$Proportion)))
  }
)

# Give an adapter a clean environment (parent: globalenv) holding the helpers,
# so serialising it to another R process never drags in bulkdeconv's namespace.
bd_isolate <- function(f) {
  env <- new.env(parent = globalenv())
  fns <- c(bd_helpers, bd_native)
  for (nm in names(fns)) assign(nm, fns[[nm]], envir = env)
  for (nm in names(fns)) environment(env[[nm]]) <- env
  environment(f) <- env
  f
}

bd_get_adapter <- function(method) {
  f <- bd_adapter_fns[[method]]
  if (is.null(f)) bd_abort("No adapter is defined for method '", method, "'.")
  bd_isolate(f)
}
