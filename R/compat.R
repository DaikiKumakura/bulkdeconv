# Compatibility between a method and the request. Returns the plan: which
# single-sample strategy applies, and the validated method options.

bd_check_method <- function(row, species, expr, n_samples, opts) {
  m <- row$method
  if (row$status == "excluded") {
    bd_abort("Method '", m, "' is excluded: ", row$notes)
  }
  if (row$status == "unavailable") {
    bd_abort("Method '", m, "' cannot be run: ", row$notes)
  }
  if (!species %in% bd_split(row$species)) {
    ok <- bulkdeconv_methods(species)$method
    bd_abort("Method '", m, "' supports ", gsub("|", " and ", row$species, fixed = TRUE),
             ", not ", species, ". Methods for ", species, ": ", paste(ok, collapse = ", "), ".")
  }
  accepted <- bd_split(row$expr_accepted)
  if (!expr %in% accepted) {
    bd_abort("Method '", m, "' accepts expr = ", paste0("\"", accepted, "\"", collapse = " or "),
             ", not \"", expr, "\". ",
             if (expr == "count") "Convert counts to TPM (rsemverse::rsem_read(expr = \"tpm\")) instead." else "")
  }

  plan <- list(duplicate = FALSE, warnings = character())
  if (n_samples == 1L && row$single_sample == "duplicate") {
    plan$duplicate <- TRUE
    plan$warnings <- c(plan$warnings, paste0(
      "'", m, "' cannot run on one sample directly; the single sample was duplicated internally ",
      "and the copy was dropped from the result."))
  } else if (n_samples < row$min_samples_n) {
    bd_abort("Method '", m, "' needs at least ", row$min_samples_n, " sample(s); `data` has ",
             n_samples, ".")
  }
  if (identical(row$cohort_dependent, "yes") && n_samples < 10L) {
    plan$warnings <- c(plan$warnings, paste0(
      "'", m, "' output depends on the other samples in the run; with ", n_samples,
      " samples the values are unstable and should not be compared with other runs."))
  }
  if (nzchar(row$license) && grepl("non-commercial", row$license, fixed = TRUE)) {
    plan$warnings <- c(plan$warnings, paste0("'", m, "' licence: ", row$license, "."))
  }
  if (identical(m, "base")) {
    plan$warnings <- c(plan$warnings, "BASE runs permutations; time grows with n_permutations (default 100). The permutation seed is fixed (option `seed`).")
  }
  if (identical(m, "abis")) {
    plan$warnings <- c(plan$warnings, "ABIS signatures come from sorted blood cells, not tumor tissue.")
  }
  if (identical(m, "ssmd")) {
    plan$warnings <- c(plan$warnings, "SSMD has no declared licence; check the authors' terms before redistributing results.")
  }
  plan$opts <- bd_check_options(m, opts, species)
  plan
}

# Allowed options per method, with defaults and value checks.
bd_tcga_types <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH",
                   "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG",
                   "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "UCEC", "UCS", "UVM")

# An option has either `values` (a set of strings) or `range` (whole numbers).
# `default = NULL` makes the option required.
bd_option_spec <- function(method, species) {
  switch(method,
    timer = list(cancer_type = list(
      default = NULL,
      # TIMER's 32 types: the TCGA types without LAML, plus THYM
      values = c(setdiff(bd_tcga_types, "LAML"), "THYM"))),
    base = list(n_permutations = list(default = 100L, range = c(1L, 1000L)),
                seed = list(default = 1L, range = c(0L, 1000000L))),
    xcell2 = list(preset = list(
      default = if (species == "human") "TMECompendium" else "MouseRNAseqData",
      values = {
        p <- bd_read_registry("presets.csv")
        p$preset[p$species == species]
      })),
    consensus_tme = list(cancer_type = list(
      default = "Unfiltered",
      values = c("Unfiltered", bd_tcga_types))),
    immucellai2 = list(tissue = list(default = "tumor", values = c("tumor", "non-tumor"))),
    ssmd = list(tissue = list(default = "Inflammatory",
                              values = c("Inflammatory", "Central Nervous System",
                                         "Hematopoietic System", "Blood"))),
    list())
}

bd_check_options <- function(method, opts, species) {
  spec <- bd_option_spec(method, species)
  if (length(opts) && (is.null(names(opts)) || any(!nzchar(names(opts))))) {
    bd_abort("Method options passed through `...` must be named.")
  }
  unknown <- setdiff(names(opts), names(spec))
  if (length(unknown)) {
    bd_abort("Unknown option(s) for '", method, "': ", paste(unknown, collapse = ", "), ". ",
             if (length(spec)) paste0("Allowed: ", paste(names(spec), collapse = ", "), ".")
             else "This method takes no options.")
  }
  out <- list()
  for (nm in names(spec)) {
    val <- opts[[nm]] %||% spec[[nm]]$default
    if (is.null(val)) {
      bd_abort("Option `", nm, "` is required for '", method, "'. One of: ",
               paste(spec[[nm]]$values, collapse = ", "), ".")
    }
    if (!is.null(spec[[nm]]$range)) {
      r <- spec[[nm]]$range
      if (!is.numeric(val) || length(val) != 1L || is.na(val) || val != round(val) ||
          val < r[1L] || val > r[2L]) {
        bd_abort("Option `", nm, "` for '", method, "' must be a whole number from ", r[1L], " to ", r[2L], ".")
      }
      out[[nm]] <- as.integer(val)
      next
    }
    if (!is.character(val) || length(val) != 1L || !val %in% spec[[nm]]$values) {
      bd_abort("Option `", nm, "` for '", method, "' must be one of: ",
               paste(spec[[nm]]$values, collapse = ", "), ".")
    }
    out[[nm]] <- val
  }
  out
}
