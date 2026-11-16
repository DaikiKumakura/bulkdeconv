# Build the bundled data of three signature-based methods from their
# published sources (downloaded at pinned locations):
#   gedit       : GEDIT reference matrices (Nadel et al. 2021, MIT), presets
#                 HPCA (Human Primary Cell Atlas, 25 cell types incl. stroma)
#                 and ImmunoStates (20 immune cell types)
#   danaher     : marker genes of 14 immune cell populations (Danaher et al.
#                 2017, J Immunother Cancer, Additional file 1, Table S4)
#   charoentong : pan-cancer immune metagenes of 28 immune cell types
#                 (Charoentong et al. 2017, Cell Reports, Table S6)
#
#   Rscript data-raw/build_signatures.R
#
# Writes inst/extdata/gedit/<preset>.rds, inst/extdata/danaher/markers.rds and
# inst/extdata/charoentong/gene_sets.rds.

GEDIT_SHA <- "be9b33f89842b860821a874b6cb2586bd0d2a050"
GEDIT_PRESETS <- c(HPCA = "HPCA_Recommended.csv", ImmunoStates = "ImmunoStates.tsv")
DANAHER_URL <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs40425-017-0215-8/MediaObjects/40425_2017_215_MOESM1_ESM.xlsx"
CHAROENTONG_URL <- "https://ars.els-cdn.com/content/image/1-s2.0-S2211124716317090-mmc3.xlsx"

tmp <- tempfile(); dir.create(tmp)
get <- function(url, name) {
  f <- file.path(tmp, name)
  utils::download.file(url, f, mode = "wb", quiet = TRUE, headers = c(`User-Agent` = "Mozilla/5.0"))
  f
}

# ---- GEDIT ----
dir.create(file.path("inst", "extdata", "gedit"), recursive = TRUE, showWarnings = FALSE)
for (p in names(GEDIT_PRESETS)) {
  f <- get(sprintf("https://raw.githubusercontent.com/BNadel/GEDIT/%s/ReferenceMatrices/%s", GEDIT_SHA, GEDIT_PRESETS[[p]]),
           GEDIT_PRESETS[[p]])
  d <- utils::read.delim(f, sep = if (grepl("\\.csv$", f)) "," else "\t", check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(d[, -1]); rownames(m) <- d[[1]]; storage.mode(m) <- "double"
  stopifnot(!anyDuplicated(rownames(m)), all(m >= 0, na.rm = TRUE), !anyNA(m))
  saveRDS(m, file.path("inst", "extdata", "gedit", paste0(p, ".rds")), compress = "xz")
  cat(sprintf("gedit %s: %d genes x %d cell types\n", p, nrow(m), ncol(m)))
}

# ---- Danaher ----
dir.create(file.path("inst", "extdata", "danaher"), recursive = TRUE, showWarnings = FALSE)
s4 <- readxl::read_excel(get(DANAHER_URL, "danaher.xlsx"), sheet = "S4. Selected markers")
markers <- split(s4$Gene, s4$`Cell Type`)
stopifnot(length(markers) == 14, sum(lengths(markers)) == 60)
saveRDS(markers, file.path("inst", "extdata", "danaher", "markers.rds"), compress = "xz")
cat(sprintf("danaher: %d genes in %d populations\n", sum(lengths(markers)), length(markers)))

# ---- Charoentong ----
dir.create(file.path("inst", "extdata", "charoentong"), recursive = TRUE, showWarnings = FALSE)
s6 <- readxl::read_excel(get(CHAROENTONG_URL, "charoentong.xlsx"), skip = 2)
s6 <- s6[!is.na(s6$Metagene) & !is.na(s6$`Cell type`), ]
sets <- lapply(split(s6$Metagene, s6$`Cell type`), unique)
stopifnot(length(sets) == 28)
saveRDS(sets, file.path("inst", "extdata", "charoentong", "gene_sets.rds"), compress = "xz")
cat(sprintf("charoentong: %d genes in %d cell types\n", length(unique(unlist(sets))), length(sets)))
