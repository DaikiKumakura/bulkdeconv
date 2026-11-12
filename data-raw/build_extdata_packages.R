# Builds the parts of inst/extdata/ that come from installed R packages (their
# data objects): estimate, spatialdecon, immucellai_mouse, consensus_tme, xcell,
# mmcp_counter; and MCP-counter's marker table, downloaded from its repository at
# the pinned commit. The other parts are built by data-raw/build_extdata.R.
# Sources and licences: inst/COPYRIGHTS.
#
#   Rscript data-raw/build_extdata_packages.R <library-with-the-packages> [...]
#
# Needs: tidyestimate, SpatialDecon, ImmuCellAImouse, ConsensusTME, xCell, mMCPcounter, imsig, quantiseqr.

args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args, .libPaths()))
root <- if (file.exists("DESCRIPTION")) "." else stop("run from the package root")
out <- file.path(root, "inst", "extdata")
ld <- function(nm, pkg) { e <- new.env(); utils::data(list = nm, package = pkg, envir = e); get(nm, e) }
mk <- function(d) dir.create(file.path(out, d), showWarnings = FALSE, recursive = TRUE)
num_matrix <- function(m) { r <- rownames(m); m <- apply(m, 2, as.numeric); rownames(m) <- r; m }

mk("estimate")
saveRDS(list(common_genes = ld("common_genes", "tidyestimate"), gene_sets = ld("gene_sets", "tidyestimate")),
        file.path(out, "estimate", "estimate_data.rds"), compress = "xz")
mk("spatialdecon")
saveRDS(ld("safeTME", "SpatialDecon"), file.path(out, "spatialdecon", "safeTME.rds"), compress = "xz")
mk("immucellai_mouse")
saveRDS(list(marker_exp = ld("marker_exp", "ImmuCellAImouse"),
             l1_marker = ld("l1_marker", "ImmuCellAImouse"), l2_marker = ld("l2_marker", "ImmuCellAImouse"),
             l3_marker = ld("l3_marker", "ImmuCellAImouse"),
             l1_correction = num_matrix(ld("l1_cell_correction_matrix_new", "ImmuCellAImouse")),
             l2_correction = num_matrix(ld("l2_cell_correction_matrix_new", "ImmuCellAImouse")),
             l3_correction = num_matrix(ld("l3_cell_correction_matrix_new", "ImmuCellAImouse"))),
        file.path(out, "immucellai_mouse", "reference.rds"), compress = "xz")

mk("consensus_tme")
saveRDS(ld("consensusGeneSets", "ConsensusTME"), file.path(out, "consensus_tme", "gene_sets.rds"), compress = "xz")

mk("xcell")
xd <- ld("xCell.data", "xCell")
sigs <- xd$signatures
sigs <- if (is.list(sigs) && !inherits(sigs, "GeneSetCollection")) sigs else {
  suppressMessages(library(GSEABase)); setNames(lapply(sigs, GSEABase::geneIds), vapply(sigs, GSEABase::setName, ""))
}
saveRDS(list(signatures = sigs, genes = xd$genes, fit_values = xd$spill$fv, spill = xd$spill$K),
        file.path(out, "xcell", "xcell_data.rds"), compress = "xz")

mk("mcp_counter")
f <- tempfile()
utils::download.file("https://raw.githubusercontent.com/ebecht/MCPcounter/b6eac73e91c246fcff0bb1a5c68a816cd588fc48/Signatures/genes.txt", f, quiet = TRUE, mode = "wb")
gn <- utils::read.table(f, sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)
saveRDS(gn, file.path(out, "mcp_counter", "markers.rds"), compress = "xz")

mk("imsig")
saveRDS(ld("sig", "imsig"), file.path(out, "imsig", "sig.rds"), compress = "xz")

mk("quantiseq")
qd <- system.file("extdata", package = "quantiseqr")
rd <- function(f, ...) utils::read.table(file.path(qd, f), header = FALSE, sep = "	", stringsAsFactors = FALSE, ...)
sg <- utils::read.table(file.path(qd, "TIL10_signature.txt"), header = TRUE, sep = "	", row.names = 1)
sc <- rd("TIL10_mRNA_scaling.txt")
hgnc <- utils::read.csv(file.path(qd, "HGNC_genenames_20170418.txt.gz"), header = TRUE, sep = "	", stringsAsFactors = FALSE)
saveRDS(list(signature = as.matrix(sg), scaling = setNames(sc[[2]], sc[[1]]), abgenes = rd("TIL10_TCGA_aberrant_immune_genes.txt")[[1]],
             rmgenes = rd("TIL10_rmgenes.txt")[[1]],
             hgnc = hgnc[, c("ApprovedSymbol", "Status", "ApprovedName", "PreviousSymbols")]),
        file.path(out, "quantiseq", "til10.rds"), compress = "xz")

mk("mmcp_counter")
saveRDS(ld("mMCPcounter_signatures_GCRm39", "mMCPcounter"), file.path(out, "mmcp_counter", "signatures_GCRm39.rds"), compress = "xz")

for (f in list.files(out, recursive = TRUE, full.names = TRUE)) cat(sprintf("%-55s %6.2f MB\n", sub(paste0(out, "/"), "", f, fixed = TRUE), file.size(f) / 1e6))
