# Run from the package root, with an optional output directory as the first argument.
stage <- tempfile("bulkdeconv-pkgdown-source-")
dir.create(stage)
source_files <- c("_pkgdown.yml", "DESCRIPTION", "NAMESPACE", "README.md", "NEWS.md",
                  "LICENSE.md", "R", "man")
stopifnot(all(file.copy(source_files, stage, recursive = TRUE)))
args <- commandArgs(trailingOnly = TRUE)
destination <- if (length(args)) args[[1L]] else "_site"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
destination <- normalizePath(destination, winslash = "/")
stopifnot(!length(list.files(destination, all.files = TRUE, no.. = TRUE)))

pkgdown::build_site(pkg = stage, preview = FALSE, install = TRUE, new_process = FALSE,
                    override = list(destination = destination))
pages <- c("index.html", "reference/index.html", "reference/bulkdeconv.html", "news/index.html")
stopifnot(all(file.exists(file.path(destination, pages))),
          !file.exists(file.path(destination, "AGENTS.html")))
cat("Verified pkgdown preview:", destination, "\n")
