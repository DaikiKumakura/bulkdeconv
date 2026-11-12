#!/usr/bin/env bash
# Run on Ubuntu 24.04 with R 4.5.x. Pass an existing check library to reuse it.
set -euo pipefail
package_dir=$(cd "$(dirname "$0")/.." && pwd)
stage=$(mktemp -d /tmp/bulkdeconv-linux.XXXXXXXX)
printf 'Linux check directory: %s\n' "$stage"
lib=${1:-"$stage/lib"}
mkdir -p "$lib"
lib=$(cd "$lib" && pwd)
cd "$stage"
cp -a "$package_dir" bulkdeconv
export R_LIBS="$lib" R_LIBS_USER="$lib" R_LIBS_SITE="$lib"
Rscript --vanilla -e '
  stopifnot(getRversion() >= "4.5.0", getRversion() < "4.6.0")
  options(timeout = 300, Ncpus = 4L,
          HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
            paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))
  deps <- c("callr", "testthat", "remotes", "BiocManager", "quadprog", "pracma")
  missing <- setdiff(deps, rownames(installed.packages()))
  if (length(missing)) {
    install.packages(missing, repos = "https://packagemanager.posit.co/cran/__linux__/noble/latest")
  }
  stopifnot(all(vapply(deps, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1L))))
'
R CMD build --no-build-vignettes bulkdeconv
export _R_CHECK_FORCE_SUGGESTS_=false _R_CHECK_CRAN_INCOMING_REMOTE_=false
R CMD check --no-manual bulkdeconv_0.1.0.tar.gz
