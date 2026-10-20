#!/usr/bin/env bash
# Build and verify an R/Bioconductor package.
#
# Usage, from the directory CONTAINING the package directory:
#   bash build_check.sh <PackageDir> [--fast] [--lib <dir>]
#
#   --fast          roxygen -> install -> testthat only
#   --no-vignettes  skip vignette building (INTERIM ONLY - a working
#                   vignette is required for Bioconductor submission;
#                   use this to verify everything else while fixing pandoc)
#   --lib DIR   library directory to install into (default: ./Rlib)
#
# Writes check.log and bioccheck.log beside the package and prints the
# status line of each. Exits non-zero if R CMD check reports an ERROR or
# BiocCheck reports an ERROR; environmental warnings and notes do not
# fail the run, because they are not package defects (see SKILL.md).

set -uo pipefail

PKG="${1:-}"

# Accept being run from INSIDE the package directory, which is where an
# RStudio project puts you, as well as from the parent with the package
# named. R CMD build/check must run from the parent, so chdir there.
if [ -z "$PKG" ] || [ ! -d "$PKG" ]; then
    if [ -f "DESCRIPTION" ] && [ -d "R" ]; then
        PKG="$(basename "$PWD")"
        cd .. || exit 2
        echo "(running from inside the package; switched to $PWD)"
    else
        echo "usage: bash build_check.sh [PackageDir] [--fast] [--no-vignettes] [--lib <dir>]" >&2
        echo "  run it from inside the package directory, or from the" >&2
        echo "  parent with the package directory named." >&2
        exit 2
    fi
else
    shift
fi

FAST=0
NOVIG=0
LIB="$PWD/Rlib"
while [ $# -gt 0 ]; do
    case "$1" in
        --fast) FAST=1; shift ;;
        --no-vignettes) NOVIG=1; shift ;;
        --lib)  LIB="$2"; shift 2 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

mkdir -p "$LIB"
# R CMD subprocesses read R_LIBS; Rscript code below also sets .libPaths().
export R_LIBS="$LIB:$(Rscript -e 'cat(.libPaths()[1])')"
# Suggests that cannot be installed must not fail the check.
export _R_CHECK_FORCE_SUGGESTS_=false

PKGNAME=$(awk '/^Package:/ {print $2; exit}' "$PKG/DESCRIPTION")
VERSION=$(awk '/^Version:/ {print $2; exit}' "$PKG/DESCRIPTION")
TARBALL="${PKGNAME}_${VERSION}.tar.gz"
echo "== $PKGNAME $VERSION (library: $LIB)"

# A source tree reconstructed by untarring a BUILT tarball carries fields that
# R CMD build generated (Author, Maintainer, Packaged). BiocCheck errors when
# Author/Maintainer sit alongside Authors@R, so catch it before the long run.
if grep -q "^Authors@R:" "$PKG/DESCRIPTION" && \
   grep -qE "^(Author|Maintainer|Packaged): " "$PKG/DESCRIPTION"; then
    echo "   DESCRIPTION carries build-generated fields alongside Authors@R:"
    grep -nE "^(Author|Maintainer|Packaged): " "$PKG/DESCRIPTION" | sed 's/^/     /'
    echo "   Remove them (BiocCheck reports them as ERRORs):"
    echo "     sed -i.bak -E '/^(Author|Maintainer|Packaged): /d' $PKG/DESCRIPTION"
    exit 1
fi

echo "== environment"
# Fail early and legibly: an absent roxygen2 otherwise dies with an R
# traceback, and an R older than the package requires wastes a long install.
PKG="$PKG" Rscript - <<'RCHECK' || exit 1
pkg <- Sys.getenv("PKG")
desc <- read.dcf(file.path(pkg, "DESCRIPTION"))
dep <- if ("Depends" %in% colnames(desc)) desc[1, "Depends"] else ""
need <- sub(".*R \\(>= *([0-9.]+)\\).*", "\\1", gsub("[\r\n]+", " ", dep))
cat("   R ", as.character(getRversion()), sep = "")
if (grepl("^[0-9.]+$", need)) {
    cat("  (package requires >= ", need, ")", sep = "")
    if (getRversion() < need) {
        cat("\n   BLOCKER: this R is older than the package requires.\n")
        cat("   Upgrade R before installing anything else.\n")
        quit(status = 1)
    }
}
cat("\n")
installed <- rownames(installed.packages())

## The package's own hard dependencies, taken from DESCRIPTION so this
## stays correct as the package evolves. Base packages ship with R.
fields <- intersect(c("Depends", "Imports"), colnames(desc))
decl <- unlist(strsplit(paste(desc[1, fields], collapse = ","), ","))
decl <- trimws(sub("\\(.*", "", decl))
base <- c("R", "methods", "stats", "utils", "grDevices", "graphics",
          "tools", "parallel", "grid", "splines", "stats4")
deps <- setdiff(decl[nzchar(decl)], base)

tools <- c("roxygen2", "testthat")
missTools <- setdiff(tools, installed)
missDeps <- setdiff(deps, installed)

if (length(missTools) > 0L || length(missDeps) > 0L) {
    if (length(missTools) > 0L) {
        cat("   missing verification tool(s): ",
            paste(missTools, collapse = ", "), "\n", sep = "")
    }
    if (length(missDeps) > 0L) {
        cat("   missing package dependencies: ",
            paste(missDeps, collapse = ", "), "\n", sep = "")
    }
    cat("   install everything at once with:\n\n")
    cat('     if (!requireNamespace("BiocManager", quietly = TRUE))\n')
    cat('         install.packages("BiocManager")\n')
    cat("     BiocManager::install(c(",
        paste0('"', c(missTools, missDeps), '"', collapse = ", "),
        "))\n\n", sep = "")
    quit(status = 1)
}
cat("   toolchain and ", length(deps), " dependencies present\n", sep = "")
RCHECK

echo "== roxygen2"
Rscript -e "suppressMessages(roxygen2::roxygenise('$PKG'))" || exit 1

echo "== R CMD INSTALL"
# Install rather than pkgload::load_all(): load_all can corrupt S4
# dispatch for classes extending Bioconductor containers.
R CMD INSTALL --no-docs -l "$LIB" "$PKG" > install.log 2>&1 || {
    echo "   INSTALL FAILED - see install.log"; tail -20 install.log; exit 1; }
echo "   ok"

if [ -d "$PKG/tests/testthat" ]; then
    echo "== testthat"
    Rscript -e "
      .libPaths(c('$LIB', .libPaths()))
      suppressMessages({library(testthat); library($PKGNAME)})
      res <- test_dir('$PKG/tests/testthat', package='$PKGNAME',
                      reporter='summary', stop_on_failure=FALSE)
      df <- as.data.frame(res)
      cat('   failures:', sum(df\$failed), ' errors:', sum(df\$error),
          ' warnings:', sum(df\$warning), '\n')
    " 2>&1 | grep -E "^(   |[a-z0-9_-]+: )" | tail -20
fi

if [ "$FAST" -eq 1 ]; then
    echo "== fast mode: stopping before build/check"
    exit 0
fi

BUILD_ARGS=""
CHECK_ARGS=""
if [ "$NOVIG" -eq 1 ]; then
    BUILD_ARGS="--no-build-vignettes"
    CHECK_ARGS="--no-vignettes --no-build-vignettes"
    echo "   NOTE: skipping vignettes. R CMD check will report two extra"
    echo "         WARNINGs ('files in vignettes', 'package vignettes')."
    echo "         This is interim only - Bioconductor requires a vignette"
    echo "         that builds, so a clean run must not use this flag."
fi

echo "== R CMD build"
# Stale *.BiocCheck / *.Rcheck folders inside the package are themselves
# a BiocCheck ERROR, and would otherwise be committed to the repository.
rm -rf "${PKGNAME}.Rcheck" "${PKGNAME}.BiocCheck" "$TARBALL" \
       "$PKG"/*.Rcheck "$PKG"/*.BiocCheck "$PKG/build"
R CMD build $BUILD_ARGS "$PKG" > build.log 2>&1 || {
    echo "   BUILD FAILED - see build.log"; tail -20 build.log; exit 1; }
echo "   $TARBALL"

echo "== R CMD check --as-cran"
R CMD check --no-manual --as-cran $CHECK_ARGS "$TARBALL" > check.log 2>&1
grep "^Status" check.log || echo "   (no status line - see check.log)"
# Match wrapped lines too: a long check name pushes WARNING onto its own line.
grep -B1 -E "(WARNING|ERROR)$" check.log | grep -E "checking|WARNING|ERROR" \
    | sed 's/^/   /' | head -20

RC=0
if grep -q "^Status.*ERROR" check.log; then RC=1; fi

if Rscript -e ".libPaths(c('$LIB', .libPaths())); quit(status = !requireNamespace('BiocCheck', quietly=TRUE))"; then
    echo "== BiocCheck (on the package directory)"
    Rscript -e "
      .libPaths(c('$LIB', .libPaths()))
      # The directory is the target Bioconductor CI and reviewers use; it
      # catches DESCRIPTION issues a tarball run does not.
      invisible(BiocCheck::BiocCheck('$PKG',
          \`no-check-R-ver\`=TRUE, \`no-check-bioc-help\`=TRUE,
          \`no-check-remotes\`=TRUE))
    " > bioccheck.log 2>&1
    if grep -qE "ERRORS \|" bioccheck.log; then
        grep -E "ERRORS \|" bioccheck.log | sed 's/^/   /'
        grep -E "^! WARNING" bioccheck.log | sed 's/^/   /'
        if grep -qE "x [1-9][0-9]* ERRORS" bioccheck.log; then RC=1; fi
    else
        # BiocCheck produced no verdict line, so it died before finishing.
        # Surface that instead of printing nothing at all.
        echo "   BiocCheck did not complete - no summary line produced."
        echo "   Last lines of bioccheck.log:"
        tail -15 bioccheck.log | sed 's/^/     /'
        RC=1
    fi
else
    echo "== BiocCheck not installed - skipped"
    echo "   install with: conda install -c bioconda bioconductor-bioccheck"
fi

echo "== logs: check.log bioccheck.log build.log install.log"
exit $RC
