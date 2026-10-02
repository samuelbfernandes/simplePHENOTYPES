#!/usr/bin/env bash
# test-installed.sh — run the test suite against an INSTALLED copy of the package,
# the way R CMD check / CI sees it (no R/, docs/, vignettes/ or benchmarks/ next to
# tests/, nothing from the source tree on the search path).
#
#   bash dev/test-installed.sh               # build, install, run every test file
#   FILTER=grammar bash dev/test-installed.sh  # only test files matching a regex
#   KEEP=1 bash dev/test-installed.sh        # keep the temporary directory
#
# What it does
#   1. R CMD build (no vignette build, no manual) of the working tree into a temp
#      directory: only what R CMD build would ship (.Rbuildignore is honoured), and
#      the source tree is left untouched (no .o/.so/target written into it);
#   2. R CMD INSTALL --no-docs into a temporary library (compiles the Rust core:
#      allow several minutes the first time);
#   3. copies tests/ from the built tarball (data files next to tests/testthat
#      included) to a temporary directory away from the sources;
#   4. runs testthat::test_dir(package = "simplePHENOTYPES", load_package =
#      "installed") with NOT_CRAN=true, so the NOT_CRAN-gated tests run;
#   5. exits non-zero on any failed test or error (skips are fine).
# All temporary files live under mktemp -d in ${TMPDIR:-/tmp}.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PKG="simplePHENOTYPES"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/sp-installed.XXXXXX")"
cleanup() {
  if [ "${KEEP:-0}" = "1" ]; then
    echo "kept: $WORK"
  else
    rm -rf "$WORK"
  fi
}
trap cleanup EXIT

LIB="$WORK/lib"
BUILD="$WORK/build"
SRCX="$WORK/extracted"
TESTS="$WORK/tests"
mkdir -p "$LIB" "$BUILD" "$SRCX"

# All tools below must see the same temporary directory and nothing else.
export TMPDIR="$WORK"
export NOT_CRAN=true
unset R_TESTS || true

echo "== [1/4] R CMD build ($ROOT)"
( cd "$BUILD" && R CMD build --no-build-vignettes --no-manual "$ROOT" )
TARBALL="$(ls "$BUILD"/${PKG}_*.tar.gz | head -n 1)"
[ -f "$TARBALL" ] || { echo "no tarball was built" >&2; exit 2; }

echo "== [2/4] R CMD INSTALL --no-docs into $LIB"
R CMD INSTALL --no-docs --library="$LIB" "$TARBALL"

echo "== [3/4] copy tests/ away from the sources"
tar -xzf "$TARBALL" -C "$SRCX" "$PKG/tests"
cp -R "$SRCX/$PKG/tests" "$TESTS"
[ -d "$TESTS/testthat" ] || { echo "tests/testthat is missing from the tarball" >&2; exit 2; }

echo "== [4/4] testthat::test_dir(load_package = \"installed\")"
cd "$TESTS/testthat"     # relative paths now resolve inside the temporary copy
SP_LIB="$LIB" SP_FILTER="${FILTER:-}" Rscript --no-save --no-restore -e '
lib <- Sys.getenv("SP_LIB")
.libPaths(c(lib, .libPaths()))
loc <- normalizePath(find.package("simplePHENOTYPES"))
if (!startsWith(loc, normalizePath(lib))) {
  stop("the package under test is not the temporary installation: ", loc)
}
cat("package under test:", loc, "\n")
cat("working directory :", getwd(), "\n")
filt <- Sys.getenv("SP_FILTER")
res <- testthat::test_dir(".", package = "simplePHENOTYPES",
                          load_package = "installed",
                          filter = if (nzchar(filt)) filt else NULL,
                          reporter = "summary", stop_on_failure = FALSE,
                          stop_on_warning = FALSE)
df <- as.data.frame(res)
cat(sprintf("\nfiles: %d  tests: %d  expectations: %d  failed: %d  errors: %d  skipped: %d  warnings: %d\n",
            length(unique(df$file)), nrow(df), sum(df$nb), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
bad <- df[df$failed > 0 | df$error, c("file", "test")]
if (nrow(bad) > 0L) {
  cat("FAILED in installed mode:\n")
  cat(sprintf("  %s :: %s\n", bad$file, bad$test), sep = "")
  quit(status = 1L)
}
cat("installed-package tests: OK\n")
'
