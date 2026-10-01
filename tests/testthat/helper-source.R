# Helpers for tests that read package sources or use optional packages.
# An installed-package check (R CMD check) has no R/, docs/ or benchmarks/ next
# to tests/, so such tests skip there instead of failing.
skip_if_no_source <- function(...) {
  p <- testthat::test_path("..", "..", ...)
  testthat::skip_if_not(all(file.exists(p)),
                        "package sources are not available (installed-package check)")
}

# an optional package's function, looked up at run time so that R CMD check does
# not require the package to be declared
optional_fun <- function(pkg, fun) {
  testthat::skip_if_not_installed(pkg)
  getExportedValue(pkg, fun)
}
