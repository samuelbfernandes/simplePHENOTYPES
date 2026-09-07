# setup.R — testthat 3e setup/teardown for the package test suite.
#
# create_phenotypes() (frozen legacy) writes intermediate GDS files into the
# working directory during format conversion; under devtools::test() / R CMD
# check that directory is tests/testthat/. Remove any such artifacts after the
# whole suite runs so the source tree and CRAN bundle stay clean (BUGS.md).
withr::defer(
  unlink(list.files(getwd(), pattern = "^file.*\\.gds$", full.names = TRUE)),
  teardown_env()
)
