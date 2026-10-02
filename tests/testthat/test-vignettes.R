# Pre-PR gate: every vignette's R code runs. Each vignettes/*.Rmd is tangled
# with knitr::purl() into a temporary file and evaluated expression by
# expression in a fresh environment (parent: the global environment, where the
# vignette's own library(simplePHENOTYPES) puts the package, as for a reader).
# Visible values are printed, as knitting would, so print methods are exercised
# too. An error fails the test, naming the vignette, the chunk and the
# expression, except in chunks declared `error = TRUE` (their errors are part of
# the narrative). Chunks with `eval = FALSE` are commented out by purl() and
# are therefore not run.
#
# Nothing is written outside tempdir(): the working directory is moved into a
# temporary folder for the run, graphics go to a null device, and options and
# the RNG state are restored afterwards.

run_vignette_code <- function(rmd, max_seconds = 120) {
  vname <- basename(rmd)
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  utils::capture.output(
    suppressMessages(knitr::purl(rmd, output = script, quiet = TRUE,
                                 documentation = 1)),
    file = nullfile())
  lines <- readLines(script, warn = FALSE)
  exprs <- parse(text = lines, keep.source = TRUE)
  refs <- attr(exprs, "srcref")
  # chunk header lines look like "## ----label, option = value------"
  hdr <- grep("^## ----", lines)
  chunk_of <- function(line) {
    h <- hdr[hdr <= line]
    if (length(h) == 0L) return(list(label = "(top)", error_ok = FALSE))
    txt <- sub("-{3,}$", "", sub("^## ----", "", lines[max(h)]))
    list(label = txt, error_ok = grepl("error\\s*=\\s*TRUE", txt))
  }

  workdir <- tempfile("vignette-run-")
  dir.create(workdir)
  old_wd <- setwd(workdir)
  old_opts <- options()
  old_seed <- .Random.seed_safe()
  dev_before <- grDevices::dev.list()
  grDevices::pdf(NULL)
  on.exit({
    # close every device opened by the run, then restore the session
    while (!is.null(d <- grDevices::dev.list()) &&
           length(setdiff(d, dev_before)) > 0L) {
      grDevices::dev.off()
    }
    .restore_seed(old_seed)
    options(old_opts)
    setwd(old_wd)
    unlink(workdir, recursive = TRUE)
  }, add = TRUE)

  env <- new.env(parent = globalenv())
  failure <- NULL
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_along(exprs)) {
    ch <- chunk_of(refs[[i]][1L])
    res <- tryCatch({
      utils::capture.output({
        v <- withCallingHandlers(
          withVisible(eval(exprs[[i]], env)),
          message = function(m) invokeRestart("muffleMessage"),
          warning = function(w) invokeRestart("muffleWarning"))
        if (v$visible) print(v$value)
      }, file = nullfile())
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(res) && !ch$error_ok) {
      failure <- sprintf("vignette %s, chunk '%s', line %d of the tangled code: %s\n  > %s",
                         vname, ch$label, refs[[i]][1L], res,
                         paste(deparse(exprs[[i]], width.cutoff = 80L)[1L]))
      break
    }
  }
  list(vignette = vname, failure = failure,
       seconds = proc.time()[["elapsed"]] - t0, n_expr = length(exprs))
}

test_that("every vignette's R code runs without error", {
  skip_on_cran()
  skip_if_no_source("vignettes")
  skip_if_not_installed("knitr")
  rmds <- sort(list.files(testthat::test_path("..", "..", "vignettes"),
                          pattern = "\\.Rmd$", full.names = TRUE))
  expect_gt(length(rmds), 0L)
  for (rmd in rmds) {
    r <- run_vignette_code(rmd)
    expect_gt(r$n_expr, 0L, label = paste(r$vignette, "has code"))
    if (!is.null(r$failure)) {
      fail(r$failure)
    } else {
      succeed()
    }
    expect_lt(r$seconds, 120, label = paste(r$vignette, "run time (s)"))
  }
})

test_that("running the vignettes leaves no files in the working directory", {
  skip_on_cran()
  skip_if_no_source("vignettes")
  skip_if_not_installed("knitr")
  before <- list.files(getwd(), all.files = TRUE, recursive = TRUE)
  r <- run_vignette_code(testthat::test_path("..", "..", "vignettes",
                                             "genetic-maps.Rmd"))
  expect_null(r$failure)
  expect_identical(list.files(getwd(), all.files = TRUE, recursive = TRUE),
                   before)
})
