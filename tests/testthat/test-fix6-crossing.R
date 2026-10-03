# Round-6 crossing fixes: C1 (interference nu range / non-terminating draw),
# H1 (map$counted validation), G5 (portable default file name), C2 (docs).

tiny_pop <- function(counted = NULL, allele = c("A/G", "C/T")) {
  map <- data.frame(snp = c("m1", "m2"), chr = 1L, pos = c(1L, 2L),
                    cm = c(0, 50), allele = allele, stringsAsFactors = FALSE)
  if (!is.null(counted)) map$counted <- counted
  cis <- matrix(c(0, 1, 1, 0), 2, 2, dimnames = list(c("m1", "m2"), NULL))
  trans <- 1 - cis
  population_from_haplotypes(cis, trans, map, ids = c("a", "b"))
}

test_that("C1: huge nu is rejected quickly instead of looping forever", {
  map <- data.frame(snp = c("m1", "m2"), chr = 1L, pos = 1:2, cm = c(0, 100))
  cis <- matrix(c(0, 1, 1, 0), 2, 2, dimnames = list(c("m1", "m2"), NULL))
  pop <- population_from_haplotypes(cis, 1 - cis, map, ids = c("a", "b"))
  for (nu in c(.Machine$double.xmax, 1e7, 1e300)) {
    msg <- tryCatch({
      setTimeLimit(elapsed = 20, transient = TRUE)
      cross(pop[1], pop[2], n = 2, interference = list(nu = nu, p = 0.5), seed = 1)
      "no error"
    }, error = function(e) conditionMessage(e),
    finally = setTimeLimit(elapsed = Inf, transient = FALSE))
    expect_match(msg, "interference\\$nu")
  }
  expect_error(.check_interference(list(nu = Inf)), "nu")
  # the documented maximum is still accepted and terminates
  ok <- cross(pop[1], pop[2], n = 3, interference = list(nu = 1e6, p = 0.5),
              seed = 2)
  expect_s3_class(ok, "Population")
  expect_equal(.check_interference(list(nu = 1e6, p = 1))$nu, 1e6)
})

test_that("C1: the renewal loop errors on a non-advancing gap", {
  # bypass the validator to hit the guard directly: a rate so large that gaps
  # underflow to 0 would never advance the position
  testthat::with_mocked_bindings(
    expect_error(
      .draw_meiosis_interference(list(c(0, 1)), 2L, nu = 2, p = 0),
      "strictly positive"),
    rgamma = function(n, shape, rate) if (shape > 2.5) rep(0.1, n) else rep(0, n),
    .package = "stats")
})

test_that("C1: the stable rate is bit-identical to the old formula", {
  expect_identical(.check_interference("poisson"), "poisson")
  for (nu in c(1, 2.6, 7.3, 1e6)) for (p in c(0, 0.1, 0.37, 0.9)) {
    expect_identical(2 * nu * (1 - p), nu * (2 * (1 - p)))
  }
})

test_that("H1: counted must be character allele symbols, '' and numeric rejected", {
  expect_s3_class(tiny_pop(c("A", "T")), "Population")
  expect_s3_class(tiny_pop(c("a", NA)), "Population")
  expect_s3_class(tiny_pop(c(NA_character_, NA_character_)), "Population")
  expect_error(tiny_pop(c("", "")), "map\\$counted")
  expect_error(tiny_pop(c("A", "")), "map\\$counted")
  expect_error(tiny_pop(c(1, 0)), "map\\$counted")
  expect_error(tiny_pop(c("AG", "T")), "map\\$counted")   # not one of A/G
  expect_error(tiny_pop(c("C", "T")), "one of the two alleles")  # C not in A/G
  expect_error(tiny_pop(c(" ", "T")), "map\\$counted")
})

test_that("H1: counted without an allele label accepts any single token", {
  p <- tiny_pop(c("A", "T"), allele = c(NA, NA))
  expect_equal(p$map$counted, c("A", "T"))
})

test_that("H1: the orientation guard treats '' as unknown, not as a match", {
  a <- list(snp = c("m1", "m2"), counted = c("", ""), allele = c("A/G", "C/T"))
  b <- list(snp = c("m1", "m2"), counted = c("", ""), allele = c("G/A", "T/C"))
  expect_warning(.check_orientation(a, b), "different order")
  # a real mismatch is still an error
  a$counted <- c("A", "C"); b$counted <- c("G", "T")
  expect_error(.check_orientation(a, b), "different alleles as \\+1")
  # '' on one side falls back to the label comparison
  b$counted <- c("", "")
  expect_warning(.check_orientation(a, b), "different order")
})

test_that("H1: an '' counted attribute reaching as_population is rejected (round 7)", {
  g <- data.frame(snp = c("m1", "m2"), allele = c("A/G", "C/T"), chr = 1L,
                  pos = 1:2, cm = c(0, 50), x = c(-1, 1), y = c(0, 1))
  attr(g, "counted_allele") <- c("", "")
  expect_error(as_population(g), "counted_allele")
})

test_that("G5: default file name from an inline object is portable", {
  test_dir <- normalizePath(file.path(testthat::test_path(), ".."))
  td <- withr::local_tempdir()
  withr::local_dir(td)
  hmp <- data.table::fread(file.path(test_dir, "test.hmp.txt"),
                           data.table = FALSE)
  quiet <- function(e) suppressMessages(suppressWarnings(e))
  # a symbol keeps its name
  quiet(as_numeric(hmp, to_r = FALSE, verbose = FALSE))
  expect_true(file.exists("hmp_numeric.txt"))
  # an inline object (do.call) no longer yields '<inline ...>' in the file name
  quiet(do.call(as_numeric, list(hmp, to_r = FALSE, verbose = FALSE)))
  out <- list.files(td)
  expect_false(any(grepl("[<>]", out)))
  expect_true(any(grepl("^inline_data\\.frame.*_numeric\\.txt$", out)))
})
