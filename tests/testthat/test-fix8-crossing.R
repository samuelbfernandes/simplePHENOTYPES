# Round-7 fixes (Codex closeout section 3): default output names can coincide,
# so an existing default-named file is overwritten WITH a warning; explicit
# names are unchanged; .check_counted() rejects array-shaped input.

.fx8_path <- normalizePath(file.path(testthat::test_path(), "..",
                                     "test.hmp.txt"))
# no "default output file ..." warning (other conversion warnings are fine)
expect_no_default_warning <- function(expr) {
  w <- character()
  withCallingHandlers(suppressMessages(expr),
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
    })
  testthat::expect_false(any(grepl("default output file", w)))
}
expect_default_warning <- function(expr, pattern = "already exists and is overwritten") {
  w <- character()
  withCallingHandlers(suppressMessages(expr),
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
    })
  testthat::expect_true(any(grepl(pattern, w)))
}
.fx8_hmp <- function() data.table::fread(.fx8_path, data.table = FALSE)

test_that("hash-colliding labels: second default write warns and overwrites", {
  td <- withr::local_tempdir(); withr::local_dir(td)
  hmp <- .fx8_hmp()
  env <- new.env()
  l1 <- "a@;\\!~b"; l2 <- "a#%]@!b"
  assign(l1, hmp, envir = env)
  assign(l2, hmp[, 1:(ncol(hmp) - 1)], envir = env)
  run <- function(l) eval(bquote(as_numeric(.(as.name(l)), to_r = FALSE,
                                            verbose = FALSE)), envir = env)
  expect_no_default_warning(run(l1))
  f <- list.files(td)
  expect_length(f, 1L)
  md5a <- unname(tools::md5sum(f))
  expect_default_warning(run(l2))
  expect_identical(list.files(td), f)
  expect_false(identical(unname(tools::md5sum(f)), md5a))
})

test_that("case-variant labels and same-stem files warn on collision", {
  td <- withr::local_tempdir(); withr::local_dir(td)
  hmp <- .fx8_hmp()
  CaseLabel <- hmp; caselabel <- hmp
  expect_no_default_warning(as_numeric(CaseLabel, to_r = FALSE, verbose = FALSE))
  # on a case-insensitive file system the target already exists; on a
  # case-sensitive one it does not (and nothing collides)
  if (file.exists("caselabel_numeric.txt")) {
    expect_default_warning(as_numeric(caselabel, to_r = FALSE, verbose = FALSE))
  } else {
    expect_no_default_warning(as_numeric(caselabel, to_r = FALSE, verbose = FALSE))
  }
  # same stem from different input files
  data.table::fwrite(hmp, "same.hmp.txt", sep = "\t", quote = FALSE)
  data.table::fwrite(hmp, "same.txt", sep = "\t", quote = FALSE)
  expect_no_default_warning(as_numeric("same.hmp.txt", verbose = FALSE))
  expect_default_warning(as_numeric("same.txt", verbose = FALSE, from = "hapmap"),
                         "default output file same_numeric.txt already exists")
})

test_that("same-shape inline objects warn; explicit file names do not", {
  td <- withr::local_tempdir(); withr::local_dir(td)
  hmp <- .fx8_hmp(); hmp2 <- hmp
  hmp2[[ncol(hmp2)]] <- rev(hmp2[[ncol(hmp2)]])
  expect_no_default_warning(
    do.call(as_numeric, list(hmp, to_r = FALSE, verbose = FALSE)))
  expect_default_warning(
    do.call(as_numeric, list(hmp2, to_r = FALSE, verbose = FALSE)))
  # explicit names: no warning, overwritten exactly as before
  expect_no_default_warning(
    as_numeric(hmp, to_r = FALSE, verbose = FALSE, file_name = "out.txt"))
  expect_no_default_warning(
    as_numeric(hmp2, to_r = FALSE, verbose = FALSE, file_name = "out.txt"))
  # re-running the same default call still writes (with the warning)
  expect_no_default_warning(as_numeric(hmp, to_r = FALSE, verbose = FALSE))
  expect_default_warning(as_numeric(hmp, to_r = FALSE, verbose = FALSE))
  expect_true(file.exists("hmp_numeric.txt"))
})

test_that(".check_counted rejects array-shaped input", {
  m <- matrix(c("A", "G"), 2, 1)
  expect_error(.check_counted(m, NULL, 2, "map$counted"), "map\\$counted.*dim")
  expect_error(.check_counted(array(c("A", "G"), c(2, 1, 1)), NULL, 2), "dim")
  expect_error(.check_counted(matrix(1:2, 2, 1), NULL, 2), "dim")
  expect_identical(.check_counted(c("a", "g"), NULL, 2), c("A", "G"))
})

test_that("as_population / population_from_haplotypes reject matrix counted", {
  g <- data.frame(snp = c("m1", "m2"), allele = c("A/G", "A/G"), chr = 1L,
                  pos = 1:2, cm = c(0, 50), x = c(-1, 1), y = c(0, 1),
                  stringsAsFactors = FALSE)
  attr(g, "counted_allele") <- c("A", "G")
  expect_s3_class(as_population(g), "Population")
  attr(g, "counted_allele") <- matrix(c("A", "G"), 2, 1)
  expect_error(as_population(g), "counted_allele.*dim")
  map <- data.frame(snp = c("m1", "m2"), chr = 1L, pos = 1:2, cm = c(0, 50),
                    allele = c("A/G", "A/G"), stringsAsFactors = FALSE)
  cis <- matrix(c(0, 1, 1, 0), 2, 2); trans <- cis
  colnames(cis) <- colnames(trans) <- c("i1", "i2")
  map$counted <- c("A", "G")
  expect_s3_class(population_from_haplotypes(cis, trans, map), "Population")
  map$counted <- I(matrix(c("A", "G"), 2, 1))
  expect_error(population_from_haplotypes(cis, trans, map), "map\\$counted.*dim")
})
