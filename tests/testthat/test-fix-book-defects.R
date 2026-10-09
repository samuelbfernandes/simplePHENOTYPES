# test-fix-book-defects.R
#
# Regressions found while authoring docs/book: tibble inputs must not trip
# "Unknown or uninitialised column" warnings on optional columns, and a
# relative `home_dir` in create_phenotypes() must not be resolved twice.

test_that("mate() takes a tibble plan without a design column silently", {
  skip_if_not_installed("tibble")
  set.seed(1)
  g <- matrix(sample(c(-1L, 1L), 6 * 30, replace = TRUE), 30, 6)
  colnames(g) <- paste0("P", 1:6)
  geno <- cbind(data.frame(snp = paste0("m", 1:30), allele = "A/G",
                           chr = rep(1:3, each = 10), pos = rep(1:10, 3),
                           cm = rep(seq(0, 90, length.out = 10), 3),
                           stringsAsFactors = FALSE), as.data.frame(g))
  pop <- as_population(geno)
  plan <- tibble::tibble(mother = c("P1", "P3", "P4"),
                         father = c("P2", "P3", "P5"), n = 2)
  expect_no_warning(prog <- mate(plan, pop, seed = 9))
  df <- mate(as.data.frame(plan), pop, seed = 9)
  expect_identical(dosages(prog), dosages(df))
  expect_identical(parentage(prog), parentage(df))
  expect_equal(attr(prog, "plan")$design, c("cross", "self", "cross"))
})

test_that("population_from_haplotypes() takes a tibble map without counted/allele", {
  skip_if_not_installed("tibble")
  map <- tibble::tibble(snp = paste0("m", 1:4), chr = 1L,
                        pos = (1:4) * 1000L, cm = (0:3) * 5)
  cis <- cbind(a = c(1, 1, 0, 0), b = c(0, 1, 0, 1))
  trans <- cbind(a = c(0, 1, 1, 0), b = c(1, 1, 0, 0))
  expect_no_warning(pop <- population_from_haplotypes(cis, trans, map))
  expect_identical(pop, population_from_haplotypes(cis, trans,
                                                   as.data.frame(map)))
})

test_that("create_phenotypes() with a relative home_dir writes home_dir/output_dir", {
  e <- new.env()
  utils::data("SNP55K_maize282_maf04", envir = e)
  wd <- tempfile("relhome_")
  dir.create(file.path(wd, "results"), recursive = TRUE)
  withr::defer(unlink(wd, recursive = TRUE))
  withr::local_dir(wd)
  suppressMessages(suppressWarnings(
    create_phenotypes(geno_obj = e$SNP55K_maize282_maf04, home_dir = "results",
                      output_dir = "run1", add_QTN_num = 3, add_effect = 0.2,
                      h2 = 0.5, model = "A", rep = 1, seed = 1, to_r = TRUE,
                      verbose = FALSE)
  ))
  expect_true(dir.exists(file.path(wd, "results", "run1")))
  expect_true(length(dir(file.path(wd, "results", "run1"))) > 0)
  expect_false(dir.exists(file.path(wd, "results", "results")))
  expect_identical(normalizePath(getwd()), normalizePath(wd))
})
