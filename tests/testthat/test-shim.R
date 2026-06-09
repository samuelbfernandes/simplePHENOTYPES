# test-shim.R
#
# Smoke tests for the create_phenotypes() v1 compatibility shim.
# These verify the public API has not broken and covers the main
# model / architecture combinations used in the README.
# No phenotype values are asserted — only structure and absence of errors.

test_that("shim: single trait, model A, returns data frame", {
  data("SNP55K_maize282_maf04")
  result <- suppressMessages(
    create_phenotypes(
      geno_obj    = SNP55K_maize282_maf04,
      add_QTN_num = 3,
      add_effect  = 0.2,
      rep         = 2,
      h2          = 0.7,
      model       = "A",
      seed        = 1,
      home_dir    = tempdir()
    )
  )
  expect_null(result)   # to_r defaults to FALSE → returns invisibly NULL
})

test_that("shim: single trait to_r=TRUE returns data frame with correct rows", {
  data("SNP55K_maize282_maf04")
  n_samples <- ncol(SNP55K_maize282_maf04) - 5L   # numeric format: 5 meta cols
  result <- suppressMessages(
    create_phenotypes(
      geno_obj    = SNP55K_maize282_maf04,
      add_QTN_num = 3,
      add_effect  = 0.2,
      rep         = 3,
      h2          = 0.7,
      model       = "A",
      to_r        = TRUE,
      seed        = 1,
      home_dir    = tempdir()
    )
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), n_samples)
  expect_true(ncol(result) >= 4L)  # Trait + at least 3 rep columns
})

test_that("shim: additive + dominance model (AD)", {
  data("SNP55K_maize282_maf04")
  expect_no_error(
    suppressMessages(
      create_phenotypes(
        geno_obj    = SNP55K_maize282_maf04,
        add_QTN_num = 2,
        dom_QTN_num = 2,
        add_effect  = 0.3,
        dom_effect  = 0.2,
        rep         = 1,
        h2          = 0.5,
        model       = "AD",
        seed        = 2,
        home_dir    = tempdir()
      )
    )
  )
})

test_that("shim: pleiotropic architecture, 3 traits", {
  data("SNP55K_maize282_maf04")
  expect_no_error(
    suppressMessages(
      create_phenotypes(
        geno_obj     = SNP55K_maize282_maf04,
        add_QTN_num  = 3,
        add_effect   = c(0.04, 0.2, 0.1),
        ntraits      = 3,
        rep          = 1,
        h2           = c(0.2, 0.4, 0.4),
        architecture = "pleiotropic",
        model        = "A",
        vary_QTN     = FALSE,
        seed         = 10,
        home_dir     = tempdir()
      )
    )
  )
})

test_that("shim: HapMap file input (geno_file)", {
  hmp_file <- normalizePath(file.path(testthat::test_path(), "..", "test.hmp.txt"),
                             mustWork = FALSE)
  skip_if_not(file.exists(hmp_file))
  expect_no_error(
    suppressMessages(
      create_phenotypes(
        geno_file   = hmp_file,
        add_QTN_num = 3,
        add_effect  = 0.2,
        rep         = 1,
        h2          = 0.7,
        model       = "A",
        seed        = 1,
        home_dir    = tempdir()
      )
    )
  )
})

test_that("shim: seed reproducibility — identical results with same seed", {
  data("SNP55K_maize282_maf04")
  run <- function() suppressMessages(
    create_phenotypes(
      geno_obj    = SNP55K_maize282_maf04,
      add_QTN_num = 3,
      add_effect  = 0.2,
      rep         = 1,
      h2          = 0.7,
      model       = "A",
      to_r        = TRUE,
      seed        = 42,
      home_dir    = tempdir()
    )
  )
  expect_identical(run(), run())
})
