# Round-7 crossing fixes: H1a (as_population validates the counted_allele
# attribute; orientation guard cannot be bypassed by malformed tokens), H1b
# (all-NA non-character counted rejected), G5 (sanitized, bounded
# sanitized default file names).

mk_num <- function(allele, counted = NULL) {
  g <- data.frame(snp = c("m1", "m2"), allele = allele, chr = 1L,
                  pos = 1:2, cm = c(0, 50), x = c(-1, 1), y = c(0, 1),
                  stringsAsFactors = FALSE)
  if (!is.null(counted)) attr(g, "counted_allele") <- counted
  g
}

test_that("H1a: invalid counted_allele attribute is rejected by as_population", {
  expect_error(as_population(mk_num(c("A/G", "C/T"), c("Z", "Z"))),
               "counted_allele")
  expect_error(as_population(mk_num(c("A/G", "C/T"), c("", ""))),
               "counted_allele")
  expect_error(as_population(mk_num(c("A/G", "C/T"), c("A", "T", "G"))),
               "one entry per marker")
  # valid values are accepted and upper-cased
  p <- as_population(mk_num(c("A/G", "C/T"), c("a", "T")))
  expect_equal(p$map$counted, c("A", "T"))
  # no attribute: label-only behaviour unchanged
  expect_null(as_population(mk_num(c("A/G", "C/T")))$map$counted)
})

test_that("H1a: reversed labels with the invalid attribute no longer cross silently", {
  expect_error(
    cross(as_population(mk_num(c("A/G", "C/T"), c("Z", "Z"))),
          as_population(mk_num(c("G/A", "T/C"), c("Z", "Z")))),
    "counted_allele")
})

test_that("H1a: .check_orientation ignores tokens that are not alleles of the marker", {
  a <- list(snp = c("m1", "m2"), counted = c("Z", "Z"), allele = c("A/G", "C/T"))
  b <- list(snp = c("m1", "m2"), counted = c("Z", "Z"), allele = c("G/A", "T/C"))
  expect_warning(.check_orientation(a, b), "different order")
  # numeric counted: treated as unknown
  a$counted <- c(1, 1); b$counted <- c(1, 1)
  expect_warning(.check_orientation(a, b), "different order")
  # valid, agreeing tokens still cover the marker (no warning despite reversed labels)
  a$counted <- c("A", "C"); b$counted <- c("A", "C")
  expect_silent(.check_orientation(a, b))
  # valid, disagreeing tokens still error
  b$counted <- c("G", "T")
  expect_error(.check_orientation(a, b), "different alleles as \\+1")
  # a token not in the label of one side leaves that marker uncovered
  b$counted <- c("A", "Q")
  expect_warning(.check_orientation(a, b), "different order")
})

test_that("H1b: all-NA non-character counted columns are rejected", {
  expect_error(.check_counted(c(NA_real_, NA_real_), c("A/G", "C/T"), 2),
               "must be a character vector")
  expect_error(.check_counted(c(NA, NA), c("A/G", "C/T"), 2),
               "must be a character vector")
  expect_error(as_population(mk_num(c("A/G", "C/T"), c(NA_real_, NA_real_))),
               "counted_allele")
  expect_null(.check_counted(c(NA_character_, NA_character_), NULL, 2))
  expect_null(.check_counted(NULL))
})

test_that("G5: sanitized default file names are distinct for these labels and bounded", {
  test_dir <- normalizePath(file.path(testthat::test_path(), ".."))
  td <- withr::local_tempdir()
  withr::local_dir(td)
  hmp <- data.table::fread(file.path(test_dir, "test.hmp.txt"),
                           data.table = FALSE)
  env <- new.env()
  # run as_numeric() on the object bound to symbol `sym` in `env`
  go <- function(sym) {
    assign(sym, hmp, envir = env)
    before <- list.files(td)
    suppressMessages(suppressWarnings(
      eval(bquote(as_numeric(.(as.name(sym)), to_r = FALSE, verbose = FALSE)),
           envir = env)))
    setdiff(list.files(td), before)
  }
  f1 <- go("a+b")
  f2 <- go("a_b")
  expect_length(f1, 1L)
  expect_length(f2, 1L)
  expect_false(identical(f1, f2))
  expect_identical(f2, "a_b_numeric.txt")        # unchanged label: no hash
  expect_match(f1, "^a_b_[0-9a-f]{8}_numeric\\.txt$")
  expect_identical(.label_hash("a+b"), .label_hash("a+b"))   # deterministic
  expect_false(.label_hash("a+b") == .label_hash("a_b"))
  # a legal 300-character symbol works, with a bounded file name
  long <- paste(rep("x", 300), collapse = "")
  f3 <- go(long)
  expect_length(f3, 1L)
  expect_lte(nchar(f3), 100 + 1 + 8 + nchar("_numeric.txt"))
  expect_true(file.exists(f3))
  # two long labels that share their first 100 characters differ by the hash
  f4 <- go(paste0(long, "y"))
  expect_length(f4, 1L)
  expect_false(identical(f3, f4))
  # no windows-reserved characters in any generated name
  expect_false(any(grepl('[<>:"/\\\\|?*]', list.files(td))))
  expect_false(any(grepl('[<>:"|?*]', c(f1, f3, f4))))
})

test_that("G5: ordinary symbols keep exactly their old default name", {
  test_dir <- normalizePath(file.path(testthat::test_path(), ".."))
  td <- withr::local_tempdir()
  withr::local_dir(td)
  hmp <- data.table::fread(file.path(test_dir, "test.hmp.txt"),
                           data.table = FALSE)
  geno <- hmp
  SNP55K_maize282_maf04 <- hmp
  suppressMessages(suppressWarnings({
    as_numeric(hmp, to_r = FALSE, verbose = FALSE)
    as_numeric(geno, to_r = FALSE, verbose = FALSE)
    as_numeric(SNP55K_maize282_maf04, to_r = FALSE, verbose = FALSE)
  }))
  expect_setequal(list.files(td),
                  c("hmp_numeric.txt", "geno_numeric.txt",
                    "SNP55K_maize282_maf04_numeric.txt"))
})
