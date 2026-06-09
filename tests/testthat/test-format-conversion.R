# test-format-conversion.R
#
# Deterministic tests for the format-detection and numericalization pipeline.
# All assertions use hand-verifiable values; no random draws, no seeds needed.

test_dir <- normalizePath(file.path(testthat::test_path(), ".."))

# ---------------------------------------------------------------------------
# Helper: tiny hand-crafted 3-SNP × 4-sample HapMap-style character matrix.
#
# SNP1: AA(x2), AG(x1), GG(x1)  → AA major → coded -101: [ 1,  1,  0, -1]
# SNP2: CC(x3), TT(x1)           → CC major → coded -101: [ 1,  1,  1, -1]
# SNP3: N(x4)                    → all missing → coded:  [NA, NA, NA, NA]
# ---------------------------------------------------------------------------
.tiny_mat <- function() {
  matrix(c("AA","AA","AG","GG",
            "CC","CC","CC","TT",
            "N", "N", "N", "N"),
         nrow = 3L, byrow = TRUE)
}

# ---------------------------------------------------------------------------
# 1. detect_format() — file extensions
# ---------------------------------------------------------------------------

test_that("detect_format() recognises file extensions", {
  d <- simplePHENOTYPES:::detect_format

  expect_equal(d(file.path(test_dir, "test.hmp.txt")), "hapmap")
  expect_equal(d(file.path(test_dir, "test.vcf")),     "vcf")
  expect_equal(d(file.path(test_dir, "test.bed")),     "bed")
  expect_equal(d(file.path(test_dir, "test.ped")),     "ped")
  expect_equal(d(file.path(test_dir, "test.gds")),     "gds")
})

test_that("detect_format() recognises in-memory objects", {
  d <- simplePHENOTYPES:::detect_format

  data("SNP55K_maize282_maf04")
  expect_equal(d(SNP55K_maize282_maf04), "numeric")

  hmp_df <- data.table::fread(file.path(test_dir, "test.hmp.txt"),
                               data.table = FALSE)
  expect_equal(d(hmp_df), "hapmap")
})

# ---------------------------------------------------------------------------
# 2. parse_hapmap_chars_to_raw() — hand-verifiable raw matrix
# ---------------------------------------------------------------------------

test_that("parse_hapmap_chars_to_raw() produces correct 0/1/2/NA matrix", {
  raw <- simplePHENOTYPES:::parse_hapmap_chars_to_raw(.tiny_mat())

  # SNP1: AA→0, AA→0, AG→1, GG→2
  expect_equal(raw[1L, ], c(0L, 0L, 1L, 2L))
  # SNP2: CC→0, CC→0, CC→0, TT→2
  expect_equal(raw[2L, ], c(0L, 0L, 0L, 2L))
  # SNP3: all missing
  expect_true(all(is.na(raw[3L, ])))
})

test_that("parse_hapmap_chars_to_raw() sets non-biallelic SNP to NA", {
  tri <- matrix(c("AA","TT","CC","AA"), nrow = 1L)  # 3 distinct homozygotes
  expect_message(
    raw <- simplePHENOTYPES:::parse_hapmap_chars_to_raw(tri),
    regexp = "Non-biallelic"
  )
  expect_true(all(is.na(raw[1L, ])))
})

# ---------------------------------------------------------------------------
# 3. compute_flip() — frequency and reference methods
# ---------------------------------------------------------------------------

test_that("compute_flip() is FALSE when allele1 is major", {
  # count0 = 3 (allele1 homs), count2 = 1 (allele2 homs) → flip = FALSE
  raw <- matrix(c(0L, 0L, 0L, 2L), nrow = 1L)
  expect_false(simplePHENOTYPES:::compute_flip(raw)[[1L]])
})

test_that("compute_flip() is TRUE when allele2 is major", {
  # count0 = 1, count2 = 3 → flip = TRUE
  raw <- matrix(c(0L, 2L, 2L, 2L), nrow = 1L)
  expect_true(simplePHENOTYPES:::compute_flip(raw)[[1L]])
})

test_that("compute_flip() method=reference flips when allele1 differs from ref", {
  raw     <- matrix(c(0L, 0L, 2L, 2L), nrow = 2L)
  allele1 <- c("A", "G")
  ref     <- c("A", "A")   # SNP1: allele1==ref → no flip; SNP2: allele1≠ref → flip
  flip    <- simplePHENOTYPES:::compute_flip(raw, method = "reference",
                                              allele1 = allele1, ref = ref)
  expect_equal(flip, c(FALSE, TRUE))
})

test_that("compute_flip() method=reference errors without ref_allele", {
  raw <- matrix(c(0L, 0L), nrow = 1L)
  expect_error(
    simplePHENOTYPES:::compute_flip(raw, method = "reference"),
    regexp = "ref_allele"
  )
})

# ---------------------------------------------------------------------------
# 4. numericalize_core() — Rust kernel, hand-verifiable values
# ---------------------------------------------------------------------------

test_that("numericalize_core() -101 no-flip: major=1, het=0, minor=-1", {
  raw  <- c(0L, 0L, 1L, 2L)  # 1 SNP × 4 samples: hom1, hom1, het, hom2
  flip <- FALSE
  out  <- simplePHENOTYPES:::numericalize_core(raw, 1L, 4L, flip,
                                                "-101", "Add", "None")
  expect_equal(out, c(1L, 1L, 0L, -1L))
})

test_that("numericalize_core() -101 flip=TRUE: hom2=major=1", {
  raw  <- c(0L, 2L, 2L, 2L)
  flip <- TRUE
  out  <- simplePHENOTYPES:::numericalize_core(raw, 1L, 4L, flip,
                                                "-101", "Add", "None")
  expect_equal(out, c(-1L, 1L, 1L, 1L))
})

test_that("numericalize_core() 012 coding", {
  raw  <- c(0L, 1L, 2L, NA_integer_)
  flip <- FALSE
  out  <- simplePHENOTYPES:::numericalize_core(raw, 1L, 4L, flip,
                                                "012", "Add", "Middle")
  expect_equal(out, c(2L, 1L, 0L, 1L))  # NA imputed to het = 1
})

test_that("numericalize_core() handles mixed flip correctly (layout regression)", {
  # SNP1: [0,2,2,2] flip=TRUE  → expected: [-1, 1, 1, 1]
  # SNP2: [0,0,0,2] flip=FALSE → expected: [ 1, 1, 1,-1]
  # If the R→Rust layout were wrong (as.integer column-major instead of
  # as.vector(t()) row-major), the results would be scrambled.
  raw_mat <- matrix(c(0L, 2L, 2L, 2L,
                       0L, 0L, 0L, 2L),
                    nrow = 2L, byrow = TRUE)
  flip <- c(TRUE, FALSE)

  coded_vec <- simplePHENOTYPES:::numericalize_core(
    as.vector(t(raw_mat)),   # correct row-major
    n_snp  = 2L,
    n_samp = 4L,
    flip   = flip,
    code_as = "-101",
    model   = "Add",
    impute  = "None"
  )
  coded_mat <- matrix(coded_vec, nrow = 2L, ncol = 4L, byrow = TRUE)

  expect_equal(coded_mat[1L, ], c(-1L,  1L,  1L,  1L))  # SNP1
  expect_equal(coded_mat[2L, ], c( 1L,  1L,  1L, -1L))  # SNP2
})

# ---------------------------------------------------------------------------
# 5. End-to-end: HapMap tiny matrix → as_numeric output schema
# ---------------------------------------------------------------------------

test_that("as_numeric() tiny HapMap: correct -101 values", {
  hmp_names <- c("rs#","alleles","chrom","pos","strand","assembly#",
                 "center","protLSID","assayLSID","panelLSID","QCcode")
  meta <- as.data.frame(matrix(
    c("SNP1","A/G","1","100",rep("?",7),
      "SNP2","C/T","1","200",rep("?",7)),
    nrow = 2L, byrow = TRUE, dimnames = list(NULL, hmp_names)
  ), stringsAsFactors = FALSE)
  geno <- data.frame(s1 = c("AA","CC"), s2 = c("AA","CC"),
                     s3 = c("AG","CC"), s4 = c("GG","TT"),
                     stringsAsFactors = FALSE)
  df <- cbind(meta, geno)

  res <- suppressMessages(as_numeric(df, to_r = TRUE, verbose = FALSE))

  expect_equal(nrow(res), 2L)
  expect_equal(ncol(res), 9L)   # 5 meta + 4 samples
  expect_equal(as.integer(res[1L, 6:9]), c( 1L,  1L,  0L, -1L))  # SNP1
  expect_equal(as.integer(res[2L, 6:9]), c( 1L,  1L,  1L, -1L))  # SNP2
})

# ---------------------------------------------------------------------------
# 6. Output schema: always 5 metadata columns for every input format
# ---------------------------------------------------------------------------

test_that("all input formats produce 5-column metadata schema", {
  .meta_cols <- function(res) colnames(res)[1:5]
  expected   <- c("snp", "allele", "chr", "pos", "cm")

  hmp <- suppressMessages(as_numeric(
    file.path(test_dir, "test.hmp.txt"), to_r = TRUE, verbose = FALSE))
  vcf <- suppressMessages(as_numeric(
    file.path(test_dir, "test.vcf"), to_r = TRUE, verbose = FALSE))
  gds <- suppressMessages(as_numeric(
    file.path(test_dir, "test.gds"), to_r = TRUE, verbose = FALSE))

  expect_equal(.meta_cols(hmp), expected)
  expect_equal(.meta_cols(vcf), expected)
  expect_equal(.meta_cols(gds), expected)
})
