test_dir  <- normalizePath(file.path(testthat::test_path(), ".."))
hmp_file  <- file.path(test_dir, "test.hmp.txt")
vcf_file  <- file.path(test_dir, "test.vcf")
bed_file  <- file.path(test_dir, "test.bed")
ped_file  <- file.path(test_dir, "test.ped")
gds_file  <- file.path(test_dir, "test.gds")
snp55k_hmp <- file.path(test_dir, "SNP55K_maize282_AGPv2_20100513_1.hmp.txt")
snp55k_ref <- file.path(test_dir, "SNP55K_maize282_AGPv2_20100513_NUM.txt")

N_SNP  <- 1077L
N_SAMP <- 659L

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

.quiet <- function(expr) suppressMessages(suppressWarnings(expr))

# All formats now return the same 5-column metadata schema:
#   snp, allele, chr, pos, cm  +  sample columns
META_COLS <- 5L

.geno_mat <- function(result) {
  as.matrix(result[, (META_COLS + 1L):ncol(result)])
}

# ---------------------------------------------------------------------------
# 1. HapMap file path — 5-col schema
# ---------------------------------------------------------------------------

test_that("as_numeric() handles HapMap file path", {
  res <- .quiet(as_numeric(hmp_file, to_r = TRUE, verbose = FALSE))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), N_SNP)
  expect_equal(ncol(res) - META_COLS, N_SAMP)
  mat <- .geno_mat(res)
  expect_true(all(mat[!is.na(mat)] %in% c(-1L, 0L, 1L)))
})

# ---------------------------------------------------------------------------
# 2. HapMap in-memory data.frame
# ---------------------------------------------------------------------------

test_that("as_numeric() handles HapMap in-memory data.frame", {
  hmp_df <- data.table::fread(hmp_file, data.table = FALSE)
  res <- .quiet(as_numeric(hmp_df, to_r = TRUE, verbose = FALSE))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), N_SNP)
  expect_equal(ncol(res) - META_COLS, N_SAMP)
})

# ---------------------------------------------------------------------------
# 3. In-memory result matches file-path result
# ---------------------------------------------------------------------------

test_that("in-memory and file-path HapMap produce identical numeric matrices", {
  hmp_df  <- data.table::fread(hmp_file, data.table = FALSE)
  res_obj  <- .quiet(as_numeric(hmp_df,   to_r = TRUE, verbose = FALSE))
  res_file <- .quiet(as_numeric(hmp_file, to_r = TRUE, verbose = FALSE))
  expect_equal(.geno_mat(res_obj), .geno_mat(res_file))
})

# ---------------------------------------------------------------------------
# 4. VCF file path
# ---------------------------------------------------------------------------

test_that("as_numeric() handles VCF file path", {
  res <- .quiet(as_numeric(vcf_file, to_r = TRUE, verbose = FALSE))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), N_SNP)
  expect_equal(ncol(res) - META_COLS, N_SAMP)
  mat <- .geno_mat(res)
  expect_true(all(mat[!is.na(mat)] %in% c(-1L, 0L, 1L)))
})

# ---------------------------------------------------------------------------
# 5. BED file path
# ---------------------------------------------------------------------------

test_that("as_numeric() handles BED file path", {
  res <- .quiet(as_numeric(bed_file, to_r = TRUE, verbose = FALSE))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), N_SNP)
  expect_equal(ncol(res) - META_COLS, N_SAMP)
  mat <- .geno_mat(res)
  expect_true(all(mat[!is.na(mat)] %in% c(-1L, 0L, 1L)))
})

# ---------------------------------------------------------------------------
# 6. PED file path
# ---------------------------------------------------------------------------

test_that("as_numeric() handles PED file path", {
  res <- .quiet(as_numeric(ped_file, to_r = TRUE, verbose = FALSE))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), N_SNP)
  expect_equal(ncol(res) - META_COLS, N_SAMP)
  mat <- .geno_mat(res)
  expect_true(all(mat[!is.na(mat)] %in% c(-1L, 0L, 1L)))
})

# ---------------------------------------------------------------------------
# 7. GDS file path
# ---------------------------------------------------------------------------

test_that("as_numeric() handles GDS file path", {
  res <- .quiet(as_numeric(gds_file, to_r = TRUE, verbose = FALSE))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), N_SNP)
  expect_equal(ncol(res) - META_COLS, N_SAMP)
  mat <- .geno_mat(res)
  expect_true(all(mat[!is.na(mat)] %in% c(-1L, 0L, 1L)))
})

# ---------------------------------------------------------------------------
# 8. 012 coding
# ---------------------------------------------------------------------------

test_that("code_as='012' produces values in {0, 1, 2}", {
  res <- .quiet(as_numeric(hmp_file, to_r = TRUE, code_as = "012", verbose = FALSE))
  mat <- .geno_mat(res)
  expect_true(all(mat[!is.na(mat)] %in% c(0L, 1L, 2L)))
})

# ---------------------------------------------------------------------------
# 9. Missing-data imputation
# ---------------------------------------------------------------------------

test_that("impute='Middle' removes all NAs from HapMap result", {
  res <- .quiet(as_numeric(hmp_file, to_r = TRUE, impute = "Middle", verbose = FALSE))
  mat <- .geno_mat(res)
  expect_equal(sum(is.na(mat)), 0L)
})

# ---------------------------------------------------------------------------
# 10. Cross-format: SNP IDs identical
# ---------------------------------------------------------------------------

test_that("cross-format: SNP IDs are identical across all input formats", {
  hmp <- .quiet(as_numeric(hmp_file, to_r = TRUE, verbose = FALSE))
  vcf <- .quiet(as_numeric(vcf_file, to_r = TRUE, verbose = FALSE))
  bed <- .quiet(as_numeric(bed_file, to_r = TRUE, verbose = FALSE))
  ped <- .quiet(as_numeric(ped_file, to_r = TRUE, verbose = FALSE))
  gds <- .quiet(as_numeric(gds_file, to_r = TRUE, verbose = FALSE))
  ref_ids <- hmp[[1]]
  expect_equal(vcf[[1]], ref_ids)
  expect_equal(bed[[1]], ref_ids)
  expect_equal(ped[[1]], ref_ids)
  expect_equal(gds[[1]], ref_ids)
})

# ---------------------------------------------------------------------------
# 11. Cross-format: sample names consistent
# ---------------------------------------------------------------------------

test_that("cross-format: sample names are consistent across all input formats", {
  hmp <- .quiet(as_numeric(hmp_file, to_r = TRUE, verbose = FALSE))
  vcf <- .quiet(as_numeric(vcf_file, to_r = TRUE, verbose = FALSE))
  bed <- .quiet(as_numeric(bed_file, to_r = TRUE, verbose = FALSE))
  ped <- .quiet(as_numeric(ped_file, to_r = TRUE, verbose = FALSE))
  gds <- .quiet(as_numeric(gds_file, to_r = TRUE, verbose = FALSE))
  # All formats use META_COLS = 5 metadata columns
  ref_samp <- sort(colnames(hmp)[(META_COLS + 1L):ncol(hmp)])
  expect_equal(sort(colnames(vcf)[(META_COLS + 1L):ncol(vcf)]), ref_samp)
  expect_equal(sort(colnames(bed)[(META_COLS + 1L):ncol(bed)]), ref_samp)
  expect_equal(sort(colnames(ped)[(META_COLS + 1L):ncol(ped)]), ref_samp)
  expect_equal(sort(colnames(gds)[(META_COLS + 1L):ncol(gds)]), ref_samp)
})

# ---------------------------------------------------------------------------
# 12. Cross-format: heterozygosity rates equal (orientation-independent)
# ---------------------------------------------------------------------------

test_that("cross-format: heterozygosity rates are equal across all formats", {
  het_rate <- function(mat) mean(mat == 0L, na.rm = TRUE)
  hmp <- .quiet(as_numeric(hmp_file, to_r = TRUE, verbose = FALSE))
  vcf <- .quiet(as_numeric(vcf_file, to_r = TRUE, verbose = FALSE))
  bed <- .quiet(as_numeric(bed_file, to_r = TRUE, verbose = FALSE))
  ped <- .quiet(as_numeric(ped_file, to_r = TRUE, verbose = FALSE))
  gds <- .quiet(as_numeric(gds_file, to_r = TRUE, verbose = FALSE))
  h_ref <- het_rate(as.matrix(hmp[, (META_COLS + 1L):ncol(hmp)]))
  expect_equal(het_rate(as.matrix(vcf[, (META_COLS + 1L):ncol(vcf)])), h_ref, tolerance = 1e-6)
  expect_equal(het_rate(as.matrix(bed[, (META_COLS + 1L):ncol(bed)])), h_ref, tolerance = 1e-6)
  expect_equal(het_rate(as.matrix(ped[, (META_COLS + 1L):ncol(ped)])), h_ref, tolerance = 1e-6)
  expect_equal(het_rate(as.matrix(gds[, (META_COLS + 1L):ncol(gds)])), h_ref, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 13. Cross-format: per-SNP orientation consistency
#     Differences between formats must be exactly ±2 or 0 (valid flips).
# ---------------------------------------------------------------------------

test_that("cross-format: per-SNP allele coding is orientation-consistent", {
  hmp <- .quiet(as_numeric(hmp_file, to_r = TRUE, verbose = FALSE))
  vcf <- .quiet(as_numeric(vcf_file, to_r = TRUE, verbose = FALSE))
  bed <- .quiet(as_numeric(bed_file, to_r = TRUE, verbose = FALSE))
  ped <- .quiet(as_numeric(ped_file, to_r = TRUE, verbose = FALSE))
  gds <- .quiet(as_numeric(gds_file, to_r = TRUE, verbose = FALSE))

  hmp_samp <- colnames(hmp)[(META_COLS + 1L):ncol(hmp)]
  hmp_mat  <- as.matrix(hmp[, (META_COLS + 1L):ncol(hmp)])

  check_orientation <- function(other, label) {
    other_mat <- as.matrix(other[, (META_COLS + 1L):ncol(other)])[, hmp_samp]
    diffs <- as.vector(hmp_mat - other_mat)
    diffs <- diffs[!is.na(diffs)]
    expect_true(all(diffs %in% c(-2L, 0L, 2L)),
                info = paste("orientation check failed for", label))
  }
  check_orientation(vcf, "VCF")
  check_orientation(bed, "BED")
  check_orientation(ped, "PED")
  check_orientation(gds, "GDS")
})

# ---------------------------------------------------------------------------
# 14. SNP55K: output structure matches NUM.txt reference
# ---------------------------------------------------------------------------

test_that("SNP55K HapMap matches NUM.txt reference structure", {
  skip_if_not(file.exists(snp55k_hmp), "SNP55K HapMap file not available")
  skip_if_not(file.exists(snp55k_ref), "SNP55K NUM reference file not available")

  ref <- data.table::fread(snp55k_ref, data.table = FALSE)
  res <- .quiet(as_numeric(snp55k_hmp, to_r = TRUE,
                            code_as = "012", impute = "Middle", verbose = FALSE))

  # Dimensions: same SNP count, same sample count
  expect_equal(nrow(res), nrow(ref))
  # res has 5 meta cols, ref has 5 meta cols
  expect_equal(ncol(res) - META_COLS, ncol(ref) - META_COLS)

  # SNP identifiers
  expect_equal(res[[1]], ref[[1]])
  expect_equal(sort(colnames(res)[(META_COLS + 1L):ncol(res)]),
               sort(colnames(ref)[(META_COLS + 1L):ncol(ref)]))

  # Encoding: values in {0, 1, 2} with no NAs after middle imputation
  mat <- as.matrix(res[, (META_COLS + 1L):ncol(res)])
  expect_true(all(mat %in% c(0L, 1L, 2L)))
  expect_equal(sum(is.na(mat)), 0L)
})

# ---------------------------------------------------------------------------
# 15. Error handling: non-existent file
# ---------------------------------------------------------------------------

test_that("as_numeric() errors on a non-existent file", {
  expect_error(
    as_numeric("no_such_file.hmp.txt", to_r = TRUE, verbose = FALSE),
    regexp = "not found"
  )
})
