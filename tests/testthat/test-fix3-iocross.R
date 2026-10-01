# test-fix3-iocross.R
#
# Round-3 fixes (io + crossing group): R3-5 case-insensitive label fallback in
# the cross-pool orientation guard, R3-6 exactly-11-column HapMap objects,
# R3-7 whole-number double genotype columns become integer, R3-12 exact
# heterosis retention criteria (numeric check of the documented wording).

fx3_hmp <- function(calls, alleles = "A/G") {
  m <- nrow(calls)
  calls <- as.matrix(calls)
  colnames(calls) <- paste0("L", seq_len(ncol(calls)))
  meta <- data.frame(
    `rs#` = paste0("snp", seq_len(m)), alleles = rep_len(alleles, m),
    chrom = 1L, pos = seq_len(m) * 100L, strand = "+", `assembly#` = NA,
    center = NA, protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA,
    check.names = FALSE, stringsAsFactors = FALSE)
  cbind(meta, as.data.frame(calls, stringsAsFactors = FALSE))
}

# ---------------------------------------------------------------------------
# R3-5
# ---------------------------------------------------------------------------

test_that("R3-5: legacy label fallback ignores allele-label case", {
  mk <- function(allele) {
    n <- suppressMessages(as_numeric(fx3_hmp(matrix("AA", 3, 4)), to_r = TRUE,
                                     verbose = FALSE))
    attr(n, "counted_allele") <- NULL          # legacy input: labels only
    n$allele <- allele
    n$cm <- synthetic_map(n$chr, n$pos)
    as_population(n, pool = NA)
  }
  pa <- mk(rep("a/g", 3))
  pb <- mk(rep("A/G", 3))
  expect_no_error(suppressWarnings(simplePHENOTYPES:::.check_orientation(pa$map, pb$map)))
  expect_no_warning(simplePHENOTYPES:::.check_orientation(pa$map, pb$map))
  # opposite order still warns, disjoint alleles still error, regardless of case
  expect_warning(simplePHENOTYPES:::.check_orientation(mk(rep("g/a", 3))$map, pb$map),
                 "different order")
  expect_error(simplePHENOTYPES:::.check_orientation(mk(rep("c/t", 3))$map, pb$map),
               "different alleles")
})

# ---------------------------------------------------------------------------
# R3-6
# ---------------------------------------------------------------------------

test_that("R3-6: 11 HapMap metadata columns and no sample is not HapMap", {
  h <- fx3_hmp(matrix("AA", 3, 1))
  h0 <- h[, 1:11]
  expect_false(simplePHENOTYPES:::.hmp_header_match(names(h0)))
  expect_true(simplePHENOTYPES:::.hmp_header_match(names(h)))
  res <- tryCatch(as_numeric(h0, to_r = TRUE, verbose = FALSE),
                  error = function(e) e)
  expect_s3_class(res, "error")
  expect_false(grepl("subscript out of bounds", conditionMessage(res)))
  # a real one-sample HapMap table still converts
  expect_equal(ncol(suppressMessages(as_numeric(h, to_r = TRUE, verbose = FALSE))), 6L)
})

# ---------------------------------------------------------------------------
# R3-7
# ---------------------------------------------------------------------------

test_that("R3-7: whole-number double genotype columns become integer", {
  df <- data.frame(snp = c("a", "b"), allele = "A/G", chr = "1", pos = c(1, 2),
                   cm = c(0, 1), s_double = c(1, -1), s_zero = c(0, 0),
                   s_na = c(NA, 1), s_logical = c(TRUE, FALSE),
                   stringsAsFactors = FALSE)
  out <- suppressMessages(as_numeric(df, to_r = TRUE, verbose = FALSE))
  for (j in 6:9) expect_type(out[[j]], "integer")
  expect_equal(out$s_double, c(1, -1))
  expect_equal(out$s_na, c(NA, 1))
  d012 <- df[, 1:7]; d012$s_double <- c(2, 0)
  o2 <- suppressMessages(as_numeric(d012, to_r = TRUE, verbose = FALSE,
                                    code_as = "012"))
  expect_type(o2$s_double, "integer")
  expect_equal(o2$s_double, c(2, 0))
  # non-whole values are left as they are
  nz <- simplePHENOTYPES:::.normalize_numeric_schema(
    data.frame(snp = "a", allele = "A/G", chr = "1", pos = 1, cm = 0,
               s = 0.5, stringsAsFactors = FALSE))
  expect_type(nz$s, "double")
  expect_equal(nz$s, 0.5)
})

# ---------------------------------------------------------------------------
# R3-12
# ---------------------------------------------------------------------------

test_that("R3-12: exact retention criteria hold numerically", {
  ret <- function(pA, pB, hA, hB) {
    hAB <- pA * (1 - pB) + (1 - pA) * pB
    H1 <- hAB - (hA + hB) / 2
    g <- (pA + pB) / 2
    hBC <- g * (1 - pA) + (1 - g) * pA        # F1 x breed A (gametes independent)
    hF2 <- 2 * g * (1 - g)
    c(bc = (hBC - 0.75 * hA - 0.25 * hB) / H1,
      f2 = (hF2 - 0.5 * hA - 0.5 * hB) / H1)
  }
  pA <- 0.3; pB <- 0.6
  hwe <- function(p) 2 * p * (1 - p)
  # deviations -0.12 / +0.12: F2 keeps 1/2 with neither breed in HWE
  r <- ret(pA, pB, hwe(pA) - 0.12, hwe(pB) + 0.12)
  expect_equal(unname(r["f2"]), 0.5)
  expect_false(isTRUE(all.equal(unname(r["bc"]), 0.5)))
  # only the recurrent breed A in HWE: backcross keeps 1/2, F2 does not
  r <- ret(pA, pB, hwe(pA), hwe(pB) + 0.1)
  expect_equal(unname(r["bc"]), 0.5)
  expect_false(isTRUE(all.equal(unname(r["f2"]), 0.5)))
  # the docstring no longer claims both breeds must be in HWE
})
