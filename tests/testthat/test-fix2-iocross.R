# test-fix2-iocross.R
#
# Round-2 fixes from the independent review (io + crossing group):
# A2 counted-allele record and the cross-pool orientation guard, B3 VCF file
# validation, B7 symmetric map identity, B8 total chromosome order, B13 short
# HapMap headers, B14 uniform numeric schema, C12 heterosis retention text.

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# A 11-metadata-column HapMap table; `calls` is markers x individuals.
fx_hmp <- function(calls, alleles = "A/G", chrom = 1L, pos = NULL) {
  m <- nrow(calls)
  pos <- if (is.null(pos)) seq_len(m) * 100L else pos
  calls <- as.matrix(calls)
  colnames(calls) <- paste0("L", seq_len(ncol(calls)))
  meta <- data.frame(
    `rs#` = paste0("snp", seq_len(m)), alleles = rep_len(alleles, m),
    chrom = rep_len(chrom, m), pos = pos, strand = "+", `assembly#` = NA,
    center = NA, protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA,
    check.names = FALSE, stringsAsFactors = FALSE)
  cbind(meta, as.data.frame(calls, stringsAsFactors = FALSE))
}

fx_num <- function(calls, ...) {
  suppressMessages(as_numeric(fx_hmp(calls, ...), to_r = TRUE, verbose = FALSE))
}

fx_pop <- function(num, pool = NA, ids = NULL) {
  num$cm <- synthetic_map(num$chr, num$pos)
  if (!is.null(ids)) names(num)[-(1:5)] <- ids
  as_population(num, pool = pool)
}

# ---------------------------------------------------------------------------
# A2: the counted (+1) allele is recorded and used by the cross-pool guard
# ---------------------------------------------------------------------------

test_that("as_numeric() records the counted allele without changing any dosage", {
  na <- fx_num(matrix("AA", 3, 4))
  ng <- fx_num(matrix("GG", 3, 4))
  # the label is the raw orientation for both (that is the reported hazard) ...
  expect_identical(na$allele, rep("A/G", 3))
  expect_identical(ng$allele, rep("A/G", 3))
  # ... both encode +1 (dosage values are unchanged) ...
  expect_true(all(na[, -(1:5)] == 1L))
  expect_true(all(ng[, -(1:5)] == 1L))
  # ... and the record tells the two panels apart
  expect_identical(attr(na, "counted_allele"), rep("A", 3))
  expect_identical(attr(ng, "counted_allele"), rep("G", 3))
})

test_that("counted allele follows frequency, reference and code_as", {
  calls <- rbind(c("AA", "AA", "AG", "GG"),    # A major
                 c("AA", "GG", "GG", "GG"))    # G major
  h <- fx_hmp(calls, alleles = c("A/G", "A/G"))
  n1 <- suppressMessages(as_numeric(h, to_r = TRUE, verbose = FALSE))
  expect_identical(attr(n1, "counted_allele"), c("A", "G"))
  n2 <- suppressMessages(as_numeric(h, to_r = TRUE, verbose = FALSE,
                                    method = "reference",
                                    ref_allele = c("G", "A")))
  expect_identical(attr(n2, "counted_allele"), c("G", "A"))
  n3 <- suppressMessages(as_numeric(h, to_r = TRUE, verbose = FALSE,
                                    code_as = "012"))
  expect_identical(attr(n3, "counted_allele"), c("A", "G"))
  # dominance counts no allele
  n4 <- suppressMessages(as_numeric(h, to_r = TRUE, verbose = FALSE,
                                    model = "Dom"))
  expect_null(attr(n4, "counted_allele"))
  # dosages are those of the (unchanged) coding
  expect_equal(as.numeric(unlist(n1[1, -(1:5)])), c(1, 1, 0, -1))
  expect_equal(as.numeric(unlist(n1[2, -(1:5)])), c(-1, 1, 1, 1))
})

test_that("crossing all-AA x all-GG panels converted separately is refused", {
  pa <- fx_pop(fx_num(matrix("AA", 3, 3)), pool = "P1", ids = paste0("a", 1:3))
  pg <- fx_pop(fx_num(matrix("GG", 3, 3)), pool = "P2", ids = paste0("g", 1:3))
  expect_error(cross(pa[1], pg[1], n = 1, seed = 1),
               "count different alleles as \\+1")
  expect_error(c(pa, pg), "count different alleles as \\+1")
  # converted jointly (or to one reference), the same panels cross correctly:
  # the F1 of an AA and a GG parent is the heterozygote, dosage 0
  joint <- suppressMessages(as_numeric(
    fx_hmp(cbind(matrix("AA", 3, 3), matrix("GG", 3, 3))),
    to_r = TRUE, verbose = FALSE, method = "reference",
    ref_allele = rep("A", 3)))
  names(joint)[-(1:5)] <- c(paste0("a", 1:3), paste0("g", 1:3))
  joint$cm <- synthetic_map(joint$chr, joint$pos)
  pj <- as_population(joint)
  f1 <- cross(pj[1], pj[4], n = 2, seed = 1)
  expect_true(all(dosages(f1) == 0))
})

test_that("panels that count the same allele cross without a warning", {
  a1 <- fx_pop(fx_num(matrix(c("AA", "AG"), 3, 4)), pool = "P1",
               ids = paste0("a", 1:4))
  a2 <- fx_pop(fx_num(matrix(c("AA", "AA", "AG"), 3, 4)), pool = "P2",
               ids = paste0("b", 1:4))
  expect_no_warning(cross(a1[1], a2[1], n = 1, seed = 3))
  expect_no_warning(c(a1, a2))
})

test_that("numeric data without the record keeps the label-only behaviour", {
  strip <- function(num) {
    attr(num, "counted_allele") <- NULL
    num
  }
  na <- fx_num(matrix("AA", 3, 2))
  ng <- fx_num(matrix("GG", 3, 2))
  pa <- fx_pop(na, ids = c("a1", "a2"))
  la <- fx_pop(strip(na), ids = c("a1", "a2"))
  lg <- fx_pop(strip(ng), ids = c("g1", "g2"))
  expect_null(la$map$counted)
  expect_false(is.null(pa$map$counted))
  # identical "A/G" labels and no record: nothing to compare, as before
  expect_no_warning(cross(la[1], lg[1], n = 1, seed = 1))
  # one side with, one without: falls back to the labels (equal -> silent)
  expect_no_warning(cross(pa[1], lg[1], n = 1, seed = 1))
  # the label check still works when the record is absent
  ng2 <- strip(ng)
  ng2$allele <- rep("G/A", 3)
  expect_warning(
    cross(la[1], fx_pop(ng2, ids = c("g1", "g2"))[1], n = 1, seed = 1),
    "different order")
})

test_that(".check_orientation() compares the record marker by marker", {
  co <- simplePHENOTYPES:::.check_orientation
  m <- function(cnt, allele = "A/G") {
    d <- data.frame(snp = paste0("m", 1:3), allele = allele,
                    stringsAsFactors = FALSE)
    if (!is.null(cnt)) d$counted <- cnt
    d
  }
  expect_true(co(m(c("A", "A", "A")), m(c("A", "A", "A"))))
  # a single opposite marker is enough, and is named
  expect_error(co(m(c("A", "A", "G")), m(c("A", "A", "A"))), "m3")
  # NA on either side: not compared
  expect_true(co(m(c("A", NA, "G")), m(c("A", "G", NA))))
  # a differing label order is not a problem when the counted alleles agree
  expect_true(co(m(c("A", "A", "A"), "A/G"), m(c("A", "A", "A"), "G/A")))
  # a marker without a record still gets the label check
  expect_warning(co(m(c("A", "A", NA), c("A/G", "A/G", "A/G")),
                    m(c("A", "A", "A"), c("A/G", "A/G", "G/A"))),
                 "different order")
})

test_that("the counted-allele record does not persist through a text file", {
  num <- fx_num(matrix(c("AA", "AG", "GG"), 2, 3), alleles = c("A/G", "A/G"))
  expect_false(is.null(attr(num, "counted_allele")))
  f <- tempfile(fileext = ".txt")
  data.table::fwrite(num, f, sep = "\t", na = NA)
  again <- suppressMessages(as_numeric(f, to_r = TRUE, verbose = FALSE))
  expect_null(attr(again, "counted_allele"))
  # ... but an in-memory numeric table passes through with it
  same <- suppressMessages(as_numeric(num, to_r = TRUE, verbose = FALSE))
  expect_identical(attr(same, "counted_allele"), attr(num, "counted_allele"))
})

test_that("VCF (file path) records the counted allele", {
  skip_if_not_installed("SNPRelate")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc",
               "1\t100\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1",
               "1\t200\tm2\tC\tT\t.\tPASS\t.\tGT\t1/1\t1/1\t0/1"), vcf)
  from_file <- suppressMessages(as_numeric(vcf, to_r = TRUE, verbose = FALSE))
  expect_identical(attr(from_file, "counted_allele"), c("A", "T"))
})

# ---------------------------------------------------------------------------
# B3: a VCF file path gets the same complete-diploid / biallelic validation
# ---------------------------------------------------------------------------

fx_vcf_lines <- function(header, rows) {
  c("##fileformat=VCFv4.2", paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                    "FILTER", "INFO", "FORMAT", header),
                                  collapse = "\t"), rows)
}

test_that("haploid and partial calls in a VCF FILE are missing, with a warning", {
  skip_if_not_installed("SNPRelate")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(fx_vcf_lines(c("dipref", "dipalt", "hapref", "hapalt", "partial"),
                          "1\t100\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t0\t1\t./1"),
             vcf)
  expect_warning(
    r <- suppressMessages(as_numeric(vcf, to_r = TRUE, verbose = FALSE)),
    "3 VCF genotype call\\(s\\) are partially missing or haploid")
  # the valid tied diploids keep the documented coding (REF homozygote = +1)
  expect_identical(r$dipref, 1L)
  expect_identical(r$dipalt, -1L)
  expect_true(all(is.na(unlist(r[, c("hapref", "hapalt", "partial")]))))
})

test_that("multiallelic calls in a VCF FILE are missing, with a counted warning", {
  skip_if_not_installed("SNPRelate")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(fx_vcf_lines(c("s1", "s2", "s3", "s4"),
                          "1\t100\tm1\tA\tG,T\t.\tPASS\t.\tGT\t0/0\t1/1\t0/2\t2/2"),
             vcf)
  expect_warning(
    r <- suppressMessages(as_numeric(vcf, to_r = TRUE, verbose = FALSE)),
    "2 VCF genotype call\\(s\\) reference alleles beyond the first ALT")
  expect_identical(c(r$s1, r$s2), c(1L, -1L))
  expect_true(is.na(r$s3) && is.na(r$s4))
})

test_that("a VCF gives the same table from a file path and from a data frame", {
  skip_if_not_installed("SNPRelate")
  rows <- c(
    "1\t100\tm1\tA\tG\t.\tPASS\t.\tGT:DP\t0/0:9\t1/1:8\t0/1:7\t./.:0\t1|1:5\t0",
    "1\t200\tm2\tC\tT\t.\tPASS\t.\tGT:DP\t1/1:9\t1/1:8\t0/0:7\t0/1:1\t0|0:5\t./1",
    "2\t300\tm3\tG\tA,C\t.\tPASS\t.\tGT:DP\t0/0:9\t0/2:8\t0/1:7\t1/1:1\t2/2:5\t0/0")
  hdr <- c("s1", "s2", "s3", "s4", "s5", "s6")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(fx_vcf_lines(hdr, rows), vcf)
  df <- utils::read.delim(vcf, skip = 1, check.names = FALSE, header = TRUE,
                          stringsAsFactors = FALSE, colClasses = "character")
  expect_warning(
    a <- suppressMessages(as_numeric(vcf, to_r = TRUE, verbose = FALSE)),
    "VCF genotype call")
  expect_warning(
    b <- suppressMessages(as_numeric(df, from = "vcf", to_r = TRUE,
                                     verbose = FALSE)),
    "VCF genotype call")
  ma <- as.matrix(a[, hdr])
  mb <- as.matrix(b[, hdr])
  dimnames(ma) <- dimnames(mb) <- NULL
  expect_identical(ma, mb)
  expect_identical(a$allele, b$allele)
  expect_identical(attr(a, "counted_allele"), attr(b, "counted_allele"))
})

test_that("a gzipped VCF file is validated too", {
  skip_if_not_installed("SNPRelate")
  vcf <- tempfile(fileext = ".vcf.gz")
  con <- gzfile(vcf, "w")
  writeLines(fx_vcf_lines(c("dipref", "dipalt", "hap"),
                          "1\t100\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t1"), con)
  close(con)
  expect_warning(
    r <- suppressMessages(as_numeric(vcf, to_r = TRUE, verbose = FALSE)),
    "partially missing or haploid")
  expect_identical(c(r$dipref, r$dipalt), c(1L, -1L))
  expect_true(is.na(r$hap))
})

test_that("clean VCF files give no warning", {
  skip_if_not_installed("SNPRelate")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(fx_vcf_lines(c("a", "b", "c"),
                          c("1\t100\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t0/1",
                            "1\t200\tm2\tC\tT\t.\tPASS\t.\tGT\t0|1\t./.\t0/0")),
             vcf)
  expect_no_warning(suppressMessages(as_numeric(vcf, to_r = TRUE,
                                                verbose = FALSE)))
})

# ---------------------------------------------------------------------------
# B7: map identity is symmetric
# ---------------------------------------------------------------------------

test_that(".same_map() gives the same answer in both argument orders", {
  sm <- simplePHENOTYPES:::.same_map
  mk <- function(pos) data.frame(snp = "m1", chr = "1", pos = pos, cm = 1)
  for (pair in list(c(1e8, 99999999), c(1e8, 99999998), c(1e8, 1e8 + 1),
                    c(1, 1 + 1e-8), c(1, 1 + 2e-8), c(0, 1e-9))) {
    expect_identical(sm(mk(pair[1]), mk(pair[2])),
                     sm(mk(pair[2]), mk(pair[1])))
  }
  expect_true(sm(mk(1e8), mk(99999999)))
  expect_true(sm(mk(99999999), mk(1e8)))
  expect_false(sm(mk(1e8), mk(99999998)))
  expect_false(sm(mk(99999998), mk(1e8)))
})

test_that("cross(A, B) and cross(B, A) agree on whether the maps match", {
  base <- fx_num(matrix(c("AA", "GG"), 2, 2), alleles = "A/G",
                 pos = c(100000000L, 200000000L))
  near <- base
  near$pos <- near$pos - c(1L, 0L)
  base$cm <- near$cm <- c(10, 20)
  pa <- as_population(base)
  pb <- as_population(near)
  r1 <- try(cross(pa[1], pb[2], n = 1, seed = 1), silent = TRUE)
  r2 <- try(cross(pb[2], pa[1], n = 1, seed = 1), silent = TRUE)
  expect_identical(inherits(r1, "try-error"), inherits(r2, "try-error"))
})

# ---------------------------------------------------------------------------
# B8: the chromosome order is a total order (row-order independent)
# ---------------------------------------------------------------------------

test_that(".chr_rank() ranks '1' and '01' the same whichever comes first", {
  cr <- simplePHENOTYPES:::.chr_rank
  expect_identical(cr(c("1", "01")), c(2L, 1L))
  expect_identical(cr(c("01", "1")), c(1L, 2L))
  expect_identical(cr(c("chr1", "chr01", "chr2")), c(2L, 1L, 3L))
  expect_identical(cr(c("chr01", "chr2", "chr1")), c(1L, 3L, 2L))
  # ordinary labels are unchanged
  expect_identical(cr(c("10", "2", "1")), c(3L, 2L, 1L))
  expect_identical(cr(c("chrX", "chr2", "chr10")), c(3L, 1L, 2L))
})

test_that("a seeded mating does not depend on the row order of tied chromosome labels", {
  build <- function(rev_rows) {
    num <- data.frame(snp = c("m1", "m01"), allele = "A/G",
                      chr = c("1", "01"), pos = c(10L, 10L), cm = c(50, 50),
                      P = c(0L, 0L), stringsAsFactors = FALSE)
    if (rev_rows) num <- num[2:1, ]
    as_population(num)
  }
  for (seed in 1:12) {
    a <- dosages(double_haploid(build(FALSE), n = 3, seed = seed))
    b <- dosages(double_haploid(build(TRUE), n = 3, seed = seed))
    expect_identical(a[c("m1", "m01"), , drop = FALSE],
                     b[c("m1", "m01"), , drop = FALSE])
  }
})

# ---------------------------------------------------------------------------
# B13: a 9-column HapMap prefix is not a HapMap object
# ---------------------------------------------------------------------------

test_that("an object with only the first 9 HapMap columns is not detected as HapMap", {
  full <- fx_hmp(matrix("AA", 2, 2))
  h <- full[, 1:9]
  expect_false(identical(detect_format(h), "hapmap"))
  expect_false(simplePHENOTYPES:::.hmp_header_match(names(h)))
  # an error that names the problem, not a subscript error
  expect_error(as_numeric(h, to_r = TRUE, verbose = FALSE),
               "format was not detected|not supported")
  expect_false(identical(detect_format(full[, 1:10]), "hapmap"))
  expect_identical(detect_format(full), "hapmap")
  expect_true(simplePHENOTYPES:::.hmp_header_match(names(full)))
})

# ---------------------------------------------------------------------------
# B14: one numeric output schema
# ---------------------------------------------------------------------------

test_that("a nucleotide table with no positions round-trips with the same types", {
  tab <- data.frame(s1 = c("AA", "GG"), s2 = c("AG", "GG"), s3 = c("AA", "AA"),
                    stringsAsFactors = FALSE)
  first <- suppressMessages(as_numeric(tab, from = "table", to_r = TRUE,
                                       verbose = FALSE))
  expect_type(first$pos, "integer")
  expect_true(all(is.na(first$pos)))
  f <- tempfile(fileext = ".txt")
  data.table::fwrite(first, f, sep = "\t", na = NA)
  again <- suppressMessages(as_numeric(f, to_r = TRUE, verbose = FALSE))
  attr(first, "counted_allele") <- NULL
  expect_identical(again, first)
  expect_type(again$pos, "integer")
  expect_type(again$cm, "double")
  expect_type(again$chr, "character")
  expect_type(again$snp, "character")
})

test_that("in-memory numeric-format input is normalized to the same schema", {
  raw <- data.frame(snp = 1:2, allele = "A/G", chr = c(1, 2), pos = c(100, 200),
                    cm = NA, s1 = c(1, -1), s2 = c(0, 1),
                    stringsAsFactors = FALSE)
  out <- suppressMessages(as_numeric(raw, from = "numeric", to_r = TRUE,
                                     verbose = FALSE))
  expect_type(out$snp, "character")
  expect_type(out$chr, "character")
  expect_type(out$pos, "integer")
  expect_type(out$cm, "double")
  expect_identical(out$snp, c("1", "2"))
  expect_identical(out$chr, c("1", "2"))
  expect_equal(out$s1, c(1, -1))
  # already-normalized input is returned unchanged
  expect_identical(suppressMessages(as_numeric(out, to_r = TRUE,
                                               verbose = FALSE)), out)
  # a position that is not a whole number stays double
  raw$pos <- c(100.5, 200)
  expect_type(suppressMessages(as_numeric(raw, from = "numeric", to_r = TRUE,
                                          verbose = FALSE))$pos, "double")
})

# ---------------------------------------------------------------------------
# C12: heterosis retention -- the corrected condition, checked numerically
# ---------------------------------------------------------------------------

test_that("retention of F1 heterosis is 1/2 exactly under the documented condition", {
  # One locus, pure dominance (a = 0, d = 1). A breed is summarised by its gamete
  # frequency p and its heterozygote frequency h. The expected mean of a cross is
  # its heterozygote frequency, from the package's own .expected_cross_means().
  ecm <- function(g1, g2) as.numeric(simplePHENOTYPES:::.expected_cross_means(
    matrix(g1), matrix(g2), a = 0, d = 1))
  retention <- function(pA, hA, pB, hB, what) {
    q <- (pA + pB) / 2                         # gamete frequency of the F1
    f1 <- ecm(pA, pB) - (hA + hB) / 2
    if (what == "bc") {
      (ecm(q, pA) - (3 * hA + hB) / 4) / f1    # backcross to A
    } else {
      (ecm(q, q) - (hA + hB) / 2) / f1         # F2
    }
  }
  hw <- function(p) 2 * p * (1 - p)
  # (1) the counterexample: A = equal AA / aa mixture (p = 1/2, h = 0), B = aa
  expect_equal(ecm(0.5, 0), 0.5)               # F1 heterozygosity
  expect_equal(ecm(0.25, 0.5), 0.5)            # backcross heterozygosity
  expect_equal(retention(0.5, 0, 0, 0, "bc"), 1)   # not 1/2
  # (2) Hardy-Weinberg breeds retain exactly 1/2 in both
  expect_equal(retention(0.7, hw(0.7), 0.2, hw(0.2), "bc"), 0.5)
  expect_equal(retention(0.7, hw(0.7), 0.2, hw(0.2), "f2"), 0.5)
  # (3) a single fixed inbred line is HW (p = 0 or 1, h = 0): 1/2
  expect_equal(retention(1, 0, 0, 0, "bc"), 0.5)
  expect_equal(retention(1, 0, 0, 0, "f2"), 0.5)
  # (4) the exact criteria stated in ?heterosis, over a grid of breeds
  g <- expand.grid(pA = c(0.1, 0.35, 0.6, 0.9), pB = c(0, 0.25, 0.5, 0.8),
                   hA = c(0, 0.2, 0.5), hB = c(0, 0.3, 0.45))
  g <- g[g$hA <= hw(g$pA) + 1e-12 & g$hB <= hw(g$pB) + 1e-12, ]
  f1 <- mapply(function(pA, hA, pB, hB) ecm(pA, pB) - (hA + hB) / 2,
               g$pA, g$hA, g$pB, g$hB)
  g <- g[abs(f1) > 1e-8, ]
  expect_gt(nrow(g), 20L)
  for (k in seq_len(nrow(g))) {
    r <- g[k, ]
    dA <- r$hA - hw(r$pA)
    dB <- r$hB - hw(r$pB)
    bc <- retention(r$pA, r$hA, r$pB, r$hB, "bc")
    f2 <- retention(r$pA, r$hA, r$pB, r$hB, "f2")
    # backcross to A: 1/2 iff h_A = 2 p_A (1 - p_A)
    expect_identical(abs(bc - 0.5) < 1e-9, abs(dA) < 1e-9)
    # F2: 1/2 iff the two deviations cancel
    expect_identical(abs(f2 - 0.5) < 1e-9, abs(dA + dB) < 1e-9)
  }
})

test_that("?heterosis no longer says that fully inbred lines give one half", {
  skip_if(!file.exists(testthat::test_path("..", "..", "R", "cross_breed.R")))
  txt <- paste(readLines(testthat::test_path("..", "..", "R", "cross_breed.R"),
                         warn = FALSE), collapse = "\n")
  expect_false(grepl("or is made of\n#' fully inbred lines", txt, fixed = TRUE))
  # round 3: the "each breed must be in HWE" wording was itself corrected to
  # the exact deviation criteria
  expect_true(grepl("Hardy-Weinberg deviations", txt, fixed = TRUE))
  expect_false(grepl("hold only when each breed", txt, fixed = TRUE))
})
