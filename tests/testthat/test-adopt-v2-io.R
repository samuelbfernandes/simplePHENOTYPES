# test-adopt-v2-io.R
#
# Remaining proposals of the independent v2 audit (group v2-io-formats: io,
# formats, as_numeric, qc / filter_geno) that the older test-audit-io.R and
# test-fix*-iocross.R files do not already cover:
#   * FinalReport metadata for the Top and AB conventions and the
#     required-column diagnostics (Fable T12/T13, Codex rows 1-2),
#   * a hand-derived end-to-end fixture for every reader (Codex row 5),
#   * filter_geno() MAF boundaries, het filters on inbred data and the
#     individuals-by-markers matrix orientation (Fable T17-T20),
#   * the PLINK-port internals: exhaustive haplotypic r^2 on every length-5
#     genotype pair with analytic checks, the two-locus log-likelihood and the
#     ML root choice, and the Haploview pair classifier (Codex rows 11-13),
#   * PLINK-parity bookkeeping invariants (Codex row 14),
#   * .qtn_var() / .tx_qtn_var() against closed-form shares (Codex row 16),
#   * write_phenotypes() wide / multi-rep / NA / separator round trips,
#   * documentation-contract checks (O7, O8, IO-F15), skipped when the package
#     sources are not available.
# All files are written to tempfile() locations.

geno_of_a <- function(res) as.matrix(res[, -(1:5), drop = FALSE])

write_fr_a <- function(rows,
                       header = "SNP Name\tSample ID\tAllele1 - AB\tAllele2 - AB\tChr\tPosition") {
  path <- tempfile(fileext = ".txt")
  writeLines(c("[Header]", "Num SNPs\t1", "[Data]", header, rows), path)
  path
}

# ---------------------------------------------------------------------------
# FinalReport: allele metadata for the Top and AB conventions (T13), required
# columns (T12, Codex row 2)
# ---------------------------------------------------------------------------

test_that("FinalReport Top and AB files give hand-derived dosages and allele metadata", {
  # Top: m1 = AA AG GG AA (A is the majority allele: 2 AA vs 1 GG -> A is +1)
  #      m2 = CC TT CT TC (tied 1 CC vs 1 TT: allele 1 is the alphabetically
  #      first allele, C; both CT and TC are heterozygotes)
  top <- write_fr_a(
    c("m1\ts1\tA\tA\t1\t100", "m1\ts2\tA\tG\t1\t100",
      "m1\ts3\tG\tG\t1\t100", "m1\ts4\tA\tA\t1\t100",
      "m2\ts1\tC\tC\t2\t50",  "m2\ts2\tT\tT\t2\t50",
      "m2\ts3\tC\tT\t2\t50",  "m2\ts4\tT\tC\t2\t50"),
    header = "SNP Name\tSample ID\tAllele1 - Top\tAllele2 - Top\tChr\tPosition")
  r <- suppressMessages(as_numeric(top, to_r = TRUE, verbose = FALSE))
  expect_identical(r$snp, c("m1", "m2"))
  expect_identical(r$allele, c("A/G", "C/T"))
  expect_false(anyNA(r$allele))
  expect_identical(r$chr, c("1", "2"))
  expect_identical(r$pos, c(100L, 50L))
  expect_true(all(is.na(r$cm)))
  expect_equal(unname(geno_of_a(r)),
               matrix(c(1, 0, -1, 1,
                        1, -1, 0, 0), 2, byrow = TRUE),
               ignore_attr = TRUE)
  expect_identical(attr(r, "counted_allele"), c("A", "C"))

  # AB: m1 = AA AB BB AA ; m2 = BB BA BB BB (B is the majority allele)
  ab <- write_fr_a(
    c("m1\ts1\tA\tA\t1\t100", "m1\ts2\tA\tB\t1\t100",
      "m1\ts3\tB\tB\t1\t100", "m1\ts4\tA\tA\t1\t100",
      "m2\ts1\tB\tB\t2\t50",  "m2\ts2\tB\tA\t2\t50",
      "m2\ts3\tB\tB\t2\t50",  "m2\ts4\tB\tB\t2\t50"))
  r2 <- suppressMessages(as_numeric(ab, to_r = TRUE, verbose = FALSE))
  expect_identical(r2$allele[1L], "A/B")
  expect_false(anyNA(r2$allele))
  expect_setequal(strsplit(r2$allele[2L], "/")[[1L]], c("A", "B"))
  expect_equal(unname(geno_of_a(r2)),
               matrix(c(1, 0, -1, 1,
                        1, 0, 1, 1), 2, byrow = TRUE),
               ignore_attr = TRUE)
})

test_that("FinalReport without Allele1, or without a sample column, names what it found", {
  no_a1 <- write_fr_a("m1\ts1\tA\t1\t100",
                      header = "SNP Name\tSample ID\tAllele2 - AB\tChr\tPosition")
  expect_error(as_numeric(no_a1, from = "finalreport", to_r = TRUE,
                          verbose = FALSE),
               "cannot find required columns.*Allele2 - AB")
  no_smp <- write_fr_a("m1\tA\tA\t1\t100",
                       header = "SNP Name\tAllele1 - AB\tAllele2 - AB\tChr\tPosition")
  expect_error(as_numeric(no_smp, from = "finalreport", to_r = TRUE,
                          verbose = FALSE),
               "cannot find required columns")
})

test_that("FinalReport without Chr / Position columns still converts, with an NA map", {
  fr <- write_fr_a(c("m1\ts1\tA\tA", "m1\ts2\tB\tB", "m1\ts3\tA\tB"),
                   header = "SNP Name\tSample ID\tAllele1 - AB\tAllele2 - AB")
  r <- suppressMessages(as_numeric(fr, to_r = TRUE, verbose = FALSE))
  expect_true(is.na(r$chr) && is.na(r$pos))
  # one AA, one BB: tied, so allele 1 (alphabetically first, A) is +1
  expect_equal(as.numeric(geno_of_a(r)), c(1, -1, 0))
})

# ---------------------------------------------------------------------------
# Every reader on one hand-derived fixture (Codex row 5)
# ---------------------------------------------------------------------------

test_that("every reader gives the same hand-derived dosages, counted alleles and map", {
  skip_if_not_installed("SNPRelate")
  # snpA (A/G): 0/0 0/0 0/1 1/1 ./.  -> A: 5 copies, G: 3 -> A is +1: 1 1 0 -1 NA
  # snpB (C/T): 1/1 1/1 1/1 0/1 0/0  -> T: 7 copies, C: 3 -> T is +1: 1 1 1 0 -1
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\ts5",
    "1\t100\tsnpA\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1\t1/1\t./.",
    "2\t50\tsnpB\tC\tT\t.\tPASS\t.\tGT\t1/1\t1/1\t1/1\t0/1\t0/0"), vcf)
  expected <- rbind(snpA = c(1, 1, 0, -1, NA), snpB = c(1, 1, 1, 0, -1))
  expected_counted <- c("A", "T")

  check <- function(res, label, check_ids = TRUE) {
    expect_equal(unname(geno_of_a(res)), unname(expected), ignore_attr = TRUE,
                 info = label)
    expect_identical(names(res)[-(1:5)], paste0("s", 1:5), info = label)
    expect_identical(res$chr, c("1", "2"), info = label)
    expect_identical(res$pos, c(100L, 50L), info = label)
    expect_identical(attr(res, "counted_allele"), expected_counted, info = label)
    expect_identical(vapply(strsplit(res$allele, "/", fixed = TRUE),
                            function(a) paste(sort(a), collapse = "/"), ""),
                     c("A/G", "C/T"), info = label)
    if (check_ids) expect_identical(res$snp, c("snpA", "snpB"), info = label)
  }

  check(suppressMessages(as_numeric(vcf, to_r = TRUE, verbose = FALSE)),
        "VCF file")
  vdf <- read.table(vcf, comment.char = "", skip = 1, header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)
  names(vdf)[1L] <- "#CHROM"
  check(suppressMessages(as_numeric(vdf, from = "vcf", to_r = TRUE,
                                    verbose = FALSE)), "VCF data frame")

  gds <- tempfile(fileext = ".gds")
  SNPRelate::snpgdsVCF2GDS(vcf, gds, method = "copy.num.of.ref",
                           verbose = FALSE)
  gf <- SNPRelate::snpgdsOpen(gds)
  bed_base <- tempfile("e2e_bed")
  ped_base <- tempfile("e2e_ped")
  SNPRelate::snpgdsGDS2BED(gf, bed_base, verbose = FALSE)
  SNPRelate::snpgdsGDS2PED(gf, ped_base, verbose = FALSE)
  SNPRelate::snpgdsClose(gf)
  check(suppressMessages(as_numeric(gds, to_r = TRUE, verbose = FALSE)), "GDS")
  # PLINK files written by SNPRelate carry numeric marker ids in the .bim
  check(suppressMessages(as_numeric(paste0(bed_base, ".bed"), to_r = TRUE,
                                    verbose = FALSE)), "BED", check_ids = FALSE)
  check(suppressMessages(as_numeric(paste0(ped_base, ".ped"), to_r = TRUE,
                                    verbose = FALSE)), "PED", check_ids = FALSE)
})

test_that("a nucleotide table and numeric data frame give hand-derived output", {
  tb <- data.frame(s1 = c("AA", "CC"), s2 = c("AA", "CC"), s3 = c("AG", "CT"),
                   s4 = c("GG", "TT"), s5 = c(NA, "CC"),
                   stringsAsFactors = FALSE)
  r <- suppressMessages(as_numeric(tb, from = "table", to_r = TRUE,
                                   verbose = FALSE))
  # row 1: 2 AA vs 1 GG -> A is +1 ; row 2: 3 CC vs 1 TT -> C is +1
  expect_equal(unname(geno_of_a(r)),
               rbind(c(1, 1, 0, -1, NA), c(1, 1, 0, -1, 1)),
               ignore_attr = TRUE)
  expect_identical(r$allele, c("A/G", "C/T"))
  expect_identical(names(r)[1:5], c("snp", "allele", "chr", "pos", "cm"))
  expect_true(all(is.na(r$chr)) && all(is.na(r$pos)) && all(is.na(r$cm)))
  expect_identical(attr(r, "counted_allele"), c("A", "C"))
  # a numeric data frame passes through unchanged
  r2 <- suppressMessages(as_numeric(as.data.frame(r), from = "numeric",
                                    to_r = TRUE, verbose = FALSE))
  expect_equal(unname(geno_of_a(r2)), unname(geno_of_a(r)), ignore_attr = TRUE)
  expect_identical(r2$snp, r$snp)
  expect_identical(r2$allele, r$allele)
})

# ---------------------------------------------------------------------------
# filter_geno(): MAF boundaries, het filters on inbred data, matrix orientation
# ---------------------------------------------------------------------------

mk_num_a <- function(mat, chr = 1L, pos = NULL) {
  m <- nrow(mat)
  data.frame(snp = paste0("m", seq_len(m)), allele = "A/G", chr = chr,
             pos = if (is.null(pos)) 1000L * seq_len(m) else pos, cm = 0,
             stats::setNames(as.data.frame(mat), paste0("i", seq_len(ncol(mat)))),
             check.names = FALSE, stringsAsFactors = FALSE)
}

test_that("MAF thresholds are inclusive at exactly representable frequencies", {
  # 40 allele copies; the minor allele is coded +1 so that p is the MAF itself
  # and 2/40, 4/40 and 16/40 are the correctly rounded doubles 0.05, 0.1, 0.4
  # (the opposite coding is the known-defect test below).
  d <- mk_num_a(rbind(
    m1 = c(rep(1, 2), rep(-1, 18)),       # MAF 0.10
    m2 = c(rep(1, 1), rep(-1, 19)),       # MAF 0.05
    m3 = c(rep(1, 8), rep(-1, 12))))      # MAF 0.40
  expect_identical(filter_geno(d, maf_above = 0.05, verbose = FALSE)$snp,
                   c("m1", "m2", "m3"))
  expect_identical(filter_geno(d, maf_above = 0.4, verbose = FALSE)$snp, "m3")
  expect_identical(filter_geno(d, maf_below = 0.05, verbose = FALSE)$snp, "m2")
  expect_identical(filter_geno(d, maf_above = 0.1, maf_below = 0.1,
                               verbose = FALSE)$snp, "m1")
  # just above / below the boundary excludes
  expect_identical(filter_geno(d, maf_above = 0.11, verbose = FALSE)$snp, "m3")
  expect_identical(filter_geno(d, maf_below = 0.09, verbose = FALSE)$snp, "m2")
})

test_that("(IO-N1) MAF inclusive boundary does not depend on which allele is +1", {
  # the same MAF 0.1 marker with the minor allele coded -1 (the usual frequency
  # coding) has p = 36/40 = 0.9 and 1 - p = 0.09999999999999998 < 0.1
  d <- mk_num_a(rbind(
    mA = c(rep(1, 2), rep(-1, 18)),       # minor = +1, p = 4/40 = 0.1
    mB = c(rep(-1, 2), rep(1, 18))))      # minor = -1, p = 36/40 = 0.9
  expect_identical(filter_geno(d, maf_above = 0.1, verbose = FALSE)$snp,
                   d$snp)
  expect_identical(filter_geno(d, maf_below = 0.1, verbose = FALSE)$snp,
                   d$snp)
})

test_that("hets = 'include' / 'remove' on a panel with no heterozygote", {
  inbred <- mk_num_a(rbind(c(-1, -1, 1, 1, 1, 1), c(-1, 1, 1, -1, 1, -1)))
  expect_error(filter_geno(inbred, hets = "include", verbose = FALSE),
               "no markers passed")
  expect_identical(filter_geno(inbred, hets = "remove", verbose = FALSE)$snp,
                   c("m1", "m2"))
  # one heterozygote anywhere is enough for 'include'
  inbred[2, "i3"] <- 0
  expect_identical(filter_geno(inbred, hets = "include", verbose = FALSE)$snp,
                   "m2")
})

test_that("a bare matrix is individuals-by-markers and keeps its orientation", {
  set.seed(7)
  M <- matrix(sample(c(-1, 0, 1), 20 * 9, replace = TRUE), nrow = 20,
              dimnames = list(paste0("ind", 1:20), paste0("mk", 1:9)))
  M[, 3] <- 1                                   # monomorphic markers (columns)
  M[, 7] <- -1
  out <- filter_geno(M, verbose = FALSE)
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(20L, 7L))
  expect_identical(colnames(out), paste0("mk", c(1, 2, 4, 5, 6, 8, 9)))
  expect_identical(rownames(out), rownames(M))
  expect_equal(out, M[, -c(3, 7)])
  # the data-frame form of the same panel (markers x individuals) agrees
  df <- mk_num_a(t(M))
  expect_identical(filter_geno(df, verbose = FALSE)$snp,
                   paste0("m", c(1, 2, 4, 5, 6, 8, 9)))
})

# ---------------------------------------------------------------------------
# The PLINK-port internals
# ---------------------------------------------------------------------------

test_that(".plink_hap_rsq() is symmetric and in [0, 1] on every pair of length-5 genotype vectors", {
  skip_on_cran()
  g <- as.matrix(expand.grid(rep(list(0:2), 5L)))          # 243 vectors
  hr <- simplePHENOTYPES:::.plink_hap_rsq
  gl <- lapply(seq_len(nrow(g)), function(i) as.integer(g[i, ]))
  worst_low <- Inf; worst_high <- -Inf; worst_asym <- 0; nonfinite <- 0L
  for (i in seq_along(gl)) for (j in i:length(gl)) {
    a <- hr(gl[[i]], gl[[j]])
    b <- hr(gl[[j]], gl[[i]])
    if (!is.finite(a)) nonfinite <- nonfinite + 1L
    else {
      worst_low <- min(worst_low, a); worst_high <- max(worst_high, a)
      worst_asym <- max(worst_asym, abs(a - b))
    }
  }
  expect_identical(nonfinite, 0L)
  expect_gte(worst_low, -1e-12)
  expect_lte(worst_high, 1 + 1e-9)
  expect_lt(worst_asym, 1e-9)
})

test_that(".plink_hap_rsq() hand tables: identical, complementary and independent vectors", {
  hr <- simplePHENOTYPES:::.plink_hap_rsq
  d <- as.integer(c(0, 0, 2, 2, 1, 1, 0, 2))
  expect_equal(hr(d, d), 1)                    # complete LD, incl. double hets
  expect_equal(hr(d, 2L - d), 1)               # relabelling one locus keeps r^2
  expect_equal(hr(as.integer(2 - d), as.integer(2 - d)), 1)
  # a full 3x3 table with equal cells has D = 0
  expect_lt(hr(rep(0:2, each = 3L), rep(0:2, times = 3L)), 1e-12)
  # a monomorphic locus returns 0 (PLINK skips the pair); an all-heterozygote
  # locus carries no linkage information either (D = 0 up to rounding)
  expect_identical(hr(rep(0L, 6), as.integer(c(0, 1, 2, 0, 1, 2))), 0)
  expect_lt(hr(rep(1L, 6), as.integer(c(0, 1, 2, 0, 1, 2))), 1e-12)
})

test_that(".plink_hap_rsq() equals the textbook haplotypic r^2 when no individual is a double heterozygote", {
  # With no double heterozygote the four haplotype counts are known exactly, so
  # r^2 = D^2 / (p1 q1 p2 q2) with D = f11 f22 - f12 f21 (Hill & Robertson 1968).
  hr <- simplePHENOTYPES:::.plink_hap_rsq
  hap_r2 <- function(a, b) {
    H <- matrix(0, 2, 2)                         # allele index 1 = dosage-0 allele
    for (i in seq_along(a)) {
      ha <- if (a[i] == 1L) c(1L, 2L) else rep(if (a[i] == 0L) 1L else 2L, 2L)
      hb <- if (b[i] == 1L) c(1L, 2L) else rep(if (b[i] == 0L) 1L else 2L, 2L)
      for (k in 1:2) H[ha[k], hb[k]] <- H[ha[k], hb[k]] + 1
    }
    f <- H / sum(H)
    p <- sum(f[1, ]); q <- sum(f[, 1])
    if (p %in% c(0, 1) || q %in% c(0, 1)) return(0)
    (f[1, 1] * f[2, 2] - f[1, 2] * f[2, 1])^2 / (p * (1 - p) * q * (1 - q))
  }
  set.seed(2)
  worst <- 0; used <- 0L
  for (it in 1:1500) {
    n <- sample(8:40, 1)
    a <- sample(0:2, n, TRUE); b <- sample(0:2, n, TRUE)
    if (any(a == 1L & b == 1L)) next
    used <- used + 1L
    worst <- max(worst, abs(hap_r2(a, b) - hr(as.integer(a), as.integer(b))))
  }
  expect_gt(used, 100L)
  expect_lt(worst, 1e-12)
})

test_that(".plink_calc_lnlike() is the explicit two-locus haplotype log-likelihood", {
  ll <- simplePHENOTYPES:::.plink_calc_lnlike
  k <- c(10, 3, 2, 7); cc <- 4
  ttr <- 1 / (sum(k) + 2 * cc)
  f <- k * ttr; hhs <- cc * ttr
  for (incr in c(0, 0.05, 0.5 * hhs, hhs)) {
    F <- c(f[1] + incr, f[2] + hhs - incr, f[3] + hhs - incr, f[4] + incr)
    expected <- cc * log(F[1] * F[4] + F[2] * F[3]) + sum(k * log(F))
    expect_equal(ll(k[1], k[2], k[3], k[4], cc, f[1], f[2], f[3], f[4], hhs, incr),
                 expected, tolerance = 1e-12)
  }
  # zero haplotype counts are skipped, not log(0)
  expect_true(is.finite(ll(0, 3, 2, 0, 4, 0, 3 * ttr, 2 * ttr, 0, hhs, 0.01)))
})

test_that(".plink_em_hethet() returns the maximum-likelihood phase split", {
  em <- simplePHENOTYPES:::.plink_em_hethet
  set.seed(5)
  worst <- 0; checked <- 0L; out_of_range <- 0L
  for (it in 1:400) {
    k <- rpois(4, sample(c(0.5, 2, 6), 1)); cc <- rpois(1, 4)
    if (cc == 0) next
    r <- em(k[1], k[2], k[3], k[4], cc)
    if (isTRUE(r$mono)) next
    ttr <- 1 / (sum(k) + 2 * cc); f <- k * ttr; hhs <- cc * ttr
    ll <- function(incr) {
      F <- c(f[1] + incr, f[2] + hhs - incr, f[3] + hhs - incr, f[4] + incr)
      F[abs(F) < 1e-15] <- 0
      q <- F[1] * F[4] + F[2] * F[3]
      if (q <= 0 || any(k > 0 & F <= 0)) return(-Inf)
      cc * log(q) + sum(ifelse(k > 0, k * log(pmax(F, 1e-300)), 0))
    }
    best <- r$f11 - f[1]                         # the chosen coupling increment
    if (best < -1e-12 || best > hhs + 1e-12) out_of_range <- out_of_range + 1L
    grid_max <- max(vapply(seq(0, hhs, length.out = 2001L), ll, 0))
    worst <- max(worst, grid_max - ll(best))
    checked <- checked + 1L
  }
  expect_gt(checked, 300L)
  expect_identical(out_of_range, 0L)
  # the exact cubic solution is never worse than a 2001-point grid search
  expect_lt(worst, 1e-9)
  # no ambiguity (no double heterozygote): the frequencies are the raw shares
  r0 <- em(6, 2, 1, 3, 0)
  expect_equal(r0$f11, 6 / 12)
  expect_equal(c(r0$f1x, r0$fx1), c(8 / 12, 7 / 12))
  # nothing observed at all is flagged monomorphic
  expect_true(em(0, 0, 0, 0, 0)$mono)
})

test_that(".plink_blocks_classify() gives PLINK's class for fixed two-locus tables", {
  cl <- function(c9) simplePHENOTYPES:::.plink_blocks_classify(c9, 82L, 52L)
  # counts9 is the 3x3 genotype table, index = genotype_a * 3 + genotype_b + 1
  # monomorphic at either locus, or nothing called: class 1 (null)
  expect_identical(cl(c(0, 0, 0, 0, 0, 0, 50, 50, 50)), 1L)
  expect_identical(cl(c(10, 0, 0, 10, 0, 0, 10, 0, 0)), 1L)
  expect_identical(cl(integer(9)), 1L)
  # complete LD with 10 or more samples per homozygote class is strong LD
  strong <- function(n) { c9 <- integer(9); c9[1] <- n; c9[9] <- n; c9 }
  for (n in c(10L, 20L, 50L, 200L)) expect_identical(cl(strong(n)), 6L)
  # the evidence builds up: once 'strong' it never reverts as the sample grows
  codes <- vapply(1:15, function(n) cl(strong(n)), 1L)
  expect_lt(codes[1L], 5L)
  expect_true(all(diff(as.integer(codes >= 5L)) >= 0L))
  # D = 0 with plenty of data: strong recombination (class 0)
  indep4 <- integer(9); indep4[c(1, 3, 7, 9)] <- 25
  expect_identical(cl(indep4), 0L)
  for (k in c(1, 5, 20)) expect_identical(cl(rep(k, 9)), 0L)
  # a single double heterozygote carries no information about phase: not strong
  expect_lt(cl(c(0, 0, 0, 0, 10, 0, 0, 0, 0)), 5L)
})

test_that(".plink_blocks_classify() is invariant to relabelling alleles and to swapping loci", {
  cl <- function(c9) simplePHENOTYPES:::.plink_blocks_classify(c9, 82L, 52L)
  relabel <- function(c9, a = FALSE, b = FALSE) {
    m <- matrix(c9, 3L, 3L, byrow = TRUE)       # rows: genotype a, cols: genotype b
    if (a) m <- m[3:1, ]
    if (b) m <- m[, 3:1]
    as.integer(t(m))
  }
  set.seed(3)
  mismatches <- 0L
  for (it in 1:150) {
    c9 <- as.integer(rmultinom(1L, sample(10:80, 1L), stats::runif(9)))
    ref <- cl(c9)
    for (a in c(FALSE, TRUE)) for (b in c(FALSE, TRUE))
      if (cl(relabel(c9, a, b)) != ref) mismatches <- mismatches + 1L
    if (cl(as.integer(t(matrix(c9, 3L, 3L, byrow = TRUE)))) != ref)
      mismatches <- mismatches + 1L
  }
  expect_identical(mismatches, 0L)
})

# ---------------------------------------------------------------------------
# PLINK parity bookkeeping (Codex row 14)
# ---------------------------------------------------------------------------

test_that("the golden PLINK results are ordered subsets of the panel and partition it", {
  fx_path <- system.file("extdata", "plink_parity", "plink19_prune.rds",
                         package = "simplePHENOTYPES")
  skip_if_not(nzchar(fx_path) && file.exists(fx_path),
              "PLINK parity fixture not installed")
  fx <- readRDS(fx_path)
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  g <- SNP55K_maize282_maf04
  expect_true(is.character(fx$plink_version) && nzchar(fx$plink_version))
  for (key in names(fx$kept)) {
    kept <- fx$kept[[key]]
    expect_false(anyDuplicated(kept) > 0L, info = key)
    expect_true(all(kept %in% g$snp), info = key)
    # PLINK writes .prune.in in map order
    expect_identical(kept, g$snp[g$snp %in% kept], info = key)
    removed <- setdiff(g$snp, kept)
    expect_identical(length(kept) + length(removed), nrow(g), info = key)
    expect_length(intersect(kept, removed), 0L)
  }
})

test_that("filter_geno() kept + removed partition the panel exactly once", {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  g <- SNP55K_maize282_maf04[SNP55K_maize282_maf04$chr %in% 1:2, ]
  kept <- filter_geno(g, indep_pairwise = c(50, 5, 0.2),
                      remove_monomorphic = FALSE, verbose = FALSE)$snp
  removed <- setdiff(g$snp, kept)
  expect_false(anyDuplicated(kept) > 0L)
  expect_identical(length(kept) + length(removed), nrow(g))
  expect_identical(kept, g$snp[g$snp %in% kept])         # input order kept
  expect_gt(length(removed), 0L)
})

# ---------------------------------------------------------------------------
# qtn_table(): closed-form marginal shares (Codex row 16)
# ---------------------------------------------------------------------------

qtn_toy <- function() {
  # 45 individuals; m3 and m4 are exactly uncorrelated (a 3x3 factorial
  # replicated 5 times), m1 and m2 are identical (perfect LD)
  gr <- expand.grid(a = -1:1, b = -1:1)[rep(1:9, 5), ]
  set.seed(1)
  M <- cbind(m1 = sample(c(-1, 0, 1), 45, TRUE), m2 = NA,
             m3 = gr$a, m4 = gr$b)
  M[, "m2"] <- M[, "m1"]
  G <- data.frame(snp = colnames(M), allele = "A/G", chr = 1L,
                  pos = 100L * seq_len(ncol(M)), cm = 0, t(M),
                  check.names = FALSE)
  names(G)[-(1:5)] <- paste0("i", 1:45)
  G
}

test_that("per-QTN shares are k^2 var(e_j G_j) / V_P: uncorrelated loci sum to prop", {
  G <- qtn_toy()
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 1) |>
    additive(prop = 0.5, qtn = c("m3", "m4"), effect = c(1, 2))
  ly <- ph$layers[[1L]]
  qe <- simplePHENOTYPES:::.layer_qtn_effect(ly, 1L, 1L)
  # both loci have the same genotype variance, so the shares split prop as
  # e^2 / sum(e^2) = 1/5 and 4/5, and they sum to prop exactly (V_P = 1 here)
  sh <- simplePHENOTYPES:::.qtn_var(ph, ly, 1L, qe$qtn, qe$effect, 1)
  expect_equal(sh, c(0.1, 0.4), tolerance = 1e-12)
  expect_equal(sum(sh), 0.5, tolerance = 1e-12)
  # the reported table is the same share on the realized-variance scale
  qt <- qtn_table(ph)
  expect_equal(qt$var_explained * stats::var(ph$pheno$value), c(0.1, 0.4),
               tolerance = 1e-12)
})

test_that("per-QTN shares do not sum to prop under linkage disequilibrium", {
  G <- qtn_toy()
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 1) |>
    additive(prop = 0.5, qtn = c("m1", "m2"), effect = c(1, 1))
  ly <- ph$layers[[1L]]
  qe <- simplePHENOTYPES:::.layer_qtn_effect(ly, 1L, 1L)
  sh <- simplePHENOTYPES:::.qtn_var(ph, ly, 1L, qe$qtn, qe$effect, 1)
  # raw = 2 G, so each locus has var(G) / var(2 G) = 1/4 of prop; the
  # covariance (the other half) is attributed to neither locus
  expect_equal(sh, c(0.125, 0.125), tolerance = 1e-12)
  expect_equal(sum(sh), 0.25, tolerance = 1e-12)
  expect_lt(sum(sh), 0.5)
  qt <- qtn_table(ph)
  expect_true(all(is.finite(qt$var_explained)))
  # opposite effects on identical loci cancel: the layer has no variance and
  # every share is 0, never NaN (additive() itself refuses such a layer, so the
  # effects are handed to .qtn_var() directly)
  expect_identical(
    simplePHENOTYPES:::.qtn_var(ph, ly, 1L, qe$qtn, c(1, -1), 1), c(0, 0))
})

test_that(".qtn_var() returns zeros, not NaN, for a zero realized V_P or a zero-prop layer", {
  G <- qtn_toy()
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 1) |>
    additive(prop = 0.5, qtn = c("m3", "m4"), effect = c(1, 2))
  ly <- ph$layers[[1L]]
  qe <- simplePHENOTYPES:::.layer_qtn_effect(ly, 1L, 1L)
  for (vp in c(0, NA_real_, Inf)) {
    expect_identical(
      simplePHENOTYPES:::.qtn_var(ph, ly, 1L, qe$qtn, qe$effect, vp), c(0, 0))
  }
  ly0 <- ly; ly0$prop <- 0
  expect_identical(
    simplePHENOTYPES:::.qtn_var(ph, ly0, 1L, qe$qtn, qe$effect, 1), c(0, 0))
})

test_that(".tx_qtn_var() gives a constant gene a zero share and never NaN", {
  skip_on_cran()
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  G <- SNP55K_maize282_maf04
  tx <- simulate_transcriptome(G, n_genes = 60, seed = 1)
  gid <- rownames(tx$expression)[c(2, 5, 9)]
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 5, transcriptome = tx) |>
    transcriptome(prop = 0.4, genes = gid, slopes = c(1, -2, 0.5))
  ly <- ph$layers[[1L]]
  var_p <- stats::var(ph$pheno$value)
  base <- simplePHENOTYPES:::.tx_qtn_var(ph, ly, 1L, 1L, var_p)
  expect_equal(unname(base / base[2L]), c(0.25, 1, 0.0625), tolerance = 1e-9)
  expect_true(all(is.finite(base)))
  # a zero V_P gives all zeros
  expect_identical(unname(simplePHENOTYPES:::.tx_qtn_var(ph, ly, 1L, 1L, 0)),
                   c(0, 0, 0))
  # make the first causal gene constant: it contributes nothing
  ph2 <- ph
  ph2$expression[gid[1L], ] <- 3
  v2 <- simplePHENOTYPES:::.tx_qtn_var(ph2, ly, 1L, 1L, var_p)
  expect_identical(unname(v2[1L]), 0)
  expect_true(all(is.finite(v2)) && all(v2[-1L] > 0))
})

# ---------------------------------------------------------------------------
# write_phenotypes(): wide layout, replications, NA, separators (Codex row 17)
# ---------------------------------------------------------------------------

test_that("write_phenotypes() round trips a wide multi-trait, multi-rep table with NA", {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, n_reps = 2,
                           h2 = 0.5, seed = 12) |> additive(n_qtn = 4)
  ph$pheno$value[c(3L, 40L)] <- NA_real_         # force a known NA pattern
  wide <- phenotypes_wide(ph)
  for (sep in c("\t", ",")) {
    f <- tempfile(fileext = ".txt")
    write_phenotypes(ph, file = f, format = "wide", sep = sep)
    back <- data.table::fread(f, sep = sep, data.table = FALSE)
    expect_identical(names(back), names(wide))
    expect_identical(as.character(back$id), as.character(wide$id))
    expect_identical(back$rep, wide$rep)
    for (tr in c("Trait_1", "Trait_2")) {
      expect_identical(is.na(back[[tr]]), is.na(wide[[tr]]), info = tr)
      expect_lt(max(abs(back[[tr]] - wide[[tr]]), na.rm = TRUE),
                1e-12 * max(1, max(abs(wide[[tr]]), na.rm = TRUE)))
    }
  }
  expect_identical(sum(is.na(wide$Trait_1)) + sum(is.na(wide$Trait_2)) > 0L, TRUE)
  # writing again to the same path replaces the content (no append)
  f <- tempfile(fileext = ".txt")
  write_phenotypes(ph, file = f, format = "wide")
  n1 <- length(readLines(f))
  write_phenotypes(ph, file = f, format = "wide")
  expect_identical(length(readLines(f)), n1)
  expect_identical(n1, nrow(wide) + 1L)
})

# ---------------------------------------------------------------------------
# Documentation contracts (skipped when the sources are not shipped)
# ---------------------------------------------------------------------------

test_that("the filter_geno() parity claim does not cite randomized differential testing as evidence", {
  skip_if_no_source("R", "qc_filter_geno.R")
  txt <- paste(readLines(testthat::test_path("..", "..", "R", "qc_filter_geno.R"),
                         warn = FALSE), collapse = "\n")
  expect_false(grepl("plus randomized differential testing", txt, fixed = TRUE))
  # the limitation is stated instead
  expect_true(grepl("randomly generated genotypes", txt, fixed = TRUE))
})

test_that("the canonical API list calls format_conversion() internal", {
  skip_if_no_source("AGENTS.md")
  ag <- readLines(testthat::test_path("..", "..", "AGENTS.md"), warn = FALSE)
  line <- grep("format_conversion", ag, value = TRUE)
  expect_true(length(line) >= 1L)
  expect_true(all(grepl("internal", line, fixed = TRUE)))
})

test_that("the PLINK fixture generator refuses a panel with missing calls", {
  skip_if_no_source("data-raw", "plink_parity_fixtures.R")
  src <- readLines(testthat::test_path("..", "..", "data-raw",
                                       "plink_parity_fixtures.R"), warn = FALSE)
  guard <- grep("stopifnot\\(!anyNA\\(Dm\\)\\)", src)
  writer <- grep("\\.tped", src)[1L]
  expect_length(guard, 1L)
  expect_lt(guard, writer)               # the guard runs before the first write
})
