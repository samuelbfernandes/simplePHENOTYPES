# Byte-exact parity of filter_geno()'s LD pruning with PLINK 1.9.
#
# The golden kept-marker sets were produced by running PLINK 1.9's
# --indep-pairwise / --indep on SNP55K_maize282_maf04 (see
# data-raw/plink_parity_fixtures.R) and committed to inst/extdata, so this test
# needs no PLINK binary. filter_geno() must reproduce each set marker-for-marker.

fixture_path <- system.file("extdata", "plink_parity", "plink19_prune.rds",
                            package = "simplePHENOTYPES")

test_that("the PLINK parity fixture is present", {
  expect_true(nzchar(fixture_path) && file.exists(fixture_path))
})

# Each verified-exact config: key in the fixture -> filter_geno() call.
run_config <- list(
  pairwise_50_5_0.2 = function(g)
    filter_geno(g, indep_pairwise = c(50, 5, 0.2),
                remove_monomorphic = FALSE, verbose = FALSE),
  pairwise_100_10_0.1 = function(g)
    filter_geno(g, indep_pairwise = c(100, 10, 0.1),
                remove_monomorphic = FALSE, verbose = FALSE),
  pairwise_20_2_0.5 = function(g)
    filter_geno(g, indep_pairwise = c(20, 2, 0.5),
                remove_monomorphic = FALSE, verbose = FALSE),
  vif_50_5_2 = function(g)
    filter_geno(g, indep = c(50, 5, 2),
                remove_monomorphic = FALSE, verbose = FALSE),
  vif_100_10_5 = function(g)
    filter_geno(g, indep = c(100, 10, 5),
                remove_monomorphic = FALSE, verbose = FALSE),
  pairwise_kb_250_5_0.2 = function(g)
    filter_geno(g, indep_pairwise = c(250, 5, 0.2), window_unit = "kb",
                remove_monomorphic = FALSE, verbose = FALSE),
  pairphase_50_5_0.2 = function(g)
    filter_geno(g, indep_pairphase = c(50, 5, 0.2),
                remove_monomorphic = FALSE, verbose = FALSE),
  pairphase_100_10_0.1 = function(g)
    filter_geno(g, indep_pairphase = c(100, 10, 0.1),
                remove_monomorphic = FALSE, verbose = FALSE),
  pairphase_20_2_0.5 = function(g)
    filter_geno(g, indep_pairphase = c(20, 2, 0.5),
                remove_monomorphic = FALSE, verbose = FALSE)
)

test_that("filter_geno() LD pruning matches PLINK 1.9 marker-for-marker", {
  skip_if_not(nzchar(fixture_path) && file.exists(fixture_path),
              "PLINK parity fixture not installed")
  fixture <- readRDS(fixture_path)
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  g <- SNP55K_maize282_maf04

  for (key in names(run_config)) {
    kept <- run_config[[key]](g)$snp
    golden <- fixture$kept[[key]]
    expect_identical(sort(kept), sort(golden),
                     info = paste("config", key, "diverged from PLINK 1.9"))
  }
})

# Small cases PLINK 1.9 was executed on directly (results hard-coded, no PLINK
# needed at test time). These pin the fixes for floating-point threshold
# boundaries and monomorphic handling -- the two ways an earlier stats::cor() /
# solve() implementation diverged from PLINK.
mk_df <- function(counts, pos = NULL) {
  m <- nrow(counts); ni <- ncol(counts)
  if (is.null(pos)) pos <- seq_len(m) * 1000L
  data.frame(snp = paste0("m", seq_len(m)), allele = "A/G", chr = 1L,
             pos = pos, cm = 0,
             stats::setNames(as.data.frame(counts - 1L), paste0("i", seq_len(ni))),
             check.names = FALSE, stringsAsFactors = FALSE)
}

test_that("r^2 / VIF at a floating-point threshold boundary match PLINK exactly", {
  # r^2(m1, m2) is exactly 49/64 = 0.765625; a threshold a hair below it must
  # prune the lower-MAF marker (m1) as PLINK does -- stats::cor()'s last-bit
  # drift used to keep both.
  d <- mk_df(rbind(c(0, 2, 2, 1, 2), c(2, 1, 1, 1, 1)))
  pw <- filter_geno(d, indep_pairwise = c(2, 1, 0.76562499999995592),
                    remove_monomorphic = FALSE, verbose = FALSE)$snp
  expect_identical(pw, "m2")
  vf <- filter_geno(d, indep = c(2, 1, 4.2666666666666559),
                    remove_monomorphic = FALSE, verbose = FALSE)$snp
  expect_identical(vf, "m2")
})

test_that("a monomorphic marker is pruned by LD even without remove_monomorphic", {
  d <- mk_df(rbind(c(0, 0, 0, 0, 0, 0), c(0, 0, 1, 1, 2, 2)))
  expect_identical(
    filter_geno(d, indep_pairwise = c(2, 1, 0.2),
                remove_monomorphic = FALSE, verbose = FALSE)$snp, "m2")
  expect_identical(
    filter_geno(d, indep = c(2, 1, 2),
                remove_monomorphic = FALSE, verbose = FALSE)$snp, "m2")
})

test_that("a marker skipped by step > window is not pruned as monomorphic", {
  # window 2, step 3: m3 (monomorphic, index 3) is never loaded into a window, so
  # PLINK keeps it -- monomorphic pruning happens at load, not globally.
  d <- mk_df(rbind(c(0, 0, 1, 1, 2, 2), c(0, 1, 0, 2, 1, 2),
                   c(0, 0, 0, 0, 0, 0), c(0, 1, 2, 0, 1, 2)))
  keep <- filter_geno(d, indep_pairwise = c(2, 3, 0.9),
                      remove_monomorphic = FALSE, verbose = FALSE)$snp
  expect_identical(sort(keep), c("m1", "m2", "m3", "m4"))
})

test_that("Gabriel --blocks definitions match PLINK 1.9 block-for-block", {
  skip_if_not(nzchar(fixture_path) && file.exists(fixture_path),
              "PLINK parity fixture not installed")
  fixture <- readRDS(fixture_path)
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  d <- SNP55K_maize282_maf04
  dose <- as.matrix(d[, -(1:5)]) + 1L
  p <- rowSums(dose) / (2 * ncol(dose))
  eps <- 0.00000000000005684341886080801486968994140625   # SMALL_EPSILON
  maf_ok <- pmin(p, 1 - p) >= 0.05 * (1 - eps)             # Haploview floor

  block_sets <- function(max_kb) {
    out <- character(0)
    for (k in unique(d$chr)) {
      idx <- which(maf_ok & d$chr == k)
      idx <- idx[order(d$pos[idx])]
      if (length(idx) < 2L) next
      bl <- simplePHENOTYPES:::.plink_blocks_chrom(
        dose[idx, , drop = FALSE], d$pos[idx], as.integer(max_kb * 1000))
      for (b in bl) out <- c(out, paste(sort(d$snp[idx][b[1]:b[2]]),
                                        collapse = "|"))
    }
    sort(out)
  }
  golden_sets <- function(bl) sort(vapply(bl, function(x) paste(sort(x),
                                                               collapse = "|"), ""))
  expect_identical(block_sets(200), golden_sets(fixture$blocks$blocks_200kb))
  expect_identical(block_sets(500), golden_sets(fixture$blocks$blocks_500kb))
})

test_that("block_max_kb -> bp uses PLINK's epsilon-nudged truncation", {
  # Two identical markers 1001 bp apart. PLINK converts 1.001 kb to 1001 bp
  # ((int)(1000 * 1.001 * (1 + SMALL_EPSILON))), so they form one block; a plain
  # 1000.999... -> 1000 truncation would miss it and keep both.
  cnt <- rbind(c(rep(0, 100), rep(2, 100)), c(rep(0, 100), rep(2, 100)))
  d <- data.frame(snp = c("m1", "m2"), allele = "A/G", chr = 1L,
                  pos = c(1L, 1002L), cm = 0,
                  stats::setNames(as.data.frame(cnt - 1L), paste0("i", 1:200)),
                  check.names = FALSE)
  kept <- filter_geno(d, blocks = TRUE, block_max_kb = 1.001,
                      remove_monomorphic = FALSE, verbose = FALSE)$snp
  expect_length(kept, 1L)                          # block found -> collapsed to a tag
  # A huge span must cap (PLINK's 2^31 - 2), not overflow the integer cast to NA.
  k2 <- suppressMessages(
    filter_geno(d, blocks = TRUE, block_max_kb = 2147484,
                remove_monomorphic = FALSE, verbose = FALSE)$snp)
  expect_length(k2, 1L)
})

test_that("block_max_kb must be a single positive number", {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  d <- SNP55K_maize282_maf04[SNP55K_maize282_maf04$chr == 1L, ][1:20, ]
  for (bad in list(numeric(0), c(0, 1.001), -1, 0, NA_real_, "x")) {
    expect_error(
      filter_geno(d, blocks = TRUE, block_max_kb = bad, verbose = FALSE),
      "single positive number")
  }
})

test_that("filter_geno(blocks = TRUE) keeps one tag per block", {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  d <- SNP55K_maize282_maf04[SNP55K_maize282_maf04$chr == 1L, ][1:200, ]
  kept <- filter_geno(d, blocks = TRUE, block_max_kb = 200, verbose = FALSE)
  expect_lt(nrow(kept), nrow(d))                  # blocks collapsed to tags
  expect_true(all(kept$snp %in% d$snp))
})

test_that("pairphase runs on missing calls without error", {
  # NA is handled per complete pair (not PLINK's missing-data path, but must not
  # crash). Regression for a mono-computation NA that errored.
  d <- data.frame(snp = c("m1", "m2"), allele = "A/G", chr = 1L,
                  pos = c(1000L, 2000L), cm = 0,
                  i1 = c(-1, -1), i2 = c(0, 0), i3 = c(1, 1),
                  i4 = c(NA, 0), i5 = c(-1, NA), i6 = c(1, 1),
                  check.names = FALSE)
  expect_silent(
    keep <- filter_geno(d, indep_pairphase = c(2, 1, 0.2),
                        remove_monomorphic = FALSE, verbose = FALSE)$snp)
  expect_true(all(keep %in% c("m1", "m2")))
})

test_that("pruning respects a prior MAF filter (still valid pruning on subset)", {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  g <- SNP55K_maize282_maf04
  # A pruned + MAF-filtered run keeps a subset of the MAF-filtered markers, and
  # never a marker the MAF filter removed.
  maf_only <- filter_geno(g, maf_above = 0.2, verbose = FALSE)$snp
  both <- filter_geno(g, maf_above = 0.2, indep_pairwise = c(50, 5, 0.2),
                      verbose = FALSE)$snp
  expect_true(all(both %in% maf_only))
  expect_lt(length(both), length(maf_only))
})
