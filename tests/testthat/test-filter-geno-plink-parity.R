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
