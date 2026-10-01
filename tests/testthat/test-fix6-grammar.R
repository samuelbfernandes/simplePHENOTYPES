# Round-6 grammar fixes: G1-G4 (reps wording verified numerically, print of
# per-trait reps). No computed numbers changed; these tests pin the corrected
# statements.

.fx6 <- function() {
  utils::data("SNP55K_maize282_maf04", package = "simplePHENOTYPES",
              envir = environment())
  get("SNP55K_maize282_maf04", envir = environment())
}

test_that("G1: realized entry-mean/record H2 carry Cov(G, e); allocation is the target", {
  geno <- .fx6()
  p <- simulate_phenotype(geno, h2 = 0.5, n_qtn = 3, seed = 6, reps = 4)
  g <- genetic_values(p)[, 1]
  y <- p$pheno$value
  e_bar <- y - g                       # stored residual = e / sqrt(reps)
  e <- e_bar * sqrt(4)                 # realized single-record residual
  vg <- stats::var(g)
  # entry mean: Var(y_bar) = V_G + V_E/r + 2 Cov(G, e)/sqrt(r)
  den <- vg + stats::var(e) / 4 + 2 * stats::cov(g, e) / sqrt(4)
  expect_equal(stats::var(y), den, tolerance = 1e-12)
  real <- simplePHENOTYPES:::.realized_h2(p)
  expect_equal(real, vg / den, tolerance = 1e-12)
  # the allocation is the target, and is NOT the realized value (Codex seed 6)
  alloc <- vg / (vg + stats::var(e_bar))
  expect_equal(alloc, 0.8, tolerance = 1e-12)
  expect_equal(real, 0.9024521, tolerance = 1e-6)
  expect_gt(abs(real - alloc), 0.05)
  # record scale: Var(y) = V_G + V_E + 2 Cov(G, e)
  rec <- simplePHENOTYPES:::.realized_h2(p, scale = "record")
  expect_equal(rec, vg / (vg + stats::var(e) + 2 * stats::cov(g, e)),
               tolerance = 1e-12)
  expect_equal(rec, 0.5826881, tolerance = 1e-6)
})

test_that("G1: print() separates target from realized (covariance named)", {
  geno <- .fx6()
  p <- simulate_phenotype(geno, h2 = 0.5, n_qtn = 3, seed = 6, reps = 4)
  out <- paste(utils::capture.output(print(p)), collapse = "\n")
  expect_match(out, "Var\\(G\\) / Var\\(y_bar\\)")
  expect_match(out, "Cov\\(G, e\\)")
  expect_match(out, "target")
  expect_match(out, "Single-record realized H.* = Var\\(G\\) / Var\\(y\\)")
  # the old covariance-free claim (realized = V_G / (V_G + V_E / reps)) is gone
  expect_false(grepl("= V_G / (V_G + V_E / reps)", out, fixed = TRUE))
})

test_that("G2: a derived transcriptome component is NOT redrawn per record", {
  geno <- .fx6()
  tx <- suppressWarnings(simulate_transcriptome(geno, n_genes = 80, h2 = 0.4,
                                                 seed = 17))
  mk <- function(rp) {
    p <- simulate_phenotype(geno, seed = 29, transcriptome = tx, reps = rp)
    transcriptome(p, prop = 0.55, n_genes = 20)
  }
  t1 <- mk(1); t5 <- mk(5)
  tt1 <- simplePHENOTYPES:::.transcriptome_matrix(t1, 1)[, 1]
  tt5 <- simplePHENOTYPES:::.transcriptome_matrix(t5, 1)[, 1]
  expect_identical(tt1, tt5)                      # persistent entry-level covariate
  txe <- function(s) {
    simplePHENOTYPES:::.transcriptome_matrix(s, 1, "total")[, 1] -
      simplePHENOTYPES:::.transcriptome_matrix(s, 1, "genetic")[, 1]
  }
  expect_identical(txe(t1), txe(t5))
  expect_gt(stats::var(txe(t5)), 0.1)
  # only the phenotype residual is divided by sqrt(reps)
  res <- function(s) {
    s$pheno$value - simplePHENOTYPES:::.genetic_matrix(s, 1)[, 1] -
      simplePHENOTYPES:::.transcriptome_matrix(s, 1)[, 1] -
      simplePHENOTYPES:::.trait_mean(s, 1)
  }
  expect_equal(res(t5), res(t1) / sqrt(5), tolerance = 1e-12)
  # realized denominator is the variance of the full stored phenotype
  g5 <- simplePHENOTYPES:::.genetic_value_matrix(t5, 1)[, 1]
  expect_equal(simplePHENOTYPES:::.realized_h2(t5),
               stats::var(g5) / stats::var(t5$pheno$value), tolerance = 1e-12)
  # print says replication is conditional on the transcriptome
  out <- paste(utils::capture.output(print(t5)), collapse = "\n")
  expect_match(out, "not redrawn per record")
})

test_that("G3: vQTL realized residual variance is [V0 + Vv + 2Cov]/reps, not nominal", {
  geno <- .fx6()
  mk <- function(rp) suppressWarnings({
    b <- simulate_phenotype(geno, seed = 78, reps = rp)
    b <- additive(b, prop = 0.4, n_qtn = 3)
    vqtl(b, prop = 0.2, n_qtn = 2)
  })
  b1 <- mk(1); b4 <- mk(4)
  e1 <- b1$pheno$value - genetic_values(b1)[, 1]
  e4 <- b4$pheno$value - genetic_values(b4)[, 1]
  expect_equal(stats::var(e1), 0.6955681, tolerance = 1e-6)
  expect_equal(stats::var(e4), stats::var(e1) / 4, tolerance = 1e-12)
  expect_equal(stats::var(e4), 0.1738920, tolerance = 1e-6)
  expect_gt(abs(stats::var(e4) - 0.6 / 4), 0.02)     # not the nominal 0.15
})

test_that("G4: print() shows the per-trait reps vector faithfully", {
  geno <- .fx6()
  line <- function(rp, nt) {
    p <- simulate_phenotype(geno, n_traits = nt, reps = rp, seed = 2)
    out <- utils::capture.output(print(p))
    grep("Entry means of", out, value = TRUE)
  }
  expect_match(line(c(1, 4, 1), 3), "reps (per trait) = [1, 4, 1]", fixed = TRUE)
  expect_match(line(c(4, 1, 4), 3), "reps (per trait) = [4, 1, 4]", fixed = TRUE)
  expect_match(line(c(2, 3), 2), "reps (per trait) = [2, 3]", fixed = TRUE)
  # scalar / all-equal vector print as before
  expect_match(line(4, 1), "reps = 4 records", fixed = TRUE)
  expect_match(line(4, 3), "reps = 4 records", fixed = TRUE)
  expect_match(line(c(4, 4), 2), "reps = 4 records", fixed = TRUE)
  # nothing printed for reps = 1
  p1 <- simulate_phenotype(geno, n_traits = 2, seed = 2)
  expect_length(grep("Entry means of", utils::capture.output(print(p1))), 0L)
})
