# observe_counts(): the optional RNA-seq count observation layer mapping the
# latent Gaussian expression of a transcriptome_sim to negative-binomial counts.

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04
tx <- simulate_transcriptome(G, n_genes = 60, seed = 1)

test_that("observe_counts draws a valid, reproducible count matrix", {
  tc <- observe_counts(tx, baseline = 4, coupling = 0.5, dispersion = 0.1, seed = 2)
  Y <- tc$counts
  expect_true(all(Y >= 0) && all(Y == round(Y)))          # non-negative integers
  expect_equal(dim(Y), c(60L, tx$n_ind))
  expect_equal(dimnames(Y), dimnames(tx$expression))
  expect_equal(
    observe_counts(tx, baseline = 4, coupling = 0.5, dispersion = 0.1, seed = 2)$counts,
    Y)                                                     # reproducible under a seed
})

test_that("the negative-binomial mean and variance match the model", {
  # coupling = 0 makes mu constant across individuals, isolating the NB
  # conditional moments: E[Y] = mu, Var[Y] = mu + phi mu^2.
  mu <- exp(4); phi <- 0.2
  tc <- observe_counts(tx, baseline = 4, coupling = 0, dispersion = phi, seed = 3)
  expect_equal(mean(tc$counts), mu, tolerance = 0.08 * mu)
  vbar <- mean(apply(tc$counts, 1L, stats::var))
  expect_equal(vbar, mu + phi * mu^2, tolerance = 0.2 * (mu + phi * mu^2))
  # Poisson limit (dispersion = 0): Var ~ mean
  tcp <- observe_counts(tx, baseline = 3, coupling = 0, dispersion = 0, seed = 4)
  expect_equal(mean(apply(tcp$counts, 1L, stats::var)) / mean(tcp$counts), 1,
               tolerance = 0.2)
})

test_that("library size scales the counts and inputs are validated", {
  a <- observe_counts(tx, baseline = 4, coupling = 0, seed = 2)
  b <- observe_counts(tx, baseline = 4, coupling = 0, library_size = 4, seed = 2)
  expect_equal(mean(b$counts) / mean(a$counts), 4, tolerance = 0.15)
  # per-gene and per-individual vectors are accepted
  tc_vec <- observe_counts(tx, baseline = rnorm(60, 4), coupling = runif(60),
                           dispersion = runif(60, 0, 0.3),
                           library_size = runif(tx$n_ind, 0.5, 2), seed = 7)
  expect_equal(dim(tc_vec$counts), c(60L, tx$n_ind))
  expect_error(observe_counts(tx, dispersion = -1), "non-negative")
  expect_error(observe_counts(tx, library_size = 0), "positive")
  expect_error(observe_counts(tx, baseline = 1:3), "length 1")
  expect_error(observe_counts(42), "transcriptome_sim")
  # an out-of-range log-mean overflows to Inf and is rejected, not returned as NA
  expect_error(observe_counts(tx, baseline = 710, coupling = 0), "overflow")
})
