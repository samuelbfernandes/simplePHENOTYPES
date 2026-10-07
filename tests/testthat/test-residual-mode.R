# Fixed mode preserves the existing stream; random mode preserves normal sampling.
test_that("fixed default preserves the residual algorithm and RNG state", {
  set.seed(33)
  rng <- .Random.seed
  actual <- .seeded_residual(73, 40, 0.6)
  expect_identical(.Random.seed, rng)
  set.seed(73)
  z <- rnorm(40, sd = sqrt(0.6))
  expected <- (z - mean(z)) / sd(z) * sqrt(0.6)
  expect_identical(actual, expected)
  expect_identical(actual, .seeded_residual(73, 40, 0.6, "fixed"))
  expect_identical(.seeded_residual(73, 40, 0.6, "random"), z)
  expect_identical(.seeded_residual(73, 40, 0, "random"), rep(0, 40))
})

test_that("random draws retain the sampling distribution of mean and variance", {
  draws <- vapply(1:2500, function(i) .seeded_residual(i, 20, 1, "random"), numeric(20))
  expect_equal(mean(colMeans(draws)), 0, tolerance = 0.015)
  expect_lt(abs(var(colMeans(draws)) - 1 / 20), 0.006)
  sample_vars <- apply(draws, 2, var)
  expect_equal(mean(sample_vars), 1, tolerance = 0.025)
  expect_lt(abs(var(sample_vars) - 2 / 19), 0.015)
})

test_that("mode changes residuals but not genetic layers, replication, or seed state", {
  data("SNP55K_maize282_maf04")
  geno <- SNP55K_maize282_maf04[1:150, 1:85]
  set.seed(123)
  rng <- .Random.seed
  fixed <- simulate_phenotype(geno, h2 = 0.4, n_qtn = 8, seed = 8, n_reps = 2)
  explicit <- simulate_phenotype(geno, h2 = 0.4, n_qtn = 8, seed = 8,
                                 n_reps = 2, residual_mode = "fixed")
  random <- simulate_phenotype(geno, h2 = 0.4, n_qtn = 8, seed = 8,
                               n_reps = 2, residual_mode = "random")
  expect_identical(fixed$pheno, explicit$pheno)
  expect_identical(genetic_values(fixed), genetic_values(random))
  expect_identical(qtn_table(fixed)[, c("snp", "effect")],
                   qtn_table(random)[, c("snp", "effect")])
  expect_identical(.Random.seed, rng)
  for (r in 1:2) {
    y <- random$pheno$value[random$pheno$rep == r]
    g <- genetic_values(random)[, 1]
    e <- .seeded_residual(.layer_seed(8, "residual_t1", r - 1), 80, 0.6, "random")
    expect_equal(unname(y - g), e, tolerance = 1e-14)
  }
  means <- simulate_phenotype(geno, h2 = 0.4, n_qtn = 8, seed = 8,
                              residual_mode = "random", reps = 4)
  g <- genetic_values(means)[, 1]
  expect_equal(means$pheno$value - g, (random$pheno$value[1:80] - g) / 2)
  expect_error(simulate_phenotype(geno, residual_mode = "unknown"), "arg")
})

test_that("random correlated errors have the requested covariance without normalization", {
  sim <- list(seed = 44, n_ind = 30000,
              resid_cor = matrix(c(1, 0.6, 0.6, 1), 2), residual_mode = "random")
  e <- .correlated_residuals(sim, 1, c(0.4, 0.8))
  target <- diag(sqrt(c(0.4, 0.8))) %*% sim$resid_cor %*% diag(sqrt(c(0.4, 0.8)))
  expect_equal(unname(cov(e)), target, tolerance = 0.025)
  expected_first <- .seeded_residual(.layer_seed(44, "residual_t1", 0), 30000, 1,
                                    "random") * sqrt(0.4)
  expect_identical(e[, 1], expected_first)
  sim$resid_cor[,] <- 1
  singular <- .correlated_residuals(sim, 1, c(0.4, 0.8))
  expect_equal(singular[, 1] / sqrt(0.4), singular[, 2] / sqrt(0.8))
  zero <- .correlated_residuals(sim, 1, c(0, 0.8))
  expect_identical(zero[, 1], rep(0, 30000))
})

test_that("complex and carried traits respect residual mode", {
  data("SNP55K_maize282_maf04")
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:40)
  s <- simulate_phenotype(pop, h2 = 0.5, n_qtn = 5, seed = 33,
                          residual_mode = "random")
  combined <- complex_phenotypes(s, s, h2 = 0.5, residual_mode = "random")
  e <- combined$pheno$value - genetic_values(combined)[, 1]
  expected <- .seeded_residual(.layer_seed(combined$seed, "complex_resid_t1", 0),
                               40, 0.5, "random")
  expect_equal(e, expected, ignore_attr = TRUE)
  top <- select_ind(s, prop = 0.5)
  inherited <- simulate_phenotype(top, seed = 54)
  expect_identical(inherited$residual_mode, "random")
  e <- inherited$pheno$value - genetic_values(inherited)[, 1]
  expect_equal(e, .seeded_residual(.layer_seed(54, "residual_t1", 0), 20,
                                  inherited$trait$var_e, "random"), ignore_attr = TRUE)
  override <- simulate_phenotype(top, seed = 54, residual_mode = "fixed")
  expect_equal(var(override$pheno$value - genetic_values(override)[, 1]),
               inherited$trait$var_e)
  expect_identical(genetic_values(override), genetic_values(inherited))
  top$trait$residual_mode <- NULL
  expect_identical(simulate_phenotype(top, seed = 54)$residual_mode, "fixed")
})

test_that("random vQTL errors use conditional rather than sampled variance scaling", {
  data("SNP55K_maize282_maf04")
  s <- simulate_phenotype(SNP55K_maize282_maf04, seed = 31, residual_mode = "random") |>
    additive(prop = 0.4, n_qtn = 5) |>
    vqtl(prop = 0.2)
  layer <- s$layers[[2]]
  qe <- .layer_qtn_effect(layer, 1, 1)
  raw <- as.numeric(.geno_cols(s, qe$qtn) %*% qe$effect)
  loading <- sqrt(0.2) * (raw - mean(raw)) / sd(raw)
  weights <- exp(0.5 * (loading - max(loading)))
  conditional_sd <- weights / sqrt(mean(weights^2)) * sqrt(0.2)
  expect_equal(mean(conditional_sd^2), 0.2)
  e0 <- .seeded_residual(.layer_seed(31, "residual_t1", 0), s$n_ind, 0.4, "random")
  z <- .seeded_residual(.layer_seed(31, "vqtl_residual_t1", 0), s$n_ind, 1, "random")
  expect_equal(s$pheno$value - genetic_values(s)[, 1], e0 + conditional_sd * z,
               ignore_attr = TRUE, tolerance = 1e-14)
})
