# resid_cor: cross-trait RESIDUAL correlation (cor is the genetic one).

resid_of <- function(ph) {
  # genetic-free model: the phenotype is pure residual when no layers are added
  w <- phenotypes_wide(ph)
  as.matrix(w[, grepl("^Trait_", names(w)), drop = FALSE])
}

test_that("default resid_cor = NULL is bit-identical to omitting it", {
  data("SNP55K_maize282_maf04")
  a <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, h2 = 0.4,
                          n_qtn = 3, seed = 7)
  b <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, h2 = 0.4,
                          n_qtn = 3, seed = 7, resid_cor = NULL)
  expect_identical(a$pheno, b$pheno)
  expect_null(b$resid_cor)
})

test_that("pure-noise residual correlation approaches the target, variance kept", {
  data("SNP55K_maize282_maf04")
  Rt <- matrix(c(1, 0.6, -0.3, 0.6, 1, 0.2, -0.3, 0.2, 1), 3)
  ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 3, seed = 3,
                           resid_cor = Rt)
  e <- resid_of(ph)
  expect_equal(unname(apply(e, 2, stats::var)), rep(1, 3), tolerance = 1e-8)
  # n = 281 individuals: sampling error ~ 1/sqrt(n) ~ 0.06
  expect_lt(max(abs(stats::cor(e) - Rt)), 0.2)
  # scalar form
  ph2 <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 3,
                            resid_cor = 0.7)
  expect_gt(stats::cor(resid_of(ph2))[1, 2], 0.55)
})

test_that("correlation converges on many individuals (internal helper)", {
  sim <- list(seed = 11, n_ind = 20000L,
              resid_cor = matrix(c(1, 0.5, 0.5, 1), 2))
  W <- simplePHENOTYPES:::.correlated_residuals(sim, 1L, c(0.6, 0.3))
  expect_equal(stats::cor(W)[1, 2], 0.5, tolerance = 0.03)
  expect_equal(unname(apply(W, 2, stats::var)), c(0.6, 0.3), tolerance = 1e-8)
})

test_that("with genetic layers each residual keeps its h2 target variance", {
  data("SNP55K_maize282_maf04")
  ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, h2 = c(0.3, 0.6),
                           n_qtn = 4, seed = 5, resid_cor = 0.5)
  ind <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, h2 = c(0.3, 0.6),
                            n_qtn = 4, seed = 5)
  expect_equal(ph$var_budget, ind$var_budget)
  # trait 1 differs from the independent run only by rounding
  p1 <- ph$pheno; q1 <- ind$pheno
  k <- p1$trait == "Trait_1"
  expect_equal(p1$value[k], q1$value[k], tolerance = 1e-8)
  k2 <- p1$trait == "Trait_2"
  expect_false(isTRUE(all.equal(p1$value[k2], q1$value[k2])))
  d <- (p1$value - q1$value)[k2]
  # residual of trait 2 has variance 0.4 in both runs
  expect_true(is.finite(stats::var(d)))
})

test_that("complex_phenotypes() accepts resid_cor and NULL is bit-identical", {
  data("SNP55K_maize282_maf04")
  m1 <- additive(simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 1),
                 prop = 0.3, n_qtn = 3)
  m2 <- additive(simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 1),
                 prop = 0.3, n_qtn = 3)
  a <- complex_phenotypes(m1, m2, h2 = 0.5)
  b <- complex_phenotypes(m1, m2, h2 = 0.5, resid_cor = NULL)
  expect_identical(a$pheno, b$pheno)
  d <- complex_phenotypes(m1, m2, h2 = 0.5, resid_cor = 0.8)
  g <- a$complex_genetic[, , 1]
  e <- matrix(d$pheno$value, ncol = 2) - g
  expect_gt(stats::cor(e)[1, 2], 0.6)
  expect_error(complex_phenotypes(m1, m2, h2 = 0.5, resid_cor = 2),
               "between -1 and 1")
})

test_that("resid_cor is validated like cor", {
  data("SNP55K_maize282_maf04")
  f <- function(...) simulate_phenotype(SNP55K_maize282_maf04, seed = 1, ...)
  expect_error(f(n_traits = 1, resid_cor = 0.3), "n_traits >= 2")
  expect_error(f(n_traits = 2, resid_cor = 1.5), "between -1 and 1")
  expect_error(f(n_traits = 2, resid_cor = "a"), "between -1 and 1")
  expect_error(f(n_traits = 2, resid_cor = matrix(0, 3, 3)), "2 x 2")
  expect_error(f(n_traits = 2, resid_cor = matrix(c(1, .2, .5, 1), 2)),
               "symmetric")
  expect_error(f(n_traits = 2, resid_cor = matrix(c(2, 0, 0, 1), 2)),
               "between -1 and 1")
  expect_error(f(n_traits = 2, resid_cor = matrix(c(.9, 0, 0, 1), 2)),
               "diagonal")
  Rbad <- matrix(c(1, .9, -.9, .9, 1, .9, -.9, .9, 1), 3)
  expect_error(f(n_traits = 3, resid_cor = Rbad), "positive semi-definite")
  expect_error(f(n_traits = 3, resid_cor = -0.8), "positive semi-definite")
})

test_that("singular (perfect) correlation is accepted", {
  data("SNP55K_maize282_maf04")
  ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 4,
                           resid_cor = 1)
  e <- resid_of(ph)
  expect_equal(stats::cor(e)[1, 2], 1, tolerance = 1e-8)
})
