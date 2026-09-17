# test-ocs.R
#
# Genomic relationship matrix (VanRaden), optimum contribution selection
# (Meuwissen 1997) and its dependency-free optimizer, and the cross-usefulness
# predictor (Zhong & Jannink / Lehermeier).

data("SNP55K_maize282_maf04")

.f2 <- function(n = 60, seed = 2) {
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:40)
  f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
  suppressMessages(selfcross(f1, n = n, seed = seed))
}
.ph <- function(f2) {
  suppressMessages(simulate_phenotype(f2, h2 = 0.5, seed = 3) |>
                     additive(n_qtn = 30))
}

test_that("g_matrix has the VanRaden scaling and id dimnames", {
  f2 <- .f2(40)
  G <- g_matrix(f2)
  expect_equal(dim(G), c(40L, 40L))
  expect_equal(rownames(G), f2$ids)
  expect_true(isSymmetric(G))
  # VanRaden G has an average diagonal near 1 (1 + mean genomic inbreeding)
  expect_lt(abs(mean(diag(G)) - 1), 0.15)
  expect_error(g_matrix(matrix(1, 2, 2)), "column names")
})

test_that("optimum_contribution gives valid, monotone frontier contributions", {
  ph <- .ph(.f2(60))
  hi <- optimum_contribution(ph, merit = "gv", lambda = 1)
  lo <- optimum_contribution(ph, merit = "gv", lambda = 50)
  expect_s3_class(hi, "ocs")
  # contributions are a valid simplex point
  expect_equal(sum(hi$contributions), 1, tolerance = 1e-6)
  expect_true(all(hi$contributions >= -1e-9))
  # more penalty -> lower coancestry and lower (or equal) gain
  expect_lt(lo$coancestry, hi$coancestry)
  expect_lte(lo$merit, hi$merit + 1e-8)
  # more penalty spreads contributions over more parents
  expect_gte(lo$n_parents, hi$n_parents)
})

test_that("optimum_contribution hits a requested coancestry target", {
  ph <- .ph(.f2(60))
  ref <- optimum_contribution(ph, merit = "gv", lambda = 1)
  tgt <- ref$coancestry * 0.8
  got <- optimum_contribution(ph, merit = "gv", target_coancestry = tgt)
  expect_lt(abs(got$coancestry - tgt), 0.02 * tgt + 1e-3)
})

test_that("optimum_contribution validates the trade-off argument", {
  ph <- .ph(.f2(40))
  expect_error(optimum_contribution(ph, merit = "gv"), "exactly one")
  expect_error(
    optimum_contribution(ph, merit = "gv", lambda = 1, target_coancestry = 0.1),
    "exactly one")
  expect_error(optimum_contribution(ph, merit = "gv", lambda = -1), ">= 0")
})

test_that("a below-minimum coancestry target warns and returns the minimum", {
  ph <- .ph(.f2(40))
  expect_warning(
    optimum_contribution(ph, merit = "gv", target_coancestry = -1),
    "minimum attainable")
})

test_that("sample_parents draws the requested number weighted by contribution", {
  f2 <- .f2(40)
  ph <- .ph(f2)
  oc <- optimum_contribution(ph, merit = "gv", lambda = 5)
  mates <- sample_parents(oc, f2, n = 12, seed = 4)
  expect_s3_class(mates, "Population")
  expect_equal(n_individuals(mates), 12L)
})

test_that("optimum_contribution defaults to the transmissible breeding value", {
  ph <- .ph(.f2(40))
  # The default merit is now "bv" (additive average effect), the transmissible
  # merit Meuwissen's OCS targets -- not the total genotypic value "gv".
  d_default <- optimum_contribution(ph, lambda = 5)
  d_bv      <- optimum_contribution(ph, merit = "bv", lambda = 5)
  expect_equal(d_default$contributions, d_bv$contributions)   # default IS bv
  # For a purely additive model the transmissible breeding value equals the total
  # genetic value (analytic average effects, exact under the F2's LD); they diverge
  # only under dominance/epistasis.
  bv_vec <- simplePHENOTYPES:::.breeding_value_matrix(ph, 1L)[, 1]
  gv_vec <- genetic_values(ph)[, 1]
  expect_gt(stats::cor(bv_vec, gv_vec), 0.999)
})

test_that("sample_parents gives unique ids when a parent is drawn repeatedly", {
  f2 <- .f2(40)
  # Concentrate all contribution on two parents so n > 2 forces repeated draws.
  contr <- stats::setNames(rep(0, n_individuals(f2)), f2$ids)
  contr[1:2] <- 0.5
  oc <- structure(
    list(contributions = contr, parents = f2$ids[1:2], n_parents = 2L,
         merit = 0, coancestry = 0, lambda = 0), class = "ocs")
  mates <- sample_parents(oc, f2, n = 12, seed = 3)
  expect_equal(n_individuals(mates), 12L)
  expect_false(anyDuplicated(mates$ids) > 0)         # disambiguated, not dropped
  expect_identical(colnames(dosages(mates)), mates$ids)
  # The documented downstream workflow (which requires unique names) now runs.
  expect_s3_class(
    suppressMessages(simulate_phenotype(mates, h2 = 0.5, seed = 1) |>
                       additive(n_qtn = 10)),
    "phenotype_sim")
})

test_that("g_matrix supports a fixed allele-frequency base", {
  # One marker, gene contents (2,2,2,0): current base p = .75, fixed base p = .1.
  d <- matrix(c(1, 1, 1, -1), nrow = 1)          # dosage -1/0/1 => gene content 2,2,2,0
  colnames(d) <- paste0("i", 1:4)
  Fcur  <- diag(g_matrix(d)) - 1                  # moving base (current p = .75)
  Fbase <- diag(g_matrix(d, base_freq = 0.1)) - 1 # fixed base p = .1
  expect_equal(unname(Fcur),  c(-1/3, -1/3, -1/3, 5), tolerance = 1e-6)
  expect_equal(unname(Fbase), c(17, 17, 17, -7/9), tolerance = 1e-3)
  expect_error(g_matrix(d, base_freq = c(0.1, 0.2)), "one allele frequency")
})

test_that("optimum_contribution rejects a non-PSD relationship matrix", {
  f2   <- .f2(40)
  pop2 <- f2[1:2]
  badG <- matrix(c(1, 2, 2, 1), 2, 2, dimnames = list(pop2$ids, pop2$ids))
  merit <- stats::setNames(c(1, 2), pop2$ids)
  expect_error(
    optimum_contribution(pop2, merit = merit, G = badG, lambda = 1),
    "positive semi-definite")
})

test_that("g_matrix enforces the -1/0/1 dosage domain", {
  d <- matrix(c(2, -1), nrow = 1, dimnames = list(NULL, c("a", "b")))
  expect_error(g_matrix(d), "coded -1/0/1")
})

test_that("optimum_contribution rejects a non-finite trade-off argument", {
  ph <- .ph(.f2(40))
  expect_error(optimum_contribution(ph, target_coancestry = Inf), "finite")
  expect_error(optimum_contribution(ph, target_coancestry = "bad"), "finite")
  expect_error(optimum_contribution(ph, lambda = Inf), "finite")
})

test_that("cross_usefulness captures between-cross mean and variance", {
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:8)
  sim <- suppressMessages(simulate_phenotype(pop, h2 = 0.5, seed = 1) |>
                            additive(n_qtn = 30))
  u <- suppressMessages(
    cross_usefulness(sim, scheme = "dh", n_progeny = 40, seed = 2))
  expect_equal(nrow(u), choose(8L, 2L))
  expect_named(u, c("parent1", "parent2", "mean", "sd", "usefulness"))
  # families differ (fixed-effect scoring is not re-standardized per family)
  expect_gt(stats::sd(u$sd), 0)
  expect_gt(stats::sd(u$mean), 0)
  # sorted best-first, and U = mean + i*sd
  expect_equal(u$usefulness, u$mean + attr(u, "intensity") * u$sd,
               tolerance = 1e-8)
  expect_false(is.unsorted(rev(u$usefulness)))
})

test_that("cross_usefulness accepts explicit pairs and a low direction", {
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:6)
  sim <- suppressMessages(simulate_phenotype(pop, h2 = 0.5, seed = 1) |>
                            additive(n_qtn = 20))
  pr <- rbind(c(1, 2), c(3, 4))
  u <- suppressMessages(
    cross_usefulness(sim, pairs = pr, scheme = "dh", n_progeny = 20,
                     direction = "low", seed = 2))
  expect_equal(nrow(u), 2L)
  expect_equal(u$usefulness, u$mean - attr(u, "intensity") * u$sd,
               tolerance = 1e-8)
})
