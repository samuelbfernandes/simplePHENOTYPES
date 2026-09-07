# test-grammar.R
#
# Statistical / structural validation of the v2 grammar (DECISION-009).
# The grammar owes v1 no bit-for-bit parity; it is checked against properties:
# variance-partition identity, seed-threading invariance, realized genetic
# correlation for pleiotropy, QTN-count structure, and the documented errors
# and warnings.

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04

# realized genetic value matrix (no residual), via the internal helper
gen_mat <- function(sim) simplePHENOTYPES:::.genetic_matrix(sim)

# ---------------------------------------------------------------------------
# Foundation, one-call detection
# ---------------------------------------------------------------------------
test_that("foundation with no layers is pure noise (h2 = 0)", {
  ph <- simulate_phenotype(G, seed = 1)
  expect_s3_class(ph, "phenotype_sim")
  expect_length(ph$layers, 0)
  expect_equal(stats::var(ph$pheno$value), 1, tolerance = 0.05)
})

test_that("one-call equals the equivalent piped build", {
  oc <- simulate_phenotype(G, h2 = 0.5, n_qtn = 3, seed = 1)
  pp <- additive(simulate_phenotype(G, seed = 1), prop = 0.5, n_qtn = 3)
  expect_equal(oc$pheno$value, pp$pheno$value)
  expect_length(oc$layers, 1)
})

# ---------------------------------------------------------------------------
# Variance-partition identity
# ---------------------------------------------------------------------------
test_that("realized h2 tracks the sum of genetic proportions", {
  ph <- additive(simulate_phenotype(G, seed = 7), prop = 0.6, n_qtn = 5)
  g <- gen_mat(ph)[, 1]
  y <- ph$pheno$value
  expect_equal(stats::var(g) / stats::var(y), 0.6, tolerance = 0.08)

  ph2 <- simulate_phenotype(G, seed = 7)
  ph2 <- additive(ph2, prop = 0.4, n_qtn = 4)
  ph2 <- dominance(ph2, prop = 0.1, same_as_add = TRUE, degree = 0.5)
  expect_equal(sum(vapply(ph2$layers, function(l) l$prop, 0)), 0.5)
  expect_equal(stats::var(gen_mat(ph2)[, 1]) / stats::var(ph2$pheno$value),
               0.5, tolerance = 0.1)
})

# ---------------------------------------------------------------------------
# Seed-threading invariance (SPEC §6)
# ---------------------------------------------------------------------------
test_that("reordering layers of different types does not change their QTNs", {
  a1 <- simulate_phenotype(G, seed = 5, n_qtn = 3)
  a1 <- additive(a1, 0.3)
  a1 <- epistasis(a1, 0.2, n_pairs = 2)
  a2 <- simulate_phenotype(G, seed = 5, n_qtn = 3)
  a2 <- epistasis(a2, 0.2, n_pairs = 2)
  a2 <- additive(a2, 0.3)
  add1 <- Filter(function(l) l$type == "additive", a1$layers)[[1]]$qtn[[1]]
  add2 <- Filter(function(l) l$type == "additive", a2$layers)[[1]]$qtn[[1]]
  expect_identical(add1, add2)
})

test_that("same seed reproduces identical phenotypes", {
  p1 <- additive(simulate_phenotype(G, seed = 99), prop = 0.5, n_qtn = 4)
  p2 <- additive(simulate_phenotype(G, seed = 99), prop = 0.5, n_qtn = 4)
  expect_equal(p1$pheno$value, p2$pheno$value)
})

# ---------------------------------------------------------------------------
# Pleiotropy: realized genetic correlation (PleioArch, 2 traits)
# ---------------------------------------------------------------------------
test_that("2-trait pleiotropy realizes the target genetic correlation", {
  geno_cor <- function(target, seed) {
    ph <- simulate_phenotype(G, architecture = "pleiotropy", n_traits = 2,
                             seed = seed, cor = target, prop_var_major = 0)
    ph <- additive(ph, prop = c(0.5, 0.5), n_qtn = 60)
    stats::cor(gen_mat(ph))[1, 2]
  }
  for (target in c(0.0, 0.3, 0.6)) {
    m <- mean(vapply(1:30, function(s) geno_cor(target, s), 0))
    expect_equal(m, target, tolerance = 0.1)
  }
})

test_that("pleiotropy enforces cor^2 <= pi_target * pi_secondary", {
  ph <- simulate_phenotype(G, architecture = "pleiotropy", n_traits = 2,
                           seed = 1, cor = 0.9, pi_target = 0.5,
                           pi_secondary = 0.5)
  expect_error(additive(ph, prop = c(0.5, 0.5), n_qtn = 10),
               "Biological constraint")
})

test_that("pleiotropy with n_traits > 2 warns and uses the Cholesky fallback", {
  ph <- simulate_phenotype(G, architecture = "pleiotropy", n_traits = 3,
                           seed = 10, cor = 0.4)
  expect_warning(additive(ph, prop = c(0.4, 0.4, 0.4), n_qtn = 20),
                 "Cholesky")
})

test_that("pleiotropy with n_traits = 1 warns and proceeds as independent", {
  expect_warning(
    ph <- simulate_phenotype(G, architecture = "pleiotropy", n_traits = 1,
                             seed = 1),
    "requires n_traits > 1")
  expect_equal(ph$architecture, "independent")
})

# ---------------------------------------------------------------------------
# QTN-count structure
# ---------------------------------------------------------------------------
test_that("independent architecture draws distinct QTNs per trait", {
  ph <- additive(simulate_phenotype(G, n_traits = 2, seed = 3), prop = 0.4,
                 n_qtn = 5)
  q <- ph$layers[[1]]$qtn
  expect_length(q[[1]], 5)
  expect_length(q[[2]], 5)
  expect_false(identical(q[[1]], q[[2]]))
})

test_that("same_as_add reuses the additive QTNs", {
  ph <- simulate_phenotype(G, seed = 3)
  ph <- additive(ph, prop = 0.4, n_qtn = 5)
  ph <- dominance(ph, prop = 0.1, same_as_add = TRUE)
  expect_identical(ph$layers[[1]]$qtn, ph$layers[[2]]$qtn)
})

# ---------------------------------------------------------------------------
# complex_phenotypes()
# ---------------------------------------------------------------------------
test_that("complex_phenotypes combines models and warns on differing seeds", {
  pleio <- additive(simulate_phenotype(G, architecture = "pleiotropy",
                                       n_traits = 2, seed = 10, cor = 0.5),
                    prop = 0.4, n_qtn = 20)
  indep <- additive(simulate_phenotype(G, n_traits = 2, seed = 11),
                    prop = 0.3, n_qtn = 20)
  expect_warning(both <- complex_phenotypes(pleio, indep, h2 = 0.5),
                 "different seeds")
  expect_equal(both$architecture, "complex")
  expect_equal(both$seed, 10)
  w <- phenotypes_wide(both)
  expect_equal(stats::var(w[[3]]), 1, tolerance = 0.1)
})

# ---------------------------------------------------------------------------
# Validation: errors and warnings
# ---------------------------------------------------------------------------
test_that("genetic proportions summing above 1 error", {
  ph <- additive(simulate_phenotype(G, seed = 1), prop = 0.7, n_qtn = 3)
  expect_error(dominance(ph, prop = 0.5, same_as_add = TRUE), "above 1")
})

test_that("per-layer n_qtn overriding the baseline warns", {
  ph <- simulate_phenotype(G, seed = 1, n_qtn = 3)
  expect_warning(additive(ph, prop = 0.5, n_qtn = 7), "overrides the baseline")
})

test_that("prop length must be 1 or n_traits", {
  ph <- simulate_phenotype(G, n_traits = 3, seed = 1)
  expect_error(additive(ph, prop = c(0.2, 0.3), n_qtn = 3),
               "length 1 or n_traits")
})

# ---------------------------------------------------------------------------
# Output exporters
# ---------------------------------------------------------------------------
test_that("long and wide exporters round-trip the same values", {
  ph <- additive(simulate_phenotype(G, n_traits = 2, seed = 2), prop = 0.4,
                 n_qtn = 3)
  lng <- phenotypes_long(ph)
  wid <- phenotypes_wide(ph)
  expect_named(lng, c("id", "trait", "rep", "value"))
  expect_equal(nrow(wid), length(unique(lng$id)))
  expect_true(all(c("Trait_1", "Trait_2") %in% names(wid)))
})
