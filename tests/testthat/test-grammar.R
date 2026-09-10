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
gen_mat <- function(sim) genetic_values(sim)

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
  ph2 <- dominance(ph2, prop = 0.1, same_as_add = TRUE)
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

test_that("pleiotropy realizes the target cor for more than two traits", {
  # The exact (PleioArch) engine now covers any number of traits, so this must
  # not warn and must hit the target on every pair.
  expect_no_warning(
    ph <- simulate_phenotype(G, architecture = "pleiotropy", n_traits = 4,
                             seed = 10, cor = 0.5) |>
      additive(prop = 0.5, n_qtn = 200)
  )
  m <- rowMeans(vapply(1:12, function(s) {
    g <- genetic_values(
      simulate_phenotype(G, architecture = "pleiotropy", n_traits = 4,
                         seed = s, cor = 0.5) |>
        additive(prop = 0.5, n_qtn = 200))
    cor(g)[upper.tri(diag(4))]
  }, numeric(6)))
  expect_equal(mean(m), 0.5, tolerance = 0.1)
})

test_that("pleiotropy accepts a full correlation matrix, negatives included", {
  R <- matrix(c(1, 0.7, -0.3, 0.7, 1, -0.2, -0.3, -0.2, 1), 3, 3)
  m <- Reduce(`+`, lapply(1:12, function(s) {
    cor(genetic_values(
      simulate_phenotype(G, architecture = "pleiotropy", n_traits = 3,
                         seed = s, cor = R) |>
        additive(prop = 0.5, n_qtn = 200)))
  })) / 12
  expect_equal(m[1, 2], R[1, 2], tolerance = 0.1)
  expect_equal(m[1, 3], R[1, 3], tolerance = 0.1)
  expect_equal(m[2, 3], R[2, 3], tolerance = 0.1)
})

test_that("unattainable multi-trait correlations are rejected", {
  bad <- matrix(c(1, -0.9, -0.9, -0.9, 1, -0.9, -0.9, -0.9, 1), 3, 3)
  ph <- simulate_phenotype(G, architecture = "pleiotropy", n_traits = 3,
                           seed = 1, cor = bad)
  expect_error(additive(ph, prop = 0.5, n_qtn = 20), "not attainable")
})

test_that("pleiotropic share can vary by trait", {
  ph <- simulate_phenotype(G, architecture = "pleiotropy", n_traits = 3,
                           seed = 2, cor = 0.5, pi = c(1, 0.8, 0.7)) |>
    additive(prop = 0.5, n_qtn = 100)
  q <- ph$layers[[1]]$qtn
  expect_length(q, 3)
  # Every trait keeps n_qtn loci; a shared core carries the covariance.
  expect_true(all(vapply(q, length, 0L) == 100))
  expect_gt(length(Reduce(intersect, q)), 0)
})

test_that("pleiotropy with n_traits = 1 errors", {
  expect_error(
    simulate_phenotype(G, architecture = "pleiotropy", n_traits = 1,
                       seed = 1),
    "requires n_traits > 1")
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

# ---------------------------------------------------------------------------
# PleioArch citation notice (once per session, only when cor is used)
# ---------------------------------------------------------------------------
test_that("controlling cor emits the citation once per session", {
  id <- "simplePHENOTYPES_pleioarch_citation"
  rlang::reset_message_verbosity(id)
  on.exit(rlang::reset_message_verbosity(id), add = TRUE)

  expect_message(
    simulate_phenotype(G, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.5, seed = 1) |> additive(prop = 0.5, n_qtn = 5),
    "please cite Prado et al\\."   # year left loose: the reference is not final yet
  )
  # Second use in the same session stays quiet.
  expect_no_message(
    simulate_phenotype(G, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.5, seed = 2) |> additive(prop = 0.5, n_qtn = 5)
  )
})

test_that("the citation is not emitted without correlation control", {
  id <- "simplePHENOTYPES_pleioarch_citation"
  rlang::reset_message_verbosity(id)
  on.exit(rlang::reset_message_verbosity(id), add = TRUE)

  # Pleiotropy but no `cor`: no correlation is being controlled.
  expect_no_message(
    simulate_phenotype(G, architecture = "pleiotropy", n_traits = 2,
                       seed = 3) |> additive(prop = 0.5, n_qtn = 5)
  )
  expect_no_message(
    simulate_phenotype(G, seed = 4) |> additive(prop = 0.5, n_qtn = 5)
  )
})

# ---------------------------------------------------------------------------
# Reporting accessors: genetic_values(), qtn_table(), realized h2
# ---------------------------------------------------------------------------
test_that("genetic_values() returns the genetic component with dimnames", {
  ph <- simulate_phenotype(G, n_traits = 2, h2 = 0.5, seed = 1) |>
    additive(n_qtn = 5)
  g <- genetic_values(ph)
  expect_true(is.matrix(g))
  expect_identical(dim(g), c(length(ph$ids), 2L))
  expect_identical(rownames(g), ph$ids)
  expect_identical(colnames(g), c("Trait_1", "Trait_2"))
  # matches the internal matrix it reports on
  expect_equal(unname(g), unname(simplePHENOTYPES:::.genetic_matrix(ph)))
})

test_that("realized h2 tracks the requested h2 for mean-effect layers", {
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 2) |> additive(n_qtn = 20)
  g <- genetic_values(ph)[, 1]
  y <- phenotypes_wide(ph)$Trait_1
  expect_equal(stats::var(g) / stats::var(y), 0.5, tolerance = 0.1)
  expect_equal(simplePHENOTYPES:::.realized_h2(ph), 0.5, tolerance = 0.1)
})

test_that("vqtl uses phenotypic variance but is not part of broad h2", {
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 3) |>
    additive(prop = 0.3, n_qtn = 4) |> vqtl(prop = 0.2)
  # realized reflects the additive layer only, not the 0.5 nominal total
  expect_lt(simplePHENOTYPES:::.realized_h2(ph), 0.45)
  expect_output(print(ph), "vqtl is a residual-heterogeneity")
})

test_that("qtn_table() names the loci and effects, expanding epistatic sets", {
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 4) |>
    additive(prop = 0.3, n_qtn = 3) |>
    epistasis(prop = 0.2, n_pairs = 2, interaction = 2)
  tab <- qtn_table(ph)
  expect_named(tab, c("trait", "layer", "set", "snp", "chr", "pos", "maf",
                      "effect", "var_explained", "QTN_t1", "QTN_t2", "ld_r2"))
  # 3 additive rows + 2 pairs x 2 members
  expect_identical(nrow(tab), 7L)
  expect_identical(sum(tab$layer == "epistasis"), 4L)
  expect_true(all(is.na(tab$set[tab$layer == "additive"])))
  # each epistatic set shares one effect across its members
  eff_by_set <- tapply(tab$effect[tab$layer == "epistasis"],
                       tab$set[tab$layer == "epistasis"], function(x) length(unique(x)))
  expect_true(all(eff_by_set == 1))
  # marker names resolve against the map
  expect_true(all(tab$snp %in% ph$map$snp))
})

test_that("qtn_table() is empty but well-formed with no layers", {
  tab <- qtn_table(simulate_phenotype(G, seed = 5))
  expect_identical(nrow(tab), 0L)
  expect_named(tab, c("trait", "layer", "set", "snp", "chr", "pos", "maf",
                      "effect", "var_explained", "QTN_t1", "QTN_t2", "ld_r2"))
})
