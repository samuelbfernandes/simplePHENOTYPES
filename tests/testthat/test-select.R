# test-select.R
#
# Selection engine (select_ind) and the named breeding-scheme wrappers
# (single_seed_descent, bulk, pedigree, recurrent_selection) plus c.Population.
#
# Checks are behavioral: truncation keeps the right count, response goes in the
# requested direction, the family/index methods run, schemes preserve size and
# reproduce under a seed, and inbreeding schemes drive heterozygosity down.

data("SNP55K_maize282_maf04")

# A segregating F2 (F1 of two inbreds is uniform, so we self it once) gives the
# genetic variance selection needs.
.f2 <- function(n = 60, seed = 2) {
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
  f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
  suppressMessages(selfcross(f1, n = n, seed = seed))
}

.ph <- function(f2, ...) {
  suppressMessages(
    simulate_phenotype(f2, h2 = 0.5, seed = 3, ...) |> additive(n_qtn = 20)
  )
}

test_that("select_ind keeps the requested count and returns a Population", {
  ph <- .ph(.f2(60))
  by_prop <- select_ind(ph, prop = 0.25, on = "pheno")
  expect_s3_class(by_prop, "Population")
  expect_equal(n_individuals(by_prop), 15L)

  by_n <- select_ind(ph, n = 7, on = "pheno")
  expect_equal(n_individuals(by_n), 7L)

  by_i <- select_ind(ph, intensity = 1.4, on = "pheno")
  expect_gt(n_individuals(by_i), 0L)
  # realized intensity is close to the request
  expect_lt(abs(attr(by_i, "intensity") - 1.4), 0.4)
})

test_that("selection differential is positive for high, negative for low", {
  ph <- .ph(.f2(60))
  hi <- select_ind(ph, prop = 0.2, on = "gv", direction = "high")
  lo <- select_ind(ph, prop = 0.2, on = "gv", direction = "low")
  expect_gt(attr(hi, "differential"), 0)
  expect_lt(attr(lo, "differential"), 0)
})

test_that("gv selection separates genetic values more than the mean", {
  ph <- .ph(.f2(60))
  sel <- select_ind(ph, prop = 0.2, on = "gv")
  gv <- genetic_values(ph)[, 1]
  chosen <- attr(sel, "selected")
  expect_gt(mean(gv[chosen]), mean(gv))
})

test_that("a custom criterion (vector and function) is accepted", {
  ph <- .ph(.f2(40))
  score <- stats::setNames(seq_len(ph$n_ind), ph$ids)  # id 1 worst, last best
  top <- select_ind(ph, n = 5, on = score)
  expect_setequal(attr(top, "selected"), utils::tail(ph$ids, 5))

  same <- select_ind(ph, n = 5, on = function(s) score[s$ids])
  expect_setequal(attr(same, "selected"), attr(top, "selected"))
})

test_that("select_ind validates its inputs", {
  ph <- .ph(.f2(40))
  expect_error(select_ind(ph, n = 5, prop = 0.2), "exactly one")
  expect_error(select_ind(ph), "exactly one")
  expect_error(select_ind(ph, prop = 1.5), "in \\(0, 1\\]")
  expect_error(select_ind(ph, prop = 0.2, method = "within_family"),
               "needs a `family`")
  expect_error(select_ind(ph, n = 5, method = "index"), "needs `weights`")
})

test_that("family methods and the Smith-Hazel index run", {
  f2 <- .f2(60)
  ph <- .ph(f2)
  fam <- rep(1:6, each = 10)
  within <- select_ind(ph, prop = 0.3, method = "within_family", family = fam)
  among  <- select_ind(ph, prop = 0.3, method = "among_family", family = fam)
  expect_s3_class(within, "Population")
  expect_s3_class(among, "Population")
  # among-family keeps whole families, so counts land in family-sized steps
  expect_equal(n_individuals(among) %% 10, 0L)

  ph2 <- suppressMessages(
    simulate_phenotype(f2, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.3, h2 = 0.5, seed = 5) |> additive(n_qtn = 20)
  )
  idx <- select_ind(ph2, n = 10, method = "index", weights = c(1, 0.5))
  expect_equal(n_individuals(idx), 10L)
})

test_that("quadratic_index (QGSI) runs and differs from the linear index", {
  f2 <- .f2(60)
  ph2 <- suppressMessages(
    simulate_phenotype(f2, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.2, h2 = 0.5, seed = 5) |> additive(n_qtn = 25)
  )
  expect_error(select_ind(ph2, n = 10, method = "quadratic_index"),
               "needs `weights`")
  expect_error(
    select_ind(ph2, n = 10, method = "quadratic_index", weights = c(1, 0.5),
               quad_weights = matrix(0, 3, 3)),
    "n_traits x n_traits")
  # zero quad_weights reduces to a plain linear index on the same weights
  W0 <- select_ind(ph2, n = 10, method = "quadratic_index", weights = c(1, 0.5))
  expect_equal(n_individuals(W0), 10L)
  # a non-zero W (interaction + intermediate optimum) changes the ranking
  qi <- select_ind(ph2, n = 10, method = "quadratic_index", weights = c(1, 0.5),
                   quad_weights = matrix(c(0, 0.3, 0.3, -0.5), 2, 2))
  expect_s3_class(qi, "Population")
  expect_false(setequal(attr(W0, "selected"), attr(qi, "selected")))
})

test_that("Smith-Hazel index tolerates a singular phenotypic covariance", {
  # A rank-1 P (perfectly collinear traits, e.g. pleiotropy cor = 1) is a valid
  # model but makes solve() error; the index falls back to a pseudo-inverse.
  P <- matrix(c(1, 1, 1, 1), 2, 2)      # rank 1
  r <- c(1, 1)
  expect_warning(b <- simplePHENOTYPES:::.index_weights(P, r), "singular")
  # r lies in the column space, so the min-norm solution reproduces it: P b = r.
  expect_equal(as.numeric(P %*% b), r, tolerance = 1e-8)
  # A well-posed P matches solve() exactly and does not warn.
  P2 <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  expect_silent(b2 <- simplePHENOTYPES:::.index_weights(P2, r))
  expect_equal(as.numeric(b2), as.numeric(solve(P2, r)), tolerance = 1e-10)
})

test_that("index scorers reject a non-phenotype_sim with a clear error", {
  bad <- structure(list(), class = "bad")
  expect_error(simplePHENOTYPES:::.index_score(bad, 1, 1L), "phenotype_sim")
  expect_error(simplePHENOTYPES:::.quadratic_index_score(bad, 1, NULL, 1L),
               "phenotype_sim")
})

test_that("index methods reject non-finite weights instead of ranking on NA", {
  f2 <- .f2(60)
  ph2 <- suppressMessages(
    simulate_phenotype(f2, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.3, h2 = 0.5, seed = 5) |> additive(n_qtn = 20))
  expect_error(
    select_ind(ph2, n = 10, method = "index", weights = c(NA, 1)),
    "finite")
  expect_error(
    select_ind(ph2, n = 10, method = "quadratic_index", weights = c(Inf, 1)),
    "finite")
  expect_error(
    select_ind(ph2, n = 10, method = "quadratic_index", weights = c(1, 1),
               quad_weights = matrix(c(0, NA, NA, 0), 2, 2)),
    "finite")
})

test_that("index methods warn that a supplied `on` is ignored", {
  f2 <- .f2(60)
  ph2 <- suppressMessages(
    simulate_phenotype(f2, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.3, h2 = 0.5, seed = 5) |> additive(n_qtn = 20))
  expect_warning(
    select_ind(ph2, n = 10, method = "index", weights = c(1, 0.5),
               on = function(s) stats::runif(s$n_ind)),
    "ignores `on`")
  # the default on = "pheno" does not trigger the warning
  expect_silent(
    select_ind(ph2, n = 10, method = "index", weights = c(1, 0.5)))
})

test_that("select_ind can rank on the transmissible breeding value", {
  ph <- .ph(.f2(60))
  bv <- select_ind(ph, prop = 0.2, on = "bv")
  expect_s3_class(bv, "Population")
  # For an additive model the transmissible breeding value tracks the total
  # genetic value, so selecting on bv also raises the mean genetic value.
  gv_all <- genetic_values(ph)[, 1]
  expect_gt(mean(gv_all[attr(bv, "selected")]), mean(gv_all))
})

test_that(".avg_effect gives the classical transmissible average effect (Codex O1)", {
  # Codex O1 counterexample: dominance-only locus a = 0, d = sqrt(2), counted-allele
  # frequency p = 3/8. The classical average effect is
  # alpha = a + d(q - p) = sqrt(2) * (1 - 2*0.375) = sqrt(2)/4, giving the Mendelian
  # expected-offspring breeding values below -- not the sample-regression slope
  # (which is ~64% smaller because the sample is not in Hardy-Weinberg).
  alpha <- simplePHENOTYPES:::.avg_effect(a = 0, d = sqrt(2), p = 0.375)
  expect_equal(alpha, sqrt(2) / 4, tolerance = 1e-8)
  xg <- c(0, 0, 1, 2)
  expect_equal(alpha * (xg - 2 * 0.375),
               c(-0.265165, -0.265165, 0.088388, 0.441942), tolerance = 1e-6)
})

test_that("breeding value / index are refused under an epistasis layer", {
  # An epistatic term has no per-locus a/d, so its induced average effects cannot
  # be reconstructed; the breeding value would be incomplete, so it errors rather
  # than return a partial (misranking) value.
  f2 <- .f2(60)
  ph_epi <- suppressMessages(
    simulate_phenotype(f2, h2 = 0.5, seed = 3) |>
      additive(n_qtn = 15, prop = 0.3) |> epistasis(prop = 0.2, n_pairs = 3))
  expect_error(select_ind(ph_epi, prop = 0.2, on = "bv"), "epistasis")
  expect_error(select_ind(ph_epi, n = 5, method = "index", weights = 1),
               "epistasis")
})

test_that("combined selection is restricted to phenotypic records", {
  f2 <- .f2(60)
  ph  <- .ph(f2)
  fam <- rep(1:6, each = 10)
  # The Lush index assumes Var(own)=V_P, Cov(A,own)=V_A; applying it to a bv/gv/
  # custom score misweights, so those are rejected.
  expect_error(
    select_ind(ph, prop = 0.3, method = "combined", family = fam, h2 = 0.5,
               on = "gv"),
    "phenotypic records")
  expect_s3_class(
    select_ind(ph, prop = 0.3, method = "combined", family = fam, h2 = 0.5),
    "Population")            # default on = "pheno" still works
})

test_that("breeding value equals total genetic value for an additive F2 (LD-robust)", {
  # An F2 is in strong genome-wide LD yet HWE per locus; the analytic average-effect
  # breeding value returns the exact additive genetic value there (a sample
  # projection would too, but a per-locus sample regression would not).
  ph <- .ph(.f2(60))
  bv <- simplePHENOTYPES:::.breeding_value_matrix(ph, 1L)[, 1]
  gv <- genetic_values(ph)[, 1]
  expect_gt(stats::cor(bv, gv), 0.999)
})

test_that("combined selection needs h2 and produces sane weights", {
  f2 <- .f2(60)
  ph <- .ph(f2)
  fam <- rep(1:6, each = 10)
  expect_error(select_ind(ph, prop = 0.3, method = "combined", family = fam),
               "needs `h2`")
  cmb <- select_ind(ph, prop = 0.3, method = "combined", family = fam,
                    h2 = 0.5, family_relationship = 0.5)
  expect_s3_class(cmb, "Population")
  expect_equal(n_individuals(cmb), 18L)

  # as h2 -> 1 the family mean drops out, so combined collapses onto mass
  mass <- select_ind(ph, n = 18, on = "pheno")
  near1 <- select_ind(ph, n = 18, method = "combined", family = fam,
                      h2 = 0.999, family_relationship = 0.5)
  expect_setequal(attr(near1, "selected"), attr(mass, "selected"))
})

test_that("additive_value scores on a fixed cross-generational scale", {
  f2  <- .f2(40)
  qtn <- c(1L, 5L, 9L)
  eff <- c(0.5, -1, 2)
  av  <- additive_value(f2, qtn = qtn, effect = eff)
  expect_length(av, n_individuals(f2))
  expect_named(av, f2$ids)
  # exactly sum_j dosage_ij * effect_j
  d <- dosages(f2)
  expect_equal(unname(av), unname(colSums(d[qtn, ] * eff)), tolerance = 1e-12)
  # fixed scale: scoring a subset does NOT change the values (the whole point --
  # genetic_values() would re-centre/re-scale per population)
  sub <- f2[1:10]
  expect_equal(additive_value(sub, qtn, eff), av[sub$ids], tolerance = 1e-12)
  # loci by name resolve identically
  expect_equal(additive_value(f2, qtn = f2$map$snp[qtn], effect = eff), av,
               tolerance = 1e-12)
  # validation
  expect_error(additive_value(f2, qtn = qtn, effect = c(1, 2)),
               "one value per locus")
  expect_error(additive_value(f2, qtn = c(1, 999999L), effect = c(1, 2)),
               "marker indices")
})

test_that("c.Population pools populations and rejects mismatched maps", {
  f2 <- .f2(20)
  a <- f2[1:5]
  b <- f2[6:10]
  pooled <- c(a, b)
  expect_s3_class(pooled, "Population")
  expect_equal(n_individuals(pooled), 10L)
  expect_false(anyDuplicated(pooled$ids) > 0)

  other <- as_population(SNP55K_maize282_maf04, individuals = 1:3)
  small <- other[1]
  small$map <- small$map[1:10, ]
  small$cis <- small$cis[1:10, , drop = FALSE]
  small$trans <- small$trans[1:10, , drop = FALSE]
  expect_error(c(f2[1], small), "different marker maps")
})

test_that("single seed descent preserves line count and drives homozygosity", {
  f2 <- .f2(30)
  het0 <- mean(dosages(f2) == 0)
  ril <- suppressMessages(single_seed_descent(f2, generations = 4, seed = 2))
  expect_equal(n_individuals(ril), 30L)
  expect_lt(mean(dosages(ril) == 0), het0)
})

test_that("schemes are reproducible under a seed", {
  f2 <- .f2(30)
  a <- suppressMessages(single_seed_descent(f2, generations = 3, seed = 11))
  b <- suppressMessages(single_seed_descent(f2, generations = 3, seed = 11))
  expect_identical(dosages(a), dosages(b))
})

test_that("bulk carries the requested number forward", {
  f2 <- .f2(40)
  bk <- suppressMessages(bulk(f2, generations = 3, n = 25, seed = 2))
  expect_equal(n_individuals(bk), 25L)
})

test_that("pedigree selection maintains size and records history", {
  f2 <- .f2(60)
  pheno <- function(p) {
    suppressMessages(simulate_phenotype(p, h2 = 0.5, seed = 7) |>
                       additive(n_qtn = 30))
  }
  out <- suppressMessages(
    pedigree(f2, pheno, generations = 3, prop = 0.3, seed = 2))
  expect_equal(n_individuals(out), 60L)
  h <- attr(out, "history")
  expect_equal(nrow(h), 3L)
  expect_true(all(h$differential > 0))
})

test_that("recurrent selection intercrosses and records history", {
  f2 <- .f2(60)
  pheno <- function(p) {
    suppressMessages(simulate_phenotype(p, h2 = 0.5, seed = 7) |>
                       additive(n_qtn = 30))
  }
  out <- suppressMessages(
    recurrent_selection(f2, pheno, cycles = 2, n_parents = 8,
                        progeny_per_cross = 6, seed = 2))
  expect_s3_class(out, "Population")
  h <- attr(out, "history")
  expect_equal(nrow(h), 2L)
  expect_true(all(h$n_parents == 8L))
})
