# test-combining.R -- DECISION-026: combining ability on a frozen A + D
# architecture. Known answers from docs/SPEC-block3b.md item 1.

.ca_geno <- function(n, m = 60, seed = 1, p = NULL, inbred = FALSE) {
  set.seed(seed)
  if (is.null(p)) p <- stats::runif(m, 0.2, 0.8)
  g <- if (inbred) {
    t(vapply(p, function(pp) 2L * stats::rbinom(n, 1, pp) - 1L, integer(n)))
  } else {
    t(vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1L, integer(n)))
  }
  colnames(g) <- paste0("I", seq_len(n))
  cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                   chr = rep(1:6, each = m / 6), pos = rep(seq_len(m / 6), 6),
                   cm = rep(seq(0, 90, length.out = m / 6), 6),
                   stringsAsFactors = FALSE), as.data.frame(g))
}

test_that("single-locus known answer: inbred tester aa", {
  g <- data.frame(snp = "q", allele = "A/G", chr = 1, pos = 1, cm = 0,
                  AA = 1L, aa = -1L, T = -1L, stringsAsFactors = FALSE)
  p <- as_population(g)
  ca <- combining_ability(p[1:2], p[3], qtn = "q", a = 1, d = 0.5)
  expect_equal(unname(ca$cross_means[, 1]), c(0.5, -1))    # d and -a
  expect_equal(unname(ca$gca[1] - ca$gca[2]), 1.5)          # a + d (p_T = 0)
})

test_that("testers = candidates (HWE panel): GCA is half the breeding value", {
  G <- .ca_geno(80)
  pop <- as_population(G)
  sim <- simulate_phenotype(pop, h2 = 0.6, seed = 2) |>
    additive(prop = 0.4, n_qtn = 20) |> dominance(prop = 0.2)
  te <- template_effects(sim)
  ca <- combining_ability(pop, pop, qtn = te$qtn, a = te$a, d = te$d)
  bv <- .breeding_value_matrix(sim)[, 1]
  expect_equal(unname(ca$gca), unname(bv / 2 - mean(bv / 2)), tolerance = 1e-10)
})

test_that("template_effects reproduces the simulation's A + D value", {
  G <- .ca_geno(60, seed = 4)
  sim <- simulate_phenotype(as_population(G), h2 = 0.5, seed = 5) |>
    additive(prop = 0.3, n_qtn = 10) |> dominance(prop = 0.2)
  te <- template_effects(sim)
  gv <- genotypic_value(sim$geno, te$qtn, te$a, te$d)
  g_sim <- genetic_values(sim)[, 1]
  expect_equal(unname(gv - mean(gv)), unname(g_sim - mean(g_sim)),
               tolerance = 1e-10)
  epi <- simulate_phenotype(as_population(G), h2 = 0.5, seed = 5) |>
    additive(prop = 0.3, n_qtn = 10) |> epistasis(prop = 0.1, n_pairs = 3)
  expect_error(template_effects(epi), "epistasis")
})

test_that("additive architecture: SCA is zero and GCA ranks as additive_value", {
  pop <- as_population(.ca_geno(30, seed = 6))
  tes <- as_population(.ca_geno(5, seed = 7))
  q <- 1:40; a <- stats::rnorm(40)
  ca <- combining_ability(pop, tes, qtn = q, a = a, d = 0)
  expect_equal(max(abs(ca$sca)), 0, tolerance = 1e-12)
  expect_equal(rank(ca$gca), rank(additive_value(pop, q, a)))
})

test_that("design identities: GCA sums to zero, SCA rows/columns sum to zero", {
  pop <- as_population(.ca_geno(8, seed = 8, inbred = TRUE))
  set.seed(1); a <- stats::rnorm(60); d <- stats::rnorm(60)
  fac <- combining_ability(pop[1:5], pop[6:8], qtn = 1:60, a = a, d = d,
                           design = "factorial")
  expect_equal(sum(fac$gca), 0, tolerance = 1e-10)
  expect_equal(sum(fac$gca_testers), 0, tolerance = 1e-10)
  expect_equal(unname(rowSums(fac$sca)), rep(0, 5), tolerance = 1e-10)
  expect_equal(unname(colSums(fac$sca)), rep(0, 3), tolerance = 1e-10)
  dia <- combining_ability(pop, qtn = 1:60, a = a, d = d, design = "diallel")
  expect_equal(sum(dia$gca), 0, tolerance = 1e-10)
  expect_equal(unname(rowSums(dia$sca, na.rm = TRUE)), rep(0, 8), tolerance = 1e-10)
  expect_equal(dia$sca, t(dia$sca))
  # the diallel decomposition reproduces every cross mean
  fit <- dia$grand_mean + outer(dia$gca, dia$gca, "+") + dia$sca
  expect_equal(fit[upper.tri(fit)], dia$cross_means[upper.tri(fit)])
})

test_that("expected uses no random numbers; simulated is seeded", {
  pop <- as_population(.ca_geno(6, seed = 9))
  set.seed(42); before <- .Random.seed
  combining_ability(pop[1:4], pop[5:6], qtn = 1:60, a = rep(1, 60), d = 0.5)
  expect_identical(.Random.seed, before)
  s1 <- combining_ability(pop[1:4], pop[5:6], qtn = 1:60, a = rep(1, 60),
                          d = 0.5, method = "simulated", n_progeny = 5, seed = 3)
  s2 <- combining_ability(pop[1:4], pop[5:6], qtn = 1:60, a = rep(1, 60),
                          d = 0.5, method = "simulated", n_progeny = 5, seed = 3)
  expect_identical(s1$gca, s2$gca)
  expect_equal(n_individuals(s1$progeny), 4L * 2L * 5L)
  expect_equal(as.integer(table(families(s1$progeny, "full_sib"))), rep(5L, 8))
})

test_that("simulated GCA converges on the expected GCA", {
  pop <- as_population(.ca_geno(6, seed = 10))
  set.seed(2); a <- stats::rnorm(60, sd = 0.3); d <- stats::rnorm(60, sd = 0.3)
  ex <- combining_ability(pop[1:4], pop[5:6], qtn = 1:60, a = a, d = d)
  sims <- vapply(1:30, function(s) {
    combining_ability(pop[1:4], pop[5:6], qtn = 1:60, a = a, d = d,
                      method = "simulated", n_progeny = 20, seed = s)$gca
  }, numeric(4))
  m <- rowMeans(sims); se <- apply(sims, 1, stats::sd) / sqrt(30)
  expect_true(all(abs(m - ex$gca) < 3 * se + 1e-8))
})

test_that("input checks", {
  pop <- as_population(.ca_geno(6, seed = 11))
  expect_error(combining_ability(pop, qtn = 1:3, a = 1:3), "needs `testers`")
  expect_error(combining_ability(pop[1:2], qtn = 1:3, a = 1:3, design = "diallel"),
               "at least three")
  expect_error(combining_ability(pop[1:3], pop[4:6], qtn = 1:3, a = 1:2),
               "one value per locus")
  expect_error(combining_ability(pop[1:3], pop[4:6], qtn = 1:3, a = 1:3,
                                 n_progeny = 4), "apply to method")
})

test_that("the same individual twice is rejected by both methods (review A r3)", {
  p <- as_population(.ca_geno(6, m = 24))
  dup <- c(p[1:3], p[2])
  for (m in c("expected", "simulated")) {
    expect_error(combining_ability(dup, p[4], qtn = 1:2, a = c(1, 2),
                                   d = c(0.5, 0.25), method = m, n_progeny = 2,
                                   seed = 1), "same individual more than once")
  }
  expect_error(combining_ability(p[1:3], c(p[4], p[4]), qtn = 1:2, a = c(1, 2),
                                 design = "factorial"), "`testers` lists")
  expect_error(combining_ability(c(p[1:3], p[1]), qtn = 1:2, a = c(1, 2),
                                 design = "diallel"), "same individual")
})
