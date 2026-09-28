# test-progeny.R -- DECISION-027: progeny testing. Known answers from
# docs/SPEC-block3b.md item 4.

.pt_geno <- function(n, m = 300, seed = 1) {
  set.seed(seed)
  p <- stats::runif(m, 0.2, 0.8)
  g <- t(vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1L, integer(n)))
  colnames(g) <- paste0("I", seq_len(n))
  cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                   chr = rep(1:10, each = m / 10), pos = rep(seq_len(m / 10), 10),
                   cm = rep(seq(0, 100, length.out = m / 10), 10),
                   stringsAsFactors = FALSE), as.data.frame(g))
}

test_that("progeny structure: half-sib families of distinct mates", {
  pop <- as_population(.pt_geno(60, m = 60))
  pt <- progeny_test(pop[1:5], pop[6:60], qtn = 1:60, a = rep(0.1, 60),
                     n_progeny = 8, seed = 1)
  prog <- attr(pt, "progeny")
  expect_equal(n_individuals(prog), 40L)
  hs <- families(prog, "maternal_half_sib")
  expect_equal(as.integer(table(hs)), rep(8L, 5))
  fathers <- split(parentage(prog)$father, as.character(hs))
  expect_true(all(vapply(fathers, function(f) !anyDuplicated(f), logical(1))))
  expect_identical(pt, progeny_test(pop[1:5], pop[6:60], qtn = 1:60,
                                    a = rep(0.1, 60), n_progeny = 8, seed = 1))
})

test_that("with no residual and many progeny, the mean regresses on BV at 1/2", {
  pop <- as_population(.pt_geno(400))
  set.seed(3); a <- stats::rnorm(300, sd = 0.2)
  par <- pop[1:40]; mates <- pop[41:400]
  pt <- progeny_test(par, mates, qtn = 1:300, a = a, n_progeny = 150, seed = 4)
  x <- dosages(pop)[, , drop = FALSE] + 1
  p <- rowMeans(x[, 41:400]) / 2                     # mates' allele frequency
  bv <- colSums((x[, 1:40] - 2 * p) * a)
  slope <- unname(stats::coef(stats::lm(pt$progeny_mean ~ bv))[2])
  expect_equal(slope, 0.5, tolerance = 0.06 / 0.5)
  # with dominance the progeny mean still regresses at 1/2 on the breeding value,
  # now with the average effects alpha = a + d(q - p) of the mates (review A)
  set.seed(5); d <- stats::rnorm(300, sd = 0.3)
  prog <- attr(pt, "progeny")
  yd <- genotypic_value(prog, 1:300, a, d)
  md <- tapply(yd, as.character(families(prog, "maternal_half_sib")), mean)[par$ids]
  bv_d <- colSums((x[, 1:40] - 2 * p) * (a + d * (1 - 2 * p)))
  slope_d <- unname(stats::coef(stats::lm(md ~ bv_d))[2])
  expect_equal(slope_d, 0.5, tolerance = 0.06 / 0.5)
  bv_naive <- colSums((x[, 1:40] - 2 * p) * a)          # ignoring dominance
  expect_gt(stats::cor(md, bv_d), stats::cor(md, bv_naive))
})

test_that("progeny-test accuracy matches sqrt(n h2 / (4 + (n - 1) h2))", {
  pop <- as_population(.pt_geno(600, seed = 7))
  set.seed(8); a <- stats::rnorm(300, sd = 0.2)
  par <- pop[1:150]; mates <- pop[151:600]
  x <- dosages(pop) + 1
  p <- rowMeans(x) / 2
  bv <- colSums((x[, 1:150] - 2 * p) * a)
  for (cfg in list(c(n = 5, h2 = 0.2), c(n = 20, h2 = 0.2), c(n = 5, h2 = 0.5))) {
    n <- cfg[["n"]]; h2 <- cfg[["h2"]]
    r_theory <- sqrt(n * h2 / (4 + (n - 1) * h2))
    r <- vapply(1:6, function(s) {
      pt <- progeny_test(par, mates, qtn = 1:300, a = a, n_progeny = n, h2 = h2,
                         seed = 100 * n + s)
      stats::cor(pt$progeny_mean, bv)
    }, numeric(1))
    se <- (1 - r_theory^2) / sqrt(150 - 1) / sqrt(length(r))
    expect_lt(abs(mean(r) - r_theory), 3 * se + 0.02)
  }
})

test_that("input checks", {
  pop <- as_population(.pt_geno(10, m = 30))
  expect_error(progeny_test(pop[1:2], pop[3:10], qtn = 1:30, a = rep(1, 30),
                            n_progeny = 2, h2 = 0.3, var_e = 1), "at most one")
  expect_error(progeny_test(pop[1:2], pop[3:10], qtn = 1:30, a = rep(1, 29),
                            n_progeny = 2), "^progeny_test\\(\\).*one value per locus")
})

test_that("families are genuine half-sibs: no selfs, no repeated mates (review A)", {
  pop <- as_population(.pt_geno(10, m = 30))
  expect_error(progeny_test(pop[1], pop[1], qtn = 1:30, a = rep(1, 30),
                            n_progeny = 1), "0 distinct mates")
  expect_error(progeny_test(pop[1], pop[2], qtn = 1:30, a = rep(1, 30),
                            n_progeny = 4), "1 distinct mates")
  # the parent inside `mates`, and mates listed twice, count once / not at all
  mates <- c(pop[1:4], pop[2:4])
  pt <- progeny_test(pop[1], mates, qtn = 1:30, a = rep(1, 30), n_progeny = 3,
                     seed = 2)
  ped <- parentage(attr(pt, "progeny"))
  expect_equal(ped$design, rep("cross", 3))
  expect_equal(anyDuplicated(ped$father_key), 0L)
  expect_error(progeny_test(pop[1], mates, qtn = 1:30, a = rep(1, 30),
                            n_progeny = 4), "3 distinct mates")
})

test_that("a parent listed twice is rejected (review A r4)", {
  pop <- as_population(.pt_geno(10, m = 30))
  expect_error(progeny_test(c(pop[1], pop[1]), pop[3:10], qtn = 1:30,
                            a = rep(1, 30), n_progeny = 2),
               "^progeny_test\\(\\).*same individual more than once")
})

test_that("the progeny-test / own-performance crossover is n > (4 - h2)/(1 - h2)", {
  acc <- function(n, h2) sqrt(n * h2 / (4 + (n - 1) * h2))
  for (h2 in c(0.1, 0.2, 0.5)) {
    n0 <- (4 - h2) / (1 - h2)
    expect_equal(acc(n0, h2), sqrt(h2))
    expect_lt(acc(n0 - 1, h2), sqrt(h2))
    expect_gt(acc(n0 + 1, h2), sqrt(h2))
  }
})
