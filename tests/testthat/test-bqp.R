# select_ind(method = "bqp"): relatedness-penalized selection of exactly N
# individuals (Montesinos-Lopez et al. 2025). Deterministic, no RNG.

.bqp_brute <- function(cvec, G, N, lambda) {
  cm <- utils::combn(length(cvec), N)
  obj <- apply(cm, 2L, function(s) sum(cvec[s]) - lambda * sum(G[s, s]))
  list(sel = cm[, which.max(obj)], obj = max(obj))
}

.bqp_toy <- function(n, seed = 1) {
  set.seed(seed)
  Z <- matrix(rnorm(n * 40), n)
  G <- tcrossprod(Z) / 40
  list(c = rnorm(n), G = G)
}

.bqp_sim <- function(n = 30, seed = 5) {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES", envir = environment())
  pop <- as_population(SNP55K_maize282_maf04, individuals = seq_len(n))
  simulate_phenotype(pop, h2 = 0.5, seed = seed) |> additive(n_qtn = 20)
}

test_that("exact solver equals brute force on small cases", {
  for (n in c(8, 10, 12)) for (N in c(2, 4, 5)) for (lam in c(0, 0.5, 3)) {
    t <- .bqp_toy(n, seed = n + N)
    sol <- simplePHENOTYPES:::.bqp_solve(t$c, t$G, N, lam)
    bf <- .bqp_brute(t$c, t$G, N, lam)
    expect_equal(sol$method, "exact")
    expect_equal(sol$objective, bf$obj)
    expect_equal(sort(sol$sel), sort(bf$sel))
  }
})

test_that("local search reaches the enumeration optimum on a mid case", {
  t <- .bqp_toy(18, seed = 11)                       # choose(18, 6) = 18564
  ex <- simplePHENOTYPES:::.bqp_solve(t$c, t$G, 6, 2)
  ls <- simplePHENOTYPES:::.bqp_solve(t$c, t$G, 6, 2, max_enum = 10)
  expect_equal(ls$method, "local_search")
  expect_equal(ls$objective, ex$objective, tolerance = 1e-8)
  # deterministic
  ls2 <- simplePHENOTYPES:::.bqp_solve(t$c, t$G, 6, 2, max_enum = 10)
  expect_identical(ls$sel, ls2$sel)
})

test_that("min_gain constraints are respected (exact and local search)", {
  t <- .bqp_toy(14, seed = 3)
  S <- matrix(t$c, ncol = 1)
  rhs <- 4 * 60 / 100
  for (me in c(2e5, 10)) {
    sol <- simplePHENOTYPES:::.bqp_solve(-t$c, t$G, 4, 0.5, S = S, rhs = rhs,
                                         max_enum = me)
    expect_true(sol$feasible)
    expect_gte(sum(S[sol$sel, ]), rhs - 1e-8)
  }
  sol <- simplePHENOTYPES:::.bqp_solve(t$c, t$G, 4, 0.5, S = S, rhs = 1e6)
  expect_false(sol$feasible)
})

test_that("select_ind(method = 'bqp'): lambda = 0 is top-N; relatedness drops with lambda", {
  ph <- .bqp_sim()
  top <- select_ind(ph, n = 5, on = "pheno", method = "mass")
  b0 <- suppressMessages(select_ind(ph, n = 5, on = "pheno", method = "bqp",
                                    lambda = 0))
  expect_setequal(attr(b0, "selected"), attr(top, "selected"))
  expect_equal(n_individuals(b0), 5L)
  G <- g_matrix(ph)
  rel <- function(o) { s <- attr(o, "selected"); mean(G[s, s]) }
  b5 <- suppressMessages(select_ind(ph, n = 5, on = "pheno", method = "bqp",
                                    lambda = 50))
  expect_lte(rel(b5), rel(b0) + 1e-12)
  expect_identical(attr(b5, "method"), "bqp")
  expect_equal(attr(b5, "bqp")$solver, "exact")
  b5b <- suppressMessages(select_ind(ph, n = 5, on = "pheno", method = "bqp",
                                     lambda = 50))
  expect_identical(attr(b5, "selected"), attr(b5b, "selected"))
  # direction = "low" with lambda = 0 is bottom-N
  low <- suppressMessages(select_ind(ph, n = 5, method = "bqp", lambda = 0,
                                     direction = "low"))
  bot <- select_ind(ph, n = 5, method = "mass", direction = "low")
  expect_setequal(attr(low, "selected"), attr(bot, "selected"))
})

test_that("the citation notice is emitted", {
  ph <- .bqp_sim()
  # once per session: the first call in the session may already have shown it
  expect_no_error(suppressMessages(select_ind(ph, n = 4, method = "bqp")))
  expect_match(simplePHENOTYPES:::.cite_main("x"), "Please also cite")
})

test_that("bqp argument validation", {
  ph <- .bqp_sim()
  f <- function(...) suppressMessages(select_ind(ph, n = 4, method = "bqp", ...))
  expect_error(f(lambda = -1), "lambda")
  expect_error(f(lambda = c(1, 2)), "lambda")
  expect_error(f(lambda = NA_real_), "lambda")
  expect_error(f(min_gain = "a"), "min_gain")
  expect_error(f(min_gain = c(1, 2, 3)), "min_gain")
  expect_error(f(min_gain = 1e6), "no set")
  expect_error(select_ind(ph, n = 4, method = "mass", lambda = 2), "bqp")
  expect_error(select_ind(ph, n = 4, method = "mass", min_gain = 1), "bqp")
  expect_error(f(weights = c(1, 2)), "weights")
  expect_error(f(on = function(s) rnorm(30), weights = 1), "named criterion")
})
