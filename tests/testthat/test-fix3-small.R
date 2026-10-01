# Round-3 small fixes: R3-3 (supplied-K symmetrization is the correctly
# rounded mean) and R3-8 (seed-overflow message states the inclusive bound).

.f3_s <- 4.94065645841247e-324   # smallest positive subnormal

test_that("R3-3: triangles s and 2s average to 2s (correlation 0.5), not s", {
  K <- matrix(c(4 * .f3_s, .f3_s, 2 * .f3_s, 4 * .f3_s), 2, 2)
  Ks <- .check_relationship(K)
  expect_true(isSymmetric(Ks, tol = 0))
  expect_identical(Ks[1, 2], 2 * .f3_s)
  expect_identical(Ks[2, 1], 2 * .f3_s)
  expect_equal(Ks[1, 2] / Ks[1, 1], 0.5)
})

test_that("R3-3: exactly symmetric subnormal K stays bit-identical; 1.7e308 rank-one stays finite", {
  K <- diag(c(.f3_s, 3 * .f3_s))
  expect_identical(.check_relationship(K), K)
  K2 <- matrix(c(4 * .f3_s, 3 * .f3_s, 3 * .f3_s, 4 * .f3_s), 2, 2)
  expect_identical(.check_relationship(K2), K2)
  Kb <- matrix(1.7e308, 2, 2)
  expect_identical(.check_relationship(Kb), Kb)
  Ko <- matrix(c(1.7e308, 1.7e308, 1.7e308 * (1 - 1e-12), 1.7e308), 2, 2)
  Kos <- .check_relationship(Ko)
  expect_true(all(is.finite(Kos)))
  expect_true(isSymmetric(Kos, tol = 0))
  # opposite signs near the overflow limit: the sum cannot overflow
  Kn <- matrix(c(1.7e308, 1.7e308, 1.7e308 * (1 - 1e-12), 1.7e308), 2, 2)
  Kn[1, 2] <- -1.7e308
  Kn[2, 1] <- -1.7e308 * (1 - 1e-12)
  expect_true(all(is.finite(.check_relationship(Kn))))
})

test_that("R3-3: randomized property test against an exact high-precision mean", {
  skip_if_not_installed("Rmpfr")
  set.seed(303)
  n <- 400
  # off-diagonal scales from subnormal to 1e300; diag 1 keeps the matrix PSD
  # only while |K12| <= 1, so scale the diagonal with the off-diagonal
  mkpair <- function(i) {
    e <- runif(1, -323.5, 300)
    x <- runif(1, 0.1, 1) * 10^e * sample(c(-1, 1), 1)
    rel <- sample(c(1e-16, 1e-12, 1e-9), 1)
    y <- x * (1 + rel * runif(1, -1, 1))
    if (i %% 4 == 0) { # adjacent subnormals / tiny lattice neighbours
      k <- sample(1:40, 1)
      x <- k * .f3_s; y <- (k + sample(c(-1, 1, 2), 1)) * .f3_s
      if (y <= 0) y <- x
    }
    c(x, y)
  }
  bad <- 0L
  for (i in seq_len(n)) {
    p <- mkpair(i)
    x <- p[1]; y <- p[2]
    if (x == y) next
    d <- max(abs(x), abs(y)) * 4
    if (!is.finite(d) || d == 0) next
    K <- matrix(c(d, x, y, d), 2, 2)
    Ks <- .check_relationship(K)
    ref <- as.numeric(
      (Rmpfr::mpfr(x, 2300) + Rmpfr::mpfr(y, 2300)) / 2)
    if (!identical(Ks[1, 2], ref) || !identical(Ks[2, 1], ref)) bad <- bad + 1L
  }
  expect_identical(bad, 0L)
})

test_that("R3-3: on moderate values the mean equals (x + y) / 2 (round-2 formula)", {
  set.seed(304)
  for (i in 1:50) {
    x <- runif(1, 0.1, 0.9)
    y <- x * (1 + 1e-10 * runif(1, -1, 1))
    K <- matrix(c(1, x, y, 1), 2, 2)
    Ks <- .check_relationship(K)
    expect_identical(Ks[1, 2], (x + y) / 2)
    expect_identical(Ks[2, 1], (x + y) / 2)
  }
})

test_that("R3-8: the seed-overflow message states the inclusive magnitude bound", {
  v <- simplePHENOTYPES:::.v1_validate_seed_arith
  h2 <- matrix(.5)
  # mult = round(10 * .5) = 5: (seed + rep) * 5 <= integer.max
  # default (non-LD) call with h2 = .5, rep = 1: derive the real bound
  msg <- function(seed, ...) tryCatch(v(seed, 1, h2, FALSE, FALSE, ...),
                                      error = function(e) conditionMessage(e))
  m <- msg(.Machine$integer.max - 1)
  expect_match(m, "abs\\(seed\\) <= [0-9]+ for this call")
  lim <- as.numeric(sub(".*abs\\(seed\\) <= ([0-9]+) for this call.*", "\\1", m))
  # the stated bound is accepted for both signs; the next integer is rejected
  expect_silent(v(lim, 1, h2, FALSE, FALSE))
  expect_silent(v(-lim, 1, h2, FALSE, FALSE))
  expect_error(v(lim + 1, 1, h2, FALSE, FALSE), "too large")
  # a rejected negative seed gets the same (magnitude) statement
  mneg <- msg(-(.Machine$integer.max - 1))
  expect_match(mneg, "abs\\(seed\\) <= ")
  limn <- as.numeric(sub(".*abs\\(seed\\) <= ([0-9]+) for this call.*", "\\1", mneg))
  expect_silent(v(-limn, 1, h2, FALSE, FALSE))
  expect_false(grepl("below", mneg))
})

test_that("R3-8: LD bound is exactly 214748364 (accepted) / 214748365 (rejected) and is printed", {
  v <- simplePHENOTYPES:::.v1_validate_seed_arith
  h2 <- matrix(.5)
  expect_silent(v(214748364, 1, h2, FALSE, FALSE, ld = TRUE, n_qtn = 0))
  # with rep = 1 and n_qtn = 0 the LD extra is 2, so 10 * 214748364 + 2 fits
  e <- tryCatch(v(214748365, 1, h2, FALSE, FALSE, ld = TRUE, n_qtn = 0),
                error = function(e) conditionMessage(e))
  expect_match(e, "abs\\(seed\\) <= 214748364 for this call")
  e2 <- tryCatch(v(-214748365, 1, h2, FALSE, FALSE, ld = TRUE, n_qtn = 0),
                 error = function(e) conditionMessage(e))
  expect_match(e2, "abs\\(seed\\) <= 214748364 for this call")
})
