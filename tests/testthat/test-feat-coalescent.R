# Coalescent founders, phase (a): the SMC' sampler (docs/SPEC-coalescent.md,
# DECISION-048). Analytic gates under the standard neutral coalescent, time in
# 4 N0 units (theta = 4 N0 mu per chromosome).

a1 <- function(n) sum(1 / seq_len(n - 1))
a2 <- function(n) sum(1 / seq_len(n - 1)^2)

test_that("E[S] = theta * sum 1/i (Watterson), without and with recombination", {
  n <- 10
  theta <- 5
  reps <- 1500
  set.seed(101)
  for (rho in c(0, 20)) {
    s <- vapply(seq_len(reps), function(i)
      .coalescent_chromosome(n, theta, rho)$n_total, numeric(1))
    # Var(S) = theta a1 + theta^2 a2 without recombination (an upper bound with it)
    se <- sqrt((theta * a1(n) + theta^2 * a2(n)) / reps)
    expect_lt(abs(mean(s) - theta * a1(n)), 5 * se)
  }
})

test_that("the site frequency spectrum has E[xi_i] = theta / i", {
  n <- 8
  theta <- 10
  reps <- 1500
  set.seed(102)
  counts <- matrix(0, reps, n - 1)
  for (r in seq_len(reps)) {
    h <- .coalescent_chromosome(n, theta, rho = 5)$hap
    k <- rowSums(h)
    counts[r, ] <- tabulate(k, nbins = n - 1)
  }
  expected <- theta / seq_len(n - 1)
  se <- apply(counts, 2, sd) / sqrt(reps)
  expect_true(all(abs(colMeans(counts) - expected) < 5 * se))
})

test_that("pair TMRCA follows a size history (ms -eN)", {
  # size 1 before t1, l2 after: E[T] = (1 - e^{-2 t1}) / 2 + e^{-2 t1} l2 / 2
  t1 <- 0.3
  l2 <- 5
  hist <- data.frame(time = t1, size = l2)
  set.seed(103)
  tm <- vapply(seq_len(4000), function(i)
    .coalescent_chromosome(2, 0, 0, history = hist)$tmrca, numeric(1))
  expected <- (1 - exp(-2 * t1)) / 2 + exp(-2 * t1) * l2 / 2
  expect_lt(abs(mean(tm) - expected), 5 * sd(tm) / sqrt(length(tm)))
})

test_that("set.seed() and seed = reproduce the haplotypes", {
  set.seed(7); a <- .coalescent_chromosome(12, 30, 20)
  set.seed(7); b <- .coalescent_chromosome(12, 30, 20)
  expect_identical(a$hap, b$hap)
  expect_identical(a$pos, b$pos)
  c1 <- .coalescent_chromosome(12, 30, 20, seed = 5)
  c2 <- .coalescent_chromosome(12, 30, 20, seed = 5)
  expect_identical(c1$hap, c2$hap)
  # a seed = call leaves R's stream alone (no draw)
  set.seed(9); before <- runif(1)
  set.seed(9); .coalescent_chromosome(12, 30, 20, seed = 5); expect_identical(runif(1), before)
})

test_that("output layout: sites x haplotypes 0/1, sorted positions, seg_sites subset", {
  set.seed(8)
  x <- .coalescent_chromosome(15, 80, 40, seg_sites = 25)
  expect_equal(dim(x$hap), c(25L, 15L))
  expect_true(all(x$hap %in% 0:1))
  expect_false(is.unsorted(x$pos))
  expect_true(all(x$pos >= 0 & x$pos < 1))
  expect_true(all(rowSums(x$hap) %in% 1:14))           # every site segregates
  expect_gt(x$n_total, 25)
  all_sites <- .coalescent_chromosome(15, 80, 40, seed = 3)
  expect_equal(nrow(all_sites$hap), all_sites$n_total)
})

test_that("inputs are validated (R and Rust errors are ordinary R errors)", {
  expect_error(.coalescent_chromosome(1, 1, 1), "n_hap")
  expect_error(.coalescent_chromosome(5, -1, 1), "theta")
  expect_error(.coalescent_chromosome(5, 1, NA), "rho")
  expect_error(.coalescent_chromosome(5, 1, 1, history = list(time = 1)), "data frame")
  expect_error(.coalescent_chromosome(5, 1, 1, history = data.frame(time = c(2, 1), size = 1)),
               "strictly increasing")
  expect_error(.coalescent_chromosome(5, 1, 1, history = data.frame(time = 1, size = 0)),
               "sizes in")
})

test_that("the recorded seed of any run replays it (full 32-bit range)", {
  set.seed(42)
  x <- .coalescent_chromosome(6, 10, 5)
  expect_gt(x$seed, .Machine$integer.max)    # this stream draws a seed above 2^31 - 1
  y <- .coalescent_chromosome(6, 10, 5, seed = x$seed)
  expect_identical(y$hap, x$hap)
  expect_identical(y$pos, x$pos)
  expect_error(.coalescent_chromosome(6, 10, 5, seed = 2^32), "2\\^32")
  expect_error(.coalescent_chromosome(6, 10, 5, seed = -1), "2\\^32")
})

test_that("a factor-valued history is refused, not read as level codes", {
  h <- data.frame(time = factor("0.3"), size = factor("5"))
  expect_error(.coalescent_chromosome(2, 0, 0, history = h), "numeric columns")
})

test_that("classed numbers reach the core by value (bit64::integer64)", {
  skip_if_not_installed("bit64")   # not a dependency: looked up indirectly
  i64 <- getExportedValue("bit64", "as.integer64")
  a <- .coalescent_chromosome(8, 20, 20, seed = 17)
  b <- .coalescent_chromosome(8, i64(20), i64(20), seed = i64(17))
  expect_identical(b$hap, a$hap)
  expect_identical(b$pos, a$pos)
})

test_that("rates that overflow or explode are errors, never hangs", {
  expect_error(.coalescent_chromosome(2, 1e308, 1e308, seed = 1), "finite")
  expect_error(.coalescent_chromosome(2, 1e15, 1e15, seed = 1), "too large")
})

test_that("an extreme history is refused, never aborts R or collapses times (Codex probes)", {
  # size 1e308 overflowed to an index panic (seed 41); size 1e-308 made the
  # coalescence rate Inf, so every waiting time was 0
  expect_error(.coalescent_chromosome(2, 0, 0, history = data.frame(time = 5e-324, size = 1e308),
                                      seed = 41), "sizes in")
  expect_error(.coalescent_chromosome(2, 0, 0, history = data.frame(time = 5e-324, size = 1e-308),
                                      seed = 1), "sizes in")
  expect_error(.coalescent_chromosome(2, 0, 0, history = data.frame(time = 1e13, size = 1),
                                      seed = 1), "times must be in")
})

test_that("two-locus TMRCA correlation is SMC' (near the ARG), not SMC", {
  # pair TMRCAs at the two ends of a sequence of scaled length rho: the exact ARG
  # gives (rho + 18) / (rho^2 + 13 rho + 18), the SMC 1 / (1 + rho); the SMC' lies
  # just below the ARG value (Wilton, Carmi & Hobolth 2015, Genetics 200:343-355)
  rho <- 1
  set.seed(104)
  x <- t(vapply(seq_len(20000), function(i) {
    z <- .coalescent_chromosome(2, 0, rho)
    c(z$tmrca, z$tmrca_end)
  }, numeric(2)))
  r <- cor(x[, 1], x[, 2])
  se <- (1 - r^2) / sqrt(nrow(x))
  expect_gt(r, 1 / (1 + rho) + 4 * se)
  expect_lt(r, (rho + 18) / (rho^2 + 13 * rho + 18) + 3 * se)
})
