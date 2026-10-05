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

# ---------------------------------------------------------------------------
# founders_coalescent(): the exported founder generator (phase c)
# ---------------------------------------------------------------------------

test_that("founders_coalescent returns a crossable Population with a rebased linear map", {
  f <- founders_coalescent(30, n_chr = 3, seg_sites = 50, seed = 1)
  expect_s3_class(f, "Population")
  expect_identical(n_individuals(f), 30L)
  expect_identical(nrow(f$map), 150L)
  expect_identical(f$map$snp[1:2], c("1_1", "1_2"))
  for (k in 1:3) {
    cm <- f$map$cm[f$map$chr == k]
    expect_identical(cm[1], 0)                       # rebased to the first site
    expect_false(is.unsorted(cm))
    expect_lte(max(cm), 100)                         # GENERIC: 1 Morgan
  }
  expect_true(all(dosages(f) %in% c(-1, 0, 1)))
  h <- haplotypes(f)
  carriers <- rowSums(h$cis + h$trans)                # derived-allele copies of 60
  expect_true(all(carriers >= 1 & carriers <= 59))    # every kept site segregates
  kid <- cross(f[1], f[2], n = 3, seed = 2)          # usable by the crossing engine
  expect_identical(n_individuals(kid), 3L)
  meta <- attr(f, "coalescent")
  expect_identical(meta$species, "GENERIC")
  expect_length(meta$seeds, 3)
  expect_true(all(meta$n_total >= 50))
})

test_that("inbred founders are fully homozygous and use n_ind haplotypes", {
  f <- founders_coalescent(12, seg_sites = 40, inbred = TRUE, seed = 3)
  expect_true(all(dosages(f) %in% c(-1, 1)))
  expect_identical(attr(f, "coalescent")$n_hap, 12L)
  h <- haplotypes(f)
  expect_identical(h$cis, h$trans)
})

test_that("seed reproduces and restores the RNG; set.seed also reproduces", {
  a <- founders_coalescent(10, n_chr = 2, seg_sites = 30, seed = 9)
  b <- founders_coalescent(10, n_chr = 2, seg_sites = 30, seed = 9)
  expect_identical(dosages(a), dosages(b))
  expect_identical(a$map, b$map)
  set.seed(4); before <- runif(1)
  set.seed(4); founders_coalescent(10, seg_sites = 30, seed = 9); expect_identical(runif(1), before)
  set.seed(5); c1 <- founders_coalescent(10, seg_sites = 30)
  set.seed(5); c2 <- founders_coalescent(10, seg_sites = 30)
  expect_identical(dosages(c1), dosages(c2))
})

test_that("every runMacs() preset runs and reports its parameters", {
  for (sp in c("GENERIC", "MAIZE", "WHEAT", "CATTLE")) {
    f <- founders_coalescent(8, seg_sites = 20, species = sp, seed = 6)
    expect_identical(nrow(f$map), 20L, info = sp)
    meta <- attr(f, "coalescent")
    expect_identical(meta$species, sp)
  }
  # the presets' per-chromosome theta / rho / map length, from the runMacs() source
  pre <- simplePHENOTYPES:::.coalescent_preset
  expect_equal(c(pre("GENERIC")$theta, pre("GENERIC")$rho, pre("GENERIC")$morgans), c(1000, 400, 1))
  expect_equal(c(pre("MAIZE")$theta, pre("MAIZE")$rho, pre("MAIZE")$morgans), c(1000, 800, 2))
  expect_equal(c(pre("WHEAT")$theta, pre("WHEAT")$rho, pre("WHEAT")$morgans), c(320, 288, 1.43))
  expect_equal(nrow(pre("WHEAT")$history), 26L)
  expect_equal(nrow(pre("MAIZE")$history), 15L)
  ct <- pre("CATTLE")
  expect_equal(ct$theta, 2.8e9 / 30 * 9.4e-9 * 360)
  expect_equal(ct$morgans, 9.26e-9 * 2.8e9 / 30)
  expect_equal(ct$history$time[1], 3 / 360)
})

test_that("overrides replace the preset; too few sites is an error (as runMacs)", {
  f <- founders_coalescent(10, seg_sites = 10, theta = 20, rho = 5, morgans = 0.5,
                           history = data.frame(time = 1, size = 2), seed = 7)
  meta <- attr(f, "coalescent")
  expect_equal(c(meta$theta, meta$rho, meta$morgans), c(20, 5, 0.5))
  expect_lte(max(f$map$cm), 50)
  expect_error(founders_coalescent(5, seg_sites = 1e6, theta = 1, rho = 1, seed = 1),
               "fewer than the 1000000 requested")
  expect_error(founders_coalescent(5, seg_sites = 0), "seg_sites")
  expect_error(founders_coalescent(1, inbred = TRUE), "two haplotypes")
  expect_error(founders_coalescent(5, split = -1), "split")
  # an odd n_ind would put one individual's two haplotypes in different demes (Codex)
  expect_error(founders_coalescent(3, split = 10), "even `n_ind`")
  expect_error(founders_coalescent(3, inbred = TRUE, split = 10), "even `n_ind`")
  expect_error(founders_coalescent(5, theta = c(1, 2), n_chr = 3), "theta")
})

test_that("split: between- minus within-subpopulation diversity is 2 theta J", {
  # pairwise differences = theta x pair tree length = 2 theta T (4 N0 units);
  # E[T] is 1/2 within a subpopulation and J + 1/2 between (constant size), so
  # the difference of means is 2 theta J. split = 200 generations with the
  # GENERIC Ne = 100 gives J = 200 / 400 + 1e-6.
  theta <- 30
  # constant size (no change before 1e11) is required for the closed form
  f <- founders_coalescent(4, n_chr = 300, theta = theta, rho = 2, split = 200,
                           history = data.frame(time = 1e11, size = 1), seed = 8)
  h <- haplotypes(f)
  hap <- cbind(h$cis[, 1], h$trans[, 1], h$cis[, 2], h$trans[, 2],
               h$cis[, 3], h$trans[, 3], h$cis[, 4], h$trans[, 4])
  chr <- f$map$chr
  diffs <- function(i, j) tapply(hap[, i] != hap[, j], chr, sum)
  within <- c(diffs(1, 2), diffs(1, 3), diffs(5, 6), diffs(7, 8))
  between <- c(diffs(1, 5), diffs(2, 6), diffs(3, 7), diffs(4, 8))
  J <- 200 / 400 + 1e-6
  est <- mean(between) - mean(within)
  se <- sqrt(var(between) / length(between) + var(within) / length(within))
  expect_lt(abs(est - 2 * theta * J), 5 * se)
  expect_lt(abs(mean(within) - theta), 5 * sd(within) / sqrt(length(within)))
})
