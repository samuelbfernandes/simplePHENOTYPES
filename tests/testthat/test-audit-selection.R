# test-audit-selection.R
#
# Regression tests for the v2 selection audit (select_ind() and the scheme
# wrappers): input contracts (vector trait, trait range, NA family, direction
# vectors, prop + n_select), name-based alignment of the index scores, honest
# statistics for method = "random", a genuinely random bulk(), explicit failure
# modes for degenerate schemes, RNG hygiene, and equation-level pins.

data("SNP55K_maize282_maf04")

.af2 <- function(n = 60, seed = 2, parents = 1:2) {
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
  f1 <- cross(pop[parents[1]], pop[parents[2]], n = 1, seed = 1)
  suppressMessages(selfcross(f1, n = n, seed = seed))
}

.aph <- function(f2, n_traits = 1, ...) {
  suppressMessages(
    simulate_phenotype(f2, h2 = 0.5, n_traits = n_traits, seed = 3, ...) |>
      additive(n_qtn = 20)
  )
}

.apheno <- function(p) {
  suppressMessages(simulate_phenotype(p, h2 = 0.5, seed = 7) |>
                     additive(n_qtn = 30))
}
.apheno2 <- function(p) {
  suppressMessages(simulate_phenotype(p, h2 = 0.5, n_traits = 2, seed = 7) |>
                     additive(n_qtn = 30))
}

# ---- O1 / F8: `trait` must be one in-range index ---------------------------

test_that("a vector `trait` is an error, never recycled (O1)", {
  ph2 <- .aph(.af2(60), n_traits = 2)
  for (on in c("pheno", "gv", "bv")) {
    expect_error(select_ind(ph2, n = 5, on = on, trait = c(1, 2)),
                 "one trait index")
  }
  # a valid scalar still works and differs by trait
  s1 <- attr(select_ind(ph2, n = 5, trait = 1), "selected")
  s2 <- attr(select_ind(ph2, n = 5, trait = 2), "selected")
  expect_false(setequal(s1, s2))
  # the vector form remains valid where it is defined (culling)
  expect_no_error(suppressMessages(
    select_ind(ph2, method = "culling", culling = c(0.6, 0.6), trait = c(1, 2))))
})

test_that("a trait out of range gets an accurate error, incl. tandem (F8)", {
  ph <- .aph(.af2(60))
  expect_error(select_ind(ph, n = 5, trait = 3), "trait index in 1\\.\\.1")
  expect_error(select_ind(ph, n = 5, trait = 0), "trait index in 1\\.\\.1")
  f2 <- .af2(60)
  expect_error(
    suppressMessages(pedigree(f2, .apheno, generations = 2, prop = 0.3,
                              trait = c(1, 2), seed = 1)),
    "returned 1 trait")
  expect_error(
    suppressMessages(recurrent_selection(f2, .apheno, cycles = 2, n_parents = 5,
                                         progeny_per_cross = 5, trait = c(1, 2),
                                         seed = 1)),
    "returned 1 trait")
})

# ---- F1: index scores are joined to individuals by NAME --------------------

test_that("index / Smith-Hazel / QGSI / Lush selection ignores sim$pheno row order (F1)", {
  f2 <- .af2(60)
  ph2 <- .aph(f2, n_traits = 2)
  fam <- rep(1:6, each = 10)
  run <- function(s) {
    list(
      idx  = attr(select_ind(s, n = 10, method = "index", weights = c(1, 0.5)),
                  "selected"),
      qgsi = attr(select_ind(s, n = 10, method = "quadratic_index",
                             weights = c(1, 0.5),
                             quad_weights = matrix(c(0.2, 0.1, 0.1, 0.3), 2)),
                  "selected"),
      comb = attr(select_ind(s, n = 10, method = "combined", family = fam,
                             h2 = 0.5), "selected"),
      mass = attr(select_ind(s, n = 10), "selected"))
  }
  base <- run(ph2)
  set.seed(99)
  shuf <- ph2
  shuf$pheno <- ph2$pheno[sample(nrow(ph2$pheno)), ]
  expect_false(identical(shuf$pheno$id, ph2$pheno$id))
  after <- run(shuf)
  for (m in names(base)) expect_setequal(after[[m]], base[[m]])
})

test_that("the Smith-Hazel score equals b = P^-1 G a computed by hand (T2)", {
  ph2 <- .aph(.af2(60), n_traits = 2)
  w <- c(1, 0.5)
  G <- .breeding_value_matrix(ph2, 1)
  P <- cbind(.criterion_values(ph2, "pheno", 1, 1),
             .criterion_values(ph2, "pheno", 2, 1))
  b <- solve(stats::cov(P), stats::cov(G) %*% w)
  hand <- ph2$ids[order(as.numeric(P %*% b), decreasing = TRUE)[1:10]]
  sel <- select_ind(ph2, n = 10, method = "index", weights = w)
  expect_setequal(attr(sel, "selected"), hand)
})

# ---- F2 / O10: method = "random" -------------------------------------------

test_that("method = 'random' reports honest S and i and validates `on` (F2/O10)", {
  ph <- .aph(.af2(60))
  set.seed(1)
  s <- select_ind(ph, n = 5, method = "random")
  y <- .criterion_values(ph, "pheno", 1, 1)
  sel <- attr(s, "selected")
  # S and i are measured on the phenotype actually reported, not on the runif
  # score that drove the draw
  expect_equal(attr(s, "differential"), mean(y[sel]) - mean(y), tolerance = 1e-12)
  expect_equal(attr(s, "intensity"), (mean(y[sel]) - mean(y)) / stats::sd(y),
               tolerance = 1e-12)
  expect_equal(attr(s, "criterion"), "pheno")
  # unbiased: the average differential over many draws is near zero, far from
  # the uniform-score i(p) ~ 1.6 that was reported before
  ds <- vapply(1:300, function(k) {
    attr(select_ind(ph, n = 5, method = "random"), "differential")
  }, numeric(1))
  expect_lt(abs(mean(ds)), 0.35 * stats::sd(y))
  expect_error(select_ind(ph, n = 5, method = "random", on = "foo"),
               "`on` must be")
  # the RNG draw itself is unchanged by the reporting fix
  set.seed(1)
  a <- attr(select_ind(ph, n = 5, method = "random"), "selected")
  set.seed(1)
  expect_identical(a, ph$ids[order(stats::runif(60), decreasing = TRUE)[1:5]])
})

# ---- F3: bulk() is a random bulk -------------------------------------------

.bulk_contrib <- function(f2, size, seed, gens = 1L) {
  bk <- suppressMessages(bulk(f2, generations = gens, n = size, seed = seed))
  pp <- parentage(bk)
  as.integer(table(factor(pp$mother_key, levels = unique(parentage(f2)$key))))
}

test_that("bulk() draws its progeny at random from the pool (F3)", {
  f2 <- .af2(60)[1:10]
  # n = k * N was fully deterministic (every parent contributes exactly k)
  ctr <- lapply(1:12, function(s) .bulk_contrib(f2, 30, s))
  expect_true(all(vapply(ctr, sum, numeric(1)) == 30))
  expect_gt(length(unique(lapply(ctr, sort))), 1L)
  v <- vapply(ctr, stats::var, numeric(1))
  expect_gt(mean(v), 0.5 * 3.0)          # multinomial(30, 1/10): E[var] = 3.0
  # n = N: contributions are no longer capped at 2 per parent
  ctr2 <- lapply(1:60, function(s) .bulk_contrib(f2, 10, s))
  expect_true(any(vapply(ctr2, max, numeric(1)) > 2))
  expect_gt(mean(vapply(ctr2, stats::var, numeric(1))), 0.75)   # multinomial ~ 1
  # sizes smaller than the current population still work
  bk <- suppressMessages(bulk(f2, generations = 2, n = 4, seed = 1))
  expect_equal(n_individuals(bk), 4L)
  expect_false(anyDuplicated(bk$ids) > 0)
})

# ---- F5: recurrent_selection with n_parents >= N ---------------------------

test_that("recurrent_selection() refuses n_parents >= the phenotyped population (F5)", {
  f2 <- .af2(60)
  expect_error(
    suppressMessages(recurrent_selection(f2[1:5], .apheno, cycles = 1,
                                         n_parents = 10, seed = 1)),
    "n_parents")
  expect_error(
    suppressMessages(recurrent_selection(f2[1:5], .apheno, cycles = 1,
                                         n_parents = 5, seed = 1)),
    "n_parents")
  # a later cycle that would be too small is caught up front, before any work
  expect_error(
    suppressMessages(recurrent_selection(f2, .apheno, cycles = 2, n_parents = 8,
                                         n_crosses = 2, progeny_per_cross = 3,
                                         seed = 1)),
    "n_crosses \\* progeny_per_cross")
})

# ---- F6: among_family count semantics --------------------------------------

test_that("among_family keeps whole families, so n is a floor (F6)", {
  ph <- .aph(.af2(60))
  fam <- rep(1:6, each = 10)
  expect_equal(n_individuals(select_ind(ph, n = 3, method = "among_family",
                                        family = fam)), 10L)
  expect_equal(n_individuals(select_ind(ph, n = 10, method = "among_family",
                                        family = fam)), 10L)
  expect_equal(n_individuals(select_ind(ph, n = 15, method = "among_family",
                                        family = fam)), 20L)
  # within_family, by contrast, returns exactly n
  expect_equal(n_individuals(select_ind(ph, n = 7, method = "within_family",
                                        family = fam)), 7L)
  # singleton families with keep = 1 never overshoot (largest-remainder)
  expect_equal(n_individuals(select_ind(ph, n = 1, method = "within_family",
                                        family = seq_len(60))), 1L)
})

# ---- F7 / O6: fixed-scale caveat -------------------------------------------

test_that("simulate_phenotype callbacks re-standardise; phenotype_value() does not (F7/O6)", {
  f2 <- .af2(80)
  d <- dosages(f2)
  poly <- which(apply(d, 1, stats::sd) > 0)
  q <- poly[round(seq(5, length(poly) - 5, length.out = 5))]
  eff <- c(0.5, -0.3, 0.4, 0.2, -0.6)
  mk <- function(p) suppressMessages(
    simulate_phenotype(p, h2 = 0.5, seed = 3) |> additive(qtn = q, effect = eff))
  base <- mk(f2)
  top <- select_ind(base, prop = 0.25, on = "gv")
  sub <- mk(top)
  # the sim's own genetic layer is re-scaled to `prop` on every population ...
  expect_equal(stats::var(genetic_values(sub)[, 1]),
               stats::var(genetic_values(base)[, 1]), tolerance = 1e-6)
  # ... whereas on the fixed scale the additive variance really falls
  expect_lt(stats::var(additive_value(top, q, eff)),
            0.9 * stats::var(additive_value(f2, q, eff)))
  # the documented fixed-scale recipe: rank on phenotype_value() via `on =`
  ve <- attr(phenotype_value(f2, q, eff, h2 = 0.5, seed = 1), "var_e")
  fixed <- function(s) phenotype_value(s$geno, q, eff, var_e = ve)
  out <- suppressMessages(
    pedigree(f2, mk, generations = 2, prop = 0.25, on = fixed, seed = 4))
  h <- attr(out, "history")
  expect_equal(nrow(h), 2L)
  expect_true(all(is.finite(h$differential)) && all(h$differential > 0))
})

# ---- F9 / O5: NA in `family` -----------------------------------------------

test_that("NA in `family` is an error for every family method (F9/O5)", {
  ph <- .aph(.af2(60))
  fam <- c(NA, rep(1:2, each = 30)[-1])
  for (m in c("within_family", "among_family", "combined")) {
    expect_error(select_ind(ph, n = 5, method = m, family = fam, h2 = 0.5),
                 "NA")
  }
})

# ---- F10 / F11: argument contracts of the schemes --------------------------

test_that("a direction vector is an error, not silently truncated (F10)", {
  f2 <- .af2(60)
  ph <- .aph(f2)
  expect_error(select_ind(ph, n = 5, direction = c("high", "low")), "direction")
  expect_error(suppressMessages(pedigree(f2, .apheno, generations = 1, prop = 0.2,
                                         direction = c("high", "low"))),
               "direction")
  expect_error(suppressMessages(recurrent_selection(f2, .apheno, cycles = 1,
                                                    n_parents = 5,
                                                    direction = c("high", "low"))),
               "direction")
  expect_error(suppressMessages(pedigree(f2, .apheno, generations = 1, prop = 0.2,
                                         direction = "sideways")), "direction")
  # "low" is honoured
  out <- suppressMessages(pedigree(f2, .apheno, generations = 1, prop = 0.2,
                                   direction = "low", seed = 1))
  expect_lt(attr(out, "history")$differential, 0)
})

test_that("pedigree() rejects prop and n_select together (F11/O8)", {
  f2 <- .af2(60)
  expect_error(suppressMessages(pedigree(f2, .apheno, generations = 1,
                                         prop = 0.9, n_select = 2)),
               "not both")
  out <- suppressMessages(pedigree(f2, .apheno, generations = 1, n_select = 2,
                                   seed = 1))
  expect_equal(attr(out, "history")$n_selected, 2L)
})

# ---- O2: callbacks that do not advance the scheme --------------------------

test_that("a callback not backed by the population it was given is an error (O2)", {
  fA <- .af2(40, parents = 1:2)
  fB <- .af2(40, parents = 3:4)
  simB <- .apheno(fB)
  expect_error(suppressMessages(pedigree(fA, function(p) simB, generations = 2,
                                         prop = 0.3, seed = 1)),
               "not backed by the population")
  expect_error(suppressMessages(recurrent_selection(fA, function(p) simB,
                                                    cycles = 2, n_parents = 5,
                                                    progeny_per_cross = 4,
                                                    seed = 1)),
               "not backed by the population")
  # the realistic typo: the closure captures the base population, so the second
  # generation would silently re-select from it
  expect_error(suppressMessages(pedigree(fA, function(p) .apheno(fA),
                                         generations = 3, prop = 0.3, seed = 1)),
               "not backed by the population")
  # a callback that uses its argument is fine
  expect_no_error(suppressMessages(pedigree(fA, .apheno, generations = 2,
                                            prop = 0.3, seed = 1)))
})

# ---- O3: overflow guard ----------------------------------------------------

test_that("overflowing scores / SD give a clear error, never NaN statistics (O3)", {
  ph2 <- .aph(.af2(60), n_traits = 2)
  expect_error(select_ind(ph2, n = 5, method = "quadratic_index",
                          weights = c(1, 1), quad_weights = diag(1e308, 2)),
               "non-finite|overflow")
  ph <- .aph(.af2(60))
  crit <- rep(c(1e308, -1e308), length.out = 60)
  expect_error(select_ind(ph, n = 5, on = crit), "overflow|non-finite")
  # a finite large-scale criterion still works
  ok <- select_ind(ph, n = 5, on = crit / 1e300)
  expect_true(is.finite(attr(ok, "differential")))
})

test_that(".index_weights validates its covariance inputs (Codex)", {
  expect_error(.index_weights(matrix(c(1, NA, NA, 1), 2), c(1, 1)),
               "non-finite|missing")
  expect_error(.index_weights(diag(2), c(1, NA)), "non-finite|missing")
})

# ---- O7 / O9 / F13: c.Population tolerance, seeds, RNG hygiene --------------

test_that("c.Population map tolerance is explicit and pinned (O7)", {
  f2 <- .af2(20)
  a <- f2[1:3]
  b <- f2[4:6]
  b1 <- b
  b1$map$cm[1] <- b1$map$cm[1] + 1e-10       # numerical noise: pooled
  expect_equal(n_individuals(c(a, b1)), 6L)
  b2 <- b
  b2$map$cm[1] <- b2$map$cm[1] + 1e-5        # a real difference: rejected
  expect_error(c(a, b2), "different marker maps")
  b3 <- b
  b3$map$cm[500] <- b3$map$cm[500] + 0.01    # one marker moved: rejected
  expect_error(c(a, b3), "different marker maps")
})

test_that("scheme seeds are validated before any RNG use (O9)", {
  f2 <- .af2(20)
  for (bad in list(1.9, c(1, 999), -1, NA, "a")) {
    expect_error(single_seed_descent(f2, generations = 1, seed = bad), "seed")
    expect_error(bulk(f2, generations = 1, seed = bad), "seed")
    expect_error(suppressMessages(
      pedigree(f2, .apheno, generations = 1, prop = 0.3, seed = bad)), "seed")
    expect_error(suppressMessages(
      recurrent_selection(f2, .apheno, cycles = 1, n_parents = 5, seed = bad)),
      "seed")
  }
})

test_that("every wrapper reproduces under a seed and restores the caller's RNG (F13)", {
  f2 <- .af2(30)
  runs <- list(
    ssd = function(s) single_seed_descent(f2, generations = 2, seed = s),
    bulk = function(s) bulk(f2, generations = 2, n = 20, seed = s),
    ped = function(s) pedigree(f2, .apheno, generations = 2, prop = 0.3, seed = s),
    rec = function(s) recurrent_selection(f2, .apheno, cycles = 2, n_parents = 6,
                                          progeny_per_cross = 5, seed = s))
  for (nm in names(runs)) {
    a <- suppressMessages(runs[[nm]](11))
    b <- suppressMessages(runs[[nm]](11))
    expect_identical(dosages(a), dosages(b), info = nm)
    set.seed(10); expected <- stats::runif(2)
    set.seed(10)
    suppressMessages(runs[[nm]](11))
    expect_identical(stats::runif(2), expected, info = nm)
  }
})

# ---- documentation-contract pins (F12/F14/F15/F16) -------------------------

test_that("count resolution: prop rounds half-to-even with a floor of 1; intensity inverts i(p)", {
  keep <- function(p) .resolve_keep(NULL, p, NULL, 5L)
  expect_equal(vapply(c(0.1, 0.3, 0.5, 0.7, 0.9, 1), keep, numeric(1)),
               c(1, 2, 2, 4, 4, 5))
  ip <- function(k, N) stats::dnorm(stats::qnorm(1 - k / N)) / (k / N)
  for (k in c(1L, 10L, 100L, 500L, 999L)) {
    expect_equal(.resolve_keep(NULL, NULL, ip(k, 1000L), 1000L), k)
  }
})

test_that("ties, direction = 'low' and the realized intensity (sample SD) behave as documented", {
  ph <- .aph(.af2(60))
  tied <- select_ind(ph, n = 5, on = rep(1, 60))
  expect_identical(attr(tied, "selected"), ph$ids[1:5])
  expect_equal(attr(tied, "differential"), 0)
  expect_equal(attr(tied, "intensity"), 0)
  score <- stats::setNames(as.numeric(seq_len(60)), ph$ids)
  lo <- select_ind(ph, n = 5, on = score, direction = "low")
  expect_setequal(attr(lo, "selected"), ph$ids[1:5])
  expect_equal(attr(lo, "differential"), mean(1:5) - mean(1:60))
  hi <- select_ind(ph, n = 5, on = score)
  expect_equal(attr(hi, "intensity"),
               (mean(56:60) - mean(1:60)) / stats::sd(1:60), tolerance = 1e-12)
})

# ---- equation-level pins ---------------------------------------------------

test_that("the QGSI score equals w'y + y'Wy and reduces to w'y at W = 0 (Codex)", {
  ph2 <- .aph(.af2(60), n_traits = 2)
  Y <- .breeding_value_matrix(ph2, 1)
  w <- c(1, -0.5)
  W <- matrix(c(0.4, 0.15, 0.15, 0.2), 2)
  sc <- .quadratic_index_score(ph2, w, W, 1)
  expect_equal(unname(sc), as.numeric(Y %*% w + rowSums((Y %*% W) * Y)),
               tolerance = 1e-12)
  sc0 <- .quadratic_index_score(ph2, w, NULL, 1)
  expect_equal(unname(sc0), as.numeric(Y %*% w), tolerance = 1e-12)
})

test_that(".combined_score equals the direct solve of V b = c (Codex)", {
  for (n in c(2L, 3L, 10L)) {
    for (pr in list(c(0.5, 0.5), c(0.3, 0.25), c(0.8, 0.5))) {
      h2 <- pr[1]; r <- pr[2]; t <- r * h2
      set.seed(n)
      fam <- rep(c("a", "b"), each = n)
      y <- stats::rnorm(2 * n)
      sc <- .combined_score(y, fam, h2, r)
      m <- (1 + (n - 1) * t) / n
      V <- matrix(c(1, m, m, m), 2)
      cc <- h2 * c(1, (1 + (n - 1) * r) / n)
      b <- solve(V, cc)
      dev <- y - mean(y)
      fm <- tapply(dev, fam, mean)
      expect_equal(unname(sc), as.numeric(b[1] * dev + b[2] * fm[fam]),
                   tolerance = 1e-10)
    }
  }
  # pinned fixture: n = 2, h2 = r = 0.5 -> b = (1/3, 4/15)
  y <- c(1, 0, 0, 0)
  fam <- c("a", "a", "b", "b")
  sc <- .combined_score(y, fam, 0.5, 0.5)
  dev <- y - mean(y)
  fm <- tapply(dev, fam, mean)
  expect_equal(unname(sc), as.numeric(dev / 3 + 4 / 15 * fm[fam]),
               tolerance = 1e-12)
})

test_that(".sel_within_family allocates exactly keep_n by largest remainder (Codex)", {
  fam <- rep(c("s", "d", "l"), c(1, 2, 7))
  score <- c(5, 1, 9, 3, 8, 2, 7, 6, 4, 10)
  idx <- .sel_within_family(score, fam, 3L)
  expect_length(idx, 3L)
  expect_equal(as.integer(table(factor(fam[idx], levels = c("s", "d", "l")))),
               c(0L, 1L, 2L))
  # the chosen ones are the best inside their family
  expect_true(all(score[idx] %in% c(9, 10, 8)))
})

test_that("the response to mass selection follows R = Cov(A, P) / V_P * S (T1)", {
  ph <- .aph(.af2(300, seed = 5))
  sel <- select_ind(ph, prop = 0.2)
  A <- stats::setNames(.breeding_value_matrix(ph, 1)[, 1], ph$ids)
  P <- .criterion_values(ph, "pheno", 1, 1)
  S <- attr(sel, "differential")
  R_obs <- mean(A[attr(sel, "selected")]) - mean(A)
  R_hat <- stats::cov(A, P) / stats::var(P) * S
  expect_equal(R_obs / R_hat, 1, tolerance = 0.25)
})

# ---- scheme internals (Codex) ---------------------------------------------

test_that("tandem schedules and cross counts come out as documented (Codex)", {
  f2 <- .af2(60)
  out <- suppressMessages(pedigree(f2, .apheno2, generations = 3, prop = 0.2,
                                   trait = c(1, 2), seed = 1))
  expect_equal(attr(out, "history")$trait, c(1, 2, 1))
  rs <- suppressMessages(recurrent_selection(f2, .apheno2, cycles = 3,
                                             n_parents = 5, n_crosses = 3,
                                             progeny_per_cross = 4,
                                             trait = c(1, 2), seed = 1))
  expect_equal(n_individuals(rs), 12L)
  expect_equal(attr(rs, "history")$trait, c(1, 2, 1))
})

test_that(".self_each / .intermate / .relabel produce the documented counts and links (Codex)", {
  f2 <- .af2(20)
  kids <- suppressMessages(.self_each(f2[1:3], n_each = 2L, tag = "t"))
  expect_equal(n_individuals(kids), 6L)
  expect_false(anyDuplicated(kids$ids) > 0)
  # per-parent counts (used by bulk): zero-count parents are skipped
  k2 <- suppressMessages(.self_each(f2[1:3], n_each = c(2L, 0L, 3L), tag = "t"))
  expect_equal(n_individuals(k2), 5L)
  pp <- parentage(k2)
  expect_equal(sort(as.integer(table(pp$mother_key))), c(2L, 3L))
  # intermating: exact count, and each cross draws two distinct parents
  x <- suppressMessages(.intermate(f2[1:4], n_crosses = 3L,
                                   progeny_per_cross = 4L, tag = "c"))
  expect_equal(n_individuals(x), 12L)
  # relabelling to colliding display ids keeps the parent links (keys)
  keys_before <- parentage(kids)$mother_key
  rl <- .relabel(kids, rep(c("dup_1", "dup_2", "dup_3"), 2))
  expect_identical(parentage(rl)$mother_key, keys_before)
})

test_that(".as_founder_pop advances only the individuals a subset sim holds (T17)", {
  f2 <- .af2(40)
  sub <- suppressMessages(
    simulate_phenotype(f2, individuals = 1:10, h2 = 0.5, seed = 3) |>
      additive(n_qtn = 10))
  expect_equal(n_individuals(.as_founder_pop(sub)), 10L)
  out <- suppressMessages(single_seed_descent(sub, generations = 1, seed = 1))
  expect_equal(n_individuals(out), 10L)
})
