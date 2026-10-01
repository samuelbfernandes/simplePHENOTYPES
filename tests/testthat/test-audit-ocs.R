# test-audit-ocs.R
#
# Regression tests for the OCS / usefulness / marker-selection audit
# (reconciliation v2-ocs-usefulness-marker): g_matrix() on a subset simulation,
# unnamed G, above-maximum coancestry targets, lambda bracketing at extreme merit
# scales, argument validation, sample_parents() edge cases and RNG hygiene,
# cross_usefulness() input handling, plus exact hand-computed gates.

data("SNP55K_maize282_maf04")

.aud_f2 <- function(n = 60) {
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:40)
  f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
  suppressMessages(selfcross(f1, n = n, seed = 2))
}
.aud_ph <- function(f2, ...) {
  suppressMessages(simulate_phenotype(f2, h2 = 0.5, seed = 3, ...) |>
                     additive(n_qtn = 30))
}
# Two opposite-homozygous founders on `k` markers (P1 = +1, P2 = -1).
.aud_founders <- function(k = 4L) {
  g <- data.frame(snp = paste0("m", seq_len(k)), allele = "A/G",
                  chr = seq_len(k), pos = 1, cm = 0, P1 = 1L, P2 = -1L,
                  stringsAsFactors = FALSE)
  as_population(g)
}
.aud_G2 <- function() diag(2) |> `dimnames<-`(list(c("a", "b"), c("a", "b")))
.aud_ocs2 <- function(...) {
  optimum_contribution(.aud_founders(), merit = c(a = 0, b = 1), G = .aud_G2(),
                       ...)
}

# ---- OCS-F1 -----------------------------------------------------------------

test_that("g_matrix()/optimum_contribution() honour a subset phenotype_sim (OCS-F1)", {
  f2 <- .aud_f2(60)
  phs <- .aud_ph(f2, individuals = 1:20)
  G <- g_matrix(phs)
  expect_equal(dim(G), c(20L, 20L))
  expect_equal(rownames(G), f2$ids[1:20])
  expect_equal(G, g_matrix(f2[1:20]))
  oc <- optimum_contribution(phs, lambda = 1)
  expect_s3_class(oc, "ocs")
  expect_equal(length(oc$contributions), 20L)
  expect_equal(names(oc$contributions), f2$ids[1:20])
  # the full-population sim is unchanged
  expect_equal(dim(g_matrix(.aud_ph(f2))), c(60L, 60L))
})

# ---- OCS-A03 ----------------------------------------------------------------

test_that("an unnamed or inconsistently named G is rejected, never a corrupt ocs (OCS-A03)", {
  fnd <- .aud_founders()
  expect_error(
    optimum_contribution(fnd, merit = c(0, 1), G = diag(2), lambda = 1),
    "dimnames|row and column names")
  G <- diag(2); rownames(G) <- c("a", "b")
  expect_error(
    optimum_contribution(fnd, merit = c(a = 0, b = 1), G = G, lambda = 1),
    "row and column names")
  G <- diag(2); dimnames(G) <- list(c("a", "a"), c("a", "a"))
  expect_error(
    optimum_contribution(fnd, merit = c(0, 1), G = G, lambda = 1),
    "unique")
  # a valid object always has consistent contributions / parents / n_parents
  oc <- .aud_ocs2(lambda = 1)
  expect_identical(oc$n_parents, length(oc$parents))
})

# ---- OCS-F4 -----------------------------------------------------------------

test_that("a target above the unconstrained coancestry warns (OCS-F4)", {
  f2 <- .aud_f2(40)
  ph <- .aud_ph(f2)
  c0 <- optimum_contribution(ph, lambda = 0)
  expect_warning(
    optimum_contribution(ph, target_coancestry = 2 * c0$coancestry),
    "above")
  # max_coancestry is a ceiling: slack is fine and silent
  expect_no_warning(
    optimum_contribution(ph, max_coancestry = 2 * c0$coancestry))
  # a reachable target does not warn
  expect_no_warning(
    optimum_contribution(ph, target_coancestry = 0.8 * c0$coancestry))
})

# ---- OCS-A01 / closed form --------------------------------------------------

test_that("closed-form two-candidate OCS: G = I, g = (0, 1), lambda = 2 (Codex row 1)", {
  oc <- .aud_ocs2(lambda = 2)
  expect_equal(unname(oc$contributions), c(0.25, 0.75), tolerance = 1e-8)
  expect_equal(oc$merit, 0.75, tolerance = 1e-8)
  expect_equal(oc$coancestry, 0.3125, tolerance = 1e-8)
})

test_that("lambda tuning is invariant to the merit scale, with no false warning (OCS-A01)", {
  fnd <- .aud_founders()
  # exact answer for target 0.26: c = (0.4, 0.6) at every scale
  for (s in c(1, 1e6, 1e17, 1e18, 1e20, 1e100, 1e-12)) {
    expect_no_warning(
      oc <- optimum_contribution(fnd, merit = c(a = 0, b = s), G = .aud_G2(),
                                 target_coancestry = 0.26))
    expect_equal(unname(oc$contributions), c(0.4, 0.6), tolerance = 1e-6,
                 info = paste("spread", s))
    expect_equal(oc$coancestry, 0.26, tolerance = 1e-6, info = paste("spread", s))
    expect_equal(oc$lambda / s, 5, tolerance = 1e-5, info = paste("spread", s))
  }
  # a target that is genuinely below the minimum attainable still warns
  expect_warning(
    optimum_contribution(fnd, merit = c(a = 0, b = 1), G = .aud_G2(),
                         target_coancestry = 0.25 - 0.1),
    "minimum attainable")
})

# ---- OCS-A05 / A08 ----------------------------------------------------------

test_that("optimum_contribution() validates merit and control arguments (OCS-A05)", {
  fnd <- .aud_founders()
  G <- .aud_G2()
  expect_error(optimum_contribution(fnd, merit = c(a = NA, b = 1), G = G, lambda = 1),
               "merit")
  expect_error(optimum_contribution(fnd, merit = c(a = Inf, b = 1), G = G, lambda = 1),
               "merit")
  expect_error(.aud_ocs2(lambda = 1, min_contribution = NA_real_),
               "min_contribution")
  expect_error(.aud_ocs2(lambda = 1, min_contribution = -1), "min_contribution")
  expect_error(.aud_ocs2(lambda = 1, tol = NA_real_), "tol")
  expect_error(.aud_ocs2(lambda = 1, tol = 0), "tol")
  expect_error(.aud_ocs2(lambda = 1, max_iter = 0), "max_iter")
  expect_error(.aud_ocs2(lambda = 1, max_iter = 1.5), "max_iter")
  expect_error(.aud_ocs2(lambda = 1, max_iter = NA), "max_iter")
})

test_that("g_matrix() rejects a non-numeric dosage matrix (OCS-A08)", {
  d <- matrix(c("0", "1", "1", "0"), 2, dimnames = list(NULL, c("a", "b")))
  expect_error(g_matrix(d), "must be numeric")
})

test_that("marker_select(min_markers = NA) names the argument (OCS-A08)", {
  fnd <- .aud_founders()
  f2 <- suppressMessages(selfcross(cross(fnd[1], fnd[2], seed = 1), n = 20, seed = 2))
  expect_error(marker_select(f2, c("m1", "m2"), min_markers = NA_real_),
               "min_markers")
})

# ---- sample_parents ---------------------------------------------------------

test_that("sample_parents() works when one individual holds every contribution (OCS-F3)", {
  f2 <- .aud_f2(40)
  ids <- f2$ids
  for (k in c(1L, 7L, 40L)) {
    contr <- stats::setNames(rep(0, length(ids)), ids)
    contr[k] <- 1
    oc <- structure(list(contributions = contr, parents = ids[k], n_parents = 1L,
                         merit = 0, coancestry = 0, lambda = 0), class = "ocs")
    m <- sample_parents(oc, f2, n = 3, seed = 1)
    expect_equal(n_individuals(m), 3L)
    dm <- dosages(m)
    for (j in 1:3) {                       # every slot carries that individual
      expect_equal(unname(dm[, j]), unname(dosages(f2)[, k]))
    }
  }
})

test_that("sample_parents() slot frequencies follow the contributions (multinomial)", {
  fnd <- .aud_founders()
  oc <- structure(list(contributions = c(P1 = 0.8, P2 = 0.2), parents = c("P1", "P2"),
                       n_parents = 2L, merit = 0, coancestry = 0, lambda = 0),
                  class = "ocs")
  ids <- fnd$ids
  names(oc$contributions) <- ids
  n <- 20000
  m1 <- sample_parents(oc, fnd, n = n, seed = 11, method = "multinomial")
  m2 <- sample_parents(oc, fnd, n = n, seed = 11, method = "multinomial")
  expect_identical(m1$ids, m2$ids)                         # seed reproduces
  frac <- mean(startsWith(m1$ids, ids[1]))
  expect_lt(abs(frac - 0.8), 4 * sqrt(0.8 * 0.2 / n))
})

# ---- sample_parents: allocation (D5, owner decision 2026-09-29) -------------

test_that(".allocate_slots() is largest-remainder: sums to n, within one slot of n*c", {
  a <- simplePHENOTYPES:::.allocate_slots
  cc <- c(.50, .30, .10, .06, .04)
  expect_identical(a(cc, 10), c(5L, 3L, 1L, 1L, 0L))     # hand-computed
  expect_identical(a(cc, 7),  c(4L, 2L, 1L, 0L, 0L))     # fracs .5 .1 .7 .42 .28
  expect_identical(a(cc, 1),  c(1L, 0L, 0L, 0L, 0L))
  expect_identical(a(c(0, 1, 0), 4), c(0L, 4L, 0L))
  expect_identical(a(c(2, 6, 2), 10), c(2L, 6L, 2L))     # unnormalised input
  set.seed(1)
  for (k in 1:300) {                                     # property over random cases
    m <- sample(2:25, 1); n <- sample(1:120, 1)
    c0 <- stats::rexp(m); c0[sample(m, sample(0:(m - 1), 1))] <- 0
    if (sum(c0) == 0) next
    r <- a(c0, n); q <- n * c0 / sum(c0)
    expect_equal(sum(r), n)
    expect_true(all(r >= 0L))
    expect_true(all(abs(r - q) < 1 + 1e-9))
    expect_true(all(r[c0 == 0] == 0L))
  }
})

test_that("sample_parents() allocates slots to the optimized contributions (default)", {
  f2 <- .aud_f2(40)
  ph <- .aud_ph(f2)
  oc <- optimum_contribution(ph, merit = "bv", lambda = 5)
  n <- 40
  m <- sample_parents(oc, f2, n = n, seed = 4)
  expect_equal(n_individuals(m), n)
  got <- table(factor(m$keys, levels = f2$keys))
  want <- n * oc$contributions
  expect_true(all(abs(as.numeric(got) - as.numeric(want)) < 1))  # within one slot
  # counts do not depend on the seed (only the slot order does)
  m2 <- sample_parents(oc, f2, n = n, seed = 99)
  expect_identical(table(factor(m2$keys, levels = f2$keys)), got)
  expect_false(identical(m$keys, m2$keys))
  expect_false(anyDuplicated(m$ids) > 0)
})

test_that("allocation keeps the realized group coancestry at the optimum; multinomial drifts", {
  f2 <- .aud_f2(40)
  ph <- .aud_ph(f2)
  oc <- optimum_contribution(ph, merit = "bv", lambda = 5)
  Gm <- g_matrix(f2)
  co <- function(mates) {
    k <- as.numeric(table(factor(mates$keys, levels = f2$keys))) / n_individuals(mates)
    drop(0.5 * crossprod(k, Gm %*% k))
  }
  target <- 0.5 * drop(crossprod(oc$contributions, Gm %*% oc$contributions))
  n <- 30
  alloc <- co(sample_parents(oc, f2, n = n, seed = 1))
  multi <- vapply(1:40, function(s)
    co(sample_parents(oc, f2, n = n, seed = s, method = "multinomial")), numeric(1))
  expect_lt(abs(alloc - target), abs(mean(multi) - target) + 1e-12)
  expect_lt(abs(alloc - target), 0.2 * target + 0.02)
})

test_that("exact ties are broken at random, not toward the first-listed individual", {
  f2 <- .aud_f2(40)
  ids <- f2$ids[1:3]
  contr <- stats::setNames(rep(0, length(f2$ids)), f2$ids)
  contr[1:3] <- 1 / 3
  oc <- structure(list(contributions = contr, parents = ids, n_parents = 3L,
                       merit = 0, coancestry = 0, lambda = 0), class = "ocs")
  first_gets_four <- vapply(1:60, function(s) {
    m <- sample_parents(oc, f2, n = 10, seed = s)
    tab <- table(factor(m$keys, levels = f2$keys))[1:3]
    expect_identical(sort(as.integer(tab)), c(3L, 3L, 4L))
    tab[[1]] == 4L
  }, logical(1))
  expect_true(any(first_gets_four) && !all(first_gets_four))
})

test_that("method = 'multinomial' reproduces the previous draw for the same seed", {
  f2 <- .aud_f2(40)
  ph <- .aud_ph(f2)
  oc <- optimum_contribution(ph, merit = "bv", lambda = 5)
  m <- sample_parents(oc, f2, n = 12, seed = 4, method = "multinomial")
  set.seed(4)
  idx <- match(names(oc$contributions), f2$ids)
  drawn <- idx[sample.int(length(idx), size = 12, replace = TRUE,
                          prob = oc$contributions)]
  expect_identical(m$keys, f2$keys[drawn])
  expect_error(sample_parents(oc, f2, n = 4, method = "bogus"), "arg")
})

test_that("sample_parents() validates and restores the RNG (OCS-A07/F5)", {
  f2 <- .aud_f2(40)
  ph <- .aud_ph(f2)
  oc <- optimum_contribution(ph, lambda = 5)
  expect_error(sample_parents(oc, f2, n = 4, seed = 1.9), "seed")
  expect_error(sample_parents(oc, f2, n = 4, seed = NA_real_), "seed")
  set.seed(99); before <- stats::runif(1)
  set.seed(99)
  invisible(sample_parents(oc, f2, n = 4, seed = 5))
  expect_identical(stats::runif(1), before)               # ambient RNG restored
})

# ---- cross_usefulness -------------------------------------------------------

# (maize lines carry residual heterozygosity; the heterozygous-parent warning is
# tested explicitly below and muffled elsewhere)
.aud_cu <- function(...) suppressWarnings(cross_usefulness(...))
.aud_sim <- function() {
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:6)
  suppressMessages(simulate_phenotype(pop, h2 = 0.5, seed = 1) |>
                     additive(n_qtn = 20))
}

# Four fully homozygous lines on 24 unlinked-ish markers.
.aud_inbred_sim <- function() {
  set.seed(42)
  k <- 24L
  g <- data.frame(snp = paste0("m", seq_len(k)), allele = "A/G",
                  chr = rep(1:4, each = 6), pos = rep(1:6, 4) * 1e6,
                  cm = rep(seq(0, 50, length.out = 6), 4), stringsAsFactors = FALSE)
  for (l in paste0("L", 1:4)) g[[l]] <- sample(c(-1L, 1L), k, replace = TRUE)
  pop <- as_population(g)
  suppressMessages(simulate_phenotype(pop, h2 = 0.5, seed = 1) |>
                     additive(n_qtn = 12))
}

test_that("cross_usefulness() validates pairs and select_top (OCS-A06)", {
  sim <- .aud_sim()
  expect_error(cross_usefulness(sim, pairs = rbind(c(1.9, 2.1)), n_progeny = 10),
               "whole")
  expect_error(cross_usefulness(sim, pairs = matrix(numeric(0), 0, 2),
                                n_progeny = 10), "at least one")
  expect_error(cross_usefulness(sim, select_top = NA_real_, n_progeny = 10),
               "select_top")
  expect_error(cross_usefulness(sim, select_top = NaN, n_progeny = 10),
               "select_top")
})

test_that("cross_usefulness() validates the seed and restores the RNG (OCS-A07/F5)", {
  sim <- .aud_sim()
  pr <- rbind(c(1, 2))
  expect_error(cross_usefulness(sim, pairs = pr, n_progeny = 10, seed = 1.9),
               "seed")
  set.seed(99); before <- stats::runif(1)
  set.seed(99)
  invisible(.aud_cu(sim, pairs = pr, n_progeny = 10, seed = 5))
  expect_identical(stats::runif(1), before)
  a <- .aud_cu(sim, pairs = pr, n_progeny = 10, seed = 5)
  b <- .aud_cu(sim, pairs = pr, n_progeny = 10, seed = 5)
  expect_identical(a, b)
})

test_that("heterozygous parents with a dh/selfcross family warn once (USE-F1)", {
  f2 <- .aud_f2(20)                       # F2 parents are heterozygous
  sim <- .aud_ph(f2)
  pr <- rbind(c(1, 2), c(3, 4))
  expect_warning(
    cross_usefulness(sim, pairs = pr, scheme = "dh", n_progeny = 10, seed = 1),
    "heterozygous")
  w <- character()
  withCallingHandlers(
    cross_usefulness(sim, pairs = pr, scheme = "selfcross", n_progeny = 10,
                     generations = 2, seed = 1),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_equal(sum(grepl("heterozygous", w)), 1L)      # once per call, not per pair
  # inbred parents (and the single-cross scheme) do not warn
  sim0 <- .aud_inbred_sim()
  expect_no_warning(
    cross_usefulness(sim0, pairs = pr, scheme = "dh", n_progeny = 10, seed = 1))
  expect_no_warning(
    cross_usefulness(sim, pairs = pr, scheme = "cross", n_progeny = 10, seed = 1))
})

test_that("cross_usefulness() reproduces a manual family evaluation (Codex row 6)", {
  sim <- .aud_sim()
  pop <- simplePHENOTYPES:::.as_founder_pop(sim)
  model <- simplePHENOTYPES:::.additive_model(sim, 1L)
  set.seed(7)
  fam <- simplePHENOTYPES:::.make_family(pop[1], pop[2], "dh", 30L, 5L)
  gv <- simplePHENOTYPES:::.additive_gv(fam, model)
  u <- .aud_cu(sim, pairs = rbind(c(1, 2)), scheme = "dh",
               n_progeny = 30, seed = 7)
  i <- stats::dnorm(stats::qnorm(0.9)) / 0.1
  expect_equal(u$mean, mean(gv))
  expect_equal(u$sd, stats::sd(gv))
  expect_equal(u$usefulness, mean(gv) + i * stats::sd(gv))
  expect_equal(attr(u, "intensity"), i)
})

test_that(".additive_model() scores layers on the realized scale, summing repeated markers (Codex row 7)", {
  f2 <- .aud_f2(40)
  ph <- suppressMessages(
    simulate_phenotype(f2, h2 = 0.5, seed = 3) |>
      additive(n_qtn = 10, prop = 0.3) |> additive(n_qtn = 8, prop = 0.2))
  m <- simplePHENOTYPES:::.additive_model(ph, 1L)
  gv <- simplePHENOTYPES:::.additive_gv(f2, m)
  tot <- 0
  for (ly in ph$layers) {
    comp <- simplePHENOTYPES:::.component_raw(ly, ph, 1L, 1L)
    tot <- tot + comp / stats::sd(comp) * sqrt(ly$prop)
  }
  expect_equal(gv - mean(gv), tot - mean(tot), tolerance = 1e-8)
})

# ---- recurrent_parent_recovery ---------------------------------------------

test_that("recurrent_parent_recovery() warns when no marker is informative (OCS-A04)", {
  g <- data.frame(snp = paste0("m", 1:3), allele = "A/G", chr = 1:3, pos = 1,
                  cm = 0, P1 = 1L, P2 = 1L, stringsAsFactors = FALSE)
  fnd <- as_population(g)
  expect_warning(r <- recurrent_parent_recovery(fnd, fnd[1], fnd[2]),
                 "informative")
  expect_true(all(is.na(r)))
  expect_identical(attr(r, "n_markers"), 0L)
})

# ---- g_matrix exact gates ---------------------------------------------------

test_that("g_matrix() equals VanRaden's method 1 by hand (Codex row 11)", {
  # m1 = (-1, 0, 1), m2 = (1, 1, -1) over three individuals
  d <- rbind(m1 = c(-1, 0, 1), m2 = c(1, 1, -1))
  colnames(d) <- paste0("i", 1:3)
  Gm <- matrix(c(26/17, 8/17, -2, 8/17, 8/17, -16/17, -2, -16/17, 50/17), 3)
  expect_equal(unname(g_matrix(d)), Gm, tolerance = 1e-10)
  # fixed base: an all-+1 marker with base_freq .1 keeps every entry at 18
  d1 <- matrix(1, 1, 3, dimnames = list("m", paste0("i", 1:3)))
  expect_equal(unname(g_matrix(d1, base_freq = 0.1)), matrix(18, 3, 3),
               tolerance = 1e-10)
  # fully inbred / fully heterozygous individuals under base_freq = .5
  d2 <- rbind(c(1, 0, -1), c(1, 0, -1))
  colnames(d2) <- c("hom", "het", "hom2")
  Fi <- diag(g_matrix(d2, base_freq = c(0.5, 0.5))) - 1
  expect_equal(unname(Fi), c(1, -1, 1), tolerance = 1e-10)
})

# ---- usefulness analytic gates (Fable T8/T9) --------------------------------

test_that("usefulness sd/mean match the analytic family moments for inbred parents", {
  # Two inbred parents on 10 unlinked loci (one per chromosome), P1 = +1 everywhere,
  # P2 = -1 at loci 1..6 and +1 at 7..10: the six segregating loci each contribute
  # eff^2 to the DH variance (x = +/-1 with prob 1/2).
  k <- 10L
  g <- data.frame(snp = paste0("m", seq_len(k)), allele = "A/G", chr = seq_len(k),
                  pos = 1, cm = 0, P1 = 1L, P2 = c(rep(-1L, 6), rep(1L, 4)),
                  P3 = c(rep(1L, 6), rep(-1L, 4)),       # third line: >= 3 lines needed
                  stringsAsFactors = FALSE)
  pop <- as_population(g)
  sim <- suppressMessages(simulate_phenotype(pop, h2 = 0.5, seed = 1) |>
                            additive(n_qtn = k))
  mod <- simplePHENOTYPES:::.additive_model(sim, 1L)
  eff <- mod$eff[match(g$snp, mod$snp)]
  v_dh <- sum(eff[1:6]^2)
  pr <- rbind(c(1, 2))
  n <- 6000L
  u_dh <- cross_usefulness(sim, pairs = pr, scheme = "dh", n_progeny = n, seed = 3)
  expect_equal(u_dh$sd, sqrt(v_dh), tolerance = 0.05)
  # mean = mid-parent value
  mid <- sum(eff * (g$P1 + g$P2) / 2)
  expect_lt(abs(u_dh$mean - mid), 4 * sqrt(v_dh / n))
  # F2 (one selfing generation): variance is halved
  u_f2 <- cross_usefulness(sim, pairs = pr, scheme = "selfcross", generations = 1,
                           n_progeny = n, seed = 3)
  expect_equal(u_f2$sd, sqrt(v_dh / 2), tolerance = 0.05)
  # F6: variance (1 + F)/2 with F = 1 - 2^-5
  u_f6 <- cross_usefulness(sim, pairs = pr, scheme = "selfcross", generations = 5,
                           n_progeny = n, seed = 3)
  expect_equal(u_f6$sd, sqrt(v_dh * (1 + (1 - 2^-5)) / 2), tolerance = 0.05)
  # single cross of inbred parents: every progeny is the same F1
  u_f1 <- cross_usefulness(sim, pairs = pr, scheme = "cross", n_progeny = 50,
                           seed = 3)
  expect_equal(u_f1$sd, 0, tolerance = 1e-12)
})
