# test-nonadditive-cor.R
#
# DECISION-023: under architecture = "pleiotropy", dominance and epistasis
# effects are drawn with the same PleioArch covariance as the additive layer,
# so every component targets `cor` (the realized correlation converges as the
# units and individuals grow, for unlinked loci as in this HWE panel), and the total genetic correlation
# targets `cor` when the layers' per-trait `prop` profiles are proportional --
# as with the scalar `prop` and 60-120 units used here, where a few seeds suffice.
# Before this, both layers gave every trait one identical effect series, so the
# non-additive components realized a correlation of ~ +1 whatever `cor` was
# (audit 2026-09-17: target cor = 0 -> realized 1.0 for epistasis, 0.45 for AD).
#
# An outbred HWE panel is used: the bundled maize lines are ~0.4% heterozygous,
# too few to realize dominance at all on most loci.

.outbred <- function(n = 400, m = 1200, seed = 99) {
  set.seed(seed)
  p <- stats::runif(m, 0.15, 0.5)
  g <- t(vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1L, integer(n)))
  colnames(g) <- paste0("i", seq_len(n))
  cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                   chr = rep(1:10, each = m / 10),
                   pos = rep(seq_len(m / 10), 10) * 1e5, cm = NA_real_,
                   stringsAsFactors = FALSE),
        as.data.frame(g))
}
OB <- .outbred()

.pleio <- function(r, s, ...) {
  simulate_phenotype(OB, architecture = "pleiotropy", n_traits = 2, cor = r,
                     seed = s, ...)
}
.comp_cor <- function(sim, i) {
  ly <- sim$layers[[i]]
  stats::cor(.component_raw(ly, sim, 1L, 1L), .component_raw(ly, sim, 2L, 1L))
}
.total_cor <- function(sim) {
  g <- .genetic_value_matrix(sim, 1L)
  stats::cor(g[, 1], g[, 2])
}
.mean_over <- function(seeds, f) {
  mean(vapply(seeds, function(s) suppressMessages(suppressWarnings(f(s))),
              numeric(1)))
}

test_that("dominance under pleiotropy tracks cor (no longer ~ +1)", {
  for (r in c(0, 0.5)) {
    m <- .mean_over(1:4, function(s) {
      sim <- .pleio(r, s) |> additive(prop = 0.3, n_qtn = 120) |>
        dominance(prop = 0.2)
      .comp_cor(sim, 2L)
    })
    expect_lt(abs(m - r), 0.15)
  }
})

test_that("epistasis under pleiotropy tracks cor (no longer ~ +1)", {
  for (r in c(0, 0.5)) {
    m <- .mean_over(1:4, function(s) {
      sim <- .pleio(r, s) |> additive(prop = 0.3, n_qtn = 120) |>
        epistasis(prop = 0.2, n_pairs = 120)
      .comp_cor(sim, 2L)
    })
    expect_lt(abs(m - r), 0.15)
  }
})

test_that("the TOTAL genetic correlation of an A+D+E model tracks cor", {
  for (r in c(-0.5, 0.5)) {
    m <- .mean_over(1:4, function(s) {
      .total_cor(.pleio(r, s, pi = 0.7) |> additive(prop = 0.25, n_qtn = 120) |>
                   dominance(prop = 0.15) |> epistasis(prop = 0.15, n_pairs = 120))
    })
    expect_lt(abs(m - r), 0.15)
  }
})

test_that("one-call model = 'AD' under pleiotropy tracks cor", {
  m <- .mean_over(1:4, function(s) {
    .total_cor(simulate_phenotype(OB, architecture = "pleiotropy", n_traits = 2,
                                  cor = 0, h2 = 0.5, n_qtn = 120, model = "AD",
                                  seed = s))
  })
  expect_lt(abs(m), 0.15)
})

test_that("non-additive pleiotropic effects differ across traits and are seeded", {
  mk <- function() suppressMessages(
    .pleio(0, 7) |> additive(prop = 0.3, n_qtn = 20) |>
      dominance(prop = 0.1) |> epistasis(prop = 0.1, n_pairs = 10))
  a <- mk(); b <- mk()
  for (i in 2:3) {
    expect_false(isTRUE(all.equal(a$layers[[i]]$effect[[1]],
                                  a$layers[[i]]$effect[[2]])))
  }
  expect_identical(a$layers, b$layers)
})

test_that("pi partitions epistatic sets into shared + trait-specific", {
  sim <- suppressMessages(.pleio(0.3, 3, pi = 0.5) |>
    additive(prop = 0.3, n_qtn = 20) |> epistasis(prop = 0.2, n_pairs = 20))
  q <- sim$layers[[2]]$qtn
  key <- function(m) apply(m, 1, paste, collapse = "-")
  shared <- intersect(key(q[[1]]), key(q[[2]]))
  expect_equal(length(shared), 10L)                     # round(0.5 * 20)
  n_sets <- length(unique(key(rbind(q[[1]], q[[2]]))))  # 10 shared + 2 x 10
  expect_equal(n_sets, 30L)
  expect_length(unique(c(q[[1]], q[[2]])), 2L * n_sets) # loci all disjoint
})

test_that("fixed loci are validated, not rejected, under pleiotropy and ld (audit P4, DECISION-043)", {
  base_p <- .pleio(0.5, 1, pi = 0.5) |> additive(prop = 0.3, n_qtn = 10)
  # a vector is shared by every trait; with pi < 1 there are no trait-specific
  # loci to carry the specific variance (partial pleiotropy is complex_phenotypes())
  expect_error(dominance(base_p, prop = 0.1, qtn = c(1, 2)), "pi = 1")
  expect_error(epistasis(base_p, prop = 0.1, qtn = matrix(c(1, 2), 1)), "pi = 1")
  # with pi = 1 (the default) every locus is shared, which is valid
  all_shared <- .pleio(0.5, 1) |> additive(prop = 0.3, n_qtn = 10)
  expect_no_error(suppressWarnings(
    dominance(all_shared, prop = 0.1, qtn = c(1, 2), same_as_add = FALSE)))
  data("SNP55K_maize282_maf04", envir = environment())
  base_ld <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "ld",
                                n_traits = 2, seed = 1) |>
    additive(prop = 0.3, n_qtn = 3)
  expect_error(epistasis(base_ld, prop = 0.1, qtn = matrix(c(1, 2), 1)),
               "not supported under architecture = \"ld\"")
  # a vector would put the same marker on both traits of a linked pair
  expect_error(dominance(base_ld, prop = 0.1, qtn = c(1, 2)),
               "cannot be causal for both")
})

test_that("explicit effects / dist are rejected under pleiotropy", {
  base_p <- .pleio(0.5, 1) |> additive(prop = 0.3, n_qtn = 10)
  expect_error(epistasis(base_p, prop = 0.1, n_pairs = 3, effect = c(1, 1, 1)),
               "multivariate draw")
})

test_that("an unattainable cor errors on the non-additive path, naming it (O4)", {
  # cor^2 must not exceed pi_1 * pi_2 (0.6^2 > 0.3 * 0.3); no additive layer, so
  # this reaches the non-additive feasibility check, not additive()'s.
  expect_error(.pleio(0.6, 1, pi = 0.3) |> epistasis(prop = 0.2, n_pairs = 20),
               "epistasis\\(\\): Biological constraint")
})

# --- independent-review findings (2026-09-25) --------------------------------

.warnings_of <- function(expr) {
  w <- character(0)
  withCallingHandlers(suppressMessages(expr),
                      warning = function(c) {
                        w <<- c(w, conditionMessage(c))
                        invokeRestart("muffleWarning")
                      })
  w
}

test_that("non-proportional per-trait prop: total is attenuated and warned (O1)", {
  # Each component targets cor = 0.6, but the total targets
  # 0.6 * (sqrt(.4*.1) + sqrt(.1*.4)) / sqrt(.5*.5) = 0.48, not 0.6.
  mk <- function(s) .pleio(0.6, s) |> additive(prop = c(0.4, 0.1), n_qtn = 120) |>
    dominance(prop = c(0.1, 0.4))
  w <- .warnings_of(mk(1))
  expect_true(any(grepl("TOTAL genetic correlation", w) & grepl("0\\.480", w)))
  m <- .mean_over(1:4, function(s) .total_cor(mk(s)))
  expect_lt(abs(m - 0.48), 0.15)
  # scalar prop: proportional profiles, no such warning
  w2 <- .warnings_of(.pleio(0.6, 1) |> additive(prop = 0.3, n_qtn = 20) |>
                       dominance(prop = 0.2))
  expect_false(any(grepl("TOTAL genetic correlation", w2)))
})

test_that("constant units are left out of the variance allocation (O2)", {
  # 10 shared units (3 constant) + 10 trait-specific per trait, cor 0.4, pi 0.5.
  # Counting the 3 dead shared units would give 0.4*0.7/(0.5*0.7+0.5) = 0.33.
  set.seed(11)
  n <- 400
  Z <- matrix(stats::rnorm(n * 30), n, 30)
  Z[, 1:3] <- 1
  q <- list(c(1:10, 11:20), c(1:10, 21:30))
  pi_vec <- c(0.5, 0.5); vg <- c(1, 1)
  sigma <- matrix(0.4, 2, 2); diag(sigma) <- pi_vec * vg
  e1 <- .pleio_unit_effects(q, sigma, pi_vec, vg, function(j) Z[, j])
  expect_true(all(e1[[1]][1:3] == 0) && all(e1[[2]][1:3] == 0))
  r <- vapply(1:800, function(s) {
    set.seed(s)
    e <- .pleio_unit_effects(q, sigma, pi_vec, vg, function(j) Z[, j])
    stats::cor(Z[, q[[1]]] %*% e[[1]], Z[, q[[2]]] %*% e[[2]])[1]
  }, numeric(1))
  expect_lt(abs(mean(r) - 0.4), 0.04)
  # every shared unit constant: the correlation cannot be carried -> warn
  Z[, 1:10] <- 1
  expect_warning(.pleio_unit_effects(q, sigma, pi_vec, vg, function(j) Z[, j]),
                 "no longer controlled by `cor`")        # not "near 0" (round 7)
})

test_that("'ld' keeps every causal locus trait-specific (rounds 2-3, O2/O3)", {
  # Epistatic sets have no linked-distinct construction (outside the r2 window,
  # and able to reuse another layer's causal marker for the other trait), and a
  # fresh dominance draw can collide with the additive loci across layers.
  data("SNP55K_maize282_maf04", envir = environment())
  ld0 <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "ld",
                            n_traits = 2, seed = 8)
  expect_error(epistasis(ld0, prop = 0.2, n_pairs = 3),
               "not supported under architecture = \"ld\"")
  ld1 <- additive(ld0, prop = 0.4, n_qtn = 5)
  expect_error(dominance(ld1, prop = 0.1, same_as_add = FALSE, n_qtn = 3),
               "same_as_add = FALSE\\) is not supported")
  # dominance on the additive linked loci stays allowed: the same loci per
  # trait, and the two traits' loci disjoint
  d <- suppressWarnings(dominance(ld1, prop = 0.1))
  q <- d$layers[[2]]$qtn
  expect_identical(q, d$layers[[1]]$qtn)
  expect_length(intersect(q[[1]], q[[2]]), 0L)
})

test_that("one informative shared unit warns even with reused loci (round-4 O3)", {
  set.seed(3)
  n <- 200
  Z <- matrix(stats::rnorm(n * 6), n, 6)
  Z[, 1] <- 1                                   # shared unit 1 is constant
  q <- list(c(1L, 2L, 3L, 4L), c(1L, 2L, 5L, 6L))   # shared units: 1, 2
  pi_vec <- c(0.5, 0.5); vg <- c(1, 1)
  sigma <- matrix(0.2, 2, 2); diag(sigma) <- pi_vec * vg
  R <- matrix(c(1, 0.2, 0.2, 1), 2)
  ucol <- function(j) Z[, j]
  expect_warning(.pleio_unit_effects(q, sigma, pi_vec, vg, ucol, R = R,
                                     fresh = FALSE),
                 "only one informative shared")
  expect_warning(.pleio_unit_effects(q, sigma, pi_vec, vg, ucol, R = R,
                                     fresh = TRUE),       # 2 nominal -> 1 live
                 "only one informative shared")
  Z[, 1] <- stats::rnorm(n)                      # both shared units live
  expect_silent(.pleio_unit_effects(q, sigma, pi_vec, vg, ucol, R = R,
                                    fresh = FALSE))
})

test_that("single shared unit warning is accurate with and without specifics (O3)", {
  w <- .warnings_of(.pleio(0.2, 1, pi = 0.25) |> additive(prop = 0.3, n_qtn = 4))
  hit <- w[grepl("only one shared", w)]
  expect_length(hit, 1L)
  expect_match(hit, "noisy draw")
  expect_false(grepl("exactly \\+/-1", hit))
  w1 <- .warnings_of(.pleio(0.5, 1) |> additive(prop = 0.3, n_qtn = 1))
  expect_true(any(grepl("only one shared", w1) & grepl("exactly \\+/-1", w1)))
})

test_that("cor = 0 is covered by the single-shared-unit guards (round-5 O1)", {
  # one shared unit cannot realize 0 any more than 0.5: it gives +/-1
  sim <- .pleio(0, 17) |> additive(prop = 0.3, n_qtn = 60)
  w <- .warnings_of(sim |> epistasis(prop = 0.2, n_pairs = 1))
  expect_true(any(grepl("only one shared", w) & grepl("exactly \\+/-1", w)))
  w0 <- .warnings_of(.pleio(0, 1) |> additive(prop = 0.3, n_qtn = 1))
  expect_true(any(grepl("only one shared", w0)))
  # reused / degenerate path: one live shared unit, target 0
  set.seed(3)
  Z <- matrix(stats::rnorm(200 * 6), 200, 6)
  Z[, 1] <- 1
  q <- list(c(1L, 2L, 3L, 4L), c(1L, 2L, 5L, 6L))
  sigma <- diag(c(0.5, 0.5))
  expect_warning(.pleio_unit_effects(q, sigma, c(0.5, 0.5), c(1, 1),
                                     function(j) Z[, j], R = diag(2),
                                     fresh = FALSE),
                 "only one informative shared")
  # an exact +/-1 target is realizable by one unit: no warning
  expect_false(any(grepl("only one shared",
                         .warnings_of(.pleio(1, 1) |>
                                        additive(prop = 0.3, n_qtn = 1)))))
})

test_that("lost trait-specific units warn of inflation in magnitude (round-5 O3)", {
  set.seed(5)
  n <- 300
  Z <- matrix(stats::rnorm(n * 80), n, 80)
  Z[, 61:80] <- 1                               # every specific unit constant
  q <- list(c(1:60, 61:70), c(1:60, 71:80))
  pi_vec <- c(0.5, 0.5); vg <- c(1, 1)
  sigma <- matrix(-0.4, 2, 2); diag(sigma) <- pi_vec * vg
  w <- character(0)
  e <- withCallingHandlers(
    .pleio_unit_effects(q, sigma, pi_vec, vg, function(j) Z[, j]),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("inflated in magnitude", w)))
  expect_false(any(grepl("above `cor`", w)))
  r <- stats::cor(Z[, q[[1]]] %*% e[[1]], Z[, q[[2]]] %*% e[[2]])[1]
  expect_lt(r, -0.4)                            # further from 0, not above
})

test_that("'ld' allows a single additive layer (round-7 O2)", {
  base <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "ld",
                             n_traits = 2, seed = 8) |>
    additive(prop = 0.2, n_qtn = 3)
  expect_error(additive(base, prop = 0.2, n_qtn = 3),
               "supports a single additive layer")
})

test_that("a derived transcriptome layer under pleiotropy / ld warns (round-7 O1, round-9)", {
  small <- .outbred(n = 150, m = 200, seed = 5)
  sim <- suppressMessages(simulate_phenotype(small, architecture = "pleiotropy",
                                             n_traits = 2, cor = 0.3,
                                             transcriptome = TRUE, seed = 3))
  expect_warning(transcriptome(sim, prop = 0.2, n_genes = 5),
                 "not targeted at `cor`")
  ld <- suppressMessages(simulate_phenotype(small, architecture = "ld",
                                            n_traits = 2, transcriptome = TRUE,
                                            seed = 3))
  expect_warning(transcriptome(ld, prop = 0.2, n_genes = 5),
                 "no longer comes solely from distinct linked loci")
  indep <- suppressMessages(simulate_phenotype(small, n_traits = 2,
                                               transcriptome = TRUE, seed = 3))
  expect_no_warning(transcriptome(indep, prop = 0.2, n_genes = 5))
})

test_that("marker-shortage advice points pi the right way (round-8 O4)", {
  # distinct markers needed = nt * n - (nt - 1) * n_shared: raising pi lowers it
  tiny <- .outbred(n = 60, m = 20, seed = 4)
  expect_error(suppressMessages(
    simulate_phenotype(tiny, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.3, pi = 0.5, seed = 1) |>
      additive(prop = 0.3, n_qtn = 15)), "raise pi")
  expect_no_error(suppressMessages(
    simulate_phenotype(tiny, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.3, pi = 0.8, seed = 1) |>
      additive(prop = 0.3, n_qtn = 15)))
})

test_that("multi-trait feasibility does not depend on the variances (round-10 O1)", {
  R <- diag(3); R[1, 2] <- R[2, 1] <- 0.5           # 0.5^2 > 0.1 * 0.1
  pi_vec <- c(0.1, 0.1, 1)
  sig <- function(vg) {
    s <- outer(sqrt(vg), sqrt(vg)) * R; diag(s) <- pi_vec * vg; s
  }
  expect_error(.check_pleio_feasible(sig(c(1e-16, 1e-16, 0.5)), R, pi_vec),
               "not attainable")
  expect_error(.check_pleio_feasible(sig(c(0.3, 0.3, 0.3)), R, pi_vec),
               "not attainable")
  R[1, 2] <- R[2, 1] <- 0.1                          # boundary: 0.1^2 = 0.1*0.1
  expect_silent(.check_pleio_feasible(sig(c(1e-16, 1e-16, 0.5)), R, pi_vec))
})

test_that("a zero-prop transcriptome layer does not warn (round-10 O2)", {
  small <- .outbred(n = 150, m = 200, seed = 5)
  sim <- suppressMessages(simulate_phenotype(small, architecture = "pleiotropy",
                                             n_traits = 2, cor = 0.3,
                                             transcriptome = TRUE, seed = 3))
  expect_no_warning(transcriptome(sim, prop = 0, n_genes = 5))
  ld <- suppressMessages(simulate_phenotype(small, architecture = "ld",
                                            n_traits = 2, transcriptome = TRUE,
                                            seed = 3))
  expect_warning(transcriptome(ld, prop = c(0.3, 0), n_genes = 5),  # round 11
                 "no longer comes solely from distinct linked loci")
})

test_that("no inflation warning when the target correlation is 0 (round-11 O2)", {
  set.seed(5)
  Z <- matrix(stats::rnorm(300 * 80), 300, 80)
  Z[, 61:80] <- 1
  q <- list(c(1:60, 61:70), c(1:60, 71:80))
  sigma <- diag(c(0.5, 0.5))                          # cor = 0
  w <- character(0)
  withCallingHandlers(
    .pleio_unit_effects(q, sigma, c(0.5, 0.5), c(1, 1), function(j) Z[, j]),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("inflated", w)))
})

test_that("single-unit consequence is classified per trait pair (round-12 O1)", {
  R <- matrix(0.5, 3, 3); diag(R) <- 1
  txt <- .pleio_single_unit_consequence(R, c(TRUE, TRUE, FALSE))
  expect_match(txt, "traits 1-2: with no trait-specific variance .* exactly \\+/-1")
  expect_match(txt, "traits 1-3, 2-3: .* one noisy draw")
  # end to end: pi = (1, 1, 0.1) with one shared QTN
  sim3 <- simulate_phenotype(OB, architecture = "pleiotropy", n_traits = 3,
                             cor = R * 0 + diag(3) + 0.05 * (1 - diag(3)),
                             pi = c(1, 1, 0.1), seed = 2)
  w <- .warnings_of(sim3 |> additive(prop = 0.3, n_qtn = 2))
  hit <- w[grepl("only one shared", w)]
  expect_true(any(grepl("traits 1-2: .*exactly \\+/-1", hit)))
})

test_that("non-proportional total warns on relative attenuation (round-12 O2)", {
  w <- .warnings_of(.pleio(0.01, 1) |> additive(prop = c(0.49, 0.01), n_qtn = 60) |>
                      dominance(prop = c(0.01, 0.49)))
  expect_true(any(grepl("0\\.003 \\(cor 0\\.010\\)", w)))
})

test_that("zero-variance pairs get no single-unit +/-1 claim (round-13 O1)", {
  w <- .warnings_of(.pleio(0, 17, pi = 1) |> additive(prop = 0, n_qtn = 1))
  expect_false(any(grepl("only one shared", w)))
  expect_true(any(grepl("zero additive variance", w)))
  R <- matrix(0.3, 3, 3); diag(R) <- 1
  txt <- .pleio_single_unit_consequence(R, rep(TRUE, 3), vg = c(0.2, 0, 0.2))
  expect_match(txt, "^traits 1-3:")
  expect_identical(.pleio_single_unit_consequence(R, rep(TRUE, 3),
                                                  vg = c(0, 0, 0.2)), "")
})

test_that("dead shared units at cor = 0 report the variance, not the correlation (round-15 O1)", {
  set.seed(9)
  Z <- matrix(stats::rnorm(200 * 30), 200, 30)
  Z[, 1:10] <- 1                                     # every shared unit constant
  q <- list(c(1:10, 11:20), c(1:10, 21:30))
  expect_warning(.pleio_unit_effects(q, diag(c(0.5, 0.5)), c(0.5, 0.5), c(1, 1),
                                     function(j) Z[, j]),
                 "requested cor = 0 is still targeted")
  sig <- matrix(0.3, 2, 2); diag(sig) <- 0.5
  expect_warning(.pleio_unit_effects(q, sig, c(0.5, 0.5), c(1, 1),
                                     function(j) Z[, j]),
                 "no longer controlled by `cor`")
})

test_that("the total check uses each component's effective target (round-16 O2)", {
  # dead shared units: the component targets 0; dead specifics: Sigma12/(pi V)
  set.seed(9)
  Z <- matrix(stats::rnorm(200 * 30), 200, 30)
  q <- list(c(1:10, 11:20), c(1:10, 21:30))
  sig <- matrix(0.2, 2, 2); diag(sig) <- 0.5
  Zd <- Z; Zd[, 1:10] <- 1
  e <- suppressWarnings(.pleio_unit_effects(q, sig, c(0.5, 0.5), c(1, 1),
                                            function(j) Zd[, j]))
  expect_equal(attr(e, "target_cor")[1, 2], 0)
  Zs <- Z; Zs[, 11:30] <- 1
  e2 <- suppressWarnings(.pleio_unit_effects(q, sig, c(0.5, 0.5), c(1, 1),
                                             function(j) Zs[, j]))
  expect_equal(attr(e2, "target_cor")[1, 2], 0.4)
  # additive c(.4, .1) at cor .4 + a component targeting 0 with c(.1, .4):
  # total target .4 * sqrt(.04) / .5 = 0.160, not 0.320
  sim <- suppressWarnings(.pleio(0.4, 5, pi = 0.5) |>
                            additive(prop = c(0.4, 0.1), n_qtn = 60))
  w <- character(0)
  withCallingHandlers(.pleio_total_cor_check(sim, c(0.1, 0.4), diag(2)),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("0\\.160 \\(cor 0\\.400\\)", w)))
})

test_that("constant specifics turn a nominal noisy single unit into +/-1 (round-16 O3)", {
  set.seed(4)
  Z <- matrix(stats::rnorm(200 * 5), 200, 5)
  Z[, 2:5] <- 1                                   # specifics constant
  q <- list(c(1L, 2L, 3L), c(1L, 4L, 5L))         # one shared unit (fresh)
  w <- character(0)
  withCallingHandlers(
    .pleio_unit_effects(q, diag(c(0.5, 0.5)), c(0.5, 0.5), c(1, 1),
                        function(j) Z[, j], R = diag(2), fresh = TRUE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("only one informative shared", w) &
                    grepl("exactly \\+/-1", w)))
})

test_that("no transcriptome warning without genetic expression variance (round-16 O4)", {
  small <- .outbred(n = 150, m = 200, seed = 5)
  tx0 <- simulate_transcriptome(small, h2 = 0, seed = 2)
  sim <- suppressMessages(simulate_phenotype(small, architecture = "pleiotropy",
                                             n_traits = 2, cor = 0.3,
                                             transcriptome = tx0, seed = 3))
  expect_no_warning(transcriptome(sim, prop = 0.2, n_genes = 5))
})

test_that("the transcriptome check covers every vary_qtn replication (round-17 O1)", {
  small <- .outbred(n = 150, m = 200, seed = 5)
  tx <- simulate_transcriptome(small, seed = 2)
  live <- 1L                                    # only gene 1 is genome-mediated
  tx$genetic_expression[-live, ] <- 0
  mk <- function(s) suppressMessages(
    simulate_phenotype(small, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.3, transcriptome = tx, seed = s,
                       vary_qtn = TRUE, n_reps = 6))
  # find a seed whose canonical draw misses gene 1 but a replication hits it
  hit <- NULL
  for (s in 1:200) {
    ly <- suppressWarnings(transcriptome(mk(s), prop = 0.2, n_genes = 3))
    ly <- ly$layers[[length(ly$layers)]]
    canon <- unlist(ly$qtn); reps <- unlist(ly$qtn_reps)
    if (!live %in% canon && live %in% reps) { hit <- s; break }
  }
  skip_if(is.null(hit), "no seed with gene 1 only in a replication")
  expect_warning(transcriptome(mk(hit), prop = 0.2, n_genes = 3),
                 "not targeted at `cor`")
})

test_that("the total check reports replications whose target differs (round-18 O2)", {
  sim <- suppressWarnings(.pleio(0.4, 5, pi = 0.5) |>
                            additive(prop = c(0.4, 0.1), n_qtn = 60))
  R <- matrix(c(1, 0.4, 0.4, 1), 2)
  w <- character(0)
  withCallingHandlers(
    .pleio_total_cor_check(sim, c(0.1, 0.4), R, list(R, R, diag(2))),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("0\\.320 \\(cor 0\\.400\\)", w)))
  expect_true(any(grepl("0\\.160 \\(cor 0\\.400\\) in replication 3", w)))
  expect_false(any(grepl("in replication 2", w)))       # same as canonical
})

test_that("zero-prop traits' genes do not trigger the transcriptome warning (round-18 O1)", {
  small <- .outbred(n = 150, m = 200, seed = 5)
  tx <- simulate_transcriptome(small, seed = 2)
  live <- 1L
  tx$genetic_expression[-live, ] <- 0
  mk <- function(s) suppressMessages(
    simulate_phenotype(small, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.3, transcriptome = tx, seed = s))
  hit <- NULL
  for (s in 1:300) {
    ly <- suppressWarnings(transcriptome(mk(s), prop = c(0.2, 0), n_genes = 3))
    q <- ly$layers[[length(ly$layers)]]$qtn
    if (!live %in% q[[1]] && live %in% q[[2]]) { hit <- s; break }
  }
  skip_if(is.null(hit), "no seed with the live gene only on the zero-prop trait")
  expect_no_warning(transcriptome(mk(hit), prop = c(0.2, 0), n_genes = 3))
})

test_that("zero-slope genes carrying no realized genetic signal do not warn (round-19 O2)", {
  small <- .outbred(n = 150, m = 200, seed = 5)
  tx <- simulate_transcriptome(small, seed = 2)
  tx$genetic_expression[-1L, ] <- 0             # only gene 1 is genome-mediated
  sim <- suppressMessages(simulate_phenotype(small, architecture = "pleiotropy",
                                             n_traits = 2, cor = 0.3,
                                             transcriptome = tx, seed = 3))
  expect_no_warning(transcriptome(sim, prop = 0.2, genes = 1:2, slopes = c(0, 1)))
  expect_warning(transcriptome(sim, prop = 0.2, genes = 1:2, slopes = c(1, 1)),
                 "not targeted at `cor`")
})

test_that("hetless dominance loci in any vary_qtn replication warn (round-19 O1)", {
  g <- .outbred(n = 200, m = 30, seed = 6)
  ind <- 6:ncol(g)
  for (k in 1:10) {                             # markers 1-10: no heterozygotes
    v <- unlist(g[k, ind]); v[v == 0] <- 1; g[k, ind] <- v
  }
  mk <- function(s) simulate_phenotype(g, n_traits = 1, seed = s, vary_qtn = TRUE,
                                       n_reps = 5) |>
    additive(prop = 0.3, n_qtn = 4)
  hl <- function(q) any(vapply(q, function(j) all(unlist(g[j, ind]) != 0),
                               logical(1)))
  hit <- NULL
  for (s in 1:200) {
    a <- suppressWarnings(mk(s))$layers[[1]]
    if (!hl(a$qtn[[1]]) && any(vapply(a$qtn_reps, function(q) hl(q[[1]]),
                                      logical(1)))) { hit <- s; break }
  }
  skip_if(is.null(hit), "no seed with hetless loci only in a replication")
  expect_warning(dominance(suppressWarnings(mk(hit)), prop = 0.2),
                 "some \\(but not all\\) selected loci have no heterozygous")
})

test_that("tiny correlations are not rounded to 0.000 in the total warning (round-21 O3)", {
  w <- .warnings_of(.pleio(1e-4, 1) |> additive(prop = c(0.49, 0.01), n_qtn = 60) |>
                      dominance(prop = c(0.01, 0.49)))
  hit <- w[grepl("TOTAL genetic correlation", w)]
  expect_length(hit, 1L)
  expect_match(hit, "2\\.8e-05 \\(cor 0\\.0001\\)")
})
