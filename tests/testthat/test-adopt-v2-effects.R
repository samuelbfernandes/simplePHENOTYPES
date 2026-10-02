# test-adopt-v2-effects.R
#
# Remaining test-gap proposals of the 2026-09-29 audit for the v2 effects and
# architecture engines (audit group v2-effects-arch: R/effects_series.R,
# R/effects_pleioarch.R, R/arch_independent.R, R/arch_ld.R, R/qc_ld_methods.R)
# that were not adopted into test-audit-grammar.R / test-nonadditive-cor.R and
# the older files. IDs: EFFECTS-F<n> = Fable table row n, EFFECTS-C<n> = Codex
# table row n (see the handoff for the coverage table).

data("SNP55K_maize282_maf04")
G_av <- SNP55K_maize282_maf04

# outbred HWE panel (heterozygotes present, loci in linkage equilibrium); the
# bundled maize lines are ~0.4% heterozygous
.av_outbred <- function(n = 400, m = 1200, seed = 99) {
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
OB_av <- .av_outbred()

# individuals x markers -1/0/1 HWE matrix (no map)
.av_matrix <- function(n = 200, m = 300, p_lo = 0.15, p_hi = 0.5, seed = 99) {
  set.seed(seed)
  p <- stats::runif(m, p_lo, p_hi)
  g <- vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1, numeric(n))
  dimnames(g) <- list(paste0("i", seq_len(n)), paste0("m", seq_len(m)))
  g
}

.av_pleio <- function(r, s, nt = 2L, ...) {
  suppressMessages(
    simulate_phenotype(OB_av, architecture = "pleiotropy", n_traits = nt,
                       cor = r, seed = s, ...))
}
.av_comp_cor <- function(sim, i, t1 = 1L, t2 = 2L) {
  ly <- sim$layers[[i]]
  stats::cor(.component_raw(ly, sim, t1, 1L), .component_raw(ly, sim, t2, 1L))
}

# a minimal stand-in for a phenotype_sim, enough for the internal draw helpers
.av_fake_sim <- function(maf, chr = NULL, nt = 2L, arch = "independent",
                         all_het = NULL, arch_args = list()) {
  sim <- list(maf = maf, n_traits = nt, architecture = arch,
              arch_args = arch_args)
  if (!is.null(all_het)) sim$all_het <- all_het
  if (!is.null(chr)) sim$map <- data.frame(chr = chr)
  sim
}

# ---------------------------------------------------------------------------
# EFFECTS-F2: per-component realized correlation at a NEGATIVE cor
# ---------------------------------------------------------------------------
test_that("dominance and epistasis components each track a negative cor (EFFECTS-F2)", {
  skip_on_cran()
  res <- vapply(1:4, function(s) {
    suppressWarnings({
      sim <- .av_pleio(-0.5, s) |> additive(prop = 0.3, n_qtn = 120) |>
        dominance(prop = 0.2) |> epistasis(prop = 0.2, n_pairs = 120)
      c(D = .av_comp_cor(sim, 2L), E = .av_comp_cor(sim, 3L))
    })
  }, c(D = 0, E = 0))
  m <- rowMeans(res)
  expect_lt(abs(m[["D"]] - (-0.5)), 0.15)
  expect_lt(abs(m[["E"]] - (-0.5)), 0.15)
  expect_true(all(res[, 1] < 0))            # right sign, not merely close
})

# ---------------------------------------------------------------------------
# EFFECTS-F3: three traits, a signed cor matrix and per-trait pi on dominance
# ---------------------------------------------------------------------------
test_that("a 3-trait signed cor matrix with per-trait pi is realized by the dominance layer (EFFECTS-F3)", {
  skip_on_cran()
  R <- matrix(c(1, -0.4, 0.3, -0.4, 1, 0.2, 0.3, 0.2, 1), 3)
  res <- vapply(1:8, function(s) {
    suppressWarnings({
      sim <- .av_pleio(R, s, nt = 3L, pi = c(0.9, 0.7, 0.8)) |>
        additive(prop = 0.3, n_qtn = 120) |>
        dominance(prop = 0.2, same_as_add = FALSE, n_qtn = 120)
      comps <- vapply(1:3, function(t) .component_raw(sim$layers[[2]], sim, t, 1L),
                      numeric(ncol(OB_av) - 5L))
      stats::cor(comps)[upper.tri(R)]
    })
  }, numeric(3))
  m <- rowMeans(res)
  expect_lt(max(abs(m - R[upper.tri(R)])), 0.12)
  expect_identical(sign(m), sign(R[upper.tri(R)]))
})

# ---------------------------------------------------------------------------
# EFFECTS-F4: .pleio_draw reproduces the reference algorithm (hand re-implementation)
# ---------------------------------------------------------------------------
test_that(".pleio_draw equals a hand re-implementation of the PleioArch algorithm (EFFECTS-F4)", {
  sim0 <- .av_pleio(0.4, 5, pi = 0.6)
  ph <- additive(sim0, prop = 0.3, n_qtn = 40)
  ly <- ph$layers[[1]]
  # layer seeding rule: (seed, type, occurrence)
  set.seed(.layer_seed(5, "additive", 0L))
  n_qtn <- 40L; pi_t <- 0.6; V <- 0.3; rho <- 0.4
  pleio_n <- round(pi_t * n_qtn); spec_n <- n_qtn - pleio_n          # 24, 16
  cand <- .candidate_markers(sim0)
  drawn <- sample(cand, pleio_n + 2L * spec_n, replace = FALSE)
  # shared effects ~ MVN(0, Sigma / pleio_n): Sigma = [pi V, rho V; rho V, pi V],
  # via the symmetric square root of Sigma / n (eigen), not a Cholesky factor
  Sigma <- matrix(c(pi_t * V, rho * V, rho * V, pi_t * V), 2L)
  e <- eigen(Sigma / pleio_n, symmetric = TRUE)
  root <- e$vectors %*% (sqrt(e$values) * t(e$vectors))
  W <- matrix(stats::rnorm(pleio_n * 2L), ncol = 2L) %*% root
  S <- lapply(1:2, function(t) stats::rnorm(spec_n, 0, sqrt((1 - pi_t) * V / spec_n)))
  sc <- function(i) 1 / sqrt(2 * sim0$maf[i] * (1 - sim0$maf[i]))   # 1/sqrt(2pq)
  shared <- drawn[seq_len(pleio_n)]
  for (t in 1:2) {
    spec_idx <- drawn[pleio_n + (t - 1L) * spec_n + seq_len(spec_n)]
    expect_equal(ly$qtn[[t]], c(shared, spec_idx))
    expect_equal(ly$effect[[t]], c(W[, t] * sc(shared), S[[t]] * sc(spec_idx)))
  }
  expect_identical(attr(ly$qtn, "pleio_shared"), as.integer(shared))
})

# ---------------------------------------------------------------------------
# EFFECTS-C3 (+F5c): .pleio_draw moments by class: shared-major, shared-minor, specific
# ---------------------------------------------------------------------------
test_that(".pleio_draw gives each class its analytic moments on the genotype-standardized scale (EFFECTS-C3)", {
  skip_on_cran()
  vg <- c(0.3, 0.5); pi_v <- c(0.6, 0.8); rho <- -0.4
  n_qtn <- 30L; n_major <- 2L; p_major <- 0.5
  sim0 <- .av_pleio(rho, 1, pi = pi_v, n_pleio_major = n_major,
                    prop_var_major = p_major)
  pleio_n <- round(mean(pi_v) * n_qtn)                                # 21
  spec_n <- n_qtn - pleio_n                                           # 9
  stat <- vapply(1:3000, function(s) {
    d <- suppressMessages(.pleio_draw(sim0, n_qtn, vg, s))
    w <- vapply(1:2, function(t) {
      idx <- d$qtn[[t]]; d$effect[[t]] * sqrt(2 * sim0$maf[idx] * (1 - sim0$maf[idx]))
    }, numeric(n_qtn))
    maj <- seq_len(n_major); mnr <- n_major + seq_len(pleio_n - n_major)
    spc <- pleio_n + seq_len(spec_n)
    c(maj11 = sum(w[maj, 1]^2), maj12 = sum(w[maj, 1] * w[maj, 2]),
      mnr11 = sum(w[mnr, 1]^2), mnr12 = sum(w[mnr, 1] * w[mnr, 2]),
      mnr22 = sum(w[mnr, 2]^2),
      spc1 = sum(w[spc, 1]^2), spc2 = sum(w[spc, 2]^2),
      spc12 = sum(w[spc, 1] * w[spc, 2]))
  }, numeric(8))
  m <- rowMeans(stat)
  s12 <- rho * sqrt(vg[1] * vg[2])
  expect_equal(m[["maj11"]], p_major * pi_v[1] * vg[1], tolerance = 0.08)
  expect_equal(m[["maj12"]], p_major * s12, tolerance = 0.08)
  expect_equal(m[["mnr11"]], (1 - p_major) * pi_v[1] * vg[1], tolerance = 0.06)
  expect_equal(m[["mnr22"]], (1 - p_major) * pi_v[2] * vg[2], tolerance = 0.06)
  expect_equal(m[["mnr12"]], (1 - p_major) * s12, tolerance = 0.06)
  # trait-specific loci: (1 - pi_t) V_t each, independent across traits
  expect_equal(m[["spc1"]], (1 - pi_v[1]) * vg[1], tolerance = 0.06)
  expect_equal(m[["spc2"]], (1 - pi_v[2]) * vg[2], tolerance = 0.06)
  expect_lt(abs(m[["spc12"]]), 0.01)
  # total per trait is V_t and the total covariance sits in the shared loci
  expect_equal(m[["maj11"]] + m[["mnr11"]] + m[["spc1"]], vg[1], tolerance = 0.05)
})

# ---------------------------------------------------------------------------
# EFFECTS-F5: n_pleio_major guards and non-convergence
# ---------------------------------------------------------------------------
test_that("n_pleio_major / prop_var_major guards fire with their own messages (EFFECTS-F5)", {
  base <- function(...) .av_pleio(0, 1, pi = 0.1, ...)
  # (a) more major loci than shared loci (pi = 0.1, n_qtn = 20 -> 2 shared)
  expect_error(additive(base(n_pleio_major = 5, prop_var_major = 0.5),
                        prop = 0.3, n_qtn = 20), "cannot exceed")
  # (b) every shared locus is major but prop_var_major < 1 parks variance on none
  expect_error(additive(base(n_pleio_major = 2, prop_var_major = 0.5),
                        prop = 0.3, n_qtn = 20), "assigns .* to minor loci")
  # (c) the pair must be set together (rejected by the constructor)
  expect_error(base(prop_var_major = 0.5), "both be positive or both be zero")
  expect_error(base(n_pleio_major = 1, prop_var_major = 0),
               "both be positive or both be zero")
  # n_pleio_major == shared count with prop_var_major = 1 is allowed
  expect_no_error(suppressWarnings(
    additive(base(n_pleio_major = 2, prop_var_major = 1), prop = 0.3, n_qtn = 20)))
})

test_that("a major pleiotropic locus stops the realized correlation converging on cor (EFFECTS-F5)", {
  skip_on_cran()
  sd_r <- function(...) {
    stats::sd(vapply(1:20, function(s) {
      ph <- .av_pleio(0.4, s, pi = 0.6, ...) |> additive(prop = 0.5, n_qtn = 400)
      stats::cor(genetic_values(ph))[1, 2]
    }, 0))
  }
  without <- sd_r()
  with_major <- sd_r(n_pleio_major = 1, prop_var_major = 0.8)
  expect_lt(without, 0.08)         # many equal-weight loci: only sampling spread
  expect_gt(with_major, 0.2)       # 80% of the shared variance on one effect pair
  expect_gt(with_major, 3 * without)
})

# ---------------------------------------------------------------------------
# EFFECTS-F6 / C12: exactly singular feasible boundary; .draw_mvnorm
# ---------------------------------------------------------------------------
test_that("the cor^2 = pi1 pi2 boundary draws finite effects; a hair beyond it errors (EFFECTS-F6)", {
  sim <- function(pi_) .av_pleio(0.5, 1, pi = pi_)
  ph <- additive(sim(0.5), prop = 0.4, n_qtn = 40)
  ef <- unlist(ph$layers[[1]]$effect)
  expect_true(all(is.finite(ef)))
  expect_error(additive(sim(0.5 - 1e-6), prop = 0.4, n_qtn = 40),
               "Biological constraint")
  # three traits: the PSD eigenvalue test at the singular boundary (M = R with
  # diag pi is rank one here: eigenvalues 1.5, 0, 0) passes; one notch beyond fails
  R3 <- matrix(0.5, 3, 3); diag(R3) <- 1
  vg3 <- rep(0.3, 3)
  sig <- function(pi3) { s <- outer(sqrt(vg3), sqrt(vg3)) * R3; diag(s) <- pi3 * vg3; s }
  expect_no_error(.check_pleio_feasible(sig(rep(0.5, 3)), R3, rep(0.5, 3)))
  expect_error(.check_pleio_feasible(sig(rep(0.5 - 1e-6, 3)), R3, rep(0.5 - 1e-6, 3)),
               "not attainable")
})

test_that(".draw_mvnorm reproduces a singular negative covariance exactly and a regular one in moment (EFFECTS-C12)", {
  S <- matrix(c(0.5, -0.5, -0.5, 0.5), 2)                 # cor = -1 boundary
  set.seed(4)
  z <- .draw_mvnorm(25L, S)
  expect_equal(dim(z), c(25L, 2L))
  expect_lt(max(abs(rowSums(z))), 1e-12)                   # w2 = -w1 per row
  expect_equal(stats::cor(z)[1, 2], -1)
  # a nonsingular negative covariance: E[crossprod(draw)] = Sigma
  S2 <- matrix(c(0.4, -0.15, -0.15, 0.3), 2)
  set.seed(5)
  acc <- Reduce(`+`, lapply(1:3000, function(i) crossprod(.draw_mvnorm(10L, S2)))) / 3000
  expect_equal(acc, S2, tolerance = 0.05)
  # zero rows / zero diagonal: no RNG consumed, all-zero block
  set.seed(6); before <- .Random.seed
  expect_equal(.draw_mvnorm(0L, S2), matrix(0, 0, 2))
  expect_equal(.draw_mvnorm(5L, diag(0, 2)), matrix(0, 5, 2))
  expect_identical(before, .Random.seed)
})

# ---------------------------------------------------------------------------
# EFFECTS-F8: the ld window is inclusive at both ends (synthetic 4-marker panel)
# ---------------------------------------------------------------------------
.av_ld_panel <- function() {
  set.seed(11)
  n <- 60
  a <- sample(c(-1, 0, 1), n, TRUE, c(0.3, 0.4, 0.3))
  b <- a; k <- sample(n, 12); b[k] <- sample(c(-1, 0, 1), 12, TRUE)
  c3 <- sample(c(-1, 0, 1), n, TRUE)
  d4 <- c3; k <- sample(n, 25); d4[k] <- sample(c(-1, 0, 1), 25, TRUE)
  M <- rbind(a, b, c3, d4)                                  # markers x individuals
  df <- cbind(data.frame(snp = paste0("m", 1:4), allele = "A/G",
                         chr = c(1, 1, 2, 2), pos = c(1e4, 2e4, 1e4, 2e4),
                         cm = NA_real_, stringsAsFactors = FALSE),
              as.data.frame(M))
  names(df)[-(1:5)] <- paste0("i", seq_len(n))
  list(df = df, r2 = c(stats::cor(a, b)^2, stats::cor(c3, d4)^2))
}

test_that("ld r2 window is inclusive at both endpoints and exclusive just inside them (EFFECTS-F8)", {
  P <- .av_ld_panel()
  lo <- min(P$r2); hi <- max(P$r2)
  expect_gt(hi - lo, 0.05)
  mk <- function(r_lo, r_hi, nq) {
    simulate_phenotype(P$df, architecture = "ld", n_traits = 2,
                       ld_type = "direct", r2_min = r_lo, r2_max = r_hi,
                       n_qtn = nq, seed = 1) |> additive(prop = 0.4)
  }
  # both pairs sit exactly ON the endpoints: both are admissible, n_qtn = 2 works
  ph <- mk(lo, hi, 2L)
  ld <- ph$layers[[1]]$ld
  expect_setequal(c(ld$qtn_t1, ld$qtn_t2), 1:4)
  expect_equal(sort(ld$r2), sort(P$r2))
  # shrinking either end by 1e-9 (relative) removes that pair
  expect_error(mk(lo * (1 + 1e-9), hi * (1 - 1e-9), 1L), "could not find")
  expect_error(mk(lo * (1 + 1e-9), hi, 2L), "could not find")
  expect_error(mk(lo, hi * (1 - 1e-9), 2L), "could not find")
  # a window pinned to the lower pair returns only that pair
  one <- mk(lo * (1 - 1e-9), lo * (1 + 1e-9), 1L)$layers[[1]]$ld
  expect_setequal(c(one$qtn_t1, one$qtn_t2), if (P$r2[1] < P$r2[2]) 1:2 else 3:4)
  # nothing below the weakest pair
  expect_error(mk(0, lo * (1 - 1e-9), 1L), "could not find")
})

# ---------------------------------------------------------------------------
# EFFECTS-F9 / C16: recompute every ld geometry quantity from genotypes and map
# ---------------------------------------------------------------------------
test_that("ld draws recomputed from genotypes: windows, flanks, chromosomes, uniqueness (EFFECTS-F9)", {
  skip_on_cran()
  r2_of <- function(ph, i, j) {
    as.numeric(stats::cor(.geno_cols(ph, i), .geno_cols(ph, j)))^2
  }
  for (lt in c("direct", "indirect")) {
    for (s in c(2, 7)) {
      ph <- simulate_phenotype(G_av, architecture = "ld", n_traits = 2,
                               ld_type = lt, r2_min = 0.3, r2_max = 0.7,
                               n_qtn = 4, seed = s) |> additive(prop = 0.4)
      ld <- ph$layers[[1]]$ld
      chr <- ph$map$chr; pos <- ph$map$pos
      expect_identical(ph$layers[[1]]$qtn[[1]], ld$qtn_t1)
      expect_identical(ph$layers[[1]]$qtn[[2]], ld$qtn_t2)
      used <- c(ld$qtn_t1, ld$qtn_t2, ld$cause[!is.na(ld$cause)])
      expect_identical(anyDuplicated(used), 0L)            # all loci distinct
      for (i in seq_len(nrow(ld))) {
        a <- ld$qtn_t1[i]; b <- ld$qtn_t2[i]
        expect_identical(chr[a], chr[b])
        rp <- r2_of(ph, a, b)
        expect_equal(ld$r2[i], rp)                         # reported = recomputed
        if (lt == "direct") {
          expect_true(is.na(ld$cause[i]))
          expect_true(rp >= 0.3 && rp <= 0.7)
        } else {
          cz <- ld$cause[i]
          expect_identical(chr[cz], chr[a])
          expect_true(pos[a] < pos[cz] && pos[cz] < pos[b])       # t1 < cause < t2
          expect_true(all(c(r2_of(ph, cz, a), r2_of(ph, cz, b), rp) >= 0.3))
          expect_true(all(c(r2_of(ph, cz, a), r2_of(ph, cz, b), rp) <= 0.7))
        }
      }
    }
  }
})

# ---------------------------------------------------------------------------
# EFFECTS-F11 / C1: effect series boundaries and edge cases
# ---------------------------------------------------------------------------
test_that(".effect_series edge cases and exact overflow / underflow positions (EFFECTS-F11)", {
  expect_equal(.effect_series(0L), numeric(0))
  expect_equal(.effect_series(1L), 0.5)                               # default base
  expect_equal(.effect_series(1L, effect = 2), 2)                     # base 2, n = 1
  expect_equal(.effect_series(4L, effect = -0.5),
               c(-0.5, 0.25, -0.125, 0.0625))                         # alternating sign
  expect_equal(.effect_series(3L), c(0.5, 0.25, 0.125))
  expect_equal(.effect_series(3L, effect = c(2, -1, 5)), c(2, -1, 5)) # verbatim
  # wrong length / non-finite / non-numeric explicit series
  expect_error(.effect_series(3L, effect = c(1, 2)), "length n_qtn")
  expect_error(.effect_series(3L, effect = c(1, 2), count = "n_pairs",
                              arg = "a"), "`a` must be .* length n_pairs")
  expect_error(.effect_series(3L, effect = c(1, NA, 2)), "finite")
  expect_error(.effect_series(3L, effect = Inf), "finite")
  expect_error(.effect_series(3L, effect = "x"), "finite")
  expect_error(.effect_series(3L, dist = c("geometric", "x")), "one non-missing")
  # 2^1023 is the largest finite power of 2; 2^1024 overflows
  expect_true(all(is.finite(.effect_series(1023L, effect = 2))))
  expect_error(.effect_series(1024L, effect = 2), "from position 1024")
  expect_error(.effect_series(1100L, effect = -2), "from position 1024")
  # 0.5^1074 is the smallest subnormal; 0.5^1075 underflows to exactly 0
  expect_true(all(.effect_series(1074L, effect = 0.5) != 0))
  expect_error(.effect_series(1075L, effect = 0.5), "from position 1075")
  # a zero base is the documented all-zero series, never an overflow/underflow error
  expect_equal(.effect_series(1100L, effect = 0), rep(0, 1100L))
  # the failure names the series and the count argument, not "polymorphic loci"
  err <- tryCatch(.effect_series(1100L, effect = 2, count = "n_pairs"),
                  error = function(e) conditionMessage(e))
  expect_match(err, "geometric effect series", fixed = TRUE)
  expect_match(err, "n_pairs = 1100", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# EFFECTS-F13: independently drawn components add no cross-component covariance
# ---------------------------------------------------------------------------
test_that("additive trait 1 and dominance trait 2 are uncorrelated on average (EFFECTS-F13)", {
  skip_on_cran()
  r <- vapply(1:24, function(s) {
    suppressWarnings({
      sim <- .av_pleio(0.3, s) |> additive(prop = 0.3, n_qtn = 100) |>
        dominance(prop = 0.2, same_as_add = FALSE, n_qtn = 100)
      stats::cor(.component_raw(sim$layers[[1]], sim, 1L, 1L),
                 .component_raw(sim$layers[[2]], sim, 2L, 1L))
    })
  }, 0)
  expect_lt(abs(mean(r)), 0.04)
})

# ---------------------------------------------------------------------------
# EFFECTS-F14: the fixed-n sampling floor of the realized correlation
# ---------------------------------------------------------------------------
test_that("with n = 20 the realized cor has the sampling spread of a correlation, however many loci (EFFECTS-F14)", {
  skip_on_cran()
  M20 <- .av_matrix(n = 20, m = 1500)
  r <- vapply(1:60, function(s) {
    ph <- suppressMessages(
      simulate_phenotype(M20, architecture = "pleiotropy", n_traits = 2,
                         cor = 0, seed = s)) |> additive(prop = 0.5, n_qtn = 1000)
    stats::cor(genetic_values(ph))[1, 2]
  }, 0)
  # an independent-normal sample correlation at n = 20 has sd ~ 1/sqrt(19) = 0.23
  expect_gt(stats::sd(r), 0.12)
  expect_lt(stats::sd(r), 0.28)
  expect_lt(abs(mean(r)), 0.12)
})

# ---------------------------------------------------------------------------
# EFFECTS-F15 / C10: .pleio_total_cor_check
# ---------------------------------------------------------------------------
.av_warnings <- function(expr) {
  w <- character(0)
  withCallingHandlers(suppressMessages(expr),
                      warning = function(c) {
                        w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
                      })
  w
}

test_that("without cor the total-correlation check is silent even for non-proportional props (EFFECTS-F15)", {
  w <- .av_warnings(
    suppressMessages(simulate_phenotype(OB_av, architecture = "pleiotropy",
                                        n_traits = 2, seed = 1)) |>
      additive(prop = c(0.49, 0.01), n_qtn = 60) |>
      dominance(prop = c(0.01, 0.49)))
  expect_false(any(grepl("TOTAL", w)))
  sim <- simulate_phenotype(OB_av, architecture = "pleiotropy", n_traits = 2,
                            seed = 1)
  expect_null(.pleio_total_cor_check(sim, c(0.1, 0.4)))
  # the same layers with cor requested do warn (the guard is the missing cor)
  w2 <- .av_warnings(.av_pleio(0.5, 1) |>
                       additive(prop = c(0.49, 0.01), n_qtn = 60) |>
                       dominance(prop = c(0.01, 0.49)))
  expect_true(any(grepl("TOTAL", w2)))
})

test_that("the TOTAL warning equals the variance-weighted formula, signs kept, over >2 traits (EFFECTS-C10)", {
  R <- matrix(c(1, 0.5, -0.4, 0.5, 1, -0.2, -0.4, -0.2, 1), 3)
  pa <- c(0.4, 0.1, 0.3); pd <- c(0.05, 0.4, 0.1)
  w <- .av_warnings(.av_pleio(R, 3, nt = 3L) |>
                      additive(prop = pa, n_qtn = 90) |> dominance(prop = pd))
  w <- w[grepl("TOTAL", w)]
  expect_length(w, 1L)
  expect_false(grepl("NaN|NA", w))
  for (pr in list(c(1, 2), c(1, 3), c(2, 3))) {
    i <- pr[1]; j <- pr[2]
    tot <- R[i, j] * (sqrt(pa[i] * pa[j]) + sqrt(pd[i] * pd[j])) /
      sqrt((pa[i] + pd[i]) * (pa[j] + pd[j]))
    pat <- sprintf("traits %d-%d: %s \\(cor %s\\)", i, j,
                   sprintf("%.3f", tot), sprintf("%.3f", R[i, j]))
    expect_match(w, pat)
    expect_identical(sign(tot), sign(R[i, j]))
    expect_lt(abs(tot), abs(R[i, j]))                       # attenuated toward 0
  }
  # a trait with zero total variance: its pairs are skipped, never NaN
  Rz <- matrix(c(1, 0, -0.4, 0, 1, 0, -0.4, 0, 1), 3)
  w0 <- .av_warnings(.av_pleio(Rz, 3, nt = 3L) |>
                       additive(prop = c(0.4, 0, 0.3), n_qtn = 90) |>
                       dominance(prop = c(0.05, 0, 0.4)))
  tot_w <- w0[grepl("TOTAL", w0)]
  expect_length(tot_w, 1L)
  expect_false(grepl("NaN|NA|traits 1-2|traits 2-3", tot_w))
  expect_match(tot_w, "traits 1-3:")
})

# ---------------------------------------------------------------------------
# EFFECTS-C4: .pleio_check_zero_var, three-trait table
# ---------------------------------------------------------------------------
test_that(".pleio_check_zero_var errors on nonzero cor, warns once per cor = 0 pair, never reports 0 (EFFECTS-C4)", {
  Rm <- function(a12 = 0, a13 = 0, a23 = 0) {
    R <- diag(3); R[1, 2] <- R[2, 1] <- a12; R[1, 3] <- R[3, 1] <- a13
    R[2, 3] <- R[3, 2] <- a23; R
  }
  cnt <- function(vg, R) {
    n <- 0L
    withCallingHandlers(.pleio_check_zero_var(R, vg),
                        warning = function(w) { n <<- n + 1L; invokeRestart("muffleWarning") })
    n
  }
  # all variances positive: silent
  expect_silent(.pleio_check_zero_var(Rm(0.3, -0.2, 0.1), c(0.1, 0.2, 0.3)))
  # one zero trait: a nonzero correlation involving it is an error
  expect_error(.pleio_check_zero_var(Rm(a12 = 0.3), c(0.1, 0, 0.3)),
               "trait 2 has zero additive variance")
  # (pairs are visited 1-2, 1-3, 2-3: the cor = 0 pair 1-2 warns before pair 2-3 errors)
  expect_error(suppressWarnings(.pleio_check_zero_var(Rm(a23 = -0.3), c(0.1, 0, 0.3))),
               "trait 2")
  expect_error(suppressWarnings(.pleio_check_zero_var(Rm(a13 = 0.5), c(0, 0.2, 0.3))),
               "trait 1")
  # ... but a nonzero cor between the two positive-variance traits is fine
  expect_equal(cnt(c(0, 0.2, 0.3), Rm(a23 = 0.5)), 2L)    # pairs 1-2 and 1-3 only
  # one zero trait, every cor = 0: one warning per affected pair (2), none for 2-3
  expect_equal(cnt(c(0, 0.2, 0.3), Rm()), 2L)
  expect_equal(cnt(c(0.1, 0, 0.3), Rm()), 2L)
  # two zero traits: all three pairs affected
  expect_equal(cnt(c(0, 0, 0.3), Rm()), 3L)
  expect_error(.pleio_check_zero_var(Rm(a12 = 0.2), c(0, 0, 0.3)), "zero")
  # the component name is used in the message
  expect_error(.pleio_check_zero_var(Rm(a12 = 0.3), c(0.1, 0, 0.3), "dominance"),
               "zero dominance variance")
  # end to end: the zero-variance trait's genetic value is constant and its
  # correlation with the others is undefined (NA), not 0
  sim <- .av_pleio(0, 1, pi = 1)
  ph <- suppressWarnings(additive(sim, prop = c(0.3, 0), n_qtn = 20))
  gv <- genetic_values(ph)
  expect_true(all(gv[, 2] == 0))
  r <- suppressWarnings(stats::cor(gv)[1, 2])
  expect_true(is.na(r))
})

# ---------------------------------------------------------------------------
# EFFECTS-C5: .pleio_partition table
# ---------------------------------------------------------------------------
test_that(".pleio_partition follows round(mean(pi) n) with its guards (EFFECTS-C5)", {
  R <- matrix(c(1, 0.3, 0.3, 1), 2)
  part <- function(n, pi_) suppressWarnings(.pleio_partition(n, pi_, R))
  # mean(pi) = 0.5: R's round() is half-to-even (0.5 -> 0, 1.5 -> 2, 2.5 -> 2)
  pi5 <- c(0.1, 0.9)
  expect_error(part(1L, pi5), "round the shared")
  expect_equal(unlist(part(2L, pi5)), c(pleio_n = 1, spec_n = 1))
  expect_equal(unlist(part(3L, pi5)), c(pleio_n = 2, spec_n = 1))
  expect_equal(unlist(part(4L, pi5)), c(pleio_n = 2, spec_n = 2))
  expect_equal(unlist(part(5L, pi5)), c(pleio_n = 2, spec_n = 3))
  expect_equal(unlist(part(6L, pi5)), c(pleio_n = 3, spec_n = 3))
  # pi = 0.25
  for (n in 1:2) expect_error(part(n, 0.25), "round the shared")
  expect_equal(unlist(part(3L, c(0.25, 0.25))), c(pleio_n = 1, spec_n = 2))
  expect_equal(unlist(part(6L, c(0.25, 0.25))), c(pleio_n = 2, spec_n = 4))
  # pi = 1: everything shared, no trait-specific units needed
  for (n in 1:6) expect_equal(unlist(part(n, c(1, 1))), c(pleio_n = n, spec_n = 0))
  # pi = 0: no shared units, all trait-specific, no error
  expect_equal(unlist(part(4L, c(0, 0))), c(pleio_n = 0, spec_n = 4))
  # pi < 1 but no unit left for the trait-specific class
  expect_error(part(3L, c(0.9, 0.9)), "leaves no trait-specific")
  # the invariant pleio_n + spec_n = n whenever it returns
  for (n in 3:12) for (p in c(0, 0.3, 0.5, 0.7, 1)) {
    out <- tryCatch(part(n, c(p, p)), error = function(e) NULL)
    if (!is.null(out)) expect_equal(out$pleio_n + out$spec_n, n)
  }
  # exactly one shared unit, |cor| < 1 (0 included): warns; |cor| = 1: silent
  expect_warning(.pleio_partition(2L, pi5, R), "only one shared")
  expect_warning(.pleio_partition(3L, c(0.25, 0.25), diag(2)), "only one shared")
  expect_no_warning(.pleio_partition(2L, pi5, matrix(c(1, 1, 1, 1), 2)))
  # the wording arguments only affect the message
  expect_error(.pleio_partition(1L, pi5, R, arg = "n_pairs", unit = "interacting set"),
               "n_pairs = 1")
})

# ---------------------------------------------------------------------------
# EFFECTS-C6: .pleio_single_unit_consequence, one four-trait fixture
# ---------------------------------------------------------------------------
test_that(".pleio_single_unit_consequence lists exact, noisy and nothing else in one 4-trait call (EFFECTS-C6)", {
  R <- diag(4)
  R[1, 2] <- R[2, 1] <- 0.3        # both traits without specifics -> exact +/-1
  R[1, 3] <- R[3, 1] <- 1          # |cor| = 1: realizable, skipped
  R[1, 4] <- R[4, 1] <- 0.2        # trait 4 has zero variance: skipped
  R[2, 3] <- R[3, 2] <- 0.2        # trait 3 has specifics -> one noisy draw
  txt <- .pleio_single_unit_consequence(R, no_spec = c(TRUE, TRUE, FALSE, FALSE),
                                        vg = c(0.3, 0.3, 0.3, 0))
  expect_match(txt, "traits 1-2: with no trait-specific variance")
  expect_match(txt, "exactly \\+/-1")
  expect_match(txt, "traits 2-3: the whole cross-trait covariance rests on that single locus")
  # nothing else is listed: no 1-3 (|cor| = 1), no pair with the zero-variance trait 4
  expect_false(grepl("1-3|1-4|2-4|3-4", txt))
  # the exact clause comes first, the noisy one second
  expect_lt(regexpr("traits 1-2", txt), regexpr("traits 2-3", txt))
  # nothing to report -> "" (no pair left)
  expect_identical(.pleio_single_unit_consequence(diag(3), c(TRUE, TRUE, TRUE),
                                                  vg = c(0, 0, 0)), "")
  expect_identical(.pleio_single_unit_consequence(matrix(1, 2, 2), c(TRUE, TRUE)), "")
  # `one` names the unit in the noisy clause
  expect_match(.pleio_single_unit_consequence(R, c(TRUE, TRUE, FALSE, FALSE),
                                              one = "set", vg = c(0.3, 0.3, 0.3, 0)),
               "single set")
})

# ---------------------------------------------------------------------------
# EFFECTS-C7 / C9: non-additive raw-component moments (V_t, pi V, Sigma_ij)
# ---------------------------------------------------------------------------
test_that(".pleio_unit_effects: total variance -> V, shared-only -> pi V, covariance -> Sigma_ij, effective target -> cor (EFFECTS-C7/C9)", {
  # exactly orthogonal, centered, unit-variance design columns: then
  # Var(c_t) = sum_u e_tu^2 and Cov(c_1, c_2) = sum over shared units e_1u e_2u
  n <- 200L; n_col <- 30L
  set.seed(12)
  Q <- qr.Q(qr(cbind(1, matrix(stats::rnorm(n * n_col), n))))[, -1L]
  Z <- Q * sqrt(n - 1)
  expect_equal(unname(apply(Z, 2, stats::var)), rep(1, n_col))
  expect_lt(max(abs(crossprod(Z)[upper.tri(diag(n_col))])), 1e-8)
  V <- 1; pi_v <- c(0.5, 0.5); rho <- 0.3
  sigma <- matrix(c(pi_v[1] * V, rho * V, rho * V, pi_v[2] * V), 2)
  q <- list(c(1:10, 11:20), c(1:10, 21:30))           # 10 shared + 10 specific each
  set.seed(13)
  draws <- lapply(1:400, function(i) {
    e <- .pleio_unit_effects(q, sigma, pi_v, c(V, V), function(j) Z[, j])
    list(e = e, tgt = attr(e, "target_cor"))
  })
  stat <- vapply(draws, function(d) {
    e <- d$e
    c(v1 = sum(e[[1]]^2), v2 = sum(e[[2]]^2),
      sh1 = sum(e[[1]][1:10]^2), sh2 = sum(e[[2]][1:10]^2),
      cv = sum(e[[1]][1:10] * e[[2]][1:10]))
  }, numeric(5))
  m <- rowMeans(stat)
  expect_equal(m[["v1"]], V, tolerance = 0.08)         # V_t, NOT Sigma_tt = 0.5
  expect_equal(m[["v2"]], V, tolerance = 0.08)
  expect_equal(m[["sh1"]], pi_v[1] * V, tolerance = 0.12)
  expect_equal(m[["sh2"]], pi_v[2] * V, tolerance = 0.12)
  expect_equal(m[["cv"]], rho * V, tolerance = 0.12)
  # each unit is divided by its design sd (here 1): the effects are the raw draws,
  # and the effective target the unit-effect routine reports is cor
  expect_equal(draws[[1]]$tgt, matrix(c(1, rho, rho, 1), 2))
  # specific units of different traits are independent: no shared covariance
  expect_lt(abs(mean(vapply(draws, function(d) sum(d$e[[1]][11:20] * d$e[[2]][11:20]), 0))), 0.05)
})

test_that("dominance under pleiotropy: raw component variance averages V_t, covariance Sigma_12 (EFFECTS-C7)", {
  skip_on_cran()
  V <- 0.2; pi_ <- 0.5; rho <- 0.3
  st <- vapply(1:20, function(s) {
    suppressWarnings({
      sim <- .av_pleio(rho, s, pi = pi_) |>
        additive(prop = 0.3, n_qtn = 60) |>
        dominance(prop = V, same_as_add = FALSE, n_qtn = 60)
      ly <- sim$layers[[2]]
      c1 <- .component_raw(ly, sim, 1L, 1L); c2 <- .component_raw(ly, sim, 2L, 1L)
      c(stats::var(c1), stats::var(c2), stats::cov(c1, c2))
    })
  }, numeric(3))
  m <- rowMeans(st)
  expect_equal(m[1], V, tolerance = 0.2)                # V_t (0.2), not pi V (0.1)
  expect_equal(m[2], V, tolerance = 0.2)
  expect_equal(m[3], rho * V, tolerance = 0.35)          # Sigma_12 = cor V = 0.06
  expect_gt(m[1], 1.5 * pi_ * V)
})

# ---------------------------------------------------------------------------
# EFFECTS-C8: .pleio_units disjointness / ordering
# ---------------------------------------------------------------------------
test_that(".pleio_units: shared units first and common, everything else globally distinct (EFFECTS-C8)", {
  pi_v <- c(0.8, 0.5, 0.3); R <- matrix(0.2, 3, 3); diag(R) <- 1
  sim <- .av_pleio(R, 1, nt = 3L, pi = pi_v)
  pleio_n <- round(mean(pi_v) * 10); spec_n <- 10 - pleio_n            # 5 + 5
  # epistatic sets, interaction = 3
  set.seed(1)
  q <- .pleio_units(sim, 10L, 3L, pi_v, R, "n_pairs")
  expect_length(q, 3L)
  for (t in 1:3) expect_equal(dim(q[[t]]), c(10L, 3L))
  for (t in 2:3) expect_identical(q[[t]][seq_len(pleio_n), ], q[[1]][seq_len(pleio_n), ])
  all_m <- c(q[[1]][seq_len(pleio_n), ],
             unlist(lapply(q, function(m) m[pleio_n + seq_len(spec_n), ])))
  expect_identical(anyDuplicated(all_m), 0L)
  expect_length(all_m, (pleio_n + 3 * spec_n) * 3L)                    # 60 markers
  expect_true(all(all_m %in% .candidate_markers(sim)))
  # dominance-style single locus per unit: vectors, shared first
  set.seed(2)
  q1 <- .pleio_units(sim, 10L, 1L, pi_v, R, "n_qtn")
  for (t in 1:3) expect_length(q1[[t]], 10L)
  for (t in 2:3) expect_identical(q1[[t]][seq_len(pleio_n)], q1[[1]][seq_len(pleio_n)])
  expect_identical(anyDuplicated(c(q1[[1]], q1[[2]][-seq_len(pleio_n)],
                                   q1[[3]][-seq_len(pleio_n)])), 0L)
  # marker shortage is an error that explains which way pi moves the need
  expect_error(.pleio_units(sim, 600L, 3L, pi_v, R, "n_pairs"), "distinct markers")
})

# ---------------------------------------------------------------------------
# EFFECTS-C11: .pleio_pi_vector / .pleio_cor_matrix
# ---------------------------------------------------------------------------
test_that(".pleio_pi_vector validates and expands every spelling of pi (EFFECTS-C11)", {
  s <- function(nt, ...) list(n_traits = nt, arch_args = list(...))
  expect_equal(.pleio_pi_vector(s(3L)), rep(1, 3))                      # default pi = 1
  expect_equal(.pleio_pi_vector(s(3L, pi = 0.4)), rep(0.4, 3))
  expect_equal(.pleio_pi_vector(s(3L, pi = c(0.2, 0.5, 1))), c(0.2, 0.5, 1))
  expect_equal(.pleio_pi_vector(s(2L, pi_target = 0.7, pi_secondary = 0.4)), c(0.7, 0.4))
  expect_error(.pleio_pi_vector(s(3L, pi = c(0.5, 0.5))), "length 1 or n_traits")
  expect_error(.pleio_pi_vector(s(2L, pi = 0.5, pi_target = 0.5)), "either `pi`")
  expect_error(.pleio_pi_vector(s(2L, pi = 0.5, pi_secondary = 0.5)), "either `pi`")
  expect_error(.pleio_pi_vector(s(3L, pi_target = 0.5, pi_secondary = 0.5)),
               "two-trait interface")
  expect_error(.pleio_pi_vector(s(2L, pi = 1.2)), "between 0 and 1")
  expect_error(.pleio_pi_vector(s(2L, pi = -0.1)), "between 0 and 1")
  expect_error(.pleio_pi_vector(s(2L, pi = c(0.5, NA))), "between 0 and 1")
  expect_error(.pleio_pi_vector(s(2L, pi = Inf)), "between 0 and 1")
  expect_error(.pleio_pi_vector(s(2L, pi = "a")), "between 0 and 1")
  expect_equal(.pleio_pi_vector(s(2L, pi = c(0, 1))), c(0, 1))          # closed bounds
})

test_that("pi_target alone / pi_secondary alone are accepted (EFFECTS-C11, EFFECTS-D1)", {
  s <- function(nt, ...) list(n_traits = nt, arch_args = list(...))
  expect_equal(.pleio_pi_vector(s(2L, pi_target = 0.7)), c(0.7, 0.7))
  expect_equal(.pleio_pi_vector(s(2L, pi_secondary = 0.4)), c(1, 0.4))
  expect_error(.pleio_pi_vector(s(3L, pi_target = 0.5)), "two-trait interface")
  expect_error(.pleio_pi_vector(s(3L, pi_secondary = 0.5)), "two-trait interface")
  # through the public constructor
  ph <- .av_pleio(0.3, 1, pi_target = 0.7) |> additive(prop = 0.3, n_qtn = 30)
  expect_length(ph$layers[[1]]$effect[[1]], 30L)
})

test_that(".pleio_cor_matrix expands a scalar and rejects every malformed matrix (EFFECTS-C11)", {
  s <- function(nt, cor = NULL) list(n_traits = nt, arch_args = if (is.null(cor)) list() else list(cor = cor))
  expect_equal(.pleio_cor_matrix(s(3L)), diag(3))                       # absent -> 0
  expect_equal(.pleio_cor_matrix(s(3L, -0.25)),
               matrix(c(1, -0.25, -0.25, -0.25, 1, -0.25, -0.25, -0.25, 1), 3))
  expect_equal(.pleio_cor_matrix(s(2L, 1)), matrix(1, 2, 2))
  expect_equal(.pleio_cor_matrix(s(2L, -1)), matrix(c(1, -1, -1, 1), 2))
  Rok <- matrix(c(1, 0.3, -0.2, 0.3, 1, 0.1, -0.2, 0.1, 1), 3)
  expect_identical(.pleio_cor_matrix(s(3L, Rok)), Rok)
  expect_error(.pleio_cor_matrix(s(3L, diag(2))), "3 x 3")               # wrong size
  bad <- Rok; bad[1, 2] <- NA
  expect_error(.pleio_cor_matrix(s(3L, bad)), "finite and between -1 and 1")
  bad <- Rok; bad[1, 2] <- bad[2, 1] <- Inf
  expect_error(.pleio_cor_matrix(s(3L, bad)), "finite and between -1 and 1")
  bad <- Rok; bad[1, 2] <- bad[2, 1] <- 1.2
  expect_error(.pleio_cor_matrix(s(3L, bad)), "between -1 and 1")
  bad <- Rok; bad[1, 2] <- 0.5                                            # asymmetric
  expect_error(.pleio_cor_matrix(s(3L, bad)), "symmetric")
  bad <- Rok; diag(bad) <- c(1, 0.9, 1)
  expect_error(.pleio_cor_matrix(s(3L, bad)), "diagonal")
  expect_error(.pleio_cor_matrix(s(2L, 1.5)), "between -1 and 1")
  expect_error(.pleio_cor_matrix(s(2L, NA_real_)), "between -1 and 1")
  expect_error(.pleio_cor_matrix(s(3L, c(0.2, 0.3))), "one finite value")
  expect_error(.pleio_cor_matrix(s(2L, "a")), "one finite value")
})

# ---------------------------------------------------------------------------
# EFFECTS-C13: .draw_qtn / .draw_qtn_pairs / .pleio_draw restore the caller's RNG
# ---------------------------------------------------------------------------
test_that("sub-seeded QTN draws are reproducible and leave the caller's RNG stream untouched (EFFECTS-C13)", {
  sim <- .av_pleio(0.3, 1)
  simi <- simulate_phenotype(OB_av, n_traits = 2, seed = 1)
  for (s in list(sim, simi)) {
    set.seed(10); stream <- stats::runif(3)
    set.seed(10); before <- .Random.seed
    a <- .draw_qtn(s, 7L, 4242L)
    expect_identical(.Random.seed, before)               # restored
    expect_equal(stats::runif(3), stream)                # stream continues unperturbed
    set.seed(999)                                        # caller state is irrelevant
    b <- .draw_qtn(s, 7L, 4242L)
    expect_identical(a, b)
    expect_false(identical(a, .draw_qtn(s, 7L, 4243L)))
  }
  set.seed(10); before <- .Random.seed
  p1 <- .draw_qtn_pairs(simi, 4L, 2L, 77L)
  expect_identical(.Random.seed, before)
  set.seed(5); expect_identical(p1, .draw_qtn_pairs(simi, 4L, 2L, 77L))
  set.seed(10); before <- .Random.seed
  d1 <- suppressMessages(.pleio_draw(sim, 20L, c(0.3, 0.3), 31L))
  expect_identical(.Random.seed, before)
  set.seed(3); expect_identical(d1, suppressMessages(.pleio_draw(sim, 20L, c(0.3, 0.3), 31L)))
  # without a sub-seed the draw uses (and advances) the caller's stream
  set.seed(10); a <- .draw_qtn(simi, 7L, NULL); after <- .Random.seed
  set.seed(10); expect_identical(.draw_qtn(simi, 7L, NULL), a)
  expect_false(identical(after, before))
  # architecture: pleiotropy shares one set, independent draws a set per trait
  expect_identical(.draw_qtn(sim, 7L, 5L)[[1]], .draw_qtn(sim, 7L, 5L)[[2]])
  ind <- .draw_qtn(simi, 7L, 5L)
  expect_false(identical(ind[[1]], ind[[2]]))
})

# ---------------------------------------------------------------------------
# EFFECTS-C14: .draw_qtn_distinct_chr
# ---------------------------------------------------------------------------
test_that(".draw_qtn_distinct_chr partitions chromosomes round-robin and has both shortage errors (EFFECTS-C14)", {
  # chromosome order is first occurrence: 3, 1, 2, 5, 4 -> trait 1: {3, 2, 4}, trait 2: {1, 5}
  chr <- rep(c(3, 1, 2, 5, 4), each = 10L)
  sim <- .av_fake_sim(maf = rep(0.3, 50L), chr = chr, nt = 2L)
  cand <- seq_len(50L)
  set.seed(1)
  q <- .draw_qtn_distinct_chr(sim, 6L, cand)
  expect_length(q, 2L)
  expect_true(all(chr[q[[1]]] %in% c(3, 2, 4)))
  expect_true(all(chr[q[[2]]] %in% c(1, 5)))
  expect_identical(anyDuplicated(unlist(q)), 0L)
  expect_true(all(lengths(q) == 6L))
  # three traits: chromosomes {3,5}, {1,4}, {2}
  sim3 <- .av_fake_sim(maf = rep(0.3, 50L), chr = chr, nt = 3L)
  q3 <- .draw_qtn_distinct_chr(sim3, 4L, cand)
  expect_true(all(chr[q3[[3]]] == 2))
  expect_true(all(chr[q3[[1]]] %in% c(3, 5)) && all(chr[q3[[2]]] %in% c(1, 4)))
  # fewer chromosomes than traits
  expect_error(.draw_qtn_distinct_chr(.av_fake_sim(rep(0.3, 20), rep(1:2, each = 10), 3L),
                                      2L, seq_len(20L)),
               "at least n_traits \\(3\\) chromosomes; the map has 2")
  # an assigned pool is too small: trait 1 owns only chromosome 1 (5 markers)
  chr_s <- c(rep(1, 5), rep(2, 30))
  sim_s <- .av_fake_sim(rep(0.3, 35L), chr_s, 2L)
  expect_error(.draw_qtn_distinct_chr(sim_s, 10L, seq_len(35L)),
               "trait 1 has only 5 markers")
  # the same shortage on trait 2's chromosome
  chr_t <- c(rep(1, 30), rep(2, 4))
  expect_error(.draw_qtn_distinct_chr(.av_fake_sim(rep(0.3, 34L), chr_t, 2L), 10L,
                                      seq_len(34L)), "trait 2 has only 4 markers")
  # candidates only: a monomorphic marker is not in any pool
  maf <- rep(0.3, 35L); maf[1:2] <- 0
  q_ok <- .draw_qtn_distinct_chr(.av_fake_sim(maf, chr_s, 2L), 3L,
                                 which(maf > 0))
  expect_false(any(q_ok[[1]] %in% 1:2))
  # through the public constructor (one chromosome, 2 traits)
  M1 <- .av_matrix(n = 40, m = 60)
  df1 <- cbind(data.frame(snp = colnames(M1), allele = "A/G", chr = 1,
                          pos = seq_len(60) * 1e4, cm = NA_real_,
                          stringsAsFactors = FALSE), as.data.frame(t(M1)))
  err <- tryCatch(
    simulate_phenotype(df1, n_traits = 2, distinct_chr = TRUE, seed = 1) |>
      additive(prop = 0.3, n_qtn = 3), error = function(e) conditionMessage(e))
  expect_match(err, "at least n_traits")
})

# ---------------------------------------------------------------------------
# EFFECTS-C15: .draw_qtn_pairs / .draw_qtn capacity and candidate pool
# ---------------------------------------------------------------------------
test_that(".draw_qtn_pairs samples only candidates, exactly at capacity, errors one over (EFFECTS-C15)", {
  maf <- c(0.1, NA, 0.2, 0, 0.3, 0.4, 0.25, 0.15, Inf, NaN)
  cand <- c(1L, 3L, 5L, 6L, 7L, 8L)                      # finite and > 0
  for (arch in c("independent", "pleiotropy")) {
    sim <- .av_fake_sim(maf, nt = 2L, arch = arch)
    expect_identical(.candidate_markers(sim), cand)
    set.seed(3)
    p <- .draw_qtn_pairs(sim, 3L, 2L, 11L)               # 6 = capacity
    expect_length(p, 2L)
    for (t in 1:2) {
      expect_equal(dim(p[[t]]), c(3L, 2L))
      expect_setequal(p[[t]], cand)                      # every candidate exactly once
      expect_identical(anyDuplicated(p[[t]]), 0L)
    }
    if (arch == "pleiotropy") expect_identical(p[[1]], p[[2]])
    else expect_false(identical(p[[1]], p[[2]]))
    expect_error(.draw_qtn_pairs(sim, 7L, 1L, 11L), "exceed candidate markers \\(6\\)")
    expect_error(.draw_qtn_pairs(sim, 2L, 4L, 11L), "exceed candidate markers \\(6\\)")
  }
  # .draw_qtn: n_qtn = capacity works, capacity + 1 errors before sampling
  sim <- .av_fake_sim(maf, nt = 2L, arch = "independent")
  d <- .draw_qtn(sim, 6L, 3L)
  expect_setequal(d[[1]], cand)
  set.seed(8); before <- .Random.seed
  expect_error(.draw_qtn(sim, 7L, 3L), "exceeds the number of polymorphic markers \\(6\\)")
  expect_identical(.Random.seed, before)                 # the check precedes any RNG use
  # all-heterozygous markers are excluded and named in the hint
  sim_h <- .av_fake_sim(c(0.2, 0.5, 0.3), nt = 2L, all_het = c(FALSE, TRUE, FALSE))
  expect_error(.draw_qtn(sim_h, 3L, 1L), "heterozygous in every individual")
})

# ---------------------------------------------------------------------------
# EFFECTS-C17: .gabriel_blocks tie-break, MAF floor, span, untouched markers
# ---------------------------------------------------------------------------
test_that(".gabriel_blocks keeps the first highest-MAF member, honours the MAF floor and the span (EFFECTS-C17)", {
  set.seed(3)
  n <- 120
  x <- stats::rbinom(n, 2, 0.4); y <- stats::rbinom(n, 2, 0.3)
  Dm <- rbind(x, x, x, x, y)                  # markers 1-4 identical: one perfect block
  rownames(Dm) <- NULL
  pos <- c(1000, 2000, 3000, 4000, 1e6); chr <- rep(1, 5); keep <- rep(TRUE, 5)
  gb <- function(maf, keep_ = keep, ...) .gabriel_blocks(Dm, chr, pos, keep_, maf, ...)
  # ties: equal MAF -> the FIRST member is the tag; the lone marker 5 is untouched
  expect_identical(gb(rep(0.3, 5)), c(TRUE, FALSE, FALSE, FALSE, TRUE))
  # a strictly higher MAF wins wherever it sits, and a tie among the top picks the first
  expect_identical(gb(c(0.2, 0.3, 0.35, 0.35, 0.3)), c(FALSE, FALSE, TRUE, FALSE, TRUE))
  expect_identical(gb(c(0.25, 0.3, 0.3, 0.3, 0.3)), c(FALSE, TRUE, FALSE, FALSE, TRUE))
  # MAF floor 0.05 (Haploview): at 0.05 (and within SMALL_EPSILON below it) the marker
  # takes part in the block; clearly below 0.05 it is ignored and left untouched
  expect_identical(gb(c(0.3, 0.3, 0.3, 0.05, 0.3)), c(TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(gb(c(0.3, 0.3, 0.3, 0.05 * (1 - 1e-14), 0.3)),
                   c(TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(gb(c(0.3, 0.3, 0.3, 0.05 * (1 - 1e-12), 0.3)),
                   c(TRUE, FALSE, FALSE, TRUE, TRUE))
  # markers already dropped (keep = FALSE) take no part and are not revived
  expect_identical(gb(rep(0.3, 5), c(TRUE, TRUE, TRUE, FALSE, TRUE)),
                   c(TRUE, FALSE, FALSE, FALSE, TRUE))
  # a span limit shorter than the block splits it (1.5 kb: {1,2} and {3,4})
  expect_identical(gb(rep(0.3, 5), max_kb = 1.5), c(TRUE, FALSE, TRUE, FALSE, TRUE))
  # a huge span is capped (no integer NA) and equals the unconstrained result
  expect_identical(gb(rep(0.3, 5), max_kb = 1e12), gb(rep(0.3, 5), max_kb = 500))
  # the result is a plain logical of the input length
  out <- gb(rep(0.3, 5))
  expect_type(out, "logical"); expect_length(out, 5L)
})
