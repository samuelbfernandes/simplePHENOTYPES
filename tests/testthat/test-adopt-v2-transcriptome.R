# Remaining test-gap proposals of the 2026-09 audit (group v2-transcriptome,
# reconciliation section 6) that were not adopted by test-audit-transcriptome.R,
# test-fix2/3/4-transcriptome.R or the older transcriptome tests. Small panels
# only; deterministic seeds; every tolerance is stated.

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04
.adopt_tx <- suppressWarnings(simulate_transcriptome(G, n_genes = 120, seed = 1))

# ---- TX-F5 / TX-F6: GREML calibration -------------------------------------------

test_that("TX-F5: mimic GREML targets are calibrated at the distribution level", {
  skip_on_cran()
  base <- simulate_transcriptome(G, n_genes = 150, h2 = 0.5, seed = 1)
  mm <- simulate_transcriptome(G, mimic = base$expression, seed = 2)
  gr <- mm$calibration$h2$h2_greml
  expect_equal(length(gr), 150L)
  expect_true(all(gr >= 0 & gr <= 1))
  # the MEAN of the per-gene estimates tracks the mean realized heritability of
  # the source (documented: per-gene estimates are noisy, sd ~ 0.1 at n = 280,
  # so only the distribution is calibrated); empirical gap 0.005, tolerance 0.05
  expect_lt(abs(mean(gr) - mean(base$genes$h2_realized)), 0.05)
  # and those estimates are the generator's h2 targets
  expect_equal(mm$genes$h2_target, gr)
})

test_that("TX-F6: .greml_h2 agrees with an external REML (rrBLUP::mixed.solve)", {
  skip_on_cran()
  mixed_solve <- optional_fun("rrBLUP", "mixed.solve")
  sim <- simplePHENOTYPES:::.normalize_geno(G, "geno")
  dose <- simplePHENOTYPES:::.geno_cols(sim, seq_len(sim$n_markers))
  K <- simplePHENOTYPES:::.tx_grm(sweep(dose, 2, colMeans(dose), "-"))
  n <- nrow(K)
  set.seed(3)
  eig <- eigen(K, symmetric = TRUE)
  Uh <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)))
  hs <- c(0.1, 0.3, 0.5, 0.7, 0.9, 0.4)
  Y <- t(vapply(hs, function(h) {
    as.numeric(Uh %*% rnorm(n)) * sqrt(h) + rnorm(n) * sqrt(1 - h)
  }, numeric(n)))
  Y <- t(scale(t(Y)))                       # the caller standardizes rows
  est <- simplePHENOTYPES:::.greml_h2(Y, K)
  ref <- vapply(seq_along(hs), function(i) {
    r <- mixed_solve(Y[i, ], K = K)
    r$Vu / (r$Vu + r$Ve)
  }, numeric(1))
  # observed max |difference| 3.7e-6 (optimizer tolerance); tolerance 1e-4
  expect_equal(est, ref, tolerance = 1e-4)
})

# ---- TX-F7 / TX-F8 / TX-F14: counts ----------------------------------------------

test_that("TX-F7: the NB mean/variance law holds across dispersions", {
  skip_on_cran()
  tx <- .adopt_tx
  mu <- exp(3)
  for (phi in c(0, 0.05, 0.5, 2)) {
    tc <- observe_counts(tx, baseline = 3, coupling = 0, dispersion = phi, seed = 11)
    v <- mean(apply(tc$counts, 1L, stats::var))
    expect_lt(abs(mean(tc$counts) / mu - 1), 0.03,
              label = paste("|mean ratio - 1|, phi =", phi))
    # empirical ratios 0.98-1.01 over 120 genes x 280 individuals; tolerance 5%
    expect_lt(abs(v / (mu + phi * mu^2) - 1), 0.05,
              label = paste("|variance ratio - 1|, phi =", phi))
  }
})

test_that("TX-F8: a per-individual library size scales the realized mean exactly", {
  tx <- .adopt_tx
  set.seed(2)
  L <- stats::runif(tx$n_ind, 0.2, 5)
  tc <- observe_counts(tx, baseline = 4, coupling = 0, library_size = L, seed = 3)
  # mu_gi = L_i exp(alpha) when coupling = 0 (exact)
  expect_equal(unname(tc$count_model$mu),
               unname(exp(4) * matrix(L, tx$n_genes, tx$n_ind, byrow = TRUE)),
               tolerance = 1e-10)
  # the observed per-individual mean count follows L (empirical r = 0.998)
  expect_gt(stats::cor(colMeans(tc$counts), L), 0.99)
})

test_that("TX-F14: count-scale R2 on the genetic expression is below the latent h2", {
  tx <- .adopt_tx
  tc <- observe_counts(tx, baseline = 3, coupling = 0.5, dispersion = 0.1, seed = 5)
  j <- which.max(tx$genes$h2_realized)
  r2 <- stats::cor(tc$counts[j, ], tx$genetic_expression[j, ])^2
  # empirical: latent 0.67 vs count-scale 0.29; require a clear gap (0.1)
  expect_lt(r2, tx$genes$h2_realized[j] - 0.1)
  # on average over genes the count-scale R2 does not exceed the latent h2
  r2all <- vapply(seq_len(tx$n_genes), function(g) {
    stats::cor(tc$counts[g, ], tx$genetic_expression[g, ])^2
  }, numeric(1))
  expect_lt(mean(r2all), mean(tx$genes$h2_realized))
})

test_that("C11: observe_counts stores mu = L exp(alpha + sigma z) exactly and isolates the RNG", {
  tx <- .adopt_tx
  tx$expression[1, ] <- 5                         # a constant gene: z = 0
  set.seed(1)
  L <- stats::runif(tx$n_ind, 0.5, 2)
  alpha <- seq(1, 3, length.out = tx$n_genes)
  sigma <- seq(0.2, 0.8, length.out = tx$n_genes)
  oc <- observe_counts(tx, library_size = L, baseline = alpha, coupling = sigma,
                       seed = 9)
  Z <- t(apply(tx$expression, 1L, function(r) {
    s <- stats::sd(r); if (s > 0) (r - mean(r)) / s else rep(0, length(r))
  }))
  mu <- exp(outer(alpha, rep(1, tx$n_ind)) + Z * sigma +
              matrix(log(L), tx$n_genes, tx$n_ind, byrow = TRUE))
  expect_equal(unname(oc$count_model$mu), unname(mu), tolerance = 1e-12)
  expect_true(all(Z[1, ] == 0))
  # the constant gene's mean is exp(alpha_1) L_i (sigma has no effect)
  expect_equal(unname(oc$count_model$mu[1, ]), exp(alpha[1]) * L, tolerance = 1e-12)
  # a seeded call leaves the global RNG state untouched
  set.seed(77)
  s0 <- .Random.seed
  invisible(observe_counts(tx, seed = 3))
  expect_identical(.Random.seed, s0)
})

# ---- C5: Marchenko-Pastur factor count ------------------------------------------

test_that("C5: the MP factor count recovers planted rank 0/1/3 co-expression", {
  skip_on_cran()
  mk <- function(rk, T = 300, n = 100, s = 1.5) {
    Es <- matrix(stats::rnorm(T * n), T, n)
    if (rk > 0) {
      Es <- Es + s * matrix(stats::rnorm(T * rk), T, rk) %*%
        matrix(stats::rnorm(rk * n), rk, n)
    }
    t(scale(t(Es)))
  }
  est <- function(rk, reps = 20L) {
    set.seed(1000 + rk)
    vapply(seq_len(reps), function(i) {
      simplePHENOTYPES:::.tx_estimate_factors(mk(rk))
    }, integer(1))
  }
  # the estimator floors at one factor, so rank 0 and rank 1 both give 1; the
  # false-factor rate (> 1) under rank 0 must stay small (empirical 0 / 30)
  expect_lte(mean(est(0L) > 1L), 0.1)
  expect_true(all(est(1L) == 1L))
  expect_true(all(est(3L) == 3L))
})

# ---- C6: predict() argument validation (beyond seed = 1.5 / residual = NA) ------

test_that("C6: predict() rejects non-finite, vector, string, logical and huge seeds and non-flag residuals", {
  tx <- suppressWarnings(simulate_transcriptome(G, n_genes = 20, seed = 1))
  Gn <- G[, c(1:5, 6:65)]
  for (s in list(Inf, NA, NA_integer_, "a", c(1, 2), 1e12, TRUE)) {
    expect_error(predict(tx, Gn, seed = s), "`seed` must be", info = format(s))
  }
  for (r in list(1, 0L, "TRUE", c(TRUE, FALSE), NULL)) {
    expect_error(predict(tx, Gn, residual = r), "`residual` must be TRUE or FALSE",
                 info = format(r))
  }
  # a valid call still works and honours residual = FALSE
  p <- predict(tx, Gn, residual = FALSE, seed = 3)
  expect_equal(p$expression, p$genetic_expression)
})

# ---- F10 / F11 / F12: layer seed threading and the mediation split ---------------

test_that("TX-F10: transcriptome/additive layer order does not change the phenotype", {
  tx <- .adopt_tx
  a <- simulate_phenotype(G, h2 = 0.5, seed = 3, transcriptome = tx) |>
    transcriptome(prop = 0.3, n_genes = 10) |>
    additive(prop = 0.2, n_qtn = 5)
  b <- simulate_phenotype(G, h2 = 0.5, seed = 3, transcriptome = tx) |>
    additive(prop = 0.2, n_qtn = 5) |>
    transcriptome(prop = 0.3, n_genes = 10)
  expect_identical(a$pheno$value, b$pheno$value)
  ga <- subset(qtn_table(a), layer == "transcriptome")
  gb <- subset(qtn_table(b), layer == "transcriptome")
  expect_identical(ga$snp, gb$snp)
  expect_identical(ga$effect, gb$effect)
})

test_that("TX-F11: multi-trait mediation_split reports one exact row per trait", {
  tx <- .adopt_tx
  ph <- simulate_phenotype(G, n_traits = 2, h2 = 0.5, seed = 3, transcriptome = tx) |>
    transcriptome(prop = c(0.4, 0.2), n_genes = 20)
  ms <- mediation_split(ph)
  expect_equal(nrow(ms), 2L)
  expect_identical(ms$trait, c("Trait_1", "Trait_2"))
  tot <- ms$genetic_mediated + ms$env_mediated + ms$covariance
  for (t in 1:2) {
    y <- ph$pheno$value[ph$pheno$trait == paste0("Trait_", t)]
    comp <- simplePHENOTYPES:::.transcriptome_matrix(ph, 1L)[, t]
    # the three shares sum to Var(component)/V_P, and the component is scaled to
    # exactly Var = prop_t (exact, tolerance 1e-8)
    expect_equal(tot[t], stats::var(comp) / stats::var(y), tolerance = 1e-8)
    expect_equal(stats::var(comp), c(0.4, 0.2)[t], tolerance = 1e-8)
  }
  expect_gt(tot[1], tot[2])
})

test_that("TX-F12: the mean genetic-mediated fraction tracks the causal genes' realized h2", {
  skip_on_cran()
  tx <- .adopt_tx
  res <- vapply(1:20, function(s) {
    ph <- simulate_phenotype(G, h2 = 0.5, seed = s, transcriptome = tx) |>
      transcriptome(prop = 0.5, n_genes = 10)
    m <- mediation_split(ph)
    ids <- subset(qtn_table(ph), layer == "transcriptome")$snp
    c(m$genetic_mediated / (m$genetic_mediated + m$env_mediated + m$covariance),
      mean(tx$genes$h2_realized[match(ids, tx$genes$gene_id)]))
  }, numeric(2))
  # empirical gap 0.009 over 20 seeds (per-seed sd 0.055); tolerance 0.04
  expect_lt(abs(mean(res[1, ]) - mean(res[2, ])), 0.04)
})

# ---- C10 / C15: covariance between direct marker and mediated genetic parts ------

test_that("C10: H2 = direct + mediated + 2Cov/V_P when an expression gene equals a causal marker", {
  tx <- .adopt_tx
  ph0 <- simulate_phenotype(G, h2 = 0.3, seed = 3, transcriptome = tx) |>
    additive(prop = 0.3, n_qtn = 1)
  snp <- qtn_table(ph0)$snp[1]
  ids <- colnames(tx$expression)
  dose <- unlist(G[G$snp == snp, -(1:5)])
  dose <- as.numeric(dose[match(ids, names(G)[-(1:5)])])
  dose <- dose - mean(dose)
  set.seed(5)
  tx2 <- tx
  tx2$genetic_expression[1, ] <- dose                   # gene 1's genetic part = the marker
  tx2$expression[1, ] <- dose + stats::rnorm(length(dose), sd = 0.5)
  gid <- rownames(tx2$expression)[1]
  covs <- numeric(2)
  for (i in 1:2) {
    slope <- c(1, -1)[i]
    ph <- simulate_phenotype(G, h2 = 0.3, seed = 3, transcriptome = tx2) |>
      additive(prop = 0.3, n_qtn = 1) |>
      transcriptome(prop = 0.3, genes = gid, slopes = slope)
    expect_identical(qtn_table(ph)$snp[1], snp)
    vp <- stats::var(ph$pheno$value)
    gm  <- simplePHENOTYPES:::.genetic_matrix(ph, 1L)[, 1]
    txg <- simplePHENOTYPES:::.transcriptome_matrix(ph, 1L, "genetic")[, 1]
    gv  <- genetic_values(ph)[, 1]
    expect_equal(gv, gm + txg, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(mediation_split(ph)$genetic_mediated, stats::var(txg) / vp,
                 tolerance = 1e-8)
    covs[i] <- 2 * stats::cov(gm, txg) / vp
    # exact variance identity: no test may assume H2 = direct + mediated
    expect_equal(stats::var(gv) / vp,
                 stats::var(gm) / vp + stats::var(txg) / vp + covs[i],
                 tolerance = 1e-8)
  }
  # the covariance is large and has the sign of the slope (empirical +0.36 / -1.14)
  expect_gt(covs[1], 0.1)
  expect_lt(covs[2], -0.1)
})

test_that("C15: the mediated genetic share is covariance-sensitive, not prop * mean(h2)", {
  tx <- .adopt_tx
  n <- tx$n_ind
  set.seed(8)
  g1 <- stats::rnorm(n); g1 <- g1 - mean(g1)
  Gm <- rbind(g1, -g1)                                   # two perfectly opposite genetic parts
  Em <- Gm + rbind(stats::rnorm(n), stats::rnorm(n))     # independent residuals
  dimnames(Gm) <- dimnames(Em) <- list(c("gA", "gB"), colnames(tx$expression))
  tx2 <- tx
  tx2$expression <- Em
  tx2$genetic_expression <- Gm
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 3, transcriptome = tx2) |>
    transcriptome(prop = 0.5, genes = c("gA", "gB"), slopes = c(1, 1))
  h2g <- apply(Gm, 1L, stats::var) / apply(Em, 1L, stats::var)
  expect_lt(max(abs(h2g - 0.5)), 0.05)
  # the naive prop * mean(h2) = 0.25, but the opposing genetic parts cancel in the
  # slope-weighted sum: the realized genetic-mediated share is ~0
  expect_lt(abs(0.5 * mean(h2g) - 0.25), 0.03)
  expect_lt(mediation_split(ph)$genetic_mediated, 1e-3)   # empirical 2.8e-6
})

# ---- O7: exact random-ranking baseline (benchmark 02) ----------------------------

test_that("O7: best-of-m random ranking baseline is 1 - choose(D,K)/choose(D+m,K)", {
  D <- 200L; K <- 5L
  exact <- function(m) 1 - choose(D, K) / choose(D + m, K)
  expect_equal(exact(1L), K / (D + 1L), tolerance = 1e-12)  # the old printed value, m = 1 only
  expect_gt(exact(2L), 1.9 * K / (D + 1L))                  # 4.90 % vs 2.49 %
  expect_equal(exact(3L), 0.072436, tolerance = 1e-4)
  skip_on_cran()
  set.seed(12)
  for (m in 2:3) {
    hit <- mean(vapply(seq_len(40000L), function(i) any(sample.int(D + m, K) <= m),
                       logical(1)))
    # binomial s.e. <= 0.0013 at 40000 draws; tolerance 0.008 absolute
    expect_lt(abs(hit - exact(m)), 0.008)
  }
})

test_that("benchmarks 02/04/05 carry the corrected estimands", {
  skip_if_no_source("benchmarks", "02_eqtl_recovery.R")
  skip_if_no_source("benchmarks", "04_mediation_recovery.R")
  skip_if_no_source("benchmarks", "05_twas_power.R")
  rd <- function(f) paste(readLines(testthat::test_path("..", "..", "benchmarks", f),
                                    warn = FALSE), collapse = "\n")
  b2 <- rd("02_eqtl_recovery.R"); b4 <- rd("04_mediation_recovery.R")
  b5 <- rd("05_twas_power.R")
  expect_true(grepl("choose(N_DECOY, TOPK)", b2, fixed = TRUE))
  expect_true(grepl("random_genes", b2, fixed = TRUE) && grepl("strongest", b2))
  expect_true(grepl("ONLY for independent", b4))
  expect_true(grepl("REPS", b5, fixed = TRUE) && grepl("structural", b5))
})
