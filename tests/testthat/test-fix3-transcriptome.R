# Round-3 regression tests for the transcriptome group (R3-4, R3-14, R3-15).

.fx3_panel <- function(seed = 42, n = 60L, m = 80L) {
  set.seed(seed)
  Z <- matrix(sample(-1:1, n * m, TRUE), n, m,
              dimnames = list(paste0("i", seq_len(n)), paste0("m", seq_len(m))))
  E0 <- t(vapply(seq_len(5), function(k) as.numeric(scale(Z %*% stats::rnorm(m))),
                 numeric(n)))
  dimnames(E0) <- list(paste0("g", 1:5), rownames(Z))
  list(Z = Z, E0 = E0)
}

.fx3_direct <- function(x) {
  vg <- apply(x$genetic_expression, 1L, stats::var)
  ve <- apply(x$expression, 1L, stats::var)
  vr <- apply(x$expression - x$genetic_expression, 1L, stats::var)
  list(real = unname(vg / ve), alloc = unname(vg / (vg + vr)), ve = ve)
}

# ---- R3-4: scale-free h2 fields ----------------------------------------------

test_that("R3-4: mimic at variance 1e-14 reports the direct ratios, not 0", {
  p <- .fx3_panel()
  tx <- suppressWarnings(simulate_transcriptome(p$Z, mimic = p$E0 * 1e-7, seed = 9))
  d <- .fx3_direct(tx)
  expect_true(all(d$ve < 1e-12 & d$ve > 0))
  expect_true(all(d$real > 0.9 & d$real <= 1.0001))
  expect_equal(unname(tx$genes$h2_realized), d$real, tolerance = 1e-8)
  expect_equal(unname(tx$genes$h2_var_ratio), d$real, tolerance = 1e-8)
  expect_equal(unname(tx$genes$h2_allocated), d$alloc, tolerance = 1e-8)
  expect_true(all(tx$genes$h2_allocated > 0.97))
  # scale-free: same seed at unit scale gives the same ratios
  tx1 <- suppressWarnings(simulate_transcriptome(p$Z, mimic = p$E0, seed = 9))
  expect_equal(tx$genes$h2_realized, tx1$genes$h2_realized, tolerance = 1e-6)
  expect_equal(tx$genes$h2_allocated, tx1$genes$h2_allocated, tolerance = 1e-6)
})

test_that("R3-4: predict(residual = FALSE) with positive variance reports 1", {
  p <- .fx3_panel()
  tx <- suppressWarnings(simulate_transcriptome(p$Z, mimic = p$E0 * 1e-7, seed = 9))
  p0 <- predict(tx, p$Z, residual = FALSE)
  d <- .fx3_direct(p0)
  expect_true(all(d$ve > 0))
  expect_equal(unname(p0$genes$h2_realized), rep(1, 5), tolerance = 1e-8)
  expect_equal(unname(p0$genes$h2_allocated), rep(1, 5), tolerance = 1e-8)
  # ordinary scale unchanged
  tx1 <- suppressWarnings(simulate_transcriptome(p$Z, mimic = p$E0, seed = 9))
  q0 <- predict(tx1, p$Z, residual = FALSE)
  expect_equal(unname(q0$genes$h2_realized), rep(1, 5), tolerance = 1e-8)
})

test_that("R3-4: zero-variance phenotype still reports 0", {
  tx <- suppressWarnings(simulate_transcriptome(NULL, n_ind = 20, n_genes = 3,
                                                seed = 1))
  expect_true(all(tx$genes$h2_realized == 0))
  expect_true(all(tx$genes$h2_allocated == 0))
  expect_false(anyNA(tx$genes$h2_realized))
})

# ---- R3-14: h2_realized identity ---------------------------------------------

test_that("R3-14: h2_realized = Var(G)/(Var(G)+Var(R)+gr_cov) at any scale", {
  set.seed(44)
  n <- 180L; m <- 300L
  Z <- matrix(sample(-1:1, n * m, TRUE), n, m,
              dimnames = list(paste0("i", seq_len(n)), paste0("m", seq_len(m))))
  sc <- Z[, 1:30] %*% matrix(stats::rnorm(30 * 8), 30, 8)
  E0 <- t(vapply(seq_len(8), function(k)
    as.numeric(scale(sqrt(0.4) * scale(sc[, k]) + sqrt(0.6) * stats::rnorm(n))),
    numeric(n)))
  dimnames(E0) <- list(paste0("g", 1:8), rownames(Z))
  tx <- suppressWarnings(simulate_transcriptome(Z, mimic = 3 * E0, seed = 9))
  vg <- apply(tx$genetic_expression, 1L, stats::var)
  vr <- apply(tx$expression - tx$genetic_expression, 1L, stats::var)
  gr <- tx$var_budget$gr_cov
  expect_equal(unname(tx$genes$h2_realized), unname(vg / (vg + vr + gr)),
               tolerance = 1e-8)
  # the unit-scale shortcut h2/(1 + gr_cov) is NOT the general identity here
  shortcut <- tx$genes$h2_target / (1 + gr)
  expect_gt(max(abs(tx$genes$h2_realized - shortcut)), 0.05)
  # ... but it does hold when Var(G)+Var(R) = 1 (no mimic)
  tu <- suppressWarnings(simulate_transcriptome(Z, n_genes = 10, h2 = 0.4, seed = 2))
  expect_equal(unname(tu$genes$h2_realized),
               unname(tu$genes$h2_target / (1 + tu$var_budget$gr_cov)),
               tolerance = 1e-8)
})

# ---- R3-15: marginal epistasis share, nondegenerate vs fallback --------------

.fx3_c11 <- function(seed) {
  set.seed(3)
  z <- sample(c(-1, 0, 1), 40, TRUE, prob = c(0.3, 0.4, 0.3))
  M <- cbind(m1 = z, m2 = z)
  rownames(M) <- paste0("i", seq_len(40))
  ann <- data.frame(gene_id = "g1", chr = 1, tss = 1)
  tx <- suppressWarnings(simulate_transcriptome(
    M, annotation = ann, cis_window = 0.5, h2 = 0.5, cis_fraction = 0.5,
    epistasis = 0.4, seed = seed))
  b <- tx$var_budget
  list(tx = tx,
       share = b$v_epi / (b$v_cis + b$v_trans + b$v_epi),
       s2 = 1 + b$cis_trans_cov / (b$v_cis + b$v_trans))
}

test_that("R3-15: nondegenerate blend obeys eps / (eps + (1 - eps) / s_ct^2)", {
  r <- .fx3_c11(4)                       # perfectly positively correlated scores
  expect_equal(r$s2, 2, tolerance = 1e-8)
  expect_equal(r$share, 0.4 / (0.4 + 0.6 / r$s2), tolerance = 1e-8)
  expect_equal(r$share, 0.4 / (0.4 + 0.3), tolerance = 1e-8)
})

test_that("R3-15: exact-cancellation fallback drops cis and the share is epsilon", {
  r <- .fx3_c11(1)                       # perfectly negatively correlated scores
  expect_equal(r$tx$genes$n_cis, 0L)     # cis part dropped by the fallback
  expect_equal(r$share, 0.4, tolerance = 1e-8)
  # the literal formula at the pre-fallback s_ct^2 = 0 would give 0, not 0.4
  expect_equal(0.4 / (0.4 + 0.6 / 0), 0)
})
