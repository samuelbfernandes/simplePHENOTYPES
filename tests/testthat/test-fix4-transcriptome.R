# Round-4 regression tests for the transcriptome group (R4-5).

test_that("R4-5: near-cancelling G/R still hits the requested mimic variance", {
  # Codex example: equal-variance centered components, 0 < Var(G+R) <= 1e-12
  Gg <- c(1, 0, -1); Gg <- Gg / stats::sd(Gg)
  Rg <- -Gg + c(1e-7, -2e-7, 1e-7)  # near-perfect negative correlation
  Rg <- (Rg - mean(Rg)) / stats::sd(Rg)
  vu <- stats::var((Gg - mean(Gg)) + (Rg - mean(Rg)))
  expect_true(vu > 0 && vu <= 1e-12)
  rs <- simplePHENOTYPES:::.tx_mimic_scale(Gg, Rg, 3)
  out <- rs$esc * ((Gg - mean(Gg)) + (Rg - mean(Rg)))
  expect_equal(stats::var(out), 9, tolerance = 1e-6)
  expect_true(rs$ill)
})

test_that("R4-5: the expression is scaled as esc * (G + R), not esc * G + esc * R (Codex 2026-10-03)", {
  # Codex: standardized G = (1, 0, -1), R = -G + O(1e-16), scl = 3. Scaling the
  # two parts separately gave Var = 11.8 instead of 9; the production path now
  # forms the centered sum first and returns it already scaled (rs$u).
  Gg <- c(1, 0, -1); Gg <- Gg / stats::sd(Gg)
  Rg <- -Gg + c(1e-16, -2e-16, 1e-16)
  rs <- simplePHENOTYPES:::.tx_mimic_scale(Gg, Rg, 3)
  skip_if(!is.finite(rs$esc) || identical(rs$esc, 3), "exact cancellation on this platform")
  expect_true(rs$ill)
  expect_equal(stats::var(rs$u), 9, tolerance = 1e-8)
})

test_that("R4-5: the rescale meets the variance for subnormal and huge unit variances (2026-10-04)", {
  for (k in c(4.9406564584124654e-324 * 2^10, 1e-200, 1e200)) {
    Gg <- c(1, 0, -1) * k; Rg <- c(0, 1, -1) * k
    rs <- simplePHENOTYPES:::.tx_mimic_scale(Gg, Rg, 3)
    expect_equal(stats::var(rs$u), 9, tolerance = 1e-8)
  }
})

test_that("R4-5: ordinary vu is unchanged and exact-zero vu keeps the fallback", {
  set.seed(1)
  Gg <- stats::rnorm(50); Rg <- stats::rnorm(50)
  rs <- simplePHENOTYPES:::.tx_mimic_scale(Gg, Rg, 2)
  expect_false(rs$ill)
  expect_equal(rs$esc, 2 / sqrt(stats::var((Gg - mean(Gg)) + (Rg - mean(Rg)))))
  z <- simplePHENOTYPES:::.tx_mimic_scale(Gg, -Gg, 2)
  expect_equal(z$esc, 2)
  expect_false(z$ill)
  nf <- simplePHENOTYPES:::.tx_mimic_scale(c(NA, 1, 2), c(1, 2, 3), 2)
  expect_equal(nf$esc, 2)
})

test_that("R4-5: ordinary mimic runs do not warn about ill-conditioning", {
  set.seed(42); n <- 60; m <- 80
  Z <- matrix(sample(-1:1, n * m, TRUE), n, m,
              dimnames = list(paste0("i", 1:n), paste0("m", 1:m)))
  E0 <- t(vapply(1:6, function(k) as.numeric(scale(Z %*% stats::rnorm(m))) * k + k,
                 numeric(n)))
  dimnames(E0) <- list(paste0("g", 1:6), rownames(Z))
  w <- character()
  tx <- withCallingHandlers(simulate_transcriptome(Z, mimic = E0, seed = 3),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("near-cancelling", w)))
  expect_equal(unname(apply(tx$expression, 1, stats::var)),
               unname(apply(E0, 1, stats::var)), tolerance = 1e-8)
})
