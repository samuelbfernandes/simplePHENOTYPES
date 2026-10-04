# test-fix4-ocs.R
#
# Round-4 R4-1: the above-unconstrained-target band in .tune_lambda() is purely
# relative (no absolute double.xmin floor).

.fx4_founders <- function(k = 4L) {
  g <- data.frame(snp = paste0("m", seq_len(k)), allele = "A/G",
                  chr = seq_len(k), pos = 1, cm = 0, P1 = 1L, P2 = -1L,
                  stringsAsFactors = FALSE)
  as_population(g)
}
.fx4_G <- function(s) (s * diag(2)) |>
  `dimnames<-`(list(c("a", "b"), c("a", "b")))

test_that("target = 2*double.xmin vs c0 = double.xmin warns (R4-1)", {
  fnd <- .fx4_founders()
  xm <- .Machine$double.xmin
  G <- .fx4_G(2 * xm)
  expect_warning(
    oc <- optimum_contribution(fnd, merit = c(a = 0, b = 1), G = G,
                               target_coancestry = 2 * xm),
    "above the coancestry of the unconstrained")
  expect_identical(oc$lambda, 0)
  expect_no_warning(
    oc0 <- optimum_contribution(fnd, merit = c(a = 0, b = 1), G = G,
                                target_coancestry = xm))
  expect_identical(oc0$lambda, 0)
})

test_that("subnormal G warns when target is 2x the optimum, silent when equal (R4-1)", {
  fnd <- .fx4_founders()
  G <- .fx4_G(2e-310)
  expect_warning(
    optimum_contribution(fnd, merit = c(a = 0, b = 1), G = G,
                         target_coancestry = 2e-310),
    "above the coancestry of the unconstrained")
  expect_no_warning(
    optimum_contribution(fnd, merit = c(a = 0, b = 1), G = G,
                         target_coancestry = 1e-310))
})

test_that("earlier R3-2 / R3-1 cases are unchanged (R4-1)", {
  fnd <- .fx4_founders()
  expect_warning(
    optimum_contribution(fnd, merit = c(a = 0, b = 1), G = .fx4_G(1),
                         target_coancestry = 0.5000005),
    "above the coancestry of the unconstrained")
  expect_no_warning(
    optimum_contribution(fnd, merit = c(a = 0, b = 1), G = .fx4_G(1),
                         target_coancestry = 0.5))
  Gs <- .fx4_G(2e-16)
  expect_warning(
    optimum_contribution(fnd, merit = c(a = 0, b = 1), G = Gs,
                         target_coancestry = 2e-16),
    "above the coancestry of the unconstrained")
  expect_no_warning(
    optimum_contribution(fnd, merit = c(a = 0, b = 1), G = Gs,
                         target_coancestry = 1e-16))
  oc <- optimum_contribution(fnd, merit = c(a = -1e308, b = 1e308),
                             G = .fx4_G(1), lambda = 0)
  expect_equal(unname(oc$contribution), c(0, 1))
  oc2 <- optimum_contribution(fnd, merit = c(a = 1e16, b = 1e16 + 2),
                              G = .fx4_G(1), target_coancestry = 0.26)
  expect_equal(unname(oc2$contribution), c(.4, .6), tolerance = 1e-4)
  expect_equal(oc2$lambda, 10, tolerance = 1e-3)
})

test_that("a feasible max_coancestry is met whatever the scale of G", {
  founders <- data.frame(snp = paste0("m", 1:2), allele = "A/G", chr = 1:2,
                         pos = 1, cm = 0, P1 = 1L, P2 = -1L,
                         stringsAsFactors = FALSE) |>
    as_population()
  for (s in c(1e-310, 1e-30, 1, 1e30)) {
    Gt <- diag(2) * s
    dimnames(Gt) <- list(c("a", "b"), c("a", "b"))
    oc <- expect_no_warning(optimum_contribution(founders, merit = c(a = 0, b = 1),
                                                 G = Gt, max_coancestry = 0.3 * s))
    expect_lte(oc$coancestry, 0.3 * s * (1 + 1e-8))
    expect_equal(unname(oc$contributions), c(0.5, 0.5) + c(-1, 1) / (2 * sqrt(5)),
                 tolerance = 1e-6)
    if (s >= 1e-300) expect_equal(oc$lambda * s, sqrt(5), tolerance = 1e-6)
  }
  # merit and relationships on the same extreme scale: the reported penalty stays finite
  Gt <- diag(2) * 1e308
  dimnames(Gt) <- list(c("a", "b"), c("a", "b"))
  oc <- optimum_contribution(founders, merit = c(a = 0, b = 1e308), G = Gt,
                             max_coancestry = 3e307)
  expect_equal(oc$lambda, sqrt(5), tolerance = 1e-6)
})
