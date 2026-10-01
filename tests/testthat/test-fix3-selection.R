# test-fix3-selection.R
#
# Round-3 fixes (Codex re-review of round 2): R3-1 OCS centring without
# overflow, R3-2 purely relative above-optimum band, R3-10 select_ind response
# identity needs a linear E[A | P], R3-11 re-standardised additive accuracy.

.fx3_founders <- function(k = 4L) {
  g <- data.frame(snp = paste0("m", seq_len(k)), allele = "A/G",
                  chr = seq_len(k), pos = 1, cm = 0, P1 = 1L, P2 = -1L,
                  stringsAsFactors = FALSE)
  as_population(g)
}
.fx3_G <- function(s = 1) (s * diag(2)) |>
  `dimnames<-`(list(c("a", "b"), c("a", "b")))

# ---- R3-1 -------------------------------------------------------------------

test_that("finite merits straddling +/- xmax centre without overflow (R3-1)", {
  fnd <- .fx3_founders()
  oc <- optimum_contribution(fnd, merit = c(a = -1e308, b = 1e308),
                             G = .fx3_G(), lambda = 0)
  expect_equal(unname(oc$contributions), c(0, 1))
  expect_equal(oc$merit, 1e308)            # original-scale merit reported
  expect_true(is.finite(oc$merit))
  # low direction: the minimum-merit individual is chosen
  lo <- optimum_contribution(fnd, merit = c(a = -1e308, b = 1e308),
                             G = .fx3_G(), lambda = 0, direction = "low")
  expect_equal(unname(lo$contributions), c(1, 0))
  expect_equal(lo$merit, -1e308)
  # tuned path through the same extreme range stays finite and valid
  tu <- suppressWarnings(
    optimum_contribution(fnd, merit = c(a = -1e308, b = 1e308),
                         G = .fx3_G(), target_coancestry = 0.3))
  expect_equal(sum(tu$contributions), 1)
  expect_true(all(is.finite(tu$contributions)))
})

test_that("round-2 merit-shift results are unchanged by the R3-1 centring", {
  fnd <- .fx3_founders()
  oc <- optimum_contribution(fnd, merit = c(a = 1e16, b = 1e16 + 2),
                             G = .fx3_G(), target_coancestry = 0.26)
  expect_equal(unname(oc$contributions), c(0.4, 0.6), tolerance = 1e-6)
  for (off in c(1e6, -1e16)) {
    a <- optimum_contribution(fnd, merit = c(a = off, b = off + 2),
                              G = .fx3_G(), lambda = 10)
    b <- optimum_contribution(fnd, merit = c(a = 0, b = 2),
                              G = .fx3_G(), lambda = 10)
    expect_equal(a$contributions, b$contributions, tolerance = 1e-8)
  }
})

# ---- R3-2 -------------------------------------------------------------------

test_that("the above-optimum band is purely relative (R3-2)", {
  fnd <- .fx3_founders()
  Gs <- .fx3_G(2e-16)
  # unconstrained coancestry 1e-16 < target 2e-16: must warn, lambda 0
  expect_warning(
    oc <- optimum_contribution(fnd, merit = c(a = 0, b = 1), G = Gs,
                               target_coancestry = 2e-16),
    "above the coancestry of the unconstrained")
  expect_identical(oc$lambda, 0)
  # exactly the optimum stays silent
  expect_no_warning(
    oc0 <- optimum_contribution(fnd, merit = c(a = 0, b = 1), G = Gs,
                                target_coancestry = 1e-16))
  expect_identical(oc0$lambda, 0)
  # ordinary scale behaviour is unchanged
  expect_warning(
    optimum_contribution(fnd, merit = c(a = 0, b = 1), G = .fx3_G(),
                         target_coancestry = 0.5000005),
    "above the coancestry of the unconstrained")
  expect_no_warning(
    optimum_contribution(fnd, merit = c(a = 0, b = 1), G = .fx3_G(),
                         target_coancestry = 0.5))
})

# ---- R3-10 ------------------------------------------------------------------

test_that("E[A|P] must be linear for R = i Cov(A,P)/sigma_P (R3-10)", {
  # A = (0, 1, 2) equiprobable, P = A^2, truncation of the top third.
  A <- c(0, 1, 2); P <- A^2
  mu_A <- mean(A)
  actual <- mean(A[P == max(P)]) - mu_A              # realised response
  varP <- mean(P^2) - mean(P)^2
  covAP <- mean(A * P) - mean(A) * mean(P)
  i <- (mean(P[P == max(P)]) - mean(P)) / sqrt(varP)  # standardised differential
  linear <- i * covAP / sqrt(varP)                    # i Cov(A,P)/sigma_P
  expect_equal(actual, 1)
  expect_equal(linear, 1.0769230769, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(actual, linear)))
  # the roxygen carries the assumption
  src <- readLines(test_path("..", "..", "R", "select_ind.R"), warn = FALSE)
  txt <- paste(src[seq_len(60)], collapse = " ")
  expect_match(txt, "linear", fixed = TRUE)
})

# ---- R3-11 ------------------------------------------------------------------

test_that("additive-only re-standardised accuracy is 'near', not 'held constant' (R3-11)", {
  src <- readLines(test_path("..", "..", "R", "select_schemes.R"), warn = FALSE)
  txt <- paste(src, collapse = " ")
  expect_false(grepl("also held constant", txt, fixed = TRUE))
  expect_match(txt, "stays near", fixed = TRUE)
})
