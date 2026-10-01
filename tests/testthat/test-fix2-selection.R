# test-fix2-selection.R
#
# Round-2 fixes (Codex review): B5 OCS merit-shift invariance, B6 floating-point
# target-above-optimum comparison, B17 scalar `trait` validated in every
# non-culling select_ind() branch, C3 response-identity numeric check.

.fx2_founders <- function(k = 4L) {
  g <- data.frame(snp = paste0("m", seq_len(k)), allele = "A/G",
                  chr = seq_len(k), pos = 1, cm = 0, P1 = 1L, P2 = -1L,
                  stringsAsFactors = FALSE)
  as_population(g)
}
.fx2_G2 <- function() diag(2) |> `dimnames<-`(list(c("a", "b"), c("a", "b")))

# ---- B5 ---------------------------------------------------------------------

test_that("optimum_contribution() is invariant to a constant added to every merit (B5)", {
  fnd <- .fx2_founders()
  for (off in c(0, 1e6, 1e16, -1e16)) {
    expect_no_warning(
      oc <- optimum_contribution(fnd, merit = c(a = off, b = off + 2),
                                 G = .fx2_G2(), target_coancestry = 0.26))
    expect_equal(unname(oc$contributions), c(0.4, 0.6), tolerance = 1e-8,
                 info = paste("offset", off))
    expect_equal(oc$coancestry, 0.26, tolerance = 1e-8, info = paste("offset", off))
    expect_equal(oc$lambda, 10, tolerance = 1e-6, info = paste("offset", off))
  }
  # reported merit is on the original scale: c'g = off + 1.2
  oc <- optimum_contribution(fnd, merit = c(a = 1e6, b = 1e6 + 2),
                             G = .fx2_G2(), target_coancestry = 0.26)
  expect_equal(oc$merit, 1e6 + 1.2, tolerance = 1e-9)
  # direction = "low" keeps the reported sign on the original scale
  lo <- optimum_contribution(fnd, merit = c(a = 1e6, b = 1e6 + 2), G = .fx2_G2(),
                             direction = "low", lambda = 10)
  expect_equal(unname(lo$contributions), c(0.6, 0.4), tolerance = 1e-8)
  expect_equal(lo$merit, sum(lo$contributions * c(1e6, 1e6 + 2)), tolerance = 1e-9)
  # fixed lambda: the same shift changes nothing
  a <- optimum_contribution(fnd, merit = c(a = 0, b = 2), G = .fx2_G2(), lambda = 10)
  b <- optimum_contribution(fnd, merit = c(a = 1e16, b = 1e16 + 2), G = .fx2_G2(),
                            lambda = 10)
  expect_equal(a$contributions, b$contributions, tolerance = 1e-8)
})

# ---- B6 ---------------------------------------------------------------------

test_that("a target just above the unconstrained optimum warns; equality does not (B6)", {
  fnd <- .fx2_founders()
  expect_warning(
    oc <- optimum_contribution(fnd, merit = c(a = 0, b = 1), G = .fx2_G2(),
                               target_coancestry = 0.5000005),
    "above the coancestry of the unconstrained")
  expect_identical(oc$lambda, 0)
  expect_no_warning(
    oc0 <- optimum_contribution(fnd, merit = c(a = 0, b = 1), G = .fx2_G2(),
                                target_coancestry = 0.5))
  expect_identical(oc0$lambda, 0)
  expect_equal(unname(oc0$contributions), c(0, 1))
  # max_coancestry above the optimum stays a silent slack constraint
  expect_no_warning(
    optimum_contribution(fnd, merit = c(a = 0, b = 1), G = .fx2_G2(),
                         max_coancestry = 0.5000005))
})

# ---- B17 --------------------------------------------------------------------

test_that("scalar `trait` is validated in every non-culling branch (B17)", {
  data("SNP55K_maize282_maf04")
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:40)
  f2 <- suppressMessages(selfcross(cross(pop[1], pop[2], n = 1, seed = 1),
                                   n = 30, seed = 2))
  ph <- suppressMessages(
    simulate_phenotype(f2, architecture = "pleiotropy", n_traits = 2, cor = 0.3,
                       h2 = 0.5, seed = 3) |> additive(n_qtn = 10))
  msg <- "must be one trait index"
  expect_error(select_ind(ph, n = 5, method = "index", weights = c(1, 1),
                          trait = c(1, 2)), msg)
  expect_error(select_ind(ph, n = 5, method = "index", weights = c(1, 1),
                          trait = 9), msg)
  expect_error(select_ind(ph, n = 5, method = "quadratic_index",
                          weights = c(1, 1), quad_weights = diag(2),
                          trait = c(1, 2)), msg)
  expect_error(select_ind(ph, n = 5, method = "random", on = seq_len(30),
                          trait = c(1, 2)), msg)
  expect_error(select_ind(ph, n = 5, method = "mass", on = seq_len(30),
                          trait = 0), msg)
  expect_error(select_ind(ph, n = 5, method = "mass",
                          on = function(s) seq_len(30), trait = c(1, 2)), msg)
  # valid scalar traits still work in the ignoring branches; culling is untouched
  expect_equal(n_individuals(select_ind(ph, n = 5, method = "index",
                                        weights = c(1, 1), trait = 2)), 5)
  expect_equal(n_individuals(select_ind(ph, n = 5, method = "mass",
                                        on = seq_len(30), trait = 2)), 5)
  expect_no_error(select_ind(ph, method = "culling", culling = c(0.5, 0.5),
                             trait = c(1, 2)))
})

# ---- C3 ---------------------------------------------------------------------

test_that("R = i Cov(A,P)/sigma_P differs from i V_A/sigma_P under LD epistasis (C3)", {
  # random union of haplotypes: HWE per locus (p1 = 0.3, p2 = 0.4) and LD D = 0.05
  h <- rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
  f <- c(.47, .23, .13, .17)
  g <- expand.grid(i = 1:4, j = 1:4)
  w <- f[g$i] * f[g$j]
  x1 <- h[g$i, 1] + h[g$j, 1]
  x2 <- h[g$i, 2] + h[g$j, 2]
  E <- function(v) sum(w * v)
  cv <- function(a, b) E(a * b) - E(a) * E(b)
  I <- (x1 - E(x1)) * (x2 - E(x2)); I <- I - E(I)   # centred additive-by-additive
  expect_equal(c(E(x1) / 2, E(x2) / 2), c(0.3, 0.4))
  A <- x1
  expect_equal(cv(A, A), 0.42)                       # V_A
  expect_equal(cv(A, I), 0.04)                       # Cov(A, epistasis) != 0 under LD
  expect_equal(cv(A, A + I), 0.46)                   # Cov(A, P) != V_A
})
