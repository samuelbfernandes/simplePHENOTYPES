# test-fix2-prediction.R -- round-2 regression tests for the prediction group
# (Codex review of the round-1 fixes): B4, B15, B16, C7, C8, C9, C10.

.f2_snp <- function(idx, pool = NA_character_) {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES", envir = environment())
  as_population(get("SNP55K_maize282_maf04", envir = environment()),
                individuals = idx, pool = pool)
}

.f2_q <- c("ss196442916", "ss196439337", "ss196480535")

# ---- B4: founder_f keyed by pedigree key, not display id -----------------------

.f2_ab_backcrosses <- function() {
  # two founders sharing the display id "4226" in different pools; the analysed
  # population holds only their descendants, so both are ancestry-only founders
  pa <- .f2_snp(1, pool = "A")
  pb <- .f2_snp(1, pool = "B")
  f1 <- cross(pa, pb, n = 1, seed = 1)
  bca <- cross(f1, pa, n = 1, seed = 2)
  bcb <- cross(f1, pb, n = 1, seed = 3)
  c(bca, bcb)
}

test_that("B4: same-named founders in two pools take different F when named by key", {
  pop <- .f2_ab_backcrosses()
  ped <- .ensure_pedigree(pop)$pedigree
  fnd <- ped[ped$design == "founder", ]
  expect_equal(nrow(fnd), 2L)
  expect_equal(fnd$id[1], fnd$id[2])                 # the display ids collide
  kA <- fnd$key[fnd$pool == "A"]
  kB <- fnd$key[fnd$pool == "B"]
  A <- a_matrix(pop, founder_f = stats::setNames(c(1, 0), c(kA, kB)))
  # F1 = A x B; A[F1, A] = (2 + 0) / 2 = 1 -> 1 + 1/2; A[F1, B] = (0 + 1) / 2
  expect_equal(unname(diag(A)), c(1.5, 1.25))
  # swapping the keys swaps the diagonals
  A2 <- a_matrix(pop, founder_f = stats::setNames(c(0, 1), c(kA, kB)))
  expect_equal(unname(diag(A2)), c(1.25, 1.5))
  # naming a single founder by key leaves the other non-inbred
  A3 <- a_matrix(pop, founder_f = stats::setNames(1, kA))
  expect_equal(unname(diag(A3)), c(1.5, 1.25))
  # scalar semantics unchanged
  expect_equal(unname(diag(a_matrix(pop, founder_f = 1))), c(1.5, 1.5))
  expect_equal(unname(diag(a_matrix(pop))), c(1.25, 1.25))
})

test_that("B4: an id shared by several founders is an ambiguity error, never applied to both", {
  pop <- .f2_ab_backcrosses()
  expect_error(a_matrix(pop, founder_f = c("4226" = 1)),
               "ambiguous.*2 founders share that id")
  ped <- .ensure_pedigree(pop)$pedigree
  kA <- ped$key[ped$pool %in% "A"]
  # two names resolving to the same founder (key + ambiguous-free key twice)
  expect_error(a_matrix(pop, founder_f = stats::setNames(c(1, 1), c(kA, kA))),
               "duplicated names")
  # a key of a non-founder is not a founder
  kc <- ped$key[ped$design == "cross"][1]
  expect_error(a_matrix(pop, founder_f = stats::setNames(1, kc)),
               "not founder ids or keys")
})

test_that("B4: unique display ids still work (id interface unchanged)", {
  p2 <- .f2_snp(1:2)
  f1 <- cross(p2[1], p2[2], n = 1, seed = 1)
  all <- c(p2, f1)
  An <- a_matrix(all, founder_f = c("4226" = 1))
  expect_equal(unname(diag(An)), c(2, 1, 1))
  k <- .ensure_pedigree(p2)$keys
  Ak <- a_matrix(all, founder_f = stats::setNames(1, k[1]))
  expect_equal(An, Ak)
  # an id and a key naming the same founder are one founder named twice
  expect_error(a_matrix(all, founder_f = stats::setNames(c(1, 1), c("4226", k[1]))),
               "same founder twice")
})

# ---- B15: supplied-K symmetrization ---------------------------------------------

test_that("B15: an exactly symmetric valid K is returned bit-for-bit (subnormals kept)", {
  K <- diag(c(4.9406564584124654e-324, 1))
  Ks <- .check_relationship(K)
  expect_identical(Ks, K)
  expect_gt(Ks[1, 1], 0)
  K2 <- diag(c(1e-320, 2))
  expect_identical(.check_relationship(K2), K2)
  # the overflow-safe path is kept
  Kb <- matrix(1.7e308, 2, 2)
  Kbs <- .check_relationship(Kb)
  expect_true(all(is.finite(Kbs)))
  expect_identical(Kbs, Kb)
  # rounding-level asymmetry is still removed, and only there
  Ka <- matrix(c(2, 0.5, 0.5 + 1e-12, 1), 2, 2)
  Kas <- .check_relationship(Ka)
  expect_true(isSymmetric(Kas, tol = 0))
  expect_identical(diag(Kas), c(2, 1))
  expect_equal(Kas[1, 2], 0.5 + 5e-13, tolerance = 1e-15)
  # near-overflow asymmetry uses the safe average
  Ko <- matrix(c(1.7e308, 1.7e308, 1.7e308 * (1 - 1e-12), 1.7e308), 2, 2)
  expect_true(all(is.finite(.check_relationship(Ko))))
})

# ---- B16: duplicate names whenever either vector is named -----------------------

test_that("B16: prediction_accuracy rejects a duplicated name when only one vector is named", {
  expect_error(prediction_accuracy(c(a = 1, a = 100, b = 2, c = 3), c(1, 2, 3, 4)),
               "duplicated")
  expect_error(prediction_accuracy(c(1, 2, 3, 4), c(a = 1, a = 100, b = 2, c = 3)),
               "duplicated")
  # both named (unchanged) and clean input (unchanged)
  expect_error(prediction_accuracy(c(a = 1, a = 2, b = 2), c(a = 1, b = 2, c = 3)),
               "duplicated")
  r <- prediction_accuracy(c(a = 1, b = 2, c = 3), c(1.1, 1.9, 3.2))
  expect_equal(r[["n"]], 3)
  r2 <- prediction_accuracy(c(1, 2, 3), c(1.1, 1.9, 3.2))
  expect_equal(unname(r[["accuracy"]]), unname(r2[["accuracy"]]))
})

# ---- C7: half-sib mean under dominance ------------------------------------------

test_that("C7: dominance is in the common constant, not the parent-dependent contrast at p = 1/2", {
  # locus with a = 0, d = 1; mates' counted-allele frequency p
  half_sib_mean <- function(g_i, p, a = 0, d = 1) {
    # expected progeny value of a parent with gamete frequency g_i mated at
    # random to gametes of frequency p (the package's cross-mean expansion)
    drop(.expected_cross_means(matrix(g_i, 1), matrix(p, 1), a, d))
  }
  for (p in c(0.2, 0.5)) {
    m <- vapply(c(0, 0.5, 1), half_sib_mean, numeric(1), p = p)   # aa, Aa, AA
    alpha <- 0 + 1 * (1 - 2 * p)                                  # a + d(q - p)
    expect_equal(m[3] - m[1], alpha)                              # parent contrast
    if (p == 0.2) expect_equal(m, c(0.2, 0.5, 0.8))
    if (p == 0.5) expect_equal(m, rep(0.5, 3))                    # all equal, not 0
    # the common constant: the mean over parents in HWE, a(p - q) + 2 d p q
    q <- 1 - p
    expect_equal(sum(m * c(q^2, 2 * p * q, p^2)), 0 * (p - q) + 2 * 1 * p * q)
  }
})

# ---- C8: simulated SCA, fully inbred parents ------------------------------------

test_that("C8: fully inbred parents give zero simulated SCA only without a residual", {
  pop <- .f2_snp(1:40)
  dos <- dosages(pop)[.f2_q, , drop = FALSE]
  hom <- which(colSums(dos == 0) == 0)                  # homozygous at every QTN
  expect_gte(length(hom), 5L)
  cand <- pop[hom[1:3]]; tst <- pop[hom[4:5]]
  ca0 <- combining_ability(cand, tst, .f2_q, c(1, 0.5, 0.25), d = 0,
                           method = "simulated", n_progeny = 1, seed = 1)
  expect_equal(max(abs(ca0$sca)), 0, tolerance = 1e-12)
  ca1 <- combining_ability(cand, tst, .f2_q, c(1, 0.5, 0.25), d = 0,
                           method = "simulated", n_progeny = 1, var_e = 1,
                           seed = 1)
  expect_gt(max(abs(ca1$sca)), 0.05)
})

# ---- C9: ref rescales var_a and var_e together ----------------------------------

test_that("C9: multiplying the ref variance leaves lambda, EBV and reliability unchanged", {
  pop <- .f2_snp(1:12)
  set.seed(9); y <- stats::setNames(stats::rnorm(8), pop$ids[1:8])
  ref <- stats::rnorm(20)
  e1 <- predict_ebv(pop, y, method = "pedigree", h2 = 0.4, ref = ref)
  e2 <- predict_ebv(pop, y, method = "pedigree", h2 = 0.4, ref = ref * 100)
  expect_equal(attr(e1, "lambda"), attr(e2, "lambda"))
  expect_equal(attr(e1, "lambda"), (1 - 0.4) / 0.4)
  expect_equal(as.numeric(e1), as.numeric(e2))
  expect_equal(attr(e1, "reliability"), attr(e2, "reliability"))
  # only the absolute components change
  expect_equal(attr(e2, "var_a") / attr(e1, "var_a"), 1e4)
  expect_equal(attr(e2, "var_e") / attr(e1, "var_e"), 1e4)
})

# ---- C10: one record ------------------------------------------------------------

test_that("C10: one record gives reliability 0, and NA where the prior variance is zero", {
  pop <- .f2_snp(1:3)
  K <- diag(c(1, 0, 1)); dimnames(K) <- list(pop$ids, pop$ids)
  y <- stats::setNames(1, pop$ids[1])
  e <- predict_ebv(pop, y, K = K, var_a = 1, var_e = 1)
  expect_equal(as.numeric(e), c(0, 0, 0))
  rel <- attr(e, "reliability")
  expect_equal(as.numeric(rel)[c(1, 3)], c(0, 0))
  expect_true(is.na(rel[[2]]))
})
