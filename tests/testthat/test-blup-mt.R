# test-blup-mt.R -- DECISION-032: multi-trait known-covariance BLUP through
# predict_ebv() with an individuals x traits `pheno` matrix.

.hwe_mt <- function(n, m = 200, seed = 1, n_chr = 5) {
  set.seed(seed)
  p <- stats::runif(m, 0.2, 0.8)
  g <- t(vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1L, integer(n)))
  colnames(g) <- paste0("I", seq_len(n))
  cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                   chr = rep(seq_len(n_chr), each = m / n_chr),
                   pos = rep(seq_len(m / n_chr), n_chr),
                   cm = rep(seq(0, 100, length.out = m / n_chr), n_chr),
                   stringsAsFactors = FALSE), as.data.frame(g))
}
.cov2 <- function(v1, v2, r, tr = c("t1", "t2")) {
  S <- matrix(c(v1, r * sqrt(v1 * v2), r * sqrt(v1 * v2), v2), 2)
  dimnames(S) <- list(tr, tr)
  S
}

test_that("multi-trait BLUP equals Henderson's multi-trait mixed-model equations", {
  fnd <- as_population(.hwe_mt(6, m = 30, seed = 8))
  prog <- mate(mating_design(fnd, design = "half_diallel"), fnd, seed = 9)
  pop <- c(fnd, prog)
  A <- a_matrix(pop)
  n <- nrow(A)
  set.seed(10)
  Y <- matrix(stats::rnorm(2 * n, 5), n, 2,
              dimnames = list(pop$ids, c("t1", "t2")))
  Y[1:6, ] <- NA                          # founders unrecorded
  Y[7:10, "t1"] <- NA                     # some progeny recorded on t2 only
  Y[11:13, "t2"] <- NA                    # and some on t1 only
  G0 <- .cov2(2, 1, 0.6); R0 <- .cov2(3, 2, -0.3)
  e <- predict_ebv(pop, Y, method = "pedigree", var_a = G0, var_e = R0)
  # MME over the stacked records, u ordered trait-major, (G0 (x) A)^-1
  obs <- which(!is.na(Y), arr.ind = TRUE)
  ia <- obs[, 1]; ta <- obs[, 2]; y <- Y[obs]; nr <- length(y)
  X <- outer(ta, 1:2, "==") * 1
  Z <- matrix(0, nr, 2 * n); Z[cbind(seq_len(nr), (ta - 1) * n + ia)] <- 1
  Ri <- solve(R0[ta, ta] * outer(ia, ia, "=="))
  C <- rbind(cbind(t(X) %*% Ri %*% X, t(X) %*% Ri %*% Z),
             cbind(t(Z) %*% Ri %*% X,
                   t(Z) %*% Ri %*% Z + kronecker(solve(G0), solve(A))))
  sol <- solve(C, c(t(X) %*% Ri %*% y, t(Z) %*% Ri %*% y))
  expect_equal(unname(attr(e, "mu")), sol[1:2], tolerance = 1e-9)
  expect_equal(as.numeric(e), sol[-(1:2)], tolerance = 1e-9)
  # reliability 1 - PEV / (G0_tt A_ii), PEV the diagonal of C^uu
  pev <- diag(solve(C))[-(1:2)]
  expect_equal(as.numeric(attr(e, "reliability")),
               unname(1 - pev / (rep(diag(G0), each = n) * rep(diag(A), 2))),
               tolerance = 1e-9)
  expect_identical(dimnames(e), list(pop$ids, c("t1", "t2")))
})

test_that("one trait as a matrix is single-trait BLUP; uncorrelated traits separate", {
  pop <- as_population(.hwe_mt(40, seed = 2))
  set.seed(3)
  y1 <- stats::setNames(stats::rnorm(30), pop$ids[1:30])
  y2 <- stats::setNames(stats::rnorm(30, 2), pop$ids[1:30])
  s1 <- predict_ebv(pop, y1, var_a = 2, var_e = 3)
  s2 <- predict_ebv(pop, y2, var_a = 0.5, var_e = 1)
  m1 <- predict_ebv(pop, cbind(t1 = y1), var_a = 2, var_e = 3)
  expect_equal(as.numeric(m1), as.numeric(s1), tolerance = 1e-10)
  expect_equal(as.numeric(attr(m1, "reliability")),
               as.numeric(attr(s1, "reliability")), tolerance = 1e-10)
  expect_equal(unname(attr(m1, "mu")), attr(s1, "mu"), tolerance = 1e-10)
  expect_equal(unname(attr(m1, "marker_effects")[, 1]),
               unname(attr(s1, "marker_effects")), tolerance = 1e-10)
  G0 <- diag(c(2, 0.5)); R0 <- diag(c(3, 1))
  m2 <- predict_ebv(pop, cbind(t1 = y1, t2 = y2), var_a = G0, var_e = R0)
  expect_equal(as.numeric(m2[, "t1"]), as.numeric(s1), tolerance = 1e-10)
  expect_equal(as.numeric(m2[, "t2"]), as.numeric(s2), tolerance = 1e-10)
  expect_equal(unname(attr(m2, "reliability")[, "t2"]),
               unname(attr(s2, "reliability")), tolerance = 1e-10)
})

test_that("GBLUP marker effects reproduce the multi-trait GEBVs", {
  pop <- as_population(.hwe_mt(40, seed = 4))
  set.seed(5)
  Y <- matrix(stats::rnorm(80), 40, 2, dimnames = list(pop$ids, c("a", "b")))
  Y[31:40, "a"] <- NA
  e <- predict_ebv(pop, Y, var_a = .cov2(1, 2, 0.7, c("a", "b")),
                   var_e = .cov2(1, 1, 0.2, c("a", "b")))
  M <- t(dosages(pop)) + 1
  p <- colMeans(M) / 2
  Zc <- sweep(M, 2, 2 * p, "-")
  expect_equal(unname(Zc %*% attr(e, "marker_effects")), unname(e[, ]),
               tolerance = 1e-8)
})

# u ~ N(0, G0 (x) K), e ~ N(0, R0 (x) I): the model BLUP assumes
.draw_mt <- function(K, G0, R0, seed) {
  set.seed(seed)
  n <- nrow(K)
  ek <- eigen(K, symmetric = TRUE)
  LK <- ek$vectors %*% diag(sqrt(pmax(ek$values, 0)))
  U <- LK %*% matrix(stats::rnorm(2 * n), n, 2) %*% chol(G0)
  E <- matrix(stats::rnorm(2 * n), n, 2) %*% chol(R0)
  dimnames(U) <- dimnames(E) <- list(rownames(K), colnames(G0))
  list(u = U, y = 10 + U + E)
}

test_that("a correlated recorded trait predicts an unrecorded one, unbiasedly", {
  pop <- as_population(.hwe_mt(300, m = 400, seed = 6, n_chr = 10))
  K <- g_matrix(pop)
  G0 <- .cov2(1, 1, 0.8); R0 <- .cov2(1, 1, 0.2)
  cand <- pop$ids[201:300]
  acc_mt <- acc_st <- slope <- rel <- numeric(0)
  for (r in 1:10) {
    d <- .draw_mt(K, G0, R0, seed = 100 + r)
    Y <- d$y
    Y[cand, "t1"] <- NA                   # t1 unrecorded on the candidates
    mt <- predict_ebv(pop, Y, var_a = G0, var_e = R0)
    st <- predict_ebv(pop, Y[setdiff(pop$ids, cand), "t1"], var_a = 1,
                      var_e = 1)
    a <- prediction_accuracy(mt[cand, "t1"], d$u[cand, "t1"])
    acc_mt <- c(acc_mt, a[["accuracy"]]); slope <- c(slope, a[["slope"]])
    rel <- c(rel, mean(attr(mt, "reliability")[cand, "t1"]))
    acc_st <- c(acc_st, prediction_accuracy(st[cand], d$u[cand, "t1"])[["accuracy"]])
  }
  # the candidates' own t2 records carry information on t1 (rg = 0.8)
  expect_gt(mean(acc_mt), mean(acc_st) + 0.1)
  # BLUP is unbiased under its own model: E[u | u_hat] = u_hat
  expect_lt(abs(mean(slope) - 1), 0.1)
  # the reported reliability is the expected squared accuracy
  expect_lt(abs(mean(acc_mt^2) - mean(rel)), 0.05)
})

test_that("singular G0 (genetic correlation 1) is fine; names reorder", {
  pop <- as_population(.hwe_mt(30, seed = 7))
  set.seed(8)
  Y <- matrix(stats::rnorm(60), 30, 2, dimnames = list(pop$ids, c("t1", "t2")))
  Y[1:10, "t2"] <- NA
  G1 <- .cov2(1, 4, 1)
  e <- predict_ebv(pop, Y, var_a = G1, var_e = diag(2))
  expect_true(all(is.finite(e)))
  # with rg = 1, u_t2 = 2 u_t1 exactly, so the predictions are proportional
  expect_equal(unname(e[, "t2"]), 2 * unname(e[, "t1"]), tolerance = 1e-8)
  rel <- attr(e, "reliability")
  expect_true(all(rel >= 0 & rel <= 1))
  # covariance matrices named in another order are reordered to pheno's traits
  G0 <- .cov2(1, 2, 0.5); R0 <- .cov2(2, 1, 0.1)
  a <- predict_ebv(pop, Y, var_a = G0, var_e = R0)
  b <- predict_ebv(pop, Y, var_a = G0[2:1, 2:1], var_e = R0[2:1, 2:1])
  expect_equal(a[, ], b[, ], tolerance = 1e-12)
  # a data frame of trait columns is accepted
  d <- predict_ebv(pop, as.data.frame(Y), var_a = G0, var_e = R0)
  expect_equal(d[, ], a[, ], tolerance = 1e-12)
})

test_that("multi-trait inputs are validated", {
  pop <- as_population(.hwe_mt(10, m = 40, seed = 9))
  Y <- matrix(1:20 / 3, 10, 2, dimnames = list(pop$ids, c("t1", "t2")))
  G0 <- .cov2(1, 1, 0.5); R0 <- diag(2)
  expect_error(predict_ebv(pop, Y, h2 = 0.5), "covariance matrices")
  expect_error(predict_ebv(pop, Y, var_a = G0), "`var_e`")
  expect_error(predict_ebv(pop, Y, var_a = 1, var_e = R0), "2 x 2")
  bad <- G0; dimnames(bad) <- list(c("t1", "x"), c("t1", "x"))
  expect_error(predict_ebv(pop, Y, var_a = bad, var_e = R0), "names")
  expect_error(predict_ebv(pop, Y, var_a = .cov2(1, 1, 1.5), var_e = R0),
               "positive semidefinite")
  expect_error(predict_ebv(pop, Y, var_a = diag(c(1, 0)), var_e = R0),
               "positive variance")
  Yn <- Y; Yn[1, 1] <- NaN
  expect_error(predict_ebv(pop, Yn, var_a = G0, var_e = R0), "finite")
  Yi <- Y; Yi[1, 1] <- Inf
  expect_error(predict_ebv(pop, Yi, var_a = G0, var_e = R0), "finite")
  Ye <- Y; Ye[, "t2"] <- NA
  expect_error(predict_ebv(pop, Ye, var_a = G0, var_e = R0),
               "at least one record")
  Yr <- Y; rownames(Yr)[1] <- "nobody"
  expect_error(predict_ebv(pop, Yr, var_a = G0, var_e = R0), "row names")
  Yc <- Y; colnames(Yc) <- NULL
  expect_error(predict_ebv(pop, Yc, var_a = G0, var_e = R0), "trait")
  # a residual covariance that makes the record covariance singular
  expect_error(predict_ebv(pop, Y, var_a = G0, var_e = .cov2(1, 1, 1),
                           K = { K <- matrix(0, 10, 10)
                                 dimnames(K) <- list(pop$ids, pop$ids); K }),
               "not positive definite")
})

test_that("one-marker marker effects; an overflowing indefinite matrix is rejected (review r1)", {
  g <- data.frame(snp = "m1", allele = "A/G", chr = 1, pos = 1, cm = 0,
                  stringsAsFactors = FALSE)
  for (i in 1:6) g[[paste0("I", i)]] <- c(-1L, 0L, 1L, 1L, 0L, -1L)[i]
  pop <- as_population(g)
  Y <- matrix(c(1, 2, 3, 2, 1, 0, 0, 1, 1, 2, 2, 3), 6, 2,
              dimnames = list(pop$ids, c("t1", "t2")))
  e <- predict_ebv(pop, Y, var_a = .cov2(1, 1, 0.3), var_e = diag(2))
  me <- attr(e, "marker_effects")
  expect_identical(dim(me), c(1L, 2L))
  M <- t(dosages(pop)) + 1
  expect_equal(unname((M - 2 * colMeans(M) / 2) %*% me), unname(e[, ]),
               tolerance = 1e-10)
  # eigenvalues ~ +-1e200: the correlation-scale entries overflow, which must
  # fail the PSD check rather than make its tolerance infinite
  pop4 <- as_population(.hwe_mt(4, m = 20, seed = 11))
  K <- diag(4); dimnames(K) <- list(pop4$ids, pop4$ids)
  Y4 <- matrix(c(1, 2, NA, NA, NA, NA, 3, 4), 4, 2,
               dimnames = list(pop4$ids, c("t1", "t2")))
  R0 <- matrix(c(1e-320, 1e200, 1e200, 1e-320), 2)
  expect_error(predict_ebv(pop4, Y4, var_a = diag(2), var_e = R0, K = K),
               "positive semidefinite")
  Kb <- diag(c(1e-320, 1e-320, 1, 1)); Kb[1, 2] <- Kb[2, 1] <- 1e200
  dimnames(Kb) <- list(pop4$ids, pop4$ids)
  expect_error(predict_ebv(pop4, stats::setNames(1:4 + 0.5, pop4$ids),
                           var_a = 1, var_e = 1, K = Kb),
               "positive semidefinite")
})
