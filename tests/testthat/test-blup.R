# test-blup.R -- DECISION-030: known-variance BLUP, the pedigree A matrix,
# accuracy reporting and the selection-methods manifest.

.hwe <- function(n, m = 400, seed = 1, n_chr = 10) {
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

test_that("own records with K = I: ebv = h2 (y - mean) exactly", {
  pop <- as_population(.hwe(30, m = 40))
  set.seed(2); y <- stats::setNames(stats::rnorm(20), pop$ids[1:20])
  K <- diag(30); dimnames(K) <- list(pop$ids, pop$ids)
  e <- predict_ebv(pop, y, K = K, h2 = 0.3)
  expect_equal(unname(e[1:20]), 0.3 * unname(y - mean(y)), tolerance = 1e-12)
  expect_equal(unname(e[21:30]), rep(0, 10))          # unrelated, no record
  # reliability h2 (1 - 1/n): the mean is estimated from the same n records
  expect_equal(unname(attr(e, "reliability")[1:20]), rep(0.3 * (1 - 1 / 20), 20),
               tolerance = 1e-10)
})

test_that("A-matrix known answers (tabular method; selfs and DH)", {
  fnd <- as_population(.hwe(4, m = 30))
  fs  <- cross(fnd[1], fnd[2], n = 2, seed = 1)       # full sibs
  hs  <- cross(fnd[1], fnd[3], n = 1, seed = 2)       # half sib of fs (via P1)
  s1  <- selfcross(fnd[4], n = 2, seed = 3)           # S1 sibs of a founder
  dh  <- double_haploid(fnd[4], n = 1, seed = 4)
  all <- c(fnd, fs, hs, s1, dh)
  A <- a_matrix(all)
  id <- all$ids
  expect_equal(A[id[5], id[6]], 0.5)                  # full sibs
  expect_equal(A[id[5], id[7]], 0.25)                 # half sibs
  expect_equal(A[id[1], id[5]], 0.5)                  # parent - offspring
  expect_equal(A[id[8], id[8]], 1.5)                  # S1 of a non-inbred parent
  expect_equal(A[id[8], id[9]], 1)                    # S1 sibs
  expect_equal(A[id[10], id[10]], 2)                  # doubled haploid
  expect_equal(A[id[4], id[10]], 1)                   # DH with its parent
  expect_true(isSymmetric(A))
  expect_equal(unname(diag(A)[1:4]), rep(1, 4))       # founders are the base
})

test_that("G matches A on average by relative class (fixed founder base)", {
  fnd <- as_population(.hwe(40, m = 2000, n_chr = 20, seed = 5))
  plan <- mating_design(fnd, design = "random", n_crosses = 60,
                        progeny_per_cross = 2, seed = 6)
  prog <- mate(plan, fnd, seed = 7)
  base <- rowMeans(dosages(fnd) + 1) / 2
  G <- g_matrix(prog, base_freq = base)
  A <- a_matrix(prog)
  fsib <- A == 0.5 & upper.tri(A)
  unrel <- A == 0 & upper.tri(A)
  expect_equal(mean(G[fsib]), 0.5, tolerance = 0.06 / 0.5)
  expect_lt(abs(mean(G[unrel])), 0.03)
})

test_that("the GLS solution equals Henderson's mixed-model equations", {
  fnd <- as_population(.hwe(6, m = 30, seed = 8))
  prog <- mate(mating_design(fnd, design = "half_diallel"), fnd, seed = 9)
  pop <- c(fnd, prog)
  A <- a_matrix(pop)
  set.seed(10); y <- stats::setNames(stats::rnorm(15), prog$ids)
  e <- predict_ebv(pop, y, method = "pedigree", var_a = 2, var_e = 3)
  # MME: [1'1 1'Z; Z'1 Z'Z + lambda A^-1] [mu; u] = [1'y; Z'y]
  n <- nrow(A); Z <- matrix(0, 15, n); Z[cbind(1:15, match(names(y), pop$ids))] <- 1
  X <- matrix(1, 15, 1); lam <- 3 / 2
  C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + lam * solve(A)))
  sol <- solve(C, c(crossprod(X, y), crossprod(Z, y)))
  expect_equal(attr(e, "mu"), as.numeric(sol[1]), tolerance = 1e-10)
  expect_equal(as.numeric(e), as.numeric(sol[-1]), tolerance = 1e-10)
  # reliability from the inverse coefficient matrix: 1 - C^uu_ii lambda / A_ii
  Cinv <- solve(C)[-1, -1]
  expect_equal(unname(attr(e, "reliability")),
               unname(1 - diag(Cinv) * lam / diag(A)), tolerance = 1e-10)
})

test_that("GBLUP marker effects reproduce the GEBVs", {
  pop <- as_population(.hwe(50, m = 60, seed = 11))
  set.seed(12); y <- stats::setNames(stats::rnorm(35), pop$ids[1:35])
  e <- predict_ebv(pop, y, h2 = 0.4)
  M <- t(dosages(pop)) + 1; p <- colMeans(M) / 2
  Zc <- sweep(M, 2, 2 * p, "-")
  expect_equal(as.numeric(Zc %*% attr(e, "marker_effects")), as.numeric(e),
               tolerance = 1e-8)
})

test_that("GBLUP is unbiased under its model (slope of truth on ebv ~ 1)", {
  pop <- as_population(.hwe(300, m = 300, seed = 13))
  M <- t(dosages(pop)) + 1; p <- colMeans(M) / 2
  Zc <- sweep(M, 2, 2 * p, "-")
  slopes <- vapply(1:20, function(s) {
    set.seed(100 + s)
    u <- stats::rnorm(300)                        # every marker a QTN: the model
    g <- as.numeric(Zc %*% u); ve <- stats::var(g)   # h2 = 0.5
    y <- stats::setNames(g + stats::rnorm(300, sd = sqrt(ve)), pop$ids)
    tr <- sample(300, 200)
    e <- predict_ebv(pop, y[tr], var_a = ve, var_e = ve)
    test <- setdiff(seq_len(300), tr)
    prediction_accuracy(e[test], stats::setNames(g, pop$ids)[test])[["slope"]]
  }, numeric(1))
  expect_equal(mean(slopes), 1, tolerance = 0.12)
})

test_that("accuracy rises with training size and with h2", {
  pop <- as_population(.hwe(400, m = 300, seed = 14))
  M <- t(dosages(pop)) + 1; p <- colMeans(M) / 2
  Zc <- sweep(M, 2, 2 * p, "-")
  acc <- function(n_train, h2, s) {
    set.seed(s); u <- stats::rnorm(300); g <- as.numeric(Zc %*% u)
    ve <- stats::var(g) * (1 - h2) / h2
    y <- stats::setNames(g + stats::rnorm(400, sd = sqrt(ve)), pop$ids)
    e <- predict_ebv(pop, y[1:n_train], var_a = stats::var(g), var_e = ve)
    prediction_accuracy(e[301:400], stats::setNames(g, pop$ids)[301:400])[["accuracy"]]
  }
  a_small <- mean(vapply(1:8, function(s) acc(60, 0.5, s), numeric(1)))
  a_big   <- mean(vapply(1:8, function(s) acc(300, 0.5, s), numeric(1)))
  a_loh2  <- mean(vapply(1:8, function(s) acc(300, 0.2, s), numeric(1)))
  expect_gt(a_big, a_small)
  expect_gt(a_big, a_loh2)
})

test_that("pedigree BLUP of an unphenotyped sire reaches the progeny-test accuracy", {
  fnd <- as_population(.hwe(1200, m = 300, seed = 15))
  sires <- fnd[1:120]; dams <- fnd[121:1200]
  n <- 9; h2 <- 0.25
  set.seed(16); a <- stats::rnorm(300, sd = 0.2)
  x <- dosages(fnd) + 1; pb <- rowMeans(x) / 2
  bv <- colSums((x[, 1:120] - 2 * pb) * a)
  r <- vapply(1:4, function(s) {
    pt <- progeny_test(sires, dams, qtn = 1:300, a = a, n_progeny = n, h2 = h2,
                       seed = 20 + s)
    prog <- attr(pt, "progeny")
    all <- c(sires, dams, prog)
    vp <- stats::var(attr(pt, "records"))
    e <- predict_ebv(all, attr(pt, "records"), method = "pedigree",
                     var_a = h2 * vp, var_e = (1 - h2) * vp)
    stats::cor(e[sires$ids], bv)
  }, numeric(1))
  r_pt <- sqrt(n * h2 / (4 + (n - 1) * h2))
  expect_lt(abs(mean(r) - r_pt), 0.08)
})

test_that("manifest and input checks", {
  sm <- selection_methods()
  expect_true(all(c("id", "label", "fn", "params", "source") %in% names(sm)))
  expect_equal(anyDuplicated(sm$id), 0L)
  pop <- as_population(.hwe(10, m = 20))
  y <- stats::setNames(1:5 + 0.1, pop$ids[1:5])
  expect_error(predict_ebv(pop, y), "give `h2`")
  expect_error(predict_ebv(pop, y, h2 = 0.5, var_a = 1, var_e = 1), "not both")
  expect_error(predict_ebv(pop, c(zz = 1), h2 = 0.5), "distinct ids")
  expect_error(prediction_accuracy(1:2, 1:2), "at least three")
})

test_that("a supplied K is validated; marker effects only for an internal G (review C)", {
  g <- data.frame(snp = paste0("m", 1:20), allele = "A/G", chr = 1, pos = 1:20,
                  cm = seq(0, 95, 5), stringsAsFactors = FALSE)
  set.seed(1)
  for (i in 1:3) g[[paste0("I", i)]] <- sample(c(-1L, 0L, 1L), 20, TRUE)
  pop <- as_population(g)
  y <- c(I1 = 1, I2 = 2)
  bad <- diag(3); bad[1, 2] <- bad[2, 1] <- 2              # eigenvalues 3, 1, -1
  dimnames(bad) <- list(pop$ids, pop$ids)
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 0.5, K = bad),
               "positive semidefinite")
  nf <- diag(3); nf[3, 3] <- NaN
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 0.5, K = nf), "finite")
  asym <- diag(3); asym[1, 2] <- 0.3
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 0.5, K = asym), "symmetric")
  G <- g_matrix(pop)
  own <- predict_ebv(pop, y, var_a = 1, var_e = 0.5)
  sup <- predict_ebv(pop, y, var_a = 1, var_e = 0.5, K = G)
  expect_equal(as.numeric(own), as.numeric(sup))
  expect_false(is.null(attr(own, "marker_effects")))
  expect_null(attr(sup, "marker_effects"))
})

test_that("a supplied K must name its individuals on both axes (review C r2)", {
  g <- data.frame(snp = paste0("m", 1:20), allele = "A/G", chr = 1, pos = 1:20,
                  cm = seq(0, 95, 5), stringsAsFactors = FALSE)
  set.seed(1)
  for (i in 1:3) g[[paste0("I", i)]] <- sample(c(-1L, 0L, 1L), 20, TRUE)
  pop <- as_population(g)
  y <- c(I1 = 1, I2 = 4)
  K <- diag(c(1, 2, 3)); dimnames(K) <- list(pop$ids, pop$ids)
  ok <- predict_ebv(pop, y, var_a = 1, var_e = 1, K = K)
  expect_equal(unname(as.numeric(ok)), c(-0.6, 1.2, 0), tolerance = 1e-12)
  # the same matrix listed in another order is reordered, not misread
  o <- c(3, 1, 2)
  expect_equal(predict_ebv(pop, y, var_a = 1, var_e = 1, K = K[o, o]), ok)
  wrong <- K; colnames(wrong) <- pop$ids[c(2, 1, 3)]
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 1, K = wrong), "both axes")
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 1, K = unname(K)), "both axes")
})

test_that("unrepresentable variance ratios error instead of returning NaN (review C r3)", {
  g <- data.frame(snp = paste0("m", 1:20), allele = "A/G", chr = 1, pos = 1:20,
                  cm = seq(0, 95, 5), stringsAsFactors = FALSE)
  set.seed(1)
  for (i in 1:4) g[[paste0("I", i)]] <- sample(c(-1L, 0L, 1L), 20, TRUE)
  pop <- as_population(g)
  K <- diag(4); dimnames(K) <- list(pop$ids, pop$ids)
  expect_error(predict_ebv(pop, c(I1 = 1, I2 = 2), var_a = 1e-320, var_e = 1e308,
                           K = K), "not representable")
  expect_error(predict_ebv(pop, c(I1 = 1, I2 = 2), var_a = 1e308, var_e = 1e-320,
                           K = K), "not representable")
})

test_that("K validation is local; the manifest names real arguments (review C r5)", {
  g <- data.frame(snp = paste0("m", 1:20), allele = "A/G", chr = 1, pos = 1:20,
                  cm = seq(0, 95, 5), stringsAsFactors = FALSE)
  set.seed(1)
  for (i in 1:3) g[[paste0("I", i)]] <- sample(c(-1L, 0L, 1L), 20, TRUE)
  pop <- as_population(g)
  y <- c(I1 = 1, I2 = 2)
  nm <- function(K) { dimnames(K) <- list(pop$ids, pop$ids); K }
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 0.5,
                           K = nm(diag(c(-9, 1e9, 1)))), "negative diagonal")
  asym <- diag(c(1, 1e9, 1)); asym[1, 3] <- 5
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 0.5, K = nm(asym)),
               "symmetric")
  np <- diag(c(1, 1e9, 1)); np[1, 3] <- np[3, 1] <- 2   # |corr| = 2 > 1
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 0.5, K = nm(np)),
               "positive semidefinite")
  zr <- diag(c(0, 1, 1)); zr[1, 2] <- zr[2, 1] <- 0.1    # zero variance, nonzero cov
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 0.5, K = nm(zr)),
               "positive semidefinite")
  ok <- diag(c(0, 1, 1))                                 # a zero-variance row is fine
  expect_silent(predict_ebv(pop, y, var_a = 1, var_e = 0.5, K = nm(ok)))
  # every manifest parameter is an argument of the function it names
  m <- selection_methods()
  for (i in seq_len(nrow(m))) {
    fns <- regmatches(m$fn[i], gregexpr("[a-z_]+(?=[(])", m$fn[i], perl = TRUE))[[1]]
    fns <- setdiff(fns, "c")
    ps <- trimws(unlist(strsplit(m$params[i], ",|[|]|[+]")))
    fa <- unique(unlist(lapply(fns, function(f) names(formals(get(f))))))
    expect_true(all(ps %in% fa), info = m$id[i])
  }
})

test_that("an indefinite K the ridge cannot absorb errors; tiny variances are valid (review C r6)", {
  g <- data.frame(snp = paste0("m", 1:20), allele = "A/G", chr = 1, pos = 1:20,
                  cm = seq(0, 95, 5), stringsAsFactors = FALSE)
  set.seed(1)
  for (i in 1:2) g[[paste0("I", i)]] <- sample(c(-1L, 0L, 1L), 20, TRUE)
  pop <- as_population(g)
  y <- c(I1 = 1, I2 = 2)
  nm <- function(K) { dimnames(K) <- list(pop$ids, pop$ids); K }
  e <- 5e-9
  K <- nm(matrix(c(1, 1 + e, 1 + e, 1), 2))             # eigenvalue -5e-9
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 4.99e-9, K = K),
               "positive (semi)?definite")
  # beyond floating-point rounding an indefinite K is rejected whatever the ridge
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 1, K = K),
               "positive semidefinite")
  tiny <- nm(matrix(c(1e-12, 1e-8, 1e-8, 1), 2))          # PSD, min eigenvalue ~1e-12
  expect_silent(predict_ebv(pop, y, var_a = 1, var_e = 1, K = tiny))
})

test_that("tolerated asymmetry is symmetrized; a_matrix ids are identifiers (review C r7)", {
  g <- data.frame(snp = paste0("m", 1:20), allele = "A/G", chr = 1, pos = 1:20,
                  cm = seq(0, 95, 5), stringsAsFactors = FALSE)
  set.seed(1)
  for (i in 1:3) g[[paste0("I", i)]] <- sample(c(-1L, 0L, 1L), 20, TRUE)
  pop <- as_population(g)
  y <- c(I1 = 1, I2 = 3)
  K <- diag(c(1, 1, (4.5e-9)^2)); K[3, 1] <- 9e-9
  dimnames(K) <- list(pop$ids, pop$ids)
  Ks <- (K + t(K)) / 2
  a <- predict_ebv(pop, y, var_a = 1, var_e = 0.01, K = K)
  b <- predict_ebv(pop, y, var_a = 1, var_e = 0.01, K = Ks)
  expect_identical(a, b)
  expect_lte(max(attr(a, "reliability"), na.rm = TRUE), 1)
  gp <- g; names(gp)[6:8] <- c("2", "1", "3")
  q <- as_population(gp)
  expect_error(a_matrix(q, ids = 1), "character vector")
  expect_equal(dimnames(a_matrix(q, ids = "1")), list("1", "1"))
})

test_that("near-PSD K cannot give reliability above 1 (review C r8)", {
  g <- data.frame(snp = paste0("m", 1:20), allele = "A/G", chr = 1, pos = 1:20,
                  cm = seq(0, 95, 5), stringsAsFactors = FALSE)
  set.seed(1)
  for (i in 1:3) g[[paste0("I", i)]] <- sample(c(-1L, 0L, 1L), 20, TRUE)
  pop <- as_population(g)
  y <- c(I1 = 1, I2 = 2)
  K <- matrix(c(1, 1 - 1e-12, 1e-4 / sqrt(2),
                1 - 1e-12, 1, -1e-4 / sqrt(2),
                1e-4 / sqrt(2), -1e-4 / sqrt(2), 1e-8), 3)
  dimnames(K) <- list(pop$ids, pop$ids)
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 1e-16, K = K),
               "positive semidefinite|outside \\[0, 1\\]|not positive definite")
  # a rank-deficient but valid genomic matrix is accepted
  G <- g_matrix(pop)
  e <- predict_ebv(pop, y, var_a = 1, var_e = 0.5, K = G)
  r <- attr(e, "reliability")
  expect_true(all(r[!is.na(r)] >= 0 & r[!is.na(r)] <= 1))
})

test_that("a subnormal but valid variance does not overflow the PSD check (review C r9)", {
  g <- data.frame(snp = paste0("m", 1:20), allele = "A/G", chr = 1, pos = 1:20,
                  cm = seq(0, 95, 5), stringsAsFactors = FALSE)
  set.seed(1)
  for (i in 1:2) g[[paste0("I", i)]] <- sample(c(-1L, 0L, 1L), 20, TRUE)
  pop <- as_population(g)
  K <- diag(c(1e-320, 1)); dimnames(K) <- list(pop$ids, pop$ids)
  e <- predict_ebv(pop, c(I1 = 1, I2 = 2), var_a = 1, var_e = 1, K = K)
  expect_equal(attr(e, "mu"), 4 / 3)
  expect_equal(unname(e[2]), 1 / 3)
  expect_true(abs(e[1]) < 1e-300)
})
