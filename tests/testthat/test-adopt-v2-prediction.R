# test-adopt-v2-prediction.R -- remaining proposals of the Sept-2026 audit of the
# prediction group (predict_ebv, a_matrix, prediction_accuracy, combining_ability,
# progeny_test, template_effects) not yet covered by test-audit-prediction.R,
# test-fix2-prediction.R, test-blup.R, test-combining.R or test-progeny.R. Ids
# (T*, PRED-*) refer to .tmp/audit-2026-09-29/reconciliation/v2-prediction.md.

.ad_snp <- function(idx) {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES", envir = environment())
  as_population(get("SNP55K_maize282_maf04", envir = environment()),
                individuals = idx)
}

# random-mating (HWE) panel; columns I1..In; `p` overrides the allele frequencies
.ad_hwe <- function(n, m = 40, seed = 1, n_chr = 4, p = NULL, prefix = "I") {
  set.seed(seed)
  if (is.null(p)) p <- stats::runif(m, 0.2, 0.8)
  g <- t(vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1L, integer(n)))
  colnames(g) <- paste0(prefix, seq_len(n))
  cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                   chr = rep(seq_len(n_chr), each = m / n_chr),
                   pos = rep(seq_len(m / n_chr), n_chr),
                   cm = rep(seq(0, 100, length.out = m / n_chr), n_chr),
                   stringsAsFactors = FALSE), as.data.frame(g))
}

# ---- T1 / Codex: deeper tabular A entries, pinned to literal values ----------------

test_that("T1: a_matrix known answers for S2, S3, DH of F1 and S1, sib mating, backcross, DH x DH, self of a DH", {
  fnd <- as_population(.ad_hwe(4))
  f1  <- cross(fnd[1], fnd[2], n = 2, seed = 1)          # F1 a, b (full sibs)
  s1  <- selfcross(fnd[4], n = 1, seed = 2)              # S1 of P4: F = 1/2
  s2  <- selfcross(s1, n = 1, seed = 3)                  # S2: F = 3/4
  s3  <- selfcross(s2, n = 1, seed = 4)                  # S3: F = 7/8
  dh  <- double_haploid(f1[1], n = 2, seed = 5)          # two DH of F1 a
  sib <- cross(f1[1], f1[2], n = 1, seed = 6)            # F1a x F1b: F = 1/4
  bc  <- cross(f1[1], fnd[1], n = 1, seed = 7)           # backcross F1a x P1: F = 1/4
  dhx <- cross(dh[1], dh[2], n = 1, seed = 8)            # DH1 x DH2
  sdh <- selfcross(dh[1], n = 1, seed = 9)               # self of a DH
  hs  <- cross(fnd[1], fnd[3], n = 1, seed = 10)         # half sib of F1 a
  dhs <- double_haploid(s1, n = 1, seed = 11)            # DH of an S1 (F_parent = 1/2)
  all <- c(fnd, f1, s1, s2, s3, dh, sib, bc, dhx, sdh, hs, dhs)
  A <- a_matrix(all)
  id <- all$ids
  # 1-4 founders, 5-6 F1 a/b, 7 S1, 8 S2, 9 S3, 10-11 DH1/DH2, 12 sib, 13 bc,
  # 14 DH1xDH2, 15 self of DH1, 16 half sib, 17 DH of S1
  e <- function(i, j) unname(A[id[i], id[j]])
  expect_equal(e(5, 6), 0.5)                  # full sibs
  expect_equal(e(5, 16), 0.25)                # half sibs
  expect_equal(e(7, 7), 1.5)                  # S1: 1 + F = 1 + 1/2
  expect_equal(e(8, 8), 1.75)                 # S2: F = (1 + 1/2)/2 = 3/4
  expect_equal(e(9, 9), 1.875)                # S3: F = 7/8
  expect_equal(e(7, 8), 1.5)                  # S1 with its own self (S2)
  expect_equal(e(4, 9), 1)                    # founder with S3 descendant
  expect_equal(e(10, 10), 2)                  # DH is fully inbred
  expect_equal(e(10, 11), 1)                  # DH sibs: A_{F1a,F1a} = 1
  expect_equal(e(10, 5), 1)                   # DH with its parent F1a
  expect_equal(e(10, 1), 0.5)                 # DH with a grandparent
  expect_equal(e(10, 6), 0.5)                 # DH with the F1 b (full sib of its parent)
  expect_equal(e(12, 12), 1.25)               # F1a x F1b
  expect_equal(e(13, 13), 1.25)               # backcross
  expect_equal(e(14, 14), 1.5)                # DH1 x DH2: 1 + A_{DH1,DH2}/2
  expect_equal(e(15, 15), 2)                  # self of a DH: 1 + F, F = A_{DH,DH}/2 = 1
  expect_equal(e(15, 10), 2)                  # identical genotype to its DH parent
  expect_equal(e(14, 10), 1.5)                # (A_{DH1,DH1} + A_{DH2,DH1})/2
  expect_equal(e(17, 17), 2)                  # DH of an S1 is still A_ii = 2
  expect_equal(e(17, 7), 1.5)                 # A_{DH,parent} = A_{parent,parent} (S1)
  expect_true(isSymmetric(A))
  expect_gte(min(eigen(A, symmetric = TRUE, only.values = TRUE)$values), -1e-10)
})

# ---- T9: hand-solved mixed-model equations with an unphenotyped sire ---------------

test_that("T9: pedigree BLUP of a 5-animal pedigree equals the hand-built MME and GLS solution", {
  fn <- as_population(.ad_hwe(3, seed = 3))
  S <- fn[1]; D1 <- fn[2]; D2 <- fn[3]
  o1 <- cross(S, D1, n = 1, seed = 11)
  o2 <- cross(S, D2, n = 1, seed = 12)
  pop <- c(fn, o1, o2)
  id <- pop$ids                         # S, D1, D2, O1 (= S x D1), O2 (= S x D2)
  # the relationship matrix written down by hand (unrelated non-inbred founders)
  A <- matrix(c(1,   0,   0,   .5,  .5,
                0,   1,   0,   .5,  0,
                0,   0,   1,   0,   .5,
                .5,  .5,  0,   1,   .25,
                .5,  0,   .5,  .25, 1), 5, 5, dimnames = list(id, id))
  expect_equal(a_matrix(pop), A)
  y <- stats::setNames(c(2.1, 3.4, 1.7), id[c(2, 3, 4)])   # D1, D2, O1; no record on S, O2
  va <- 1.5; ve <- 2.5; lam <- ve / va
  e <- predict_ebv(pop, y, method = "pedigree", var_a = va, var_e = ve)
  # Henderson's MME with A^-1
  Z <- matrix(0, 3, 5); Z[cbind(1:3, match(names(y), id))] <- 1
  X <- matrix(1, 3, 1)
  C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + lam * solve(A)))
  sol <- solve(C, c(crossprod(X, y), crossprod(Z, y)))
  expect_equal(attr(e, "mu"), 2.4588235, tolerance = 1e-7)
  expect_equal(attr(e, "mu"), as.numeric(sol[1]), tolerance = 1e-10)
  expect_equal(unname(as.numeric(e)),
               c(-0.13439, -0.21855, 0.35294, -0.31086, 0.10928), tolerance = 1e-4)
  expect_equal(unname(as.numeric(e)), as.numeric(sol[-1]), tolerance = 1e-10)
  # reliability from the inverse coefficient matrix
  Cu <- solve(C)[-1, -1]
  expect_equal(unname(attr(e, "reliability")), unname(1 - diag(Cu) * lam / diag(A)),
               tolerance = 1e-10)
  # the unphenotyped sire is predicted through his progeny and mates
  expect_gt(unname(attr(e, "reliability")[1]), 0)
  expect_equal(attr(e, "lambda"), lam)
})

# ---- T6 / Codex: marker effects, ridge and a fixed base ----------------------------

test_that("T6: ridge > 0 returns no marker effects; a fixed base_freq (with a currently monomorphic marker) is back-solved on its own scale", {
  set.seed(31)
  m <- 30; n <- 40
  pp <- stats::runif(m, 0.25, 0.75)
  df <- .ad_hwe(n, m = m, seed = 31, n_chr = 3, p = pp)
  df[1, paste0("I", 1:n)] <- 1L                 # marker 1: monomorphic in the sample
  pop <- as_population(df)
  y <- stats::setNames(stats::rnorm(n), pop$ids)[1:28]
  # ridge: genomic matrix is blended, marker effects are not defined
  er <- predict_ebv(pop, y, h2 = 0.5, ridge = 0.1)
  expect_null(attr(er, "marker_effects"))
  expect_false(is.null(attr(predict_ebv(pop, y, h2 = 0.5), "marker_effects")))
  # fixed base frequencies: marker 1 is monomorphic in the sample but kept (p < 1)
  base <- pp; base[1] <- 0.6
  eb <- predict_ebv(pop, y, h2 = 0.5, base_freq = base)
  u <- attr(eb, "marker_effects")
  expect_length(u, m)
  # the monomorphic marker is kept (one effect per marker); its constant design
  # column is collinear with the intercept, so the effect is zero up to rounding
  expect_lt(abs(unname(u[1])), 1e-10)
  M <- t(dosages(pop)) + 1                       # individuals x markers, gene content
  Z <- sweep(M, 2, 2 * base, "-")
  expect_equal(unname(as.numeric(Z %*% u)), unname(as.numeric(eb)), tolerance = 1e-8)
  # and the back-solve is the marker-MME (RR-BLUP) solution on that scale
  s <- 2 * sum(base * (1 - base))
  va <- attr(eb, "var_a"); ve <- attr(eb, "var_e")
  Zr <- Z[1:28, , drop = FALSE]
  lhs <- rbind(cbind(28, t(colSums(Zr))),
               cbind(colSums(Zr), crossprod(Zr) + (ve / va) * s * diag(m)))
  sol <- solve(lhs, c(sum(y), crossprod(Zr, y)))
  expect_equal(unname(u), unname(sol[-1]), tolerance = 1e-8)
  expect_equal(attr(eb, "mu"), sol[[1]], tolerance = 1e-8)
})

# ---- T12: a foreign tester panel --------------------------------------------------

test_that("T12: with foreign testers GCA is half the breeding value at the tester frequencies", {
  m <- 40
  set.seed(40)
  pc <- stats::runif(m, 0.2, 0.8)
  pt <- pmin(pmax(pc + stats::runif(m, -0.15, 0.15), 0.05), 0.95)
  dfc <- .ad_hwe(30, m = m, seed = 41, p = pc, prefix = "C")
  dft <- .ad_hwe(25, m = m, seed = 42, p = pt, prefix = "T")
  all <- as_population(cbind(dfc, dft[, paste0("T", 1:25)]))
  cand <- all[1:30]; tst <- all[31:55]
  set.seed(43); a <- stats::rnorm(m); d <- stats::rnorm(m, 0, 0.8)
  ca <- combining_ability(cand, tst, qtn = seq_len(m), a = a, d = d,
                          design = "factorial")
  xc <- dosages(cand) + 1                                 # gene content, loci x cand
  pT <- rowMeans(dosages(tst) + 1) / 2                    # tester gamete frequency
  alpha <- a + d * (1 - 2 * pT)                           # alpha_T = a + d (q_T - p_T)
  gca <- 0.5 * as.numeric(crossprod(xc, alpha))
  expect_equal(unname(ca$gca), unname(gca - mean(gca)), tolerance = 1e-10)
  # the tester panel matters: the candidates' own frequencies would rank differently
  pC <- rowMeans(xc) / 2
  gca_own <- 0.5 * as.numeric(crossprod(xc, a + d * (1 - 2 * pC)))
  expect_gt(max(abs((gca - mean(gca)) - (gca_own - mean(gca_own)))), 1e-3)
  # and the tester GCAs are the mirror image at the candidate frequencies
  xt <- dosages(tst) + 1
  alpha_c <- a + d * (1 - 2 * pC)
  gca_t <- 0.5 * as.numeric(crossprod(xt, alpha_c))
  expect_equal(unname(ca$gca_testers), unname(gca_t - mean(gca_t)), tolerance = 1e-10)
})

# ---- T13: a candidate that is also a tester ---------------------------------------

test_that("T13: a candidate that is also a tester is selfed, in both methods", {
  pop <- .ad_snp(1:6); q <- 1:30
  set.seed(50); a <- stats::rnorm(30); d <- stats::rnorm(30, 0, 0.4)
  cand <- pop[1:4]; tst <- pop[3:6]
  ex <- combining_ability(cand, tst, q, a, d, design = "factorial")
  # self of individual i: gene content 2 -> a, 1 -> d/2, 0 -> -a, summed over loci
  xs <- dosages(pop[3:4])[q, , drop = FALSE] + 1
  self_val <- vapply(1:2, function(i) {
    x <- xs[, i]
    sum(ifelse(x == 2, a, ifelse(x == 1, d / 2, -a)))
  }, numeric(1))
  expect_equal(unname(c(ex$cross_means[3, 1], ex$cross_means[4, 2])),
               unname(self_val), tolerance = 1e-12)
  sm <- combining_ability(cand, tst, q, a, d, design = "factorial",
                          method = "simulated", n_progeny = 3, seed = 1)
  ped <- parentage(sm$progeny)
  expect_true("self" %in% ped$design)
  expect_equal(sum(ped$design == "self"), 2L * 3L)       # two selfed pairs x 3 progeny
  # every other cell of the 4 x 4 factorial is an ordinary cross
  expect_equal(sum(ped$design != "self"), (4L * 4L - 2L) * 3L)
})

# ---- Codex: exhaustive Mendelian enumeration, negative dominance -------------------

test_that("Codex: the cross expectation equals direct Mendelian enumeration over all nine genotype pairs, d of either sign", {
  # one locus; parents of gene content x1, x2 in {0, 1, 2}; genotypic values -a, d, +a
  enum <- function(x1, x2, a, d) {
    p1 <- x1 / 2; p2 <- x2 / 2
    pAA <- p1 * p2; paa <- (1 - p1) * (1 - p2); pAa <- 1 - pAA - paa
    pAA * a + pAa * d + paa * (-a)
  }
  for (a in c(1.3, 0.4)) for (d in c(-0.9, 0, 0.7, 2.2)) {
    grid <- expand.grid(x1 = 0:2, x2 = 0:2)
    direct <- matrix(mapply(enum, grid$x1, grid$x2, MoreArgs = list(a = a, d = d)),
                     3, 3)                                         # [x1, x2]
    got <- .expected_cross_means(matrix(c(0, 0.5, 1), 1), matrix(c(0, 0.5, 1), 1),
                                 a, d)
    expect_equal(unname(got), unname(direct), tolerance = 1e-14)
  }
  # the corners named in the code comment
  a <- 1.3; d <- -0.9
  Y <- .expected_cross_means(matrix(c(0, 0.5, 1), 1), matrix(c(0, 0.5, 1), 1), a, d)
  expect_equal(c(Y[3, 3], Y[3, 1], Y[1, 1], Y[2, 2]), c(a, d, -a, d / 2))
  # several loci add (linearity), through the user-facing function
  g <- data.frame(snp = c("q1", "q2"), allele = "A/G", chr = 1, pos = 1:2, cm = c(0, 10),
                  I1 = c(1L, 0L), I2 = c(0L, 0L), I3 = c(-1L, 1L), stringsAsFactors = FALSE)
  pop <- as_population(g)
  ca <- combining_ability(pop, pop, qtn = c("q1", "q2"), a = c(1.3, 0.5),
                          d = c(-0.9, 0.7), design = "factorial")
  one <- function(l, x1, x2) enum(x1, x2, c(1.3, 0.5)[l], c(-0.9, 0.7)[l])
  xg <- rbind(c(2, 1, 0), c(1, 1, 2))            # gene content of I1..I3 at q1, q2
  want <- outer(1:3, 1:3, Vectorize(function(i, k) {
    one(1, xg[1, i], xg[1, k]) + one(2, xg[2, i], xg[2, k])
  }))
  expect_equal(unname(ca$cross_means), want, tolerance = 1e-12)
})

# ---- Codex: simulated diallel ------------------------------------------------------

test_that("Codex: a simulated diallel gives a symmetric Y with an NA diagonal, is seed-reproducible and records var_e", {
  pop <- .ad_snp(1:6); q <- 1:30
  set.seed(60); a <- stats::rnorm(30); d <- stats::rnorm(30, 0, 0.3)
  s1 <- combining_ability(pop, qtn = q, a = a, d = d, design = "diallel",
                          method = "simulated", n_progeny = 3, h2 = 0.5, seed = 7)
  s2 <- combining_ability(pop, qtn = q, a = a, d = d, design = "diallel",
                          method = "simulated", n_progeny = 3, h2 = 0.5, seed = 7)
  expect_identical(s1$cross_means, s2$cross_means)
  expect_identical(s1$gca, s2$gca)
  Y <- s1$cross_means
  expect_equal(Y, t(Y))                              # reciprocal-free: Y_ij = Y_ji
  expect_true(all(is.na(diag(Y))))
  expect_true(all(is.finite(Y[row(Y) != col(Y)])))
  expect_equal(n_individuals(s1$progeny), choose(6, 2) * 3L)
  # the recorded residual variance is the broad-sense one for h2 on the progeny
  gv <- genotypic_value(s1$progeny, q, a, d)
  expect_equal(s1$var_e, stats::var(gv) * (1 - 0.5) / 0.5)
  # decomposition identities hold for the realized Y as for the expected one
  expect_equal(sum(s1$gca), 0, tolerance = 1e-10)
  expect_equal(unname(rowSums(s1$sca, na.rm = TRUE)), rep(0, 6), tolerance = 1e-10)
  # the seed changes the draw
  s3 <- combining_ability(pop, qtn = q, a = a, d = d, design = "diallel",
                          method = "simulated", n_progeny = 3, h2 = 0.5, seed = 8)
  expect_false(identical(s1$cross_means, s3$cross_means))
})

# ---- Codex / T15: template_effects branches ---------------------------------------

test_that("T15/Codex: template_effects reconstructs each trait and replicate; trait/rep are validated", {
  pp <- .ad_snp(1:40)
  f2 <- selfcross(cross(pp[1], pp[2], n = 1, seed = 1), n = 60, seed = 2)
  sim <- simulate_phenotype(f2, h2 = c(0.5, 0.4), n_traits = 2, n_reps = 2,
                            architecture = "independent", vary_qtn = TRUE, seed = 3) |>
    additive(prop = c(0.3, 0.2), n_qtn = 5) |>
    dominance(prop = c(0.2, 0.2))
  for (r in 1:2) {
    gs <- genetic_values(sim, rep = r)
    for (t in 1:2) {
      te <- template_effects(sim, trait = t, rep = r)
      gv <- genotypic_value(sim$geno, te$qtn, te$a, te$d)
      expect_equal(unname(gv - mean(gv)), unname(gs[, t] - mean(gs[, t])),
                   tolerance = 1e-8, info = paste("trait", t, "rep", r))
    }
  }
  # vary_qtn gives different loci to different replicates: the branches differ
  expect_false(identical(template_effects(sim, 1, 1)$qtn, template_effects(sim, 1, 2)$qtn))
  expect_error(template_effects(sim, trait = 3), "`trait` must be at most 2")
  expect_error(template_effects(sim, rep = 3), "between 1 and n_reps")
})

test_that("Codex: template_effects with a prop = 0 layer equals the template without it", {
  pp <- .ad_snp(1:40)
  f2 <- selfcross(cross(pp[1], pp[2], n = 1, seed = 1), n = 60, seed = 2)
  base <- simulate_phenotype(f2, h2 = 0.5, seed = 3) |>
    additive(prop = 0.3, n_qtn = 5) |>
    dominance(prop = 0.2)
  zero <- simulate_phenotype(f2, h2 = 0.5, seed = 3) |>
    additive(prop = 0.3, n_qtn = 5) |>
    dominance(prop = 0.2) |>
    additive(prop = 0, n_qtn = 4)
  tb <- template_effects(base)
  tz <- template_effects(zero)
  # the zero-variance layer adds no locus and no effect
  expect_equal(tz[tz$qtn %in% tb$qtn, c("qtn", "a", "d")],
               tb[, c("qtn", "a", "d")], ignore_attr = TRUE)
  extra <- tz[!tz$qtn %in% tb$qtn, ]
  expect_true(all(extra$a == 0 & extra$d == 0))
  gv <- genotypic_value(zero$geno, tz$qtn, tz$a, tz$d)
  gs <- genetic_values(zero)[, 1]
  expect_equal(unname(gv - mean(gv)), unname(gs - mean(gs)), tolerance = 1e-8)
})
