# test-audit-prediction.R -- regression tests for the Sept-2026 audit of the
# prediction group (a_matrix, predict_ebv, prediction_accuracy, combining_ability,
# progeny_test, template_effects). Finding ids (PRED-*) refer to
# .tmp/audit-2026-09-29/reconciliation/v2-prediction.md.

.ap_snp <- function(idx) {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES", envir = environment())
  as_population(get("SNP55K_maize282_maf04", envir = environment()),
                individuals = idx)
}

.ap_q <- c("ss196442916", "ss196439337", "ss196480535")

# ---- PRED-F1: inbred founders (founder_f) ----------------------------------

test_that("PRED-F1: founder_f = 1 gives inbred founders A_ii = 2 and A[P, F1] = 1", {
  p2 <- .ap_snp(1:2)
  f1 <- cross(p2[1], p2[2], n = 2, seed = 1)
  s1 <- selfcross(f1[1], n = 1, seed = 2)
  dh <- double_haploid(f1[1], n = 1, seed = 3)
  all <- c(p2, f1, s1, dh)
  id <- all$ids
  # default (non-inbred base) is unchanged
  A0 <- a_matrix(all)
  expect_equal(A0[id[1], id[1]], 1)
  expect_equal(A0[id[1], id[3]], 0.5)
  expect_equal(A0[id[3], id[4]], 0.5)
  # inbred lines: F = 1
  A1 <- a_matrix(all, founder_f = 1)
  expect_equal(unname(diag(A1)[1:2]), c(2, 2))
  expect_equal(A1[id[1], id[2]], 0)                 # unrelated lines
  expect_equal(A1[id[1], id[3]], 1)                 # P - F1
  expect_equal(A1[id[3], id[3]], 1)                 # F1 of two inbred lines
  expect_equal(A1[id[3], id[4]], 1)                 # F1 sibs
  expect_equal(A1[id[5], id[5]], 1.5)               # S1 of the F1
  expect_equal(A1[id[6], id[6]], 2)                 # DH rule unchanged
  expect_equal(A1[id[1], id[6]], 1)                 # P - DH of the F1
  expect_true(isSymmetric(A1))
  # named vector: only the named founders are inbred
  An <- a_matrix(all, founder_f = c("4226" = 1))
  expect_equal(An[id[1], id[1]], 2)
  expect_equal(An[id[2], id[2]], 1)
  expect_equal(An[id[1], id[3]], 1)
  expect_equal(An[id[2], id[3]], 0.5)
})

test_that("PRED-F1: founder_f is validated", {
  p2 <- .ap_snp(1:2)
  expect_error(a_matrix(p2, founder_f = 2), "founder_f")
  expect_error(a_matrix(p2, founder_f = -0.1), "founder_f")
  expect_error(a_matrix(p2, founder_f = NA_real_), "founder_f")
  expect_error(a_matrix(p2, founder_f = c(1, 1)), "founder_f")           # unnamed
  expect_error(a_matrix(p2, founder_f = c(nope = 1)), "founder_f")       # unknown
  expect_error(a_matrix(p2, founder_f = c("4226" = 1, "4226" = 0)),
               "founder_f")                                              # duplicated
  expect_error(a_matrix(p2, founder_f = "1"), "founder_f")
})

test_that("PRED-F1/T2: a_matrix equals an independent coancestry recursion", {
  p2 <- .ap_snp(1:3)
  f1 <- cross(p2[1], p2[2], n = 3, seed = 1)
  f1b <- cross(p2[3], f1[1], n = 1, seed = 5)
  s1 <- selfcross(f1[1], n = 1, seed = 2)
  s2 <- selfcross(s1, n = 1, seed = 6)
  dh <- double_haploid(f1[2], n = 2, seed = 3)
  sibm <- cross(f1[1], f1[2], n = 1, seed = 7)
  all <- c(p2, f1, f1b, s1, s2, dh, sibm)
  ped <- as.data.frame(parentage(all, ancestors = TRUE))
  ped <- ped[order(ped$generation), ]
  rownames(ped) <- ped$key
  ff <- c("4226" = 1, "4722" = 0.5)   # founders: two inbred to different degrees
  memo <- new.env()
  phi <- function(a, b) {
    k <- paste(sort(c(a, b)), collapse = "|")
    if (!is.null(memo[[k]])) return(memo[[k]])
    ga <- ped[a, "generation"]; gb <- ped[b, "generation"]
    val <- if (a == b) {
      if (ped[a, "design"] == "dh") {
        1
      } else if (is.na(ped[a, "mother_key"])) {
        f <- ff[ped[a, "id"]]; (1 + if (is.na(f)) 0 else f) / 2
      } else {
        (1 + phi(ped[a, "mother_key"], ped[a, "father_key"])) / 2
      }
    } else if (ga >= gb && !is.na(ped[a, "mother_key"])) {
      (phi(ped[a, "mother_key"], b) + phi(ped[a, "father_key"], b)) / 2
    } else if (!is.na(ped[b, "mother_key"])) {
      (phi(a, ped[b, "mother_key"]) + phi(a, ped[b, "father_key"])) / 2
    } else {
      0
    }
    assign(k, val, envir = memo)
    val
  }
  keys <- .ensure_pedigree(all)$keys
  A2 <- outer(keys, keys, Vectorize(function(a, b) 2 * phi(a, b)))
  A <- a_matrix(all, founder_f = ff)
  expect_equal(unname(A), A2, tolerance = 1e-12)
})

test_that("PRED-F1: the pedigree and genomic relationships of inbred founders agree", {
  p2 <- .ap_snp(1:2)
  fs <- .ap_snp(1:2)
  expect_gt(unname(diag(g_matrix(fs))[1]), 1.5)              # genomic: about 2
  A <- a_matrix(p2, founder_f = 1)
  expect_equal(unname(diag(A)), c(2, 2))
})

# ---- PRED-01: supplied-K validation ----------------------------------------

test_that("PRED-01: a K whose correlation scale overflows is rejected", {
  pop <- .ap_snp(1:3)
  K <- matrix(0, 3, 3)
  diag(K) <- c(4.940656e-324, 1e308, 1)
  K[1, 2] <- K[2, 1] <- 1e308
  dimnames(K) <- list(pop$ids, pop$ids)
  expect_true(all(is.finite(K)))
  expect_error(.check_relationship(K), "finite|positive semidefinite")
  y <- stats::setNames(1, pop$ids[3])
  expect_error(predict_ebv(pop, y, K = K, var_a = 1, var_e = 1),
               "finite|positive semidefinite")
})

test_that("PRED-01: near-overflow entries never make the symmetrization overflow", {
  K <- matrix(c(1.7e308, 1.7e308, 1.7e308, 1.7e308), 2, 2)
  expect_error(.check_relationship(K), NA)        # valid (rank one), no overflow
  K2 <- K; K2[1, 2] <- K2[2, 1] <- -1.7e308
  expect_error(.check_relationship(K2), NA)
})

# ---- PRED-03 / PRED-07: prediction_accuracy --------------------------------

test_that("PRED-03: duplicated names in ebv or truth are an error", {
  expect_error(prediction_accuracy(c(a = 1, a = 100, b = 2, c = 3),
                                   c(a = 1, b = 2, c = 3)), "duplicated")
  expect_error(prediction_accuracy(c(a = 1, b = 2, c = 3),
                                   c(a = 1, a = 100, b = 2, c = 3)), "duplicated")
})

test_that("PRED-07/T20: a constant predictor is a documented NA with a package warning; partial overlap counts", {
  expect_warning(
    r <- prediction_accuracy(c(a = 1, b = 1, c = 1), c(a = 1, b = 2, c = 3)),
    "constant")
  expect_true(is.na(r[["accuracy"]]))
  expect_true(is.na(r[["slope"]]))
  expect_equal(r[["n"]], 3)
  expect_warning(
    r2 <- prediction_accuracy(c(a = 1, b = 2, c = 3), c(a = 5, b = 5, c = 5)),
    "constant")
  expect_true(is.na(r2[["accuracy"]]))
  # partially overlapping names: only the intersection is scored
  r3 <- prediction_accuracy(c(a = 1, b = 2, c = 3, d = 4, z = 9),
                            c(a = 1, b = 2.5, c = 3, d = 4.2, y = 0))
  expect_equal(r3[["n"]], 4)
})

# ---- PRED-06: irrelevant ref ------------------------------------------------

test_that("PRED-06: ref without h2 is an error in every entry point", {
  pop <- .ap_snp(1:6); q <- 1:30
  expect_error(
    combining_ability(pop[1:3], pop[4:6], q, rep(1, 30), method = "simulated",
                      n_progeny = 2, ref = "not a population", seed = 1),
    "ref")
  expect_error(
    combining_ability(pop[1:3], pop[4:6], q, rep(1, 30), method = "simulated",
                      n_progeny = 2, var_e = 1, ref = "x", seed = 1), "ref")
  expect_error(
    progeny_test(pop[1:2], pop[3:6], q, rep(1, 30), n_progeny = 2, ref = "x",
                 seed = 1), "ref")
  expect_error(
    progeny_test(pop[1:2], pop[3:6], q, rep(1, 30), n_progeny = 2, var_e = 1,
                 ref = "x", seed = 1), "ref")
  y <- stats::setNames(stats::rnorm(6), pop$ids)
  expect_error(predict_ebv(pop, y, var_a = 1, var_e = 1, ref = y), "ref")
  # a valid ref with h2 still works and is what sets var_a
  e <- predict_ebv(pop, y, h2 = 0.5, ref = c(1, 2, 3, 4), method = "pedigree")
  expect_equal(attr(e, "var_a"), 0.5 * stats::var(c(1, 2, 3, 4)))
  expect_error(predict_ebv(pop, y, h2 = 0.5, ref = "x", method = "pedigree"), "ref")
})

# ---- PRED-F11: duplicated individual phenotyped twice ------------------------

test_that("PRED-F11: one individual under two ids cannot be phenotyped twice", {
  pop <- .ap_snp(1:6)
  dup <- c(pop[1:3], pop[2])
  expect_equal(a_matrix(dup)[2, 4], 1)              # documented: a copy is itself
  y <- stats::setNames(c(1, 2, 3, 2.5), dup$ids)
  expect_error(predict_ebv(dup, y, method = "pedigree", var_a = 1, var_e = 1),
               "same individual")
  # one copy phenotyped is fine
  expect_error(predict_ebv(dup, y[1:3], method = "pedigree", var_a = 1,
                           var_e = 1), NA)
})

# ---- BLUP algebra adopted from the audit (T3-T8) -----------------------------

test_that("T4: variance-ratio limits", {
  pop <- .ap_snp(1:6)
  set.seed(4); y <- stats::setNames(stats::rnorm(6), pop$ids)
  K <- diag(6); dimnames(K) <- list(pop$ids, pop$ids)
  e <- predict_ebv(pop, y, K = K, var_a = 1, var_e = 1e-8)
  expect_equal(as.numeric(e), as.numeric(y - attr(e, "mu")), tolerance = 1e-6)
  e2 <- predict_ebv(pop, y, K = K, var_a = 1, var_e = 1e8)
  expect_lt(max(abs(e2)), 1e-6)
})

test_that("T7/T8: ref sets var_a; a single record; NA phenotypes are refused", {
  pop <- .ap_snp(1:6)
  set.seed(5); y <- stats::setNames(stats::rnorm(6), pop$ids)
  K <- diag(6); dimnames(K) <- list(pop$ids, pop$ids)
  ref <- stats::rnorm(20)
  e <- predict_ebv(pop, y, K = K, h2 = 0.4, ref = ref)
  expect_equal(attr(e, "var_a"), 0.4 * stats::var(ref))
  e1 <- predict_ebv(pop, y[1], K = K, var_a = 1, var_e = 1)
  expect_equal(as.numeric(e1), rep(0, 6))
  expect_equal(as.numeric(attr(e1, "reliability")), rep(0, 6))
  yn <- y; yn[2] <- NA
  expect_error(predict_ebv(pop, yn, K = K, var_a = 1, var_e = 1), "finite")
})

test_that("T5: GBLUP with a singular G equals marker-effect RR-BLUP", {
  set.seed(6)
  m <- 12; n <- 40
  pp <- stats::runif(m, 0.25, 0.75)
  g <- t(vapply(pp, function(x) stats::rbinom(n, 2, x) - 1L, integer(n)))
  colnames(g) <- paste0("I", seq_len(n))
  df <- cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                         chr = rep(1:2, each = m / 2), pos = rep(seq_len(m / 2), 2),
                         cm = rep(seq(0, 50, length.out = m / 2), 2),
                         stringsAsFactors = FALSE), as.data.frame(g))
  pop <- as_population(df)
  y <- stats::setNames(stats::rnorm(n), pop$ids)
  va <- 1.5; ve <- 2.5
  e <- predict_ebv(pop, y[1:30], var_a = va, var_e = ve)
  M <- t(dosages(pop)) + 1
  p <- colMeans(M) / 2
  Z <- sweep(M, 2, 2 * p, "-")
  s <- 2 * sum(p * (1 - p))
  lam <- ve / va
  Zr <- Z[1:30, , drop = FALSE]; yy <- y[1:30]
  # marker MME, u_m ~ N(0, I va / s), unpenalized mean
  lhs <- rbind(cbind(30, t(colSums(Zr))),
               cbind(colSums(Zr), crossprod(Zr) + lam * s * diag(m)))
  sol <- solve(lhs, c(sum(yy), crossprod(Zr, yy)))
  expect_equal(attr(e, "mu"), sol[[1]], tolerance = 1e-8)
  expect_equal(unname(as.numeric(e)), as.numeric(Z %*% sol[-1]), tolerance = 1e-8)
  expect_equal(unname(attr(e, "marker_effects")), unname(sol[-1]), tolerance = 1e-8)
})

# ---- PRED-02 / PRED-05: theory sentences pinned ------------------------------

test_that("PRED-02: dominance enters a half-sib mean through alpha = a + d(1 - 2p)", {
  m <- function(p) {
    Y <- .expected_cross_means(matrix(c(0, 0.5, 1), 1), matrix(p, 1), a = 0, d = 1)
    Y[3, 1] - Y[1, 1]                                # candidate AA minus aa
  }
  expect_equal(m(0.2), 0.6)                          # d (1 - 2p)
  expect_equal(m(0.5), 0)
})

test_that("PRED-05: simulated SCA has Mendelian sampling even with d = 0; expected is exactly 0", {
  p2 <- .ap_snp(1:6)
  f2 <- selfcross(cross(p2[1], p2[2], n = 1, seed = 1), n = 6, seed = 2)
  q <- 1:30
  ex <- combining_ability(f2[1:3], f2[4:6], q, rep(1, 30), d = 0,
                          design = "factorial")
  expect_lt(max(abs(ex$sca)), 1e-10)
  sm <- combining_ability(f2[1:3], f2[4:6], q, rep(1, 30), d = 0,
                          design = "factorial", method = "simulated",
                          n_progeny = 1, seed = 1)
  expect_gt(max(abs(sm$sca)), 1e-3)
})

# ---- PRED-F4: seed restores the ambient stream -------------------------------

test_that("PRED-F4: seed= in the simulators restores the caller's RNG state", {
  pop <- .ap_snp(1:6); q <- 1:30
  set.seed(99); before <- .Random.seed
  a1 <- combining_ability(pop[1:3], pop[4:6], q, rep(1, 30),
                          method = "simulated", n_progeny = 2, seed = 5)
  expect_identical(.Random.seed, before)
  set.seed(99); before <- .Random.seed
  pt <- progeny_test(pop[1:2], pop[3:6], q, rep(1, 30), n_progeny = 2, seed = 5)
  expect_identical(.Random.seed, before)
  # and with a residual
  set.seed(99); before <- .Random.seed
  a2 <- combining_ability(pop[1:3], pop[4:6], q, rep(1, 30),
                          method = "simulated", n_progeny = 2, h2 = 0.5, seed = 5)
  expect_identical(.Random.seed, before)
  # results are still reproducible
  a3 <- combining_ability(pop[1:3], pop[4:6], q, rep(1, 30),
                          method = "simulated", n_progeny = 2, h2 = 0.5, seed = 5)
  expect_identical(a2$cross_means, a3$cross_means)
  pt2 <- progeny_test(pop[1:2], pop[3:6], q, rep(1, 30), n_progeny = 2, seed = 5)
  expect_identical(pt$progeny_mean, pt2$progeny_mean)
})

test_that("T7 (Codex): with seed = NULL two runs from the same ambient seed agree", {
  pop <- .ap_snp(1:8); q <- 1:30
  run <- function() {
    set.seed(11)
    progeny_test(pop[1:3], pop[4:8], q, rep(1, 30), n_progeny = 2, h2 = 0.4)
  }
  r1 <- run(); r2 <- run()
  expect_identical(attr(r1, "records"), attr(r2, "records"))
  expect_identical(dosages(attr(r1, "progeny")), dosages(attr(r2, "progeny")))
  # mates -> meioses -> residuals: the residual draws do not change the genotypes
  set.seed(11)
  r3 <- progeny_test(pop[1:3], pop[4:8], q, rep(1, 30), n_progeny = 2)
  expect_identical(dosages(attr(r1, "progeny")), dosages(attr(r3, "progeny")))
})

# ---- PRED-F5: template_effects and transcriptome layers ----------------------

test_that("PRED-F5: template_effects refuses a derived transcriptome layer", {
  pp <- .ap_snp(1:60)
  sim <- simulate_phenotype(pp, h2 = 0.6, seed = 3, transcriptome = TRUE) |>
    additive(prop = 0.3, n_qtn = 5) |>
    transcriptome(prop = 0.2, n_genes = 4)
  expect_error(template_effects(sim), "transcriptome")
  # a prop = 0 transcriptome layer contributes no genetic value: allowed
  sim0 <- simulate_phenotype(pp, h2 = 0.6, seed = 3, transcriptome = TRUE) |>
    additive(prop = 0.6, n_qtn = 5) |>
    transcriptome(prop = 0, n_genes = 4)
  te <- template_effects(sim0)
  expect_equal(nrow(te), 5)
})

test_that("T15/T8 (Codex): template_effects reconstructs overlapping A/D layers", {
  pp <- .ap_snp(1:40)
  f2 <- selfcross(cross(pp[1], pp[2], n = 1, seed = 1), n = 60, seed = 2)
  sim <- simulate_phenotype(f2, h2 = 0.6, seed = 3) |>
    additive(prop = 0.2, n_qtn = 6) |>
    additive(prop = 0.2, n_qtn = 6) |>
    dominance(prop = 0.2)
  te <- template_effects(sim)
  gv <- genotypic_value(sim$geno, te$qtn, te$a, te$d)
  gs <- genetic_values(sim)[, 1]
  expect_equal(unname(gv - mean(gv)), unname(gs - mean(gs)), tolerance = 1e-8)
  # orthogonal additive layer carries its own dominance deviation
  sim2 <- simulate_phenotype(f2, h2 = 0.5, seed = 4) |>
    additive(prop = 0.5, n_qtn = 5, orthogonal = TRUE, a = 0.5, d = 0.5)
  te2 <- template_effects(sim2)
  gv2 <- genotypic_value(sim2$geno, te2$qtn, te2$a, te2$d)
  gs2 <- genetic_values(sim2)[, 1]
  expect_equal(unname(gv2 - mean(gv2)), unname(gs2 - mean(gs2)), tolerance = 1e-8)
})

# ---- Combining-ability exact forms (T10-T14, T19) ----------------------------

test_that("T10: Griffing method-4 closed forms", {
  pop <- .ap_snp(1:8); q <- 1:30
  set.seed(8); a <- stats::rnorm(30); d <- stats::rnorm(30, 0, 0.3)
  ca <- combining_ability(pop[1:6], qtn = q, a = a, d = d, design = "diallel")
  Y <- ca$cross_means; p <- nrow(Y)
  Ym <- Y; diag(Ym) <- 0
  Yi <- rowSums(Ym); Yd <- sum(Ym[upper.tri(Ym)])      # Y.. over unordered pairs
  expect_equal(unname(ca$gca), unname((p * Yi - 2 * Yd) / (p * (p - 2))),
               tolerance = 1e-10)
  S <- Ym - outer(Yi, Yi, "+") / (p - 2) + 2 * Yd / ((p - 1) * (p - 2))
  off <- row(S) != col(S)
  expect_equal(unname(ca$sca[off]), unname(S[off]), tolerance = 1e-10)
})

test_that("T11: single-locus 3 x 3 table with selfs", {
  g <- data.frame(snp = "q", allele = "A/G", chr = 1, pos = 1, cm = 0,
                  I1 = 1L, I2 = 0L, I3 = -1L)
  pop <- as_population(g)
  a <- 2; d <- 0.6
  ca <- combining_ability(pop, pop, qtn = "q", a = a, d = d,
                          design = "factorial")
  Y <- ca$cross_means
  expect_equal(unname(Y),
               unname(rbind(c(a, (a + d) / 2, d),
                            c((a + d) / 2, d / 2, (d - a) / 2),
                            c(d, (d - a) / 2, -a))))
})

test_that("T14: simulated var_e follows h2 and passes through", {
  pop <- .ap_snp(1:6); q <- 1:30
  s1 <- combining_ability(pop[1:3], pop[4:6], q, rep(1, 30), method = "simulated",
                          n_progeny = 4, h2 = 0.5, seed = 1)
  gv <- genotypic_value(s1$progeny, q, rep(1, 30), rep(0, 30))
  expect_equal(s1$var_e, stats::var(gv) * (1 - 0.5) / 0.5)
  s2 <- combining_ability(pop[1:3], pop[4:6], q, rep(1, 30), method = "simulated",
                          n_progeny = 4, var_e = 0.7, seed = 1)
  expect_equal(s2$var_e, 0.7)
})

test_that("T19: print method", {
  pop <- .ap_snp(1:6); q <- 1:30
  ca <- combining_ability(pop[1:3], pop[4:6], q, rep(1, 30))
  expect_output(print(ca), "<combining_ability>")
})

# ---- progeny_test extras (T16-T18) -------------------------------------------

test_that("T16-T18: records, vector d, var_e, n_progeny = 1, family means", {
  pop <- .ap_snp(1:10); q <- 1:30
  set.seed(9); a <- stats::rnorm(30); d <- stats::rnorm(30, 0, 0.3)
  pt <- progeny_test(pop[1:4], pop[5:10], q, a, d, n_progeny = 1, var_e = 0.5,
                     seed = 2)
  expect_equal(attr(pt, "var_e"), 0.5)
  expect_equal(pt$n, rep(1, 4))
  # no residual: records are the genotypic values of the progeny
  pt0 <- progeny_test(pop[1:4], pop[5:10], q, a, d, n_progeny = 3, seed = 2)
  prog <- attr(pt0, "progeny")
  expect_equal(unname(attr(pt0, "records")),
               unname(as.numeric(genotypic_value(prog, q, a, d))))
  fam <- families(prog, "maternal_half_sib")
  expect_equal(pt0$progeny_mean,
               as.numeric(tapply(attr(pt0, "records"), fam, mean)),
               tolerance = 1e-12)
  # the residual does not change the genotypes (mates and meioses come first)
  pth <- progeny_test(pop[1:4], pop[5:10], q, a, d, n_progeny = 3, h2 = 0.4,
                      seed = 2)
  expect_identical(dosages(attr(pth, "progeny")), dosages(prog))
})

# ---- manifest snapshot (Codex proposal 9) -------------------------------------

test_that("selection_methods(): ids, functions and sources are pinned", {
  m <- selection_methods()
  expect_equal(m$id, c("mass", "within_family", "among_family", "combined",
                       "index", "quadratic_index", "culling", "tandem", "random",
                       "ocs", "usefulness", "mabc", "mas", "combining_ability",
                       "progeny_test", "blup", "mt_blup"))
  allowed <- c("Falconer & Mackay 1996", "Lush 1947", "Smith 1936; Hazel 1943",
               "Ceron-Rojas et al. 2026", "Hazel & Lush 1942", "package",
               "Meuwissen 1997", "Zhong & Jannink 2007; Lehermeier et al. 2017",
               "Frisch & Melchinger 2001, 2005", "Lande & Thompson 1990",
               "Sprague & Tatum 1942; Griffing 1956", "package derivation",
               "Henderson 1975; VanRaden 2008", "Henderson & Quaas 1976")
  expect_true(all(m$source %in% allowed))
})
