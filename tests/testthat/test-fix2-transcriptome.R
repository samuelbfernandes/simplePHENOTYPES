# Round-2 regression tests for the transcriptome group (A1, B9, C11, C17).

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04

.m3 <- function() {
  matrix(c(-1, 0, 1, 1, 0, -1, -1, 0, 1, 1, 1, -1), 3,
         dimnames = list(paste0("i", 1:3), paste0("m", 1:4)))
}

# ---- A1: h2_realized is the realized Var(G)/Var(P) ---------------------------

test_that("A1: h2_realized = realized Var(G)/Var(E), including 2Cov(G, R)", {
  tx <- suppressWarnings(simulate_transcriptome(.m3(), n_genes = 5, h2 = 0.5,
                                                seed = 1))
  vg <- apply(tx$genetic_expression, 1L, stats::var)
  ve <- apply(tx$expression, 1L, stats::var)
  vr <- apply(tx$expression - tx$genetic_expression - rowMeans(tx$expression -
                tx$genetic_expression) + rowMeans(tx$expression -
                tx$genetic_expression), 1L, stats::var)
  # realized heritability from the realized genetic values and phenotype
  expect_equal(unname(vg / ve), tx$genes$h2_realized, tolerance = 1e-10)
  # h2_var_ratio is retained and is the same quantity
  expect_equal(tx$genes$h2_var_ratio, tx$genes$h2_realized, tolerance = 1e-12)
  # the covariance term is what separates it from the allocation quantity
  gr <- tx$var_budget$gr_cov
  expect_equal(tx$genes$h2_realized,
               unname(vg / (vg + vr + gr)), tolerance = 1e-8)
  # Codex counterexample: n = 3 realized ratio ranges far from the target 0.5
  expect_gt(max(tx$genes$h2_realized), 1)
  expect_gt(diff(range(tx$genes$h2_realized)), 0.1)
})

test_that("A1: h2_allocated is the bounded Var(G)/(Var(G)+Var(R)) allocation", {
  tx <- suppressWarnings(simulate_transcriptome(.m3(), n_genes = 5, h2 = 0.5,
                                                seed = 1))
  expect_true("h2_allocated" %in% names(tx$genes))
  expect_true(all(tx$genes$h2_allocated >= 0 & tx$genes$h2_allocated <= 1))
  # forced to the target on the reference panel
  expect_equal(tx$genes$h2_allocated, tx$genes$h2_target, tolerance = 1e-8)
  # ... so it is NOT the realized heritability
  expect_false(isTRUE(all.equal(tx$genes$h2_allocated, tx$genes$h2_realized)))
  # a larger, well-posed panel: realized ~ target, allocated == target exactly
  Gn <- G[, c(1:5, 5L + seq_len(280))]
  txn <- suppressWarnings(simulate_transcriptome(Gn, n_genes = 100, h2 = 0.6,
                                                 seed = 3))
  expect_equal(txn$genes$h2_allocated, txn$genes$h2_target, tolerance = 1e-8)
  vg <- apply(txn$genetic_expression, 1L, stats::var)
  ve <- apply(txn$expression, 1L, stats::var)
  expect_equal(unname(vg / ve), txn$genes$h2_realized, tolerance = 1e-10)
  expect_lt(abs(mean(txn$genes$h2_realized) - 0.6), 0.05)
})

test_that("A1: predict() reports the same three quantities on new individuals", {
  Gn <- G[, c(1:5, 5L + seq_len(80))]
  tx <- simulate_transcriptome(Gn[, 1:65], n_genes = 20, h2 = 0.5, seed = 4)
  pr <- predict(tx, Gn[, c(1:5, 66:85)], seed = 5)
  vg <- apply(pr$genetic_expression, 1L, stats::var)
  ve <- apply(pr$expression, 1L, stats::var)
  expect_equal(unname(vg / ve), pr$genes$h2_realized, tolerance = 1e-10)
  expect_equal(pr$genes$h2_var_ratio, pr$genes$h2_realized, tolerance = 1e-12)
  expect_true(all(pr$genes$h2_allocated >= 0 & pr$genes$h2_allocated <= 1))
  # the allocation quantity is genuinely different from the realized ratio here
  expect_false(isTRUE(all.equal(pr$genes$h2_allocated, pr$genes$h2_realized)))
  # noiseless prediction: all expression variance is genetic
  p0 <- predict(tx, Gn[, c(1:5, 66:85)], residual = FALSE)
  expect_equal(p0$genes$h2_realized[p0$genes$h2_target > 0],
               rep(1, sum(p0$genes$h2_target > 0)), tolerance = 1e-10)
})

test_that("A1: print reports the realized heritability Var(G)/Var(P)", {
  tx <- suppressWarnings(simulate_transcriptome(.m3(), n_genes = 5, h2 = 0.5,
                                                seed = 1))
  out <- paste(capture.output(print(tx)), collapse = "\n")
  expect_match(out, "Var(G)/Var(P)", fixed = TRUE)
  expect_false(grepl("Var(G)/(Var(G)+Var(R))", out, fixed = TRUE))
  expect_match(out, sprintf("%.2f", max(tx$genes$h2_realized)), fixed = TRUE)
})

test_that("A1: non-genetic genes have zero realized and allocated heritability", {
  tx <- simulate_transcriptome(NULL, n_ind = 50, n_genes = 10, seed = 1)
  expect_true(all(tx$genes$h2_realized == 0))
  expect_true(all(tx$genes$h2_allocated == 0))
})

# ---- B9: identifiability after projecting out the intercept -------------------

test_that("B9: K = a(I - 11'/n) + b 11' is flagged unidentifiable", {
  set.seed(1)
  Y <- matrix(rnorm(3 * 20), 3, 20)
  for (K in list(diag(20) - matrix(1 / 20, 20, 20),
                 3 * (diag(20) - matrix(1 / 20, 20, 20)) + 0.7,
                 diag(20) + 0.5)) {                 # I + b 11' : same family
    r <- simplePHENOTYPES:::.greml_h2(Y, K)
    expect_equal(as.numeric(r), rep(0, 3))
    expect_false(attr(r, "identifiable"))
  }
})

test_that("B9: Codex's 4-individual, 6-balanced-marker GRM is flagged", {
  Z <- vapply(utils::combn(4, 2, simplify = FALSE), function(ix) {
    z <- rep(-1, 4); z[ix] <- 1; z
  }, numeric(4))
  expect_equal(dim(Z), c(4L, 6L))
  K <- simplePHENOTYPES:::.tx_grm(sweep(Z, 2L, colMeans(Z), "-"))
  expect_equal(K, (4 / 3) * (diag(4) - matrix(1 / 4, 4, 4)),
               tolerance = 1e-12)
  set.seed(1)
  r <- simplePHENOTYPES:::.greml_h2(matrix(rnorm(5 * 4), 5, 4), K)
  expect_equal(as.numeric(r), rep(0, 5))
  expect_false(attr(r, "identifiable"))
  # public mimic path: no spurious ~1 estimates, and the user is told
  rownames(Z) <- paste0("i", 1:4); colnames(Z) <- paste0("m", 1:6)
  set.seed(2)
  E <- matrix(rnorm(5 * 4), 5, 4, dimnames = list(NULL, rownames(Z)))
  w <- character(0)
  mm <- withCallingHandlers(
    simulate_transcriptome(Z, mimic = E, seed = 1),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_true(all(mm$calibration$h2$h2_greml == 0))
  expect_true(any(grepl("identif", w)))
})

test_that("B9: identifiable K (real GRM, or spread after projection) still works", {
  set.seed(3)
  n <- 60
  Z <- matrix(rnorm(n * 200), n, 200)
  K <- simplePHENOTYPES:::.tx_grm(Z)
  Y <- matrix(rnorm(4 * n), 4, n)
  r <- simplePHENOTYPES:::.greml_h2(Y, K)
  expect_null(attr(r, "identifiable"))
  expect_true(all(r >= 0 & r <= 1))
  # a non-scalar K after projection but with a large b 11' part stays identifiable
  K2 <- K + 5
  r2 <- simplePHENOTYPES:::.greml_h2(Y, K2)
  expect_null(attr(r2, "identifiable"))
  expect_equal(as.numeric(r2), as.numeric(r), tolerance = 1e-6)
})

# ---- C11: marginal epistasis share --------------------------------------------

test_that("C11: marginal epistasis share = eps/(eps + (1-eps)/s_ct^2)", {
  set.seed(3)
  z <- sample(c(-1, 0, 1), 40, TRUE, prob = c(.3, .4, .3))
  M <- cbind(m1 = z, m2 = z)                      # perfectly correlated cis/hub
  rownames(M) <- paste0("i", 1:40)
  ann <- data.frame(gene_id = "g1", chr = 1, tss = 1)
  eps <- 0.4
  hit <- 0L; big <- 0L
  for (s in 1:8) {
    tx <- suppressWarnings(simulate_transcriptome(
      M, annotation = ann, cis_window = 0.5, h2 = 0.5, cis_fraction = 0.5,
      epistasis = eps, seed = s))
    b <- tx$var_budget
    marg <- b$v_epi / (b$v_cis + b$v_trans + b$v_epi)
    # s_ct^2 = Var(cis + trans parts) / (v_cis + v_trans) from the realized budget
    s2 <- 1 + b$cis_trans_cov / (b$v_cis + b$v_trans)
    expect_equal(marg, eps / (eps + (1 - eps) / s2), tolerance = 1e-10)
    if (tx$genes$n_cis > 0L) {
      hit <- hit + 1L
      if (s2 > 1.5) {
        big <- big + 1L
        expect_equal(marg, 4 / 7, tolerance = 1e-7)   # 0.5714286, not 0.4
        expect_gt(abs(marg - eps), 0.1)               # not "slightly"
      }
    }
  }
  expect_gt(hit, 0L); expect_gt(big, 0L)
})

test_that("C11: the roxygen text gives the formula and no longer says 'slightly'", {
  src <- readLines(testthat::test_path("..", "..", "R", "transcriptome_simulate.R"))
  expect_false(any(grepl("differs slightly otherwise", src, fixed = TRUE)))
  expect_true(any(grepl("epsilon + (1 - epsilon) / s_ct^2", src, fixed = TRUE)))
})

# ---- C17: benchmark 05 wording -------------------------------------------------

test_that("C17: benchmark 05 defines 'structural null' as no shared module/status", {
  b5 <- paste(readLines(testthat::test_path("..", "..", "benchmarks",
                                            "05_twas_power.R")), collapse = "\n")
  rd <- paste(readLines(testthat::test_path("..", "..", "benchmarks",
                                            "README.md")), collapse = "\n")
  for (txt in list(b5, rd)) {
    expect_true(grepl("no shared", txt, ignore.case = TRUE))
    expect_true(grepl("linkage disequilibrium|LD", txt))
  }
})
