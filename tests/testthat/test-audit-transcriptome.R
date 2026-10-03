# Regression tests for the 2026-09 transcriptome audit (TX-F1..F14, O1..O14).

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04

# a tiny 3 x 4 panel (audit TX-F1 falsifier)
.m3 <- function() {
  matrix(c(-1, 0, 1, 1, 0, -1, -1, 0, 1, 1, 1, -1), 3,
         dimnames = list(paste0("i", 1:3), paste0("m", 1:4)))
}

test_that("TX-F1: h2_allocated is bounded in [0, 1] for n = 3, 5, 50, 280", {
  tx3 <- suppressWarnings(simulate_transcriptome(.m3(), n_genes = 5, h2 = 0.5,
                                                 seed = 1))
  expect_true(all(tx3$genes$h2_allocated >= 0 & tx3$genes$h2_allocated <= 1))
  # the realized heritability keeps its (unbounded) finite-sample meaning
  expect_gt(max(tx3$genes$h2_var_ratio), 1)
  expect_equal(tx3$genes$h2_var_ratio,
               0.5 / (1 + tx3$var_budget$gr_cov), tolerance = 1e-6)
  for (n in c(5L, 50L, 280L)) {
    Gn <- G[, c(1:5, 5L + seq_len(n))]
    tx <- suppressWarnings(simulate_transcriptome(Gn, n_genes = 300, h2 = 0.9,
                                                  seed = 2))
    expect_true(all(tx$genes$h2_allocated >= 0 & tx$genes$h2_allocated <= 1),
                info = paste("n =", n))
    vg <- apply(tx$genetic_expression, 1L, stats::var)
    ve <- apply(tx$expression, 1L, stats::var)
    expect_equal(unname(vg / ve), tx$genes$h2_var_ratio, tolerance = 1e-10)
    expect_equal(unname(vg / ve), tx$genes$h2_realized, tolerance = 1e-10)
  }
})

test_that("TX-F1: the documented minimum n (30) is enforced by a warning", {
  expect_warning(
    simulate_transcriptome(G[, c(1:5, 6:15)], n_genes = 10, h2 = 0.5, seed = 1),
    "individuals")
  expect_warning(
    simulate_transcriptome(G[, c(1:5, 6:35)], n_genes = 10, h2 = 0.5, seed = 1),
    NA)
  # fewer than 3 stays an error; a non-genetic request never warns
  expect_error(simulate_transcriptome(.m3()[1:2, ], n_genes = 2), "at least 3")
  expect_warning(
    simulate_transcriptome(G[, c(1:5, 6:15)], n_genes = 10, h2 = 0, seed = 1), NA)
})

test_that("TX-F2: cis_fraction_realized is the realized share of Var(G), not the target", {
  tx <- simulate_transcriptome(G, n_genes = 300, h2 = 0.7, cis_fraction = 0.5,
                               seed = 3)
  vg <- apply(tx$genetic_expression, 1L, stats::var)
  expect_equal(tx$genes$cis_fraction_realized,
               unname(tx$var_budget$v_cis / vg), tolerance = 1e-9)
  both <- tx$genes$n_cis > 0 & tx$genes$trans_scale > 0
  expect_gt(sum(both), 50)
  expect_false(isTRUE(all.equal(tx$genes$cis_fraction_realized[both],
                                tx$genes$cis_fraction_target[both],
                                tolerance = 1e-6)))
  # scatters around the target but is not far from it
  expect_lt(mean(abs(tx$genes$cis_fraction_realized[both] -
                       tx$genes$cis_fraction_target[both])), 0.1)
})

test_that("TX-F9/O3: epistasis_realized exists and the marginal identity is exact", {
  tx <- simulate_transcriptome(G, n_genes = 150, epistasis = 0.4, seed = 1)
  vg <- apply(tx$genetic_expression, 1L, stats::var)
  gg <- vg > 1e-8
  expect_equal(unname(tx$genes$epistasis_realized[gg]),
               unname((tx$var_budget$v_epi / vg)[gg]), tolerance = 1e-9)
  # marginal share ~ epsilon (exactly so when the cis and trans scores are
  # uncorrelated; the cis/trans correlation shifts it slightly)
  b <- tx$var_budget
  marg <- (b$v_epi / (b$v_cis + b$v_trans + b$v_epi))[gg]
  expect_equal(mean(marg), 0.4, tolerance = 0.03)
  expect_lt(max(abs(marg - 0.4)), 0.1)
  # the realized column is the share of Var(G) (covariance rows included)
  expect_equal(mean(tx$genes$epistasis_realized[gg]), 0.4, tolerance = 0.06)
  expect_true(all(tx$genes$epistasis_target == 0.4))
})

test_that("O2: epistasis with fewer than two eligible markers warns", {
  m1 <- matrix(rep(c(-1, 0, 1), 30), ncol = 1,
               dimnames = list(paste0("i", 1:90), "m1"))
  expect_warning(
    tx <- simulate_transcriptome(m1, n_genes = 4, h2 = 0.6, cis_fraction = 1,
                                 epistasis = 1, seed = 17),
    "at least two eligible")
  expect_true(all(tx$genes$n_epi == 0))
  expect_true(all(tx$genes$epistasis_realized == 0))
})

test_that("TX-F3/O1: mimic recovers kappa and does not manufacture structure", {
  base <- simulate_transcriptome(G, n_genes = 300, h2 = 0.4,
                                 residual_module_fraction = 0.3, n_factors = 5,
                                 seed = 1)
  mm <- simulate_transcriptome(G, mimic = base$expression, seed = 2)
  expect_equal(mm$reference$n_factors, 5L)
  expect_equal(mm$reference$kappa, 0.3, tolerance = 0.1)
  # SPEC test 8 (implemented honestly): moments exact, and the STRENGTH of the
  # leading co-expression spectrum (sum of the top-Q eigenvalues of the gene-gene
  # correlation) matches within +/- 10%. Individual eigenvalues are not compared:
  # module sizes are redrawn, so they legitimately differ.
  expect_equal(unname(rowMeans(mm$expression)), unname(rowMeans(base$expression)),
               tolerance = 1e-8)
  top5 <- function(x) sum(eigen(stats::cor(t(x)), symmetric = TRUE,
                                only.values = TRUE)$values[1:5])
  expect_lt(abs(top5(mm$expression) / top5(base$expression) - 1), 0.10)
  # pure noise in, no structure out
  set.seed(5)
  noise <- matrix(rnorm(60 * 280), 60, 280, dimnames = list(NULL, names(G)[-(1:5)]))
  mn <- simulate_transcriptome(G, mimic = noise, seed = 5)
  expect_lt(mn$reference$kappa, 0.01)
  cc <- stats::cor(t(mn$expression))
  same <- outer(mn$genes$module, mn$genes$module, "==") & upper.tri(cc)
  expect_lt(abs(mean(cc[same])), 0.02)
  expect_equal(stats::median(mn$calibration$h2$h2_greml), 0)
})

test_that("TX-F3/O1: mimic does not retain loadings/signs (documented)", {
  set.seed(101)
  ids <- names(G)[-(1:5)]
  f <- rnorm(280)
  sgn <- rep(c(1, -1), each = 20)
  E <- t(vapply(sgn, function(s) s * f + rnorm(280, sd = 0.05), numeric(280)))
  dimnames(E) <- list(NULL, ids)
  mm <- simulate_transcriptome(G, mimic = E, n_factors = 1, seed = 102)
  expect_true(all(mm$loadings$loading == 1))      # sign information is not kept
  expect_equal(mm$calibration$n_factors, 1L)
  expect_true(all(c("source", "n_factors", "kappa", "h2") %in%
                    names(mm$calibration)))       # only these summarize the input
})

test_that("TX-F4: mimic warns at small n and about an ignored `h2`", {
  small <- G[, c(1:5, 6:65)]                        # n = 60
  E <- simulate_transcriptome(small, n_genes = 20, seed = 1)$expression
  expect_warning(simulate_transcriptome(small, mimic = E, seed = 2),
                 "GREML")
  Eb <- simulate_transcriptome(G, n_genes = 20, seed = 1)$expression
  expect_warning(simulate_transcriptome(G, mimic = Eb, h2 = 0.9, seed = 2),
                 "`h2` is ignored")
  expect_warning(simulate_transcriptome(G, mimic = Eb, seed = 2), NA)
})

test_that("O4: GREML on an unidentifiable K returns 0 with a flag, not noise", {
  set.seed(1)
  K0 <- matrix(0, 20, 20)
  r0 <- simplePHENOTYPES:::.greml_h2(matrix(rnorm(3 * 20), 3, 20), K0)
  expect_equal(as.numeric(r0), rep(0, 3))
  expect_false(attr(r0, "identifiable"))
  rI <- simplePHENOTYPES:::.greml_h2(matrix(rnorm(3 * 20), 3, 20), diag(20))
  expect_equal(as.numeric(rI), rep(0, 3))
})

test_that("TX-F5: observe_counts returns integer storage and guards overflow", {
  tx <- simulate_transcriptome(G, n_genes = 30, seed = 1)
  expect_identical(typeof(observe_counts(tx, dispersion = 0.1, seed = 2)$counts),
                   "integer")
  expect_identical(typeof(observe_counts(tx, dispersion = 0, seed = 2)$counts),
                   "integer")
  expect_error(observe_counts(tx, baseline = 700, coupling = 0), "integer range")
  expect_error(observe_counts(tx, baseline = 18, coupling = 0), NA)
})

test_that("TX-F8: predict() accepts a Population-backed phenotype_sim and explains n = 1", {
  tx <- simulate_transcriptome(G, n_genes = 20, seed = 1)
  pop <- as_population(G, individuals = 1:30)
  ph <- suppressMessages(
    simulate_phenotype(pop, h2 = 0.5, seed = 1) |> additive(n_qtn = 5))
  p <- predict(tx, ph, seed = 3)
  expect_equal(p$n_ind, 30L)
  expect_error(predict(tx, G[, c(1:5, 6)]), "at least two new individuals")
})

test_that("O8/O13: predict() validates seed/residual and keeps $calibration", {
  tx <- simulate_transcriptome(G, n_genes = 20, seed = 1)
  Gn <- G[, c(1:5, 6:65)]
  expect_error(predict(tx, Gn, seed = 1.5), "seed")
  expect_error(predict(tx, Gn, seed = -1), "seed")
  expect_error(predict(tx, Gn, residual = NA), "residual")
  expect_error(predict(tx, Gn, residual = "yes"), "residual")
  mm <- simulate_transcriptome(G, mimic = tx$expression, seed = 2)
  expect_false(is.null(mm$calibration))
  expect_identical(predict(mm, Gn, seed = 3)$calibration, mm$calibration)
})

test_that("O11: genotype-free mode ignores cis_fraction/epistasis consistently", {
  tx <- simulate_transcriptome(NULL, n_ind = 10, n_genes = 3, cis_fraction = 2,
                               epistasis = 0.5, seed = 1)
  expect_true(all(tx$genes$epistasis_target == 0))
  expect_true(all(tx$genes$h2_target == 0))
  # an explicit h2 = 0 gene never carries an epistasis target either
  tg <- simulate_transcriptome(G, n_genes = 10, h2 = 0, epistasis = 0.5, seed = 1)
  expect_true(all(tg$genes$epistasis_target == 0))
})

test_that("O12: transcriptome(genes =) rejects duplicate and empty gene sets", {
  tx <- simulate_transcriptome(G, n_genes = 30, seed = 1)
  ph <- simulate_phenotype(G, seed = 3, transcriptome = tx)
  g1 <- tx$genes$gene_id[1]
  expect_error(transcriptome(ph, prop = 0.2, genes = c(g1, g1)), "duplicate")
  expect_error(transcriptome(ph, prop = 0.2, genes = character(0)), "empty")
})

test_that("TX-F6/O9: no stale 'planned/deferred' status statements remain", {
  root <- testthat::test_path("..", "..")
  files <- c(file.path(root, "R", "transcriptome_layer.R"),
             file.path(root, "R", "transcriptome_simulate.R"),
             file.path(root, "docs", "SPEC-transcriptome.md"),
             file.path(root, "docs", "DECISIONS.md"))
  files <- files[file.exists(files)]
  skip_if(length(files) == 0L, "source tree not available (installed package)")
  for (f in files) {
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_false(grepl("Still-planned|remaining follow-ups|counts/GRN/tissue/epistasis deferred|separate follow-ups",
                       txt), info = f)
  }
})
