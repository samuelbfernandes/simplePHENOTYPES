# liability_threshold(), coheritability(), cor_ar1() (wish-list, 2026-10-03)

.tc_g <- function() {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES", envir = environment())
  SNP55K_maize282_maf04
}

test_that("liability_threshold() cuts the standardized liability at qnorm(cumsum(prop))", {
  g <- .tc_g()
  base <- simulate_phenotype(g, h2 = 0.5, n_qtn = 10, seed = 1)
  sim <- liability_threshold(base, prop = c(0.25, 0.5, 0.25))
  # the liability is the unthresholded phenotype, unchanged
  expect_identical(sim$liability, base$pheno)
  l <- base$pheno$value
  z <- (l - mean(l)) / stats::sd(l)
  expect_identical(sim$pheno$value,
                   as.numeric(findInterval(z, stats::qnorm(c(0.25, 0.75))) + 1L))
  expect_setequal(unique(sim$pheno$value), 1:3)
  # genetic side and the variance budget stay on the liability scale
  expect_identical(genetic_values(sim), genetic_values(base))
  expect_identical(sim$var_budget, base$var_budget)
  expect_equal(simplePHENOTYPES:::.realized_h2(sim),
               simplePHENOTYPES:::.realized_h2(base))
  # binary trait: prevalence ~ prop[2]
  b <- liability_threshold(base, prop = c(0.8, 0.2))
  expect_lt(abs(mean(b$pheno$value == 2) - 0.2), 0.06)   # ~2.5 SE at n = 282
  expect_output(print(b), "Liability-threshold")
})

test_that("liability_threshold() validates its input and keeps per-trait specs", {
  g <- .tc_g()
  two <- simulate_phenotype(g, n_traits = 2, h2 = 0.5, n_qtn = 10, seed = 2)
  expect_error(liability_threshold(two, prop = c(0.5, 0.6)), "summing to 1")
  expect_error(liability_threshold(two, prop = 1), "at least two")
  expect_error(liability_threshold(two, prop = c(0.5, 0.5), trait = 3), "trait")
  one <- liability_threshold(two, prop = c(0.5, 0.5), trait = 2)
  t1 <- one$pheno$value[one$pheno$trait == "Trait_1"]
  t2 <- one$pheno$value[one$pheno$trait == "Trait_2"]
  expect_gt(length(unique(t1)), 2L)                 # trait 1 stays continuous
  expect_setequal(unique(t2), 1:2)
  # a layer added after the threshold keeps it
  th <- simulate_phenotype(g, seed = 3) |> liability_threshold(prop = c(0.5, 0.5)) |>
    additive(prop = 0.4, n_qtn = 5)
  expect_setequal(unique(th$pheno$value), 1:2)
})

test_that("a stored trait keeps the base thresholds, so prevalence moves under selection", {
  g <- .tc_g()
  pop <- suppressMessages(as_population(g, individuals = 1:40))
  f2 <- selfcross(cross(pop[1], pop[2], n = 1, seed = 1), n = 100, seed = 2)
  s0 <- simulate_phenotype(f2, h2 = 0.6, n_qtn = 20, seed = 3) |>
    liability_threshold(prop = c(0.8, 0.2))
  top <- select_ind(s0, prop = 0.2, on = "gv")
  expect_false(is.null(top$trait$threshold_cut[[1L]]))
  f3 <- selfcross(top[1], n = 80, seed = 4)
  s3 <- simulate_phenotype(f3, seed = 5)
  expect_true(isTRUE(s3$frozen))
  cut <- top$trait$threshold_cut[[1L]]
  expect_identical(s3$pheno$value,
                   as.numeric(findInterval(s3$liability$value, cut) + 1L))
})

test_that("coheritability() is Cov(G)/sqrt(Vp Vp) with the realized H2 on the diagonal", {
  g <- .tc_g()
  sim <- suppressMessages(simulate_phenotype(g, n_traits = 2,
                                             architecture = "pleiotropy",
                                             cor = 0.6, h2 = c(0.5, 0.3),
                                             n_qtn = 20, seed = 1))
  ch <- coheritability(sim)
  G <- genetic_values(sim)
  P <- cbind(sim$pheno$value[sim$pheno$trait == "Trait_1"],
             sim$pheno$value[sim$pheno$trait == "Trait_2"])
  expect_equal(unname(ch[1, 2]), stats::cov(G[, 1], G[, 2]) /
                 (stats::sd(P[, 1]) * stats::sd(P[, 2])), tolerance = 1e-12)
  expect_equal(unname(ch[1, 2]),
               stats::cor(G)[1, 2] * sqrt(ch[1, 1] * ch[2, 2]), tolerance = 1e-12)
  expect_true(isSymmetric(ch))
  expect_identical(dimnames(ch), list(c("Trait_1", "Trait_2"), c("Trait_1", "Trait_2")))
})

test_that("cor_ar1() is rho^|i-j| and works as cor / resid_cor", {
  m <- cor_ar1(4, 0.5)
  expect_equal(unname(m[1, ]), 0.5^(0:3))
  expect_true(isSymmetric(m))
  expect_gt(min(eigen(m, only.values = TRUE)$values), 0)
  expect_error(cor_ar1(3, 1), "\\(-1, 1\\)")
  expect_error(cor_ar1(1, 0.5), "n_traits")
  g <- .tc_g()
  sim <- suppressMessages(simulate_phenotype(g, n_traits = 4,
                                             architecture = "pleiotropy",
                                             cor = cor_ar1(4, 0.9),
                                             resid_cor = cor_ar1(4, 0.5),
                                             h2 = 0.4, n_qtn = 30, seed = 1))
  expect_s3_class(sim, "phenotype_sim")
  rg <- stats::cor(genetic_values(sim))
  # adjacent genetic correlation exceeds the lag-3 one
  expect_gt(rg[1, 2], rg[1, 4])
})

test_that("variance reports stay on the liability scale after thresholding (Codex G-01)", {
  g <- .tc_g()
  # the inbred panel has few heterozygotes: dominance() notes it (expected here)
  base <- suppressWarnings(simulate_phenotype(g, h2 = 0.6, seed = 4) |>
    additive(prop = 0.4, n_qtn = 5) |> dominance(prop = 0.2))
  th <- liability_threshold(base, prop = c(0.7, 0.3))
  expect_equal(qtn_table(th)$var_explained, qtn_table(base)$var_explained,
               tolerance = 1e-12)
  expect_equal(th$ad_report, base$ad_report, tolerance = 1e-12)
})

test_that("thresholding an unseeded simulation keeps its liability (Codex G-02)", {
  g <- .tc_g()
  base <- simulate_phenotype(g, h2 = 0.5, n_qtn = 10)          # seed = NULL
  th <- liability_threshold(base, prop = c(0.5, 0.5))
  expect_identical(th$liability, base$pheno)
})
