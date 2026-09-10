data("SNP55K_maize282_maf04")
G_release <- SNP55K_maize282_maf04

test_that("simulation counts, proportions, models, and flags are validated", {
  expect_error(simulate_phenotype(G_release, n_traits = 1.5), "whole")
  expect_error(simulate_phenotype(G_release, n_reps = 0), "positive")
  expect_error(simulate_phenotype(G_release, n_qtn = -1), ">= 0")
  expect_error(simulate_phenotype(G_release, h2 = 1.1), "between 0 and 1")
  expect_error(simulate_phenotype(G_release, vary_qtn = 1), "TRUE or FALSE")
  expect_error(simulate_phenotype(G_release, model = "AX"), "arg")
  expect_error(simulate_phenotype(G_release, model = "AD"), "only used")
  expect_error(simulate_phenotype(G_release, n_qtn = 2, h2 = 0.5,
                                  model = "AA"), "arg")
})

test_that("causal loci must carry usable variation", {
  M <- cbind(monomorphic = rep(1, 100),
             polymorphic = rep(c(-1, 1), 50))
  ph <- simulate_phenotype(M, seed = 1) |>
    additive(prop = 0.5, n_qtn = 1)
  expect_identical(ph$layers[[1]]$qtn[[1]], 2L)
  expect_error(
    simulate_phenotype(M, seed = 1) |>
      additive(prop = 0.5, qtn = "monomorphic"),
    "zero usable variation"
  )
})

test_that("replicate-aware accessors report each varying architecture", {
  ph <- simulate_phenotype(G_release, n_reps = 3, vary_qtn = TRUE,
                           seed = 81) |>
    additive(prop = 0.5, n_qtn = 4)
  expect_false(isTRUE(all.equal(genetic_values(ph, rep = 1),
                                genetic_values(ph, rep = 2))))
  expect_identical(qtn_table(ph, rep = 2)$snp,
                   ph$map$snp[ph$layers[[1]]$qtn_reps[[2]][[1]]])
  expect_error(genetic_values(ph, rep = 4), "between 1 and n_reps")
  expect_error(qtn_table(ph, rep = 1.5), "whole")
})

test_that("complex models preserve replications, means, and genetic values", {
  a <- simulate_phenotype(G_release, n_traits = 2, n_reps = 2,
                          vary_qtn = TRUE, mean = c(10, 20), seed = 1) |>
    additive(prop = 0.3, n_qtn = 4)
  b <- simulate_phenotype(G_release, n_traits = 2, n_reps = 2,
                          vary_qtn = TRUE, mean = c(10, 20), seed = 1) |>
    epistasis(prop = 0.2, n_pairs = 2)
  z <- complex_phenotypes(a, b, h2 = 0.5)
  expect_identical(sort(unique(z$pheno$rep)), 1:2)
  expect_equal(nrow(z$pheno), z$n_ind * z$n_traits * z$n_reps)
  expect_equal(unname(apply(genetic_values(z, rep = 1), 2, stats::var)),
               c(0.5, 0.5), tolerance = 1e-10)
  expect_equal(as.numeric(tapply(z$pheno$value, z$pheno$trait, mean)),
               c(10, 20), tolerance = 0.25)
  expect_true(all(is.finite(simplePHENOTYPES:::.realized_h2(z))))
})

test_that("complex models reject genotype data that only share IDs", {
  M1 <- matrix(sample(c(-1, 0, 1), 80 * 20, replace = TRUE), 80, 20,
               dimnames = list(paste0("i", 1:80), paste0("m", 1:20)))
  M2 <- M1
  M2[1, 1] <- if (M2[1, 1] == 1) -1 else 1
  a <- simulate_phenotype(M1, seed = 1) |> additive(prop = 0.3, n_qtn = 3)
  b <- simulate_phenotype(M2, seed = 1) |> additive(prop = 0.3, n_qtn = 3)
  expect_error(complex_phenotypes(a, b, h2 = 0.5), "same genotype data")
})

test_that("vQTL variance is budgeted as residual heterogeneity, not H2", {
  ph <- simulate_phenotype(G_release, seed = 3) |>
    additive(prop = 0.3, n_qtn = 4) |>
    vqtl(prop = 0.2)
  expect_equal(simplePHENOTYPES:::.total_genetic_prop(ph), 0.3)
  expect_equal(simplePHENOTYPES:::.total_variance_prop(ph), 0.5)
  expect_equal(ph$var_budget$prop[ph$var_budget$component == "residual"], 0.5)
  expect_equal(stats::var(genetic_values(ph)[, 1]), 0.3, tolerance = 1e-10)
  expect_error(vqtl(simulate_phenotype(G_release, h2 = 0.5, seed = 3),
                    same_as_add = FALSE, n_qtn = 2), "explicit `prop`")
})

test_that("pleiotropy inputs define a valid correlation model", {
  expect_error(simulate_phenotype(G_release, architecture = "pleiotropy",
                                  n_traits = 2, cor = 2), "between -1 and 1")
  expect_error(simulate_phenotype(G_release, architecture = "pleiotropy",
                                  n_traits = 2, cor = matrix(c(1, .2, .3, 1), 2)),
               "symmetric")
  expect_error(simulate_phenotype(G_release, architecture = "pleiotropy",
                                  n_traits = 2, pi = 1, pi_target = 1),
               "either `pi`")
})

test_that("reference-oriented HapMap coding follows the declared allele", {
  hmp_names <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                 "center", "protLSID", "assayLSID", "panelLSID", "QCcode")
  meta <- as.data.frame(matrix(c("s1", "A/G", "1", "100", rep("?", 7)),
                               nrow = 1, dimnames = list(NULL, hmp_names)),
                        stringsAsFactors = FALSE)
  dat <- cbind(meta, data.frame(i1 = "AA", i2 = "GG", i3 = "GG", i4 = "GG"))
  a_ref <- as_numeric(dat, to_r = TRUE, method = "reference",
                      ref_allele = "A", verbose = FALSE)
  g_ref <- as_numeric(dat, to_r = TRUE, method = "reference",
                      ref_allele = "G", verbose = FALSE)
  expect_equal(as.integer(a_ref[1, 6:9]), c(1, -1, -1, -1))
  expect_equal(as.integer(g_ref[1, 6:9]), c(-1, 1, 1, 1))
})

test_that("plot dots are applied or rejected, never ignored", {
  ph <- simulate_phenotype(G_release, seed = 1) |>
    additive(prop = 0.5, n_qtn = 3)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_invisible(plot(ph, which = "hist", mar = c(3, 3, 2, 1)))
  expect_error(plot(ph, which = "hist", imaginary_parameter = 1),
               "Unsupported")
  grDevices::dev.off()
  unlink(tmp)
})
