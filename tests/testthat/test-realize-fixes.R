# Theory-review fixes O4-O8 (grammar_realize.R / select_ind.R).

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04

.rf_tx <- function() simulate_transcriptome(G, n_genes = 60, seed = 1)
.rf_sim <- function(h2, a_prop, tx_prop, tx = .rf_tx()) {
  simulate_phenotype(G, h2 = h2, seed = 2, transcriptome = tx) |>
    additive(prop = a_prop, n_qtn = 3) |>
    transcriptome(prop = tx_prop, n_genes = 10)
}

test_that("O4: marker layers must still sum to h2 when a transcriptome layer is present", {
  s <- .rf_sim(0.5, 0.1, 0.2)
  expect_error(genetic_values(s), "Incomplete h2 allocation")
  expect_error(genetic_values(.rf_sim(0.5, 0.1, 0.5)), "Incomplete h2 allocation")
  # markers fill h2 alone (transcriptome outside the budget), or markers +
  # transcriptome fill it exactly: both are accepted
  expect_error(genetic_values(.rf_sim(0.3, 0.3, 0.4)), NA)
  expect_error(genetic_values(.rf_sim(0.5, 0.2, 0.3)), NA)
})

test_that("O5: breeding value is refused for a genome-mediated transcriptome layer", {
  tx <- .rf_tx()
  s <- .rf_sim(0.3, 0.3, 0.3, tx)
  expect_gt(stats::sd(genetic_values(s)), 0)
  expect_error(.breeding_value_matrix(s), "transcriptome")
  expect_error(select_ind(s, n = 10, on = "bv"), "transcriptome")
  # a prop = 0 transcriptome layer is not genetic: breeding value still defined
  s0 <- .rf_sim(0.3, 0.3, 0, tx)
  expect_error(.breeding_value_matrix(s0), NA)
})

test_that("O6: on = 'gv' equals genetic_values() for transcriptome sims", {
  s <- .rf_sim(0.3, 0.3, 0.3)
  gv <- genetic_values(s)[, 1]
  expect_equal(unname(.criterion_values(s, "gv", 1L, 1L)), unname(gv))
  expect_gt(stats::sd(.criterion_values(s, "gv", 1L, 1L)),
            stats::sd(.genetic_matrix(s, 1L)[, 1]))
  sel <- select_ind(s, n = 10, on = "gv")
  expect_equal(sort(unname(sel)),
               sort(names(sort(gv, decreasing = TRUE))[1:10]))
})

test_that("O7: record-scale realized H2 does not depend on pheno row order", {
  s <- simulate_phenotype(G, h2 = 0.5, seed = 2, reps = 4) |>
    additive(n_qtn = 5)
  h0 <- .realized_h2(s, "record")
  s2 <- s
  set.seed(1)
  s2$pheno <- s$pheno[sample(nrow(s$pheno)), ]
  expect_equal(.realized_h2(s2, "record"), h0)
  expect_equal(.realized_h2(s2, "phenotype"), .realized_h2(s, "phenotype"))
})

test_that("O8: breeding value is twice the expected progeny deviation", {
  # one locus, a = 1, d = 0, p = 0.5, parent gene content x = 2. A parent
  # transmits one copy; the random mate contributes 2p/2 = p expected copies, so
  # the progeny mean is 1 + p against a population mean of 2p.
  p <- 0.5; x <- 2
  alpha <- .avg_effect(1, 0, p)
  A <- alpha * (x - 2 * p)
  prog_dev <- alpha * ((x / 2 + p) - 2 * p)
  expect_equal(A, 1)
  expect_equal(prog_dev, A / 2)
})
