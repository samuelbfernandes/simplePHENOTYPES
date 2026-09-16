# The transcriptome() phenotype layer: a continuous-predictor component scored on
# gene expression (real via expression=, or genome-derived via transcriptome=).
# It is a distinct *expression-mediated* variance category, NOT folded into the
# marker-based (broad-sense) heritability.

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04
tx <- simulate_transcriptome(G, n_genes = 120, seed = 1)

tx_comp <- function(ph) simplePHENOTYPES:::.transcriptome_matrix(ph, 1L)[, 1]
h2_narrow <- function(ph) {
  stats::var(genetic_values(ph)[, 1]) / stats::var(ph$pheno$value)
}

test_that("the transcriptome component realizes its target variance share", {
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 3, transcriptome = tx) |>
    transcriptome(prop = 0.5, n_genes = 25)
  # the expression-mediated component realizes prop = 0.5 of phenotypic variance
  expect_equal(stats::var(tx_comp(ph)) / stats::var(ph$pheno$value), 0.5,
               tolerance = 0.06)
})

test_that("the transcriptome layer is a separate category, not marker-based H2", {
  # additive (0.2) is marker-based genetic; transcriptome (0.4) is not
  ph <- simulate_phenotype(G, h2 = 0.6, seed = 2, transcriptome = tx) |>
    transcriptome(prop = 0.4, n_genes = 20) |>
    additive(prop = 0.2, n_qtn = 10)
  expect_equal(h2_narrow(ph), 0.2, tolerance = 0.06)          # marker-based only
  expect_equal(stats::var(tx_comp(ph)) / stats::var(ph$pheno$value), 0.4,
               tolerance = 0.06)
  # a purely non-genetic (h2 = 0) derived transcriptome must NOT inflate H2
  tx0 <- simulate_transcriptome(G, n_genes = 40, h2 = 0, seed = 8)
  ph0 <- simulate_phenotype(G, seed = 9, transcriptome = tx0) |>
    transcriptome(prop = 0.5, n_genes = 20)
  expect_lt(h2_narrow(ph0), 1e-6)                             # no genetic value
})

test_that("a real expression matrix works as a predictor alongside markers", {
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 3, expression = tx$expression) |>
    transcriptome(prop = 0.3, n_genes = 15) |>
    additive(prop = 0.2, n_qtn = 5)
  expect_s3_class(ph, "phenotype_sim")
  expect_identical(ph$expression_source, "real")
  expect_equal(h2_narrow(ph), 0.2, tolerance = 0.06)
  expect_equal(stats::var(tx_comp(ph)) / stats::var(ph$pheno$value), 0.3,
               tolerance = 0.06)
})

test_that("transcriptome = TRUE derives an expression source from the genome", {
  ph <- simulate_phenotype(G, seed = 4, transcriptome = TRUE) |>
    transcriptome(prop = 0.4, n_genes = 20)
  expect_s3_class(ph, "phenotype_sim")
  expect_identical(ph$expression_source, "derived")
  expect_equal(ncol(ph$expression), ph$n_ind)
})

test_that("the transcriptome layer is reproducible under a seed", {
  mk <- function() {
    simulate_phenotype(G, h2 = 0.6, seed = 2, transcriptome = tx) |>
      transcriptome(prop = 0.4, n_genes = 20) |>
      additive(prop = 0.2, n_qtn = 10)
  }
  expect_equal(mk()$pheno$value, mk()$pheno$value)
})

test_that("explicit causal genes and slopes are honored", {
  gid <- rownames(tx$expression)[c(3, 7, 11)]
  ph <- simulate_phenotype(G, h2 = 0.5, seed = 5, transcriptome = tx) |>
    transcriptome(prop = 0.5, genes = gid, slopes = c(1, -2, 0.5))
  ly <- ph$layers[[1]]
  expect_identical(ly$type, "transcriptome")
  expect_equal(ly$effect[[1]], c(1, -2, 0.5))
  expect_identical(rownames(tx$expression)[ly$qtn[[1]]], gid)
})

test_that("the transcriptome layer validates its inputs", {
  expect_error(
    simulate_phenotype(G, seed = 1) |> transcriptome(prop = 0.3, n_genes = 5),
    "no expression source")
  expect_error(
    simulate_phenotype(G, expression = tx$expression, transcriptome = tx),
    "only one of")
  expect_error(
    simulate_phenotype(G, transcriptome = tx) |>
      transcriptome(prop = 0.3, n_genes = 1e6),
    "exceeds")
  expect_error(
    simulate_phenotype(G, transcriptome = tx) |> transcriptome(prop = 0.3),
    "n_genes")
  expect_error(
    simulate_phenotype(G, expression = tx$expression[, 1:10]),
    "missing individual")
})

test_that("transcriptome prop is a phenotypic share, not part of the h2 budget", {
  # h2 = 0 (no marker genetics) still allows a transcriptome component
  ph <- simulate_phenotype(G, h2 = 0, seed = 10, transcriptome = tx) |>
    transcriptome(prop = 0.3, n_genes = 15)
  expect_lt(h2_narrow(ph), 1e-6)
  expect_equal(stats::var(tx_comp(ph)) / stats::var(ph$pheno$value), 0.3,
               tolerance = 0.06)
  # prop is required (not drawn from the remaining h2)
  expect_error(
    simulate_phenotype(G, transcriptome = tx) |> transcriptome(n_genes = 5),
    "explicit .prop")
})

test_that("only relative slopes matter, even for extreme finite slopes", {
  gid <- rownames(tx$expression)[1:3]
  base <- simulate_phenotype(G, h2 = 0.5, seed = 6, transcriptome = tx) |>
    transcriptome(prop = 0.4, genes = gid, slopes = c(1, -2, 0.5))
  huge <- simulate_phenotype(G, h2 = 0.5, seed = 6, transcriptome = tx) |>
    transcriptome(prop = 0.4, genes = gid, slopes = c(1e300, -2e300, 0.5e300))
  expect_equal(tx_comp(base), tx_comp(huge), tolerance = 1e-8)
})

test_that("a large valid seed does not overflow the layer seed", {
  ph <- simulate_phenotype(G, seed = 3000000, transcriptome = tx) |>
    transcriptome(prop = 0.3, n_genes = 10)
  expect_s3_class(ph, "phenotype_sim")
})

test_that("duplicate expression ids and constant causal genes are rejected", {
  # duplicate individual (column) ids -> ambiguous alignment
  E <- tx$expression
  colnames(E)[2] <- colnames(E)[1]
  expect_error(simulate_phenotype(G, expression = E), "duplicate individual")
  # duplicate gene (row) names -> ambiguous named-gene selection
  Er <- tx$expression
  rownames(Er)[2] <- rownames(Er)[1]
  expect_error(simulate_phenotype(G, expression = Er), "duplicate gene")
  # a constant-expression gene cannot be a causal gene
  Ec <- tx$expression
  Ec[1, ] <- 5                                    # gene 1 is now constant
  gid <- rownames(Ec)[1]
  expect_error(
    simulate_phenotype(G, expression = Ec) |>
      transcriptome(prop = 0.4, genes = c(gid, rownames(Ec)[2])),
    "constant expression")
})
