# test-grammar-coverage.R
#
# Coverage for grammar feature paths not exercised by test-grammar.R:
# the vqtl and epistasis layers, the same_as_add = FALSE / degree branches,
# custom and scalar effect series, multiple reps, plain-matrix input, the LD
# architecture annotation/reporting, the write_*/print exporters, and the
# documented error guards. These are property/structural checks (DECISION-009).

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04
# SNP55K is stored markers x individuals (5 metadata cols + one col per
# individual), so the number of individuals is ncol(G) - 5.
n_ind <- ncol(G) - 5L
gen_mat <- function(sim) simplePHENOTYPES:::.genetic_matrix(sim)

# ---------------------------------------------------------------------------
# epistasis layer
# ---------------------------------------------------------------------------
test_that("epistasis contributes its requested variance proportion", {
  ep <- epistasis(simulate_phenotype(G, seed = 2), prop = 0.5, n_pairs = 3,
                  interaction = 2)
  g <- gen_mat(ep)[, 1]
  expect_equal(stats::var(g) / stats::var(ep$pheno$value), 0.5, tolerance = 0.1)
})

test_that("epistasis stores n_pairs x interaction QTN matrices", {
  ep <- epistasis(simulate_phenotype(G, n_traits = 2, seed = 2), prop = 0.3,
                  n_pairs = 4, interaction = 3)
  q <- ep$layers[[1]]$qtn
  expect_equal(ep$layers[[1]]$n_pairs, 4)
  expect_equal(dim(q[[1]]), c(4L, 3L))
  # independent architecture -> distinct pairs per trait
  expect_false(identical(q[[1]], q[[2]]))
})

test_that("epistasis errors when requested markers exceed candidates", {
  M <- matrix(sample(c(-1, 0, 1), 40 * 6, replace = TRUE), nrow = 40)
  expect_error(epistasis(simulate_phenotype(M, seed = 1), prop = 0.3,
                         n_pairs = 5, interaction = 2),
               "exceed candidate markers")
})

# ---------------------------------------------------------------------------
# vqtl layer
# ---------------------------------------------------------------------------
test_that("vqtl(same_as_add = TRUE) reuses the additive QTNs and realizes", {
  ph <- additive(simulate_phenotype(G, seed = 1), prop = 0.4, n_qtn = 5)
  ph <- vqtl(ph, prop = 0.2, same_as_add = TRUE)
  expect_length(ph$layers, 2)
  expect_identical(ph$layers[[1]]$qtn, ph$layers[[2]]$qtn)
  expect_true(all(is.finite(ph$pheno$value)))
})

test_that("vqtl(same_as_add = FALSE) draws its own loci", {
  ph <- additive(simulate_phenotype(G, seed = 1, n_qtn = 5), prop = 0.4)
  expect_warning(
    ph <- vqtl(ph, prop = 0.1, same_as_add = FALSE, n_qtn = 3),
    "overrides the baseline")
  expect_length(ph$layers[[2]]$qtn[[1]], 3)
  expect_false(identical(ph$layers[[1]]$qtn, ph$layers[[2]]$qtn))
})

test_that("vqtl(same_as_add = TRUE) without a prior additive errors", {
  expect_error(vqtl(simulate_phenotype(G, seed = 1), prop = 0.1),
               "requires a prior additive")
})

# ---------------------------------------------------------------------------
# dominance: fresh-loci and degree branches
# ---------------------------------------------------------------------------
test_that("dominance(same_as_add = FALSE) draws fresh dominance loci", {
  ph <- additive(simulate_phenotype(G, seed = 3), prop = 0.4, n_qtn = 4)
  ph <- dominance(ph, prop = 0.1, same_as_add = FALSE, n_qtn = 6, degree = 0.5)
  expect_length(ph$layers[[2]]$qtn[[1]], 6)
  expect_false(identical(ph$layers[[1]]$qtn, ph$layers[[2]]$qtn))
  expect_equal(ph$layers[[2]]$degree, 0.5)
})

test_that("dominance(same_as_add = TRUE) without a prior additive errors", {
  expect_error(dominance(simulate_phenotype(G, seed = 1), prop = 0.1),
               "requires a prior additive")
})

# ---------------------------------------------------------------------------
# effect series: custom vector and scalar geometric base
# ---------------------------------------------------------------------------
test_that("a custom effect vector is used verbatim", {
  expect_equal(simplePHENOTYPES:::.effect_series(3, "geometric",
                                                 effect = c(1, 5, 2)),
               c(1, 5, 2))
})

test_that("a scalar effect is the geometric base", {
  expect_equal(simplePHENOTYPES:::.effect_series(4, "geometric", effect = 0.3),
               0.3 ^ seq_len(4))
})

test_that("an unsupported dist errors", {
  expect_error(simplePHENOTYPES:::.effect_series(3, "normal"), "geometric")
})

# ---------------------------------------------------------------------------
# multiple reps
# ---------------------------------------------------------------------------
test_that("n_reps > 1 produces one block of rows per rep", {
  rr <- additive(simulate_phenotype(G, seed = 1, n_reps = 3), prop = 0.5,
                 n_qtn = 3)
  expect_equal(sort(unique(rr$pheno$rep)), 1:3)
  expect_equal(nrow(rr$pheno), n_ind * 3)
})

# ---------------------------------------------------------------------------
# plain-matrix input and input validation
# ---------------------------------------------------------------------------
test_that("a plain -1/0/1 matrix is accepted as geno", {
  M <- matrix(sample(c(-1, 0, 1), 50 * 20, replace = TRUE), nrow = 50)
  rownames(M) <- paste0("ind", 1:50)
  ph <- additive(simulate_phenotype(M, seed = 1), prop = 0.5, n_qtn = 4)
  expect_s3_class(ph, "phenotype_sim")
  expect_equal(nrow(ph$pheno), 50)
})

test_that("a data frame without the five metadata columns errors", {
  bad <- data.frame(a = 1, b = 2, c = 3, d = 4, e = 5, f = 6)
  expect_error(simulate_phenotype(bad), "first five columns")
})

test_that("requesting more QTNs than markers errors", {
  M <- matrix(sample(c(-1, 0, 1), 40 * 10, replace = TRUE), nrow = 40)
  expect_error(additive(simulate_phenotype(M, seed = 1), prop = 0.5,
                        n_qtn = 999),
               "exceeds")
})

# ---------------------------------------------------------------------------
# LD architecture annotation and reporting
# ---------------------------------------------------------------------------
test_that("LD architecture annotates companion markers within the r2 window", {
  ld <- additive(simulate_phenotype(G, architecture = "ld", seed = 200,
                                    ld_type = "indirect", r2_min = 0.2,
                                    r2_max = 0.8),
                 prop = 0.5, n_qtn = 3)
  ann <- ld$layers[[1]]$ld[[1]]
  expect_named(ann, c("causal", "companion", "r2"))
  ok <- !is.na(ann$r2)
  expect_true(all(ann$r2[ok] >= 0.2 & ann$r2[ok] <= 0.8))
})

test_that("ld_type switches which markers are reported", {
  mk <- function(type) {
    ld <- additive(simulate_phenotype(G, architecture = "ld", seed = 200,
                                      ld_type = type), prop = 0.5, n_qtn = 3)
    simplePHENOTYPES:::.ld_reported_qtn(ld, ld$layers[[1]])[[1]]
  }
  # direct reports the causal markers themselves
  ld <- additive(simulate_phenotype(G, architecture = "ld", seed = 200,
                                    ld_type = "direct"), prop = 0.5, n_qtn = 3)
  expect_identical(mk("direct"), ld$layers[[1]]$qtn[[1]])
})

# ---------------------------------------------------------------------------
# exporters and print
# ---------------------------------------------------------------------------
test_that("write_phenotypes round-trips long and wide tables", {
  ph <- additive(simulate_phenotype(G, n_traits = 2, seed = 2), prop = 0.4,
                 n_qtn = 3)
  fl <- tempfile(fileext = ".txt")
  fw <- tempfile(fileext = ".txt")
  on.exit(unlink(c(fl, fw)), add = TRUE)

  write_phenotypes(ph, fl, format = "long")
  write_phenotypes(ph, fw, format = "wide")
  back_l <- data.table::fread(fl)
  back_w <- data.table::fread(fw)

  expect_equal(nrow(back_l), nrow(phenotypes_long(ph)))
  expect_true(all(c("Trait_1", "Trait_2") %in% names(back_w)))
  expect_equal(nrow(back_w), n_ind)
})

test_that("print.phenotype_sim reports the variance budget", {
  ph <- additive(simulate_phenotype(G, seed = 1), prop = 0.5, n_qtn = 3)
  out <- utils::capture.output(print(ph))
  expect_true(any(grepl("phenotype_sim", out)))
  expect_true(any(grepl("additive", out)))
  expect_true(any(grepl("residual", out)))
})

# ---------------------------------------------------------------------------
# complex_phenotypes error guards
# ---------------------------------------------------------------------------
test_that("complex_phenotypes rejects fewer than two inputs", {
  a <- additive(simulate_phenotype(G, n_traits = 2, seed = 1), prop = 0.3,
                n_qtn = 3)
  expect_error(complex_phenotypes(a, h2 = 0.5), "at least two")
})

test_that("complex_phenotypes rejects mismatched n_traits", {
  a <- additive(simulate_phenotype(G, n_traits = 2, seed = 1), prop = 0.3,
                n_qtn = 3)
  b <- additive(simulate_phenotype(G, n_traits = 3, seed = 1), prop = 0.3,
                n_qtn = 3)
  expect_error(complex_phenotypes(a, b, h2 = 0.5), "same n_traits")
})

test_that("complex_phenotypes rejects inputs on different genotypes", {
  M <- matrix(sample(c(-1, 0, 1), 50 * 20, replace = TRUE), nrow = 50)
  rownames(M) <- paste0("x", seq_len(50))
  a <- additive(simulate_phenotype(G, seed = 1), prop = 0.3, n_qtn = 3)
  b <- additive(simulate_phenotype(M, seed = 1), prop = 0.3, n_qtn = 3)
  expect_error(complex_phenotypes(a, b, h2 = 0.5), "same genotypes")
})
