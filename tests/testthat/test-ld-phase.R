# ld_phase: haplotype-derived coupling / repulsion for architecture = "ld".
G <- SNP55K_maize282_maf04

ld_sim <- function(phase = NULL, s = 3, type = "direct") {
  args <- list(G, architecture = "ld", n_traits = 2, ld_type = type,
               r2_min = 0.2, r2_max = 0.8, n_qtn = 6, seed = s)
  if (!is.null(phase)) args$ld_phase <- phase
  do.call(simulate_phenotype, args) |> additive(prop = 0.4)
}

pair_sign <- function(ph) {
  lay <- ph$layers[[1]]
  r <- attr(lay$ld, "r")
  sign(lay$effect[[1]] * lay$effect[[2]] * r)
}

test_that("default ld_phase is bit-identical to 'coded'", {
  a <- ld_sim(NULL)
  b <- ld_sim("coded")
  expect_identical(a$layers[[1]]$qtn, b$layers[[1]]$qtn)
  expect_identical(a$layers[[1]]$effect, b$layers[[1]]$effect)
})

test_that("coupling makes every pair sign +1 and repulsion -1, moving only trait 2", {
  base <- ld_sim("coded")
  for (type in c("direct", "indirect")) {
    cp <- ld_sim("coupling", type = type)
    rp <- ld_sim("repulsion", type = type)
    expect_true(all(pair_sign(cp) == 1))
    expect_true(all(pair_sign(rp) == -1))
  }
  cp <- ld_sim("coupling"); rp <- ld_sim("repulsion")
  expect_identical(cp$layers[[1]]$qtn, base$layers[[1]]$qtn)
  expect_identical(rp$layers[[1]]$effect[[1]], base$layers[[1]]$effect[[1]])
  expect_equal(abs(cp$layers[[1]]$effect[[2]]), abs(base$layers[[1]]$effect[[2]]))
  expect_equal(attr(cp$layers[[1]]$ld, "r"), attr(base$layers[[1]]$ld, "r"))
})

test_that("the stored signed r matches the dosage correlation and r2", {
  ph <- ld_sim("repulsion")
  ld <- ph$layers[[1]]$ld
  r <- attr(ld, "r")
  expect_equal(r^2, ld$r2, tolerance = 1e-8)
  for (i in seq_len(nrow(ld))) {
    x <- .geno_cols(ph, c(ld$qtn_t1[i], ld$qtn_t2[i]))
    expect_equal(r[i], as.numeric(stats::cor(x[, 1], x[, 2])), tolerance = 1e-8)
  }
})

test_that("ld_phase wins over additive(phase =) and is validated", {
  ph <- simulate_phenotype(G, architecture = "ld", n_traits = 2,
                           ld_phase = "repulsion", n_qtn = 6, seed = 3) |>
    additive(prop = 0.4, phase = "repulsion")
  expect_true(all(pair_sign(ph) == -1))
  expect_error(simulate_phenotype(G, architecture = "ld", n_traits = 2,
                                  ld_phase = "bogus"), "ld_phase")
  expect_error(simulate_phenotype(G, architecture = "pleiotropy", n_traits = 2,
                                  ld_phase = "repulsion"))
})
