# test-additive-effect-list.R
#
# additive(effect = list(...)): one effect specification per trait in a single
# layer (breedingDesigner SPEC-0002 re-scores frozen per-trait effects this way
# instead of one trait-masked layer per trait).

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04
Q <- c("ss196442916", "ss196439337", "ss196480535")

test_that("a per-trait effect list gives each trait its own series", {
  sim <- simulate_phenotype(G, n_traits = 2, seed = 1) |>
    additive(prop = 0.4, qtn = Q, effect = list(c(3, 2, 1), c(1, 1, 1)))
  e <- sim$layers[[1]]$effect
  expect_equal(e[[1]], c(3, 2, 1))
  expect_equal(e[[2]], c(1, 1, 1))
  # a scalar per trait is a geometric base, as for the unlisted form
  sim2 <- simulate_phenotype(G, n_traits = 2, seed = 1) |>
    additive(prop = 0.4, qtn = Q, effect = list(0.5, 0.9))
  expect_equal(sim2$layers[[1]]$effect[[2]], 0.9^(1:3))
})

test_that("each trait matches a single-trait simulation with its series", {
  e1 <- c(3, 2, 1); e2 <- c(-1, 4, 0.5)
  two <- simulate_phenotype(G, n_traits = 2, seed = 7) |>
    additive(prop = 0.4, qtn = Q, effect = list(e1, e2))
  g2 <- genetic_values(two)
  for (t in 1:2) {
    one <- simulate_phenotype(G, n_traits = 1, seed = 7) |>
      additive(prop = 0.4, qtn = Q, effect = list(e1, e2)[[t]])
    expect_equal(unname(g2[, t]), unname(genetic_values(one)[, 1]))
  }
})

test_that("re-scoring a frozen pleiotropic template reproduces its values", {
  tpl <- suppressMessages(suppressWarnings(
    simulate_phenotype(G, architecture = "pleiotropy", n_traits = 2,
                       cor = -0.6, pi = 0.8, seed = 11) |>
      additive(prop = 0.5, n_qtn = 20)))
  ly <- tpl$layers[[1]]
  frozen <- simulate_phenotype(G, n_traits = 2, seed = 11) |>
    additive(prop = 0.5, qtn = ly$qtn, effect = ly$effect)
  expect_equal(genetic_values(frozen), genetic_values(tpl))
})

test_that("invalid per-trait effect lists error clearly", {
  base <- simulate_phenotype(G, n_traits = 2, seed = 1)
  expect_error(additive(base, prop = 0.4, qtn = Q, effect = list(c(1, 2, 3))),
               "one element per trait \\(2\\); got 1")
  expect_error(additive(base, prop = 0.4, qtn = Q, effect = list(1:2, 1:3)),
               "length n_qtn")
  expect_error(additive(base, prop = 0.4, qtn = Q, effect = list(NULL, 1:3)),
               "element\\(s\\) 1 are not")
  expect_error(additive(base, prop = 0.4, qtn = Q, effect = list("a", 1:3)),
               "element\\(s\\) 1 are not")
  pl <- simulate_phenotype(G, architecture = "pleiotropy", n_traits = 2,
                           cor = 0.3, seed = 1)
  expect_error(suppressMessages(additive(pl, prop = 0.4, n_qtn = 3,
                                         effect = list(0.5, 0.5))),
               "cannot be used under")
  expect_error(additive(base, prop = 0.4, n_qtn = 3, orthogonal = TRUE,
                        a = list(1, 2)), "`a` must be a scalar")
})

test_that("epistasis() names n_pairs in its effect-length error", {
  expect_error(simulate_phenotype(G, seed = 1) |>
                 epistasis(prop = 0.2, n_pairs = 3, effect = c(1, 2)),
               "length n_pairs")
})

test_that("the orthogonal model names `a` in its effect errors", {
  pop <- as_population(G, individuals = 1:20)
  f2 <- selfcross(cross(pop[1], pop[2], n = 1, seed = 1), n = 40, seed = 2)
  expect_error(simulate_phenotype(f2, h2 = 0.6, seed = 3) |>
                 additive(orthogonal = TRUE, a = c(1, 2), n_qtn = 5),
               "`a` must be either")
})
