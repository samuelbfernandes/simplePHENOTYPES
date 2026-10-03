# breedingDesigner requests of 2026-10-03:
#   DECISION-047 -- gamma crossover interference (nu = 2.6, p = 0) is the default
#                   meiosis; interference = "poisson" keeps the isqg stream.
#   DECISION-048 -- a Population carries the trait (loci, effects, residual
#                   variance) defined in its base population; simulate_phenotype()
#                   reuses it (refit = FALSE) unless refit = TRUE.

bd_pop <- function(n = 40) {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES", envir = environment())
  suppressMessages(as_population(SNP55K_maize282_maf04, individuals = seq_len(n)))
}

# ---- DECISION-047: gamma interference default ---------------------------------

test_that("the default interference model is gamma with nu = 2.6, p = 0", {
  withr::local_options(simplePHENOTYPES.interference = NULL)
  chk <- simplePHENOTYPES:::.check_interference
  expect_identical(chk(NULL), list(nu = 2.6, p = 0))
  expect_identical(chk("poisson"), "poisson")
  expect_identical(chk("Poisson"), "poisson")
  # both resolved values are fixed points, so forwarding never re-resolves them
  expect_identical(chk(chk("poisson")), "poisson")
  expect_identical(chk(chk(NULL)), list(nu = 2.6, p = 0))
  expect_error(chk("gamma", "cross"), "cross\\(\\).*\"poisson\"")
  expect_error(chk(c("poisson", "poisson"), "cross"), "\"poisson\"")
  expect_error(chk(NA_character_, "cross"), "\"poisson\"")
})

test_that("a call without interference draws the gamma model, not Poisson", {
  withr::local_options(simplePHENOTYPES.interference = NULL)
  pop <- bd_pop(4)
  def <- suppressMessages(cross(pop[1], pop[2], n = 6, seed = 3))
  gam <- cross(pop[1], pop[2], n = 6, seed = 3,
               interference = list(nu = 2.6, p = 0))
  poi <- cross(pop[1], pop[2], n = 6, seed = 3, interference = "poisson")
  expect_identical(dosages(def), dosages(gam))
  expect_false(identical(dosages(def), dosages(poi)))
  expect_identical(dosages(selfcross(pop[1], n = 4, seed = 5)),
                   dosages(selfcross(pop[1], n = 4, seed = 5,
                                     interference = list(nu = 2.6))))
  expect_identical(dosages(double_haploid(pop[1], n = 4, seed = 5)),
                   dosages(double_haploid(pop[1], n = 4, seed = 5,
                                          interference = list(nu = 2.6))))
})

test_that("interference = \"poisson\" is the isqg stream (.draw_meiosis Poisson branch)", {
  m <- list(c(0, 0.4, 1.3), c(0.1, 0.9))
  set.seed(11)
  a <- simplePHENOTYPES:::.draw_meiosis(m, 7, "poisson")
  set.seed(11)
  b <- simplePHENOTYPES:::.draw_meiosis(m, 7, NULL)
  expect_identical(a, b)
  # the explicit argument equals the session option
  pop <- bd_pop(4)
  x <- cross(pop[1], pop[2], n = 5, seed = 8, interference = "poisson")
  withr::local_options(simplePHENOTYPES.interference = "poisson")
  y <- cross(pop[1], pop[2], n = 5, seed = 8)
  expect_identical(dosages(x), dosages(y))
})

test_that("an explicit \"poisson\" overrides a gamma session option, and schemes forward it", {
  pop <- bd_pop(4)
  withr::local_options(simplePHENOTYPES.interference = list(nu = 4))
  a <- cross(pop[1], pop[2], n = 5, seed = 8, interference = "poisson")
  withr::local_options(simplePHENOTYPES.interference = NULL)
  b <- cross(pop[1], pop[2], n = 5, seed = 8, interference = "poisson")
  expect_identical(dosages(a), dosages(b))
  # a scheme resolved to Poisson stays Poisson in its inner selfcross() calls
  f1 <- cross(pop[1], pop[2], n = 6, seed = 1, interference = "poisson")
  s1 <- single_seed_descent(f1, generations = 2, seed = 4, interference = "poisson")
  withr::local_options(simplePHENOTYPES.interference = "poisson")
  s2 <- single_seed_descent(f1, generations = 2, seed = 4)
  expect_identical(dosages(s1), dosages(s2))
})

test_that("the default keeps one crossover per Morgan per gamete", {
  withr::local_options(simplePHENOTYPES.interference = NULL)
  set.seed(1)
  itf <- simplePHENOTYPES:::.check_interference(NULL)
  d <- simplePHENOTYPES:::.draw_meiosis(list(c(0, 1.5)), 20000, itf)
  expect_equal(mean(d$counts), 1.5, tolerance = 0.03)
  # interference: under-dispersed relative to Poisson
  expect_lt(stats::var(d$counts) / mean(d$counts), 0.9)
})

# ---- DECISION-048: a population's heritability travels with it ----------------

bd_f2 <- function() {
  pop <- bd_pop(10)
  f1 <- suppressMessages(cross(pop[1], pop[2], n = 1, seed = 1))
  selfcross(f1, n = 80, seed = 9)
}
bd_pheno <- function(p, ...) {
  simulate_phenotype(p, h2 = 0.5, seed = 7, ...) |> additive(n_qtn = 30)
}
# warnings of the frozen path are throttled by rlang; make them always fire here
quiet <- function(expr) suppressWarnings(expr)

test_that("marker data and populations without a trait refit (unchanged)", {
  f2 <- bd_f2()
  expect_null(population_trait(f2))
  a <- bd_pheno(f2)
  b <- bd_pheno(f2, refit = TRUE)
  expect_identical(a$pheno, b$pheno)
  expect_false(isTRUE(a$frozen))
  expect_error(simulate_phenotype(f2, refit = FALSE), "carries none")
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES", envir = environment())
  expect_error(simulate_phenotype(SNP55K_maize282_maf04, refit = FALSE),
               "carries none")
})

test_that("select_ind() stores the trait and crossing passes it on", {
  f2 <- bd_f2()
  sim <- bd_pheno(f2)
  top <- select_ind(sim, prop = 0.25)
  pt <- population_trait(top)
  expect_equal(pt$h2, 0.5)
  expect_equal(pt$var_e, 0.5)
  expect_identical(top[1:2]$trait, top$trait)
  expect_identical(selfcross(top[1], n = 3, seed = 1)$trait, top$trait)
  expect_identical(cross(top[1], top[2], n = 3, seed = 1)$trait, top$trait)
  expect_identical(double_haploid(top[1], n = 3, seed = 1)$trait, top$trait)
  plan <- data.frame(mother = top$ids[1:2], father = top$ids[2:3], n = 2)
  expect_identical(mate(plan, top, seed = 1)$trait, top$trait)
  expect_identical(c(top[1:2], top[3:4])$trait, top$trait)
  # parents with different traits (or none) give progeny without one
  expect_null(cross(top[1], f2[1], n = 2, seed = 1)$trait)
  expect_warning(c(top[1:2], f2[1:2]), "same trait")
})

test_that("refit = FALSE reuses loci, effects and scale, and holds V_E", {
  f2 <- bd_f2()
  s0 <- bd_pheno(f2)
  top <- select_ind(s0, prop = 0.25)
  f3 <- selfcross(top[1], n = 60, seed = 3)
  s3 <- simulate_phenotype(f3, seed = 4)
  expect_true(isTRUE(s3$frozen))
  # the genetic value is the base template's genotypic value up to a constant
  tm <- template_effects(s0)
  g3 <- genotypic_value(f3, tm$qtn, tm$a, tm$d)
  gv <- genetic_values(s3)[, 1]
  expect_lt(stats::sd(gv - g3[names(gv)]), 1e-10)
  # the same template is recovered from the frozen simulation
  expect_equal(template_effects(s3)$a, tm$a)
  # residual variance is the base one (exact-variance draw), not refitted
  y <- s3$pheno$value
  expect_equal(stats::var(y - gv), 0.5, tolerance = 1e-10)
  # the realized budget closes to 1 (with its covariance row)
  vb <- s3$var_budget
  expect_true("covariance" %in% vb$component)
  expect_equal(sum(vb$prop), 1, tolerance = 1e-12)
  # the frozen trait survives selection again
  expect_identical(select_ind(s3, prop = 0.5)$trait, top$trait)
})

test_that("h2 given for a population with a trait warns; refit = TRUE refits", {
  f2 <- bd_f2()
  top <- select_ind(bd_pheno(f2), prop = 0.25)
  f3 <- selfcross(top[1], n = 40, seed = 3)
  expect_warning(simulate_phenotype(f3, h2 = 0.3, seed = 4),
                 "heritability is already set.*refit = TRUE")
  # re-specifying a layer the trait has is silent; a new layer type warns
  expect_no_warning(simulate_phenotype(f3, seed = 4) |> additive(n_qtn = 5))
  expect_warning(simulate_phenotype(f3, seed = 4) |> dominance(prop = 0.1),
                 "dominance\\(\\) is ignored")
  expect_no_warning(simulate_phenotype(f3, seed = 4))
  # the population's own h2 is not a conflict
  expect_no_warning(simulate_phenotype(f3, h2 = 0.5, seed = 4))
  expect_no_warning(simulate_phenotype(f3, h2 = 0.5, n_qtn = 30, seed = 4))
  s <- expect_no_warning(simulate_phenotype(f3, h2 = 0.3, n_qtn = 5, seed = 4,
                                            refit = TRUE))
  expect_false(isTRUE(s$frozen))
  # the refitted simulation defines a new trait for the next selection
  expect_equal(population_trait(select_ind(s, prop = 0.5))$h2, 0.3)
})

test_that("epistasis traits are stored with base-population locus centering", {
  f2 <- bd_f2()
  epi <- function(p, ...) {
    simulate_phenotype(p, h2 = 0.5, seed = 7, ...) |>
      additive(prop = 0.3, n_qtn = 10) |> epistasis(prop = 0.2, n_pairs = 3)
  }
  s0 <- epi(f2)
  top <- select_ind(s0, prop = 0.25)
  expect_false(is.null(top$trait))
  ly <- top$trait$layers[[2L]]
  expect_identical(ly$type, "epistasis")
  expect_identical(dim(ly$frozen_locus_center[[1L]]), c(3L, 2L))
  # applied to the base population itself, the stored trait reproduces it
  base <- f2
  base$trait <- top$trait
  sb <- simulate_phenotype(base, seed = 7)
  expect_equal(genetic_values(sb), genetic_values(s0), tolerance = 1e-10)
  # in a descendant the epistatic value is a fixed function of the genotype:
  # recompute it by hand from the base centers and the stored scale
  f3 <- selfcross(top[1], n = 60, seed = 3)
  s3 <- simulate_phenotype(f3, seed = 4)
  d3 <- dosages(f3)
  q <- ly$qtn[[1L]]
  ctr <- ly$frozen_locus_center[[1L]]
  raw <- rowSums(sapply(seq_len(nrow(q)), function(p) {
    (d3[q[p, 1], ] - ctr[p, 1]) * (d3[q[p, 2], ] - ctr[p, 2]) * ly$effect[[1L]][p]
  }))
  comp <- (raw - ly$frozen_center[1]) / ly$frozen_sd[1] * sqrt(0.2)
  add <- simplePHENOTYPES:::.scaled_component(s3$layers[[1L]], s3, 1L)
  expect_equal(unname(genetic_values(s3)[, 1]), unname(add + comp),
               tolerance = 1e-10)
})

test_that("traits without a fixed form (vqtl) are not stored", {
  f2 <- bd_f2()
  vq <- suppressWarnings(simulate_phenotype(f2, h2 = 0.5, seed = 7) |>
    additive(prop = 0.5, n_qtn = 10) |> vqtl(prop = 0.2))
  expect_null(select_ind(vq, prop = 0.25)$trait)
})

test_that("schemes: later generations hold V_E so h2 falls", {
  f2 <- bd_f2()
  out <- quiet(recurrent_selection(f2, bd_pheno, cycles = 5, n_parents = 10,
                                   progeny_per_cross = 8, seed = 2))
  expect_equal(population_trait(out)$var_e, 0.5)
  tm <- template_effects(bd_pheno(f2))
  vg <- function(p) stats::var(genotypic_value(p, tm$qtn, tm$a, tm$d))
  expect_lt(vg(out), vg(f2))
  # refit = TRUE in the callback is the previous rule
  r <- recurrent_selection(f2, function(p) bd_pheno(p, refit = TRUE), cycles = 2,
                           n_parents = 10, progeny_per_cross = 8, seed = 2)
  expect_s3_class(r, "Population")
})
