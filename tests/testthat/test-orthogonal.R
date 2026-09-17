# test-orthogonal.R
#
# The orthogonal genotypic model: additive(orthogonal = TRUE, a =, d =).
# Per-locus a and d give an additive/dominance partition whose variances emerge
# (Cov(A, D) = 0 in expectation under random mating -> per-locus HWE; LD is fine,
# nonrandom multilocus association is not; a finite or structured sample carries a
# real covariance, reported as add_dom_cov), the degree of dominance d/abs(a) is
# meaningful, and the breeding value picks up the dominance-induced average effect.

data("SNP55K_maize282_maf04")

# A segregating F2 (het at ~50% of the differing loci) so a dominance deviation
# is identifiable.
.f2 <- function(n = 80, seed = 2) {
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
  suppressMessages(selfcross(cross(pop[1], pop[2], n = 1, seed = 1),
                             n = n, seed = seed))
}

test_that("orthogonal with d = 0 reduces to a standard additive layer", {
  og <- suppressMessages(
    simulate_phenotype(.f2(80), h2 = 0.5, seed = 3) |>
      additive(orthogonal = TRUE, a = 0.5, d = 0, n_qtn = 10))
  vb <- og$var_budget
  expect_lt(vb$prop[vb$component == "dominance"], 1e-8)   # no dominance variance
  # a purely additive genotypic value: breeding value == total genetic value
  bv <- simplePHENOTYPES:::.breeding_value_matrix(og, 1L)[, 1]
  gv <- genetic_values(og)[, 1]
  expect_gt(stats::cor(bv, gv), 0.999)
})

test_that("orthogonal d > 0 yields an emergent dominance component, orthogonal to additive", {
  og <- suppressMessages(
    simulate_phenotype(.f2(80), h2 = 0.6, seed = 3) |>
      additive(orthogonal = TRUE, a = 0.5, d = 0.5, n_qtn = 12))
  vb  <- og$var_budget
  add <- vb$prop[vb$component == "additive"]
  dom <- vb$prop[vb$component == "dominance"]
  cov <- vb$prop[vb$component == "add_dom_cov"]
  expect_gt(dom, 0)                                # dominance variance emerges
  # additive + dominance + additive-by-dominance covariance close the layer prop
  expect_lt(abs(add + dom + cov - 0.6), 1e-6)
})

test_that("degree of dominance drives the emergent additive:dominance ratio", {
  dom_share <- function(dd) {
    og <- suppressMessages(
      simulate_phenotype(.f2(80), h2 = 0.6, seed = 3) |>
        additive(orthogonal = TRUE, a = 0.5, d = dd, n_qtn = 12))
    og$var_budget$prop[og$var_budget$component == "dominance"]
  }
  expect_lt(dom_share(0.2), dom_share(1.0))       # more dominance -> larger Vd
})

test_that("budget rows are the realized Var(A)/Var(D)/2Cov(A,D), not normalized allocations", {
  og <- suppressMessages(
    simulate_phenotype(.f2(80), h2 = 0.6, seed = 3) |>
      additive(orthogonal = TRUE, a = 0.5, d = 0.5, n_qtn = 12))
  vb   <- og$var_budget
  add  <- vb$prop[vb$component == "additive"]
  dom  <- vb$prop[vb$component == "dominance"]
  cov  <- vb$prop[vb$component == "add_dom_cov"]
  g    <- genetic_values(og)[, 1]
  bv   <- simplePHENOTYPES:::.breeding_value_matrix(og, 1L)[, 1]
  vg   <- stats::var(g); vA <- stats::var(bv); vD <- stats::var(g - bv)
  prop <- 0.6
  # each row is the *realized* variance as a fraction of Var(g), scaled to prop
  expect_lt(abs(add - prop * vA / vg), 1e-6)
  expect_lt(abs(dom - prop * vD / vg), 1e-6)
  # the finite F2 sample carries a real Cov(A, D) here (not ~ 0)
  expect_gt(abs(cov), 1e-3)
  expect_lt(abs(add + dom + cov - prop), 1e-6)     # covariance row closes the budget
})

test_that("print() shows the emergent orthogonal split, not the whole layer as additive", {
  og <- suppressMessages(
    simulate_phenotype(.f2(60), h2 = 0.6, seed = 3) |>
      additive(orthogonal = TRUE, a = 0.5, d = 0.5, n_qtn = 10))
  out <- utils::capture.output(print(og))
  expect_true(any(grepl("orthogonal a/d model", out)))   # labelled as orthogonal
  expect_true(any(grepl("add_dom_cov", out)))            # covariance row shown
  # the additive line must be the emergent share, not the full layer prop (0.6)
  add_line <- grep("^\\s*additive", out, value = TRUE)
  add_val  <- as.numeric(sub(".*additive\\s+([-0-9.]+).*", "\\1", add_line))
  expect_lt(add_val, 0.6)
})

test_that("the variance-budget plot matrix sums duplicate component rows", {
  og <- suppressMessages(
    simulate_phenotype(.f2(80), h2 = 0.5, seed = 3) |>
      additive(prop = 0.2, n_qtn = 5) |>
      additive(orthogonal = TRUE, a = 0.5, d = 0.5, prop = 0.3, n_qtn = 6))
  vb <- og$var_budget
  expect_gt(sum(vb$component == "additive"), 1)   # standard + orthogonal rows
  m  <- simplePHENOTYPES:::.var_budget_matrix(vb)
  # each trait column must reproduce that trait's full budget, not just one row
  per_trait <- tapply(vb$prop, vb$trait, sum)
  expect_equal(unname(colSums(m)), unname(per_trait[colnames(m)]),
               tolerance = 1e-9)
})

test_that("the variance plot renders an orthogonal budget with a signed covariance row", {
  og <- suppressMessages(
    simulate_phenotype(.f2(80), h2 = 0.6, seed = 3) |>
      additive(orthogonal = TRUE, a = 0.5, d = 0.5, n_qtn = 12))
  expect_true("add_dom_cov" %in% og$var_budget$component)
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)
  # data-driven ylim must accommodate the (possibly signed) covariance row
  expect_error(simplePHENOTYPES:::.plot_variance(og), NA)
})

test_that("qtn_table() carries d and counts the dominance deviation in var_explained", {
  og <- suppressMessages(
    simulate_phenotype(.f2(80), h2 = 1, seed = 3) |>
      additive(orthogonal = TRUE, a = 0.5, d = 0.4, n_qtn = 6))
  tab <- qtn_table(og)
  expect_true("d" %in% names(tab))
  expect_true(all(abs(tab$d - 0.4) < 1e-9))               # per-locus d surfaced
  expect_true(all(tab$var_explained > 0))                 # dominance loci count

  # a pure-dominance layer (a = 0) still explains variance
  pd <- suppressMessages(
    simulate_phenotype(.f2(80), h2 = 1, seed = 3) |>
      additive(orthogonal = TRUE, a = 0, d = 0.5, n_qtn = 4))
  expect_true(all(qtn_table(pd)$var_explained > 0))
})

test_that("cross_usefulness() scores an orthogonal layer on its average effect, not bare a", {
  og <- suppressMessages(
    simulate_phenotype(.f2(20), h2 = 0.6, seed = 3) |>
      additive(orthogonal = TRUE, a = 0, d = 0.6, n_qtn = 12))
  u <- suppressMessages(
    cross_usefulness(og, pairs = rbind(c(1, 2), c(3, 4)),
                     scheme = "dh", n_progeny = 30, seed = 2))
  # a = 0 but d > 0: the average effect alpha = d(1 - 2p) is non-zero away from
  # p = 0.5, so the prediction must not collapse to zero (the U2 defect, which
  # scored on bare a and returned mean = sd = usefulness = 0).
  expect_true(any(u$sd > 0))
  expect_true(any(u$usefulness != 0))
})

test_that("a non-zero d is rejected per locus, not collectively over the set", {
  f2 <- .f2(80)
  ph <- suppressMessages(simulate_phenotype(f2, h2 = 0.5, seed = 3))
  D  <- dosages(f2)                                  # markers x individuals
  hc <- rowSums(D == 0, na.rm = TRUE)
  hetless <- rownames(D)[hc == 0][1]
  segreg  <- rownames(D)[hc > 0][1]
  skip_if(is.na(hetless) || is.na(segreg))
  # d != 0 on the hetless locus must error even though the other locus segregates
  expect_error(
    additive(ph, orthogonal = TRUE, a = 0.5, qtn = c(hetless, segreg),
             d = c(0.5, 0)),
    "heterozyg")
  # d = 0 on the hetless locus (d != 0 only where there are hets) is fine
  expect_error(
    suppressMessages(
      additive(ph, orthogonal = TRUE, a = 0.5, qtn = c(hetless, segreg),
               d = c(0, 0.5))),
    NA)
})

test_that("the new orthogonal arguments do not shift the positional API", {
  f2 <- .f2(40)
  base <- function() suppressMessages(simulate_phenotype(f2, h2 = 0.5, seed = 3))
  # slots 6 and 7 must still bind phase and dist, as before the feature
  og <- suppressMessages(
    additive(base(), 0.5, 3, NULL, 0.5, "repulsion", "geometric"))
  expect_false(isTRUE(og$layers[[1]]$orthogonal))
  expect_identical(og$layers[[1]]$phase, "repulsion")
})

test_that("the orthogonal genotypic model validates its combinations", {
  f2 <- .f2(40)
  base <- function() suppressMessages(simulate_phenotype(f2, h2 = 0.5, seed = 3))

  # a/d require orthogonal = TRUE; effect is rejected under orthogonal
  expect_error(additive(base(), prop = 0.5, n_qtn = 3, a = 0.5),
               "orthogonal = TRUE")
  expect_error(additive(base(), orthogonal = TRUE, n_qtn = 3, effect = 0.5),
               "via `a`")

  # an orthogonal layer already models dominance -> no separate dominance()
  og <- suppressMessages(
    additive(base(), orthogonal = TRUE, a = 0.5, d = 0.2, n_qtn = 5))
  expect_error(dominance(og, prop = 0.1), "orthogonal")

  # correlation-controlling architectures cannot honour fixed per-locus effects
  pl <- suppressMessages(
    simulate_phenotype(f2, architecture = "pleiotropy", n_traits = 2,
                       cor = 0.3, h2 = 0.5, seed = 3))
  expect_error(additive(pl, orthogonal = TRUE, a = 0.5, n_qtn = 5),
               "cannot honour")

  # vary_qtn is unsupported (effects are fixed)
  vq <- suppressMessages(
    simulate_phenotype(f2, h2 = 0.5, seed = 3, n_reps = 2, vary_qtn = TRUE))
  expect_error(additive(vq, orthogonal = TRUE, a = 0.5, n_qtn = 5), "vary_qtn")

  # d > 0 needs heterozygotes: a hetless locus cannot carry a dominance deviation
  inb <- suppressMessages(simulate_phenotype(SNP55K_maize282_maf04, h2 = 0.5,
                                             seed = 3))
  expect_error(
    additive(inb, orthogonal = TRUE, a = 0.5, d = 0.5, qtn = "ss196419692"),
    "heterozygous")
})
