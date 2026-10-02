# Gap closure: the crossover `interference` option of cross_usefulness() is
# checked for the "dh" and "selfcross" schemes by the realized crossover-count
# dispersion (about 1 for the default Poisson process, about 1/2 at nu = 20),
# not only by counting calls to the sampler (test-fix6-propagation.R covers the
# forwarding and the "cross" scheme).

# ---- fixtures ---------------------------------------------------------------

# Two inbred parent pairs on one 2-Morgan chromosome (101 markers, 2 cM apart):
# the "ones" parents carry allele 1 and the "zeros" parents allele 0 at every
# marker, so each F1 is a phase-known heterozygote (cis = 1, trans = 0) and the
# haplotype switches along a strand of a DH or F2 individual are the crossovers
# of exactly one F1 meiosis. A short second chromosome of random homozygous
# markers (identical across parents) carries the QTN.
usefulness_parents <- function(seed = 11) {
  nm1 <- 101L
  nm2 <- 6L
  map <- data.frame(
    snp = c(paste0("a", seq_len(nm1)), paste0("b", seq_len(nm2))),
    chr = c(rep(1L, nm1), rep(2L, nm2)),
    pos = c(seq_len(nm1) * 1000L, seq_len(nm2) * 1000L),
    cm = c((seq_len(nm1) - 1) * 2, (seq_len(nm2) - 1) * 10))
  old <- .Random.seed_safe()
  on.exit(.restore_seed(old), add = TRUE)
  set.seed(seed)
  hom <- matrix(rbinom(nm2 * 4L, 1, 0.5), nm2, 4L)
  one <- rbind(matrix(1, nm1, 4L), hom)
  zero <- rbind(matrix(0, nm1, 4L), hom)
  g <- cbind(one[, 1:2], zero[, 1:2])          # parents 1, 2 = ones; 3, 4 = zeros
  ids <- paste0("p", 1:4)
  dimnames(g) <- list(map$snp, ids)
  population_from_haplotypes(g, g, map, ids = ids)
}

# (constant loci among the passed QTNs warn; irrelevant to what is tested)
parents_sim <- function(pop) {
  suppressWarnings(simulate_phenotype(pop, h2 = 0.5, seed = 7) |>
    additive(qtn = paste0("b", 1:6),
             effect = c(1, -0.5, 0.8, 0.3, -0.9, 0.6)))
}

PAIRS <- rbind(c(1, 3), c(2, 4))               # each a ones x zeros pair
IF20 <- list(nu = 20, p = 0)

# Crossover counts per strand on chromosome 1 (haplotype switches along the 101
# markers); `both = FALSE` uses the cis strand only (a DH has identical strands).
strand_counts <- function(pop, both = TRUE) {
  h <- haplotypes(pop)
  m <- if (both) cbind(h$cis, h$trans) else h$cis
  as.numeric(colSums(abs(diff(m[1:101, , drop = FALSE]))))
}

# Index of dispersion (variance / mean) and its normal-approximation SE.
disp <- function(x) stats::var(x) / mean(x)
disp_se <- function(x) sqrt(2 / (length(x) - 1))

# The progeny families cross_usefulness() builds, collected from .make_family().
grab_families <- function(sim, scheme, interference, n_progeny, seed) {
  orig <- simplePHENOTYPES:::.make_family
  fams <- list()
  testthat::local_mocked_bindings(
    .make_family = function(...) {
      f <- orig(...)
      fams[[length(fams) + 1L]] <<- f
      f
    }, .package = "simplePHENOTYPES", .env = parent.frame())
  cross_usefulness(sim, PAIRS, scheme, n_progeny = n_progeny, generations = 1,
                   seed = seed, interference = interference)
  do.call(c, fams)
}

check_dispersion <- function(base, ifx, label) {
  se <- disp_se(base)
  # Poisson (default): index of dispersion 1 within 4 SE
  expect_lt(abs(disp(base) - 1), 4 * se, label = paste(label, "default"))
  # strong interference (nu = 20): about 1/2, clearly below 1
  expect_lt(disp(ifx), 1 - 4 * disp_se(ifx), label = paste(label, "nu = 20"))
  expect_lt(disp(ifx), 0.75, label = paste(label, "nu = 20, absolute"))
  # interference redistributes crossovers but keeps the map length
  expect_lt(abs(mean(ifx) - mean(base)),
            5 * sqrt(var(base) / length(base) + var(ifx) / length(ifx)),
            label = paste(label, "mean"))
}

# ---- the tests --------------------------------------------------------------

test_that("cross_usefulness 'dh' families have under-dispersed crossovers with interference", {
  skip_on_cran()
  pop <- usefulness_parents()
  sim <- parents_sim(pop)
  # 2 pairs x 300 DH lines = 600 independent F1 gametes
  base <- strand_counts(grab_families(sim, "dh", NULL, 300, seed = 6), both = FALSE)
  ifx <- strand_counts(grab_families(sim, "dh", IF20, 300, seed = 6), both = FALSE)
  expect_length(base, 600L)
  check_dispersion(base, ifx, "dh")
})

test_that("cross_usefulness 'selfcross' families have under-dispersed crossovers with interference", {
  skip_on_cran()
  pop <- usefulness_parents()
  sim <- parents_sim(pop)
  # one selfing generation (F2): 2 pairs x 300 plants x 2 gametes = 1200 strands
  base <- strand_counts(grab_families(sim, "selfcross", NULL, 300, seed = 6))
  ifx <- strand_counts(grab_families(sim, "selfcross", IF20, 300, seed = 6))
  expect_length(base, 1200L)
  check_dispersion(base, ifx, "selfcross")
})

test_that("the interference setting moves the usefulness of dh and selfcross, reproducibly", {
  skip_on_cran()
  pop <- usefulness_parents()
  # QTN on the linked chromosome so the family variance is a linkage statistic
  sim <- suppressWarnings(simulate_phenotype(pop, h2 = 0.5, seed = 7) |>
    additive(qtn = c(paste0("a", c(10, 30, 50, 70, 90)), "b1", "b2"),
             effect = c(1, -0.5, 0.8, 0.3, -0.9, 0.2, -0.2)))
  for (sc in c("dh", "selfcross")) {
    u0 <- cross_usefulness(sim, PAIRS, sc, n_progeny = 100, generations = 1,
                           seed = 3)
    u1 <- cross_usefulness(sim, PAIRS, sc, n_progeny = 100, generations = 1,
                           seed = 3, interference = IF20)
    u2 <- cross_usefulness(sim, PAIRS, sc, n_progeny = 100, generations = 1,
                           seed = 3, interference = IF20)
    expect_identical(u1, u2, info = sc)
    expect_false(isTRUE(all.equal(u0$sd, u1$sd)), info = sc)
  }
})
