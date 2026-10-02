# Round-6 propagation: the crossover `interference` option of cross() reaches
# every function that runs meiosis -- single_seed_descent(), bulk(), pedigree(),
# recurrent_selection(), cross_usefulness(), combining_ability(method =
# "simulated"), progeny_test() -- so a whole scheme can use one meiosis model.
# `interference = NULL` (the default) must stay bit-identical to before.

# ---- fixtures ---------------------------------------------------------------

# One 2-Morgan chromosome (101 markers, 2 cM apart) on which every founder is a
# phase-known heterozygote (cis = 1, trans = 0), so each strand of a progeny of
# ONE meiosis is a gamete whose haplotype switches are its crossovers; plus a
# short chromosome of random homozygous markers that carries the QTN.
het_founders <- function(n = 8L, seed = 11) {
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
  hom <- matrix(rbinom(nm2 * n, 1, 0.5), nm2, n)
  cis <- rbind(matrix(1, nm1, n), hom)
  trans <- rbind(matrix(0, nm1, n), hom)
  ids <- paste0("f", seq_len(n))
  dimnames(cis) <- dimnames(trans) <- list(map$snp, ids)
  population_from_haplotypes(cis, trans, map, ids = ids)
}

# (the all-heterozygous F1 parents make some of these loci constant: the
# passed-QTN warning is expected and irrelevant to what is tested)
qtn_pheno <- function(p) {
  suppressWarnings(simulate_phenotype(p, h2 = 0.5, seed = 7) |>
    additive(qtn = paste0("b", 1:6),
             effect = c(1, -0.5, 0.8, 0.3, -0.9, 0.6)))
}

# Crossover counts per strand on chromosome 1: the number of haplotype switches
# along the 101 markers (2 cM apart; a double crossover inside one 2 cM interval
# is negligible under Poisson and impossible under strong interference).
strand_counts <- function(pop) {
  h <- haplotypes(pop)
  m <- cbind(h$cis, h$trans)[1:101, , drop = FALSE]
  as.numeric(colSums(abs(diff(m))))
}

# Index of dispersion (variance / mean) and its large-sample SE under a normal
# approximation: sqrt(2 / (n - 1)).
disp <- function(x) stats::var(x) / mean(x)
disp_se <- function(x) sqrt(2 / (length(x) - 1))

IF <- list(nu = 20, p = 0)

# Count every call to the interference sampler.
count_interference_draws <- function(code) {
  env <- new.env()
  env$n <- 0L
  orig <- simplePHENOTYPES:::.draw_meiosis_interference
  testthat::local_mocked_bindings(
    .draw_meiosis_interference = function(...) {
      env$n <- env$n + 1L
      orig(...)
    },
    .package = "simplePHENOTYPES")
  force(code)
  env$n
}

# ---- default is bit-identical to before -------------------------------------
# Checksums computed at HEAD (before the option existed) from the shipped maize
# panel; any change to the stream or the outputs of a default call breaks them.

ck_pop <- function(x) {
  sum(as.numeric(x$cis) * (seq_along(x$cis) %% 97)) +
    7 * sum(as.numeric(x$trans) * (seq_along(x$trans) %% 89))
}

test_that("default (interference = NULL) reproduces the pre-option results exactly", {
  skip_on_cran()
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:12)
  f1 <- cross(pop[1], pop[2], n = 6, seed = 1)
  pheno <- function(p) simulate_phenotype(p, h2 = 0.5, seed = 7) |>
    additive(n_qtn = 30)
  q <- c("ss196442916", "ss196439337", "ss196480535")

  expect_equal(ck_pop(single_seed_descent(f1, 3, seed = 2)), 11640944)
  expect_equal(ck_pop(bulk(f1, 3, n = 6, seed = 2)), 11686411)
  expect_equal(ck_pop(pedigree(f1, pheno, generations = 2, prop = 0.5,
                               seed = 2)), 11669757)
  expect_equal(ck_pop(recurrent_selection(pop, pheno, cycles = 2, n_parents = 4,
                                          progeny_per_cross = 4, seed = 2)),
               30146751)
  # ambient (unseeded) stream, and the stream is left in the same place
  set.seed(5)
  expect_equal(ck_pop(single_seed_descent(f1, 2)), 11695964)
  expect_equal(runif(1), 0.4526463, tolerance = 1e-6)
  set.seed(5)
  expect_equal(ck_pop(recurrent_selection(pop, pheno, cycles = 2, n_parents = 4,
                                          progeny_per_cross = 4)), 30710666)
  expect_equal(runif(1), 0.7214632, tolerance = 1e-6)

  sim <- simulate_phenotype(pop, h2 = 0.5, seed = 1) |> additive(n_qtn = 30)
  pr <- cbind(1:3, 4:6)
  expect_equal(cross_usefulness(sim, pr, "dh", n_progeny = 10,
                                seed = 2)$usefulness[1], 1.363007953988142)
  expect_equal(cross_usefulness(sim, pr, "selfcross", n_progeny = 10,
                                generations = 3, seed = 2)$usefulness[1],
               1.251586842836279)
  expect_equal(cross_usefulness(sim, pr, "cross", n_progeny = 10,
                                seed = 2)$usefulness[1], 0.545905465495001)
  set.seed(5)
  invisible(cross_usefulness(sim, cbind(1:2, 3:4), "selfcross", n_progeny = 5,
                             generations = 2))
  expect_equal(runif(1), 0.1826721, tolerance = 1e-6)

  ca <- combining_ability(pop[1:4], pop[5:7], qtn = q, a = c(1, .5, .25),
                          d = .3, method = "simulated", n_progeny = 5, seed = 3)
  expect_equal(ck_pop(ca$progeny), 115070440)
  expect_equal(sum(ca$cross_means), 6.54)
  cd <- combining_ability(pop[1:4], qtn = q, a = c(1, .5, .25),
                          design = "diallel", method = "simulated",
                          n_progeny = 5, seed = 3)
  expect_equal(ck_pop(cd$progeny), 59054690)
  pt <- progeny_test(pop[1:3], pop[4:12], qtn = q, a = c(1, .5, .25),
                     n_progeny = 4, h2 = .3, seed = 1)
  expect_equal(ck_pop(attr(pt, "progeny")), 22933730)
  expect_equal(sum(pt$progeny_mean), 1.20569907909044)
})

test_that("explicit interference = NULL is the same call as omitting it", {
  fx <- het_founders()
  f1 <- cross(fx[1], fx[2], n = 6, seed = 1)
  expect_identical(
    dosages(single_seed_descent(f1, 2, seed = 4)),
    dosages(single_seed_descent(f1, 2, seed = 4, interference = NULL)))
  expect_identical(
    dosages(bulk(f1, 2, seed = 4)),
    dosages(bulk(f1, 2, seed = 4, interference = NULL)))
})

# ---- every wrapper reaches the interference sampler --------------------------

test_that("interference is forwarded to every crossing call of every wrapper", {
  skip_on_cran()
  fx <- het_founders()
  f1 <- cross(fx[1], fx[2], n = 6, seed = 1)
  q <- paste0("b", 1:3)
  a <- c(1, 0.5, 0.25)

  n_with <- function(expr_fn) {
    list(none = count_interference_draws(expr_fn(NULL)),
         ifx = count_interference_draws(expr_fn(IF)))
  }
  chk <- function(r) {
    expect_identical(r$none, 0L)       # the default never enters the sampler
    expect_gt(r$ifx, 0L)
  }
  chk(n_with(function(i) single_seed_descent(f1, 2, seed = 1, interference = i)))
  chk(n_with(function(i) bulk(f1, 2, seed = 1, interference = i)))
  chk(n_with(function(i) pedigree(fx, qtn_pheno, generations = 2, prop = 0.5,
                                  seed = 1, interference = i)))
  chk(n_with(function(i) recurrent_selection(fx, qtn_pheno, cycles = 2,
                                             n_parents = 3, progeny_per_cross = 4,
                                             seed = 1, interference = i)))
  sim <- qtn_pheno(fx)
  for (sc in c("dh", "selfcross", "cross")) {
    chk(n_with(function(i) cross_usefulness(sim, cbind(1:2, 3:4), sc,
                                            n_progeny = 5, generations = 3,
                                            seed = 1, interference = i)))
  }
  chk(n_with(function(i) combining_ability(fx[1:3], fx[4:5], qtn = q, a = a,
                                           method = "simulated", n_progeny = 3,
                                           seed = 1, interference = i)))
  chk(n_with(function(i) combining_ability(fx[1:4], qtn = q, a = a,
                                           design = "diallel",
                                           method = "simulated", n_progeny = 3,
                                           seed = 1, interference = i)))
  chk(n_with(function(i) progeny_test(fx[1:2], fx[3:8], qtn = q, a = a,
                                      n_progeny = 3, seed = 1,
                                      interference = i)))
})

test_that("interference changes the stream and is seed-reproducible", {
  skip_on_cran()
  fx <- het_founders()
  f1 <- cross(fx[1], fx[2], n = 6, seed = 1)
  q <- paste0("b", 1:3)
  a <- c(1, 0.5, 0.25)
  runs <- list(
    ssd = function(i) single_seed_descent(f1, 2, seed = 3, interference = i),
    bulk = function(i) bulk(f1, 2, seed = 3, interference = i),
    ped = function(i) pedigree(fx, qtn_pheno, generations = 2, prop = 0.5,
                               seed = 3, interference = i),
    rec = function(i) recurrent_selection(fx, qtn_pheno, cycles = 2,
                                          n_parents = 4, progeny_per_cross = 10,
                                          seed = 3, interference = i),
    ca = function(i) combining_ability(fx[1:3], fx[4:5], qtn = q, a = a,
                                       method = "simulated", n_progeny = 4,
                                       seed = 3, interference = i)$progeny,
    pt = function(i) attr(progeny_test(fx[1:2], fx[3:8], qtn = q, a = a,
                                       n_progeny = 4, seed = 3,
                                       interference = i), "progeny"))
  for (nm in names(runs)) {
    base <- runs[[nm]](NULL)
    one <- runs[[nm]](IF)
    two <- runs[[nm]](IF)
    expect_identical(one$cis, two$cis, info = nm)       # reproducible
    expect_identical(one$trans, two$trans, info = nm)
    expect_false(identical(one$cis, base$cis) && identical(one$trans, base$trans),
                 info = nm)                             # a different stream
  }
  sim <- qtn_pheno(fx)
  # QTN on the heterozygous chromosome, so the family variance is a linkage
  # statistic that the meiosis model moves
  simA <- suppressWarnings(simulate_phenotype(fx, h2 = 0.5, seed = 7) |>
    additive(qtn = c(paste0("a", c(10, 30, 50, 70, 90)), "b1", "b2"),
             effect = c(1, -0.5, 0.8, 0.3, -0.9, 0.2, -0.2)))
  u0 <- cross_usefulness(simA, cbind(1:2, 3:4), "cross", n_progeny = 20, seed = 3)
  u1 <- cross_usefulness(simA, cbind(1:2, 3:4), "cross", n_progeny = 20, seed = 3,
                         interference = IF)
  u2 <- cross_usefulness(simA, cbind(1:2, 3:4), "cross", n_progeny = 20, seed = 3,
                         interference = IF)
  expect_identical(u1, u2)
  expect_false(isTRUE(all.equal(u0$sd, u1$sd)))
})

test_that("a seeded scheme with interference leaves the caller's RNG alone", {
  fx <- het_founders()
  f1 <- cross(fx[1], fx[2], n = 4, seed = 1)
  set.seed(8); a <- runif(1)
  set.seed(8)
  single_seed_descent(f1, 1, seed = 2, interference = IF)
  bulk(f1, 1, seed = 2, interference = IF)
  expect_equal(runif(1), a)
})

# ---- invalid values are rejected, naming the caller --------------------------

test_that("invalid interference is rejected up front with the caller's name", {
  fx <- het_founders()
  f1 <- cross(fx[1], fx[2], n = 4, seed = 1)
  q <- paste0("b", 1:3)
  a <- c(1, 0.5, 0.25)
  sim <- qtn_pheno(fx)
  bad <- list(nu = 0.5)
  calls <- list(
    single_seed_descent = function(i) single_seed_descent(f1, 1, interference = i),
    bulk = function(i) bulk(f1, 1, interference = i),
    pedigree = function(i) pedigree(fx, qtn_pheno, 1, interference = i),
    recurrent_selection = function(i) recurrent_selection(
      fx, qtn_pheno, 1, n_parents = 3, interference = i),
    cross_usefulness = function(i) cross_usefulness(sim, interference = i),
    combining_ability = function(i) combining_ability(
      fx[1:3], fx[4:5], qtn = q, a = a, method = "simulated", n_progeny = 2,
      interference = i),
    progeny_test = function(i) progeny_test(fx[1:2], fx[3:8], qtn = q, a = a,
                                            n_progeny = 2, interference = i))
  for (nm in names(calls)) {
    expect_error(calls[[nm]](bad), paste0("^", nm, "\\(\\): `interference\\$nu`"),
                 info = nm)
    expect_error(calls[[nm]](list(nu = 2, p = 2)), paste0("^", nm, "\\(\\)"),
                 info = nm)
    expect_error(calls[[nm]](2.6), "must be NULL or a list", info = nm)
    expect_error(calls[[nm]](list(p = 0.5)), "`nu` \\(required\\)", info = nm)
  }
  # a mistake is caught before any work (a failing callback is never reached)
  expect_error(pedigree(fx, function(p) stop("reached"), 1, interference = bad),
               "interference")
  # method = "expected" runs no meiosis: the option is refused, not ignored
  expect_error(combining_ability(fx[1:3], fx[4:5], qtn = q, a = a,
                                 interference = IF),
               "apply to method = \"simulated\" only")
})

# ---- the realized recombination reflects the interference -------------------

test_that("schemes with interference have under-dispersed crossover counts", {
  skip_on_cran()
  fx <- het_founders()
  q <- paste0("b", 1:3)
  a <- c(1, 0.5, 0.25)
  N <- 300L                              # progeny -> 600 strands (2 gametes each)
  f1 <- cross(fx[1], fx[2], n = N, seed = 1)   # selfed into one meiosis per gamete
  sim <- qtn_pheno(fx)

  check <- function(make, label) {
    base <- strand_counts(make(NULL))
    ifx <- strand_counts(make(IF))
    se <- disp_se(base)
    # Poisson: index of dispersion 1 (within 4 SE); strong interference
    # (nu = 20): about 1/2 (the 1/4-of-chromatids thinning leaves binomial noise)
    expect_lt(abs(disp(base) - 1), 4 * se, label = paste(label, "default"))
    expect_lt(disp(ifx), 1 - 4 * disp_se(ifx), label = paste(label, "nu = 20"))
    # the mean crossover count per gamete is the map length / 2 (1 per Morgan on
    # a 2 M map: one chiasma per Morgan, each kept by a gamete with prob. 1/2)
    expect_lt(abs(mean(ifx) - mean(base)),
              5 * sqrt(var(base) / length(base) + var(ifx) / length(ifx)),
              label = paste(label, "mean"))
  }

  # one selfing generation of N F1-like hets: each strand is one gamete
  het <- fx[rep(1:8, length.out = N)]
  check(function(i) single_seed_descent(het, 1, seed = 5, interference = i),
        "single_seed_descent")
  check(function(i) bulk(het, 1, n = N, seed = 5, interference = i), "bulk")
  check(function(i) pedigree(fx, qtn_pheno, generations = 1, n_select = 4,
                             pop_size = N, seed = 5, interference = i),
        "pedigree")
  check(function(i) recurrent_selection(fx, qtn_pheno, cycles = 1, n_parents = 4,
                                        n_crosses = 4,
                                        progeny_per_cross = N / 4, seed = 5,
                                        interference = i),
        "recurrent_selection")
  check(function(i) combining_ability(fx[1:4], fx[5:8], qtn = q, a = a,
                                      method = "simulated", n_progeny = 19,
                                      seed = 5, interference = i)$progeny,
        "combining_ability")
  big <- het_founders(n = 120L)
  check(function(i) attr(progeny_test(big[1:20], big[21:120], qtn = q, a = a,
                                      n_progeny = 15, seed = 5,
                                      interference = i), "progeny"),
        "progeny_test")
})

test_that("cross_usefulness simulates its families with the interference model", {
  skip_on_cran()
  fx <- het_founders()
  sim <- qtn_pheno(fx)
  orig <- simplePHENOTYPES:::.make_family
  grab <- function(i, scheme) {
    fams <- list()
    testthat::local_mocked_bindings(
      .make_family = function(...) {
        f <- orig(...)
        fams[[length(fams) + 1L]] <<- f
        f
      }, .package = "simplePHENOTYPES")
    cross_usefulness(sim, cbind(1:4, 5:8), scheme, n_progeny = 150,
                     generations = 1, seed = 6, interference = i)
    do.call(c, fams)
  }
  # "cross": every progeny strand is one meiosis of a het parent
  base <- strand_counts(grab(NULL, "cross"))
  ifx <- strand_counts(grab(IF, "cross"))
  expect_lt(abs(disp(base) - 1), 4 * disp_se(base))
  expect_lt(disp(ifx), 1 - 4 * disp_se(ifx))
})
